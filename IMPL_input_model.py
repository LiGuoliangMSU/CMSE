#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Interpolate 3D model onto mesh grids

The 3D model file: lat lon dep vp vs density Qp Qs
or
                     x  y  dep vp vs density Qp Qs
"""
import sys
import getopt
import numpy as np
from scipy import interpolate
from scipy.io import FortranFile
import pyproj
import multiprocessing
import timeit

#------ default parameters setting ------#
global flag_topo,ref_lon,ref_lat,ref_alt,ref_ellps
flag_topo=0
flag_cart=0
nslice=12

#=======END of defining control parameters========#


#--------- subrotinue for I/O of the code ------------#
def usage():
   print('''
-h or --help
-T or --topography="Parameter used to specify the name of the grided topography data; If this option is not actived, the earth topography are set as zero"
-M or --model="Parameter used to specify the name of the grided input velocity model; Data format: lat lon depth Vp Vs density Qkappa Qmu"
-B or --bin="Parameter used to specify the name of the path of the mesh bin file; A directory"
-N or --nprocessor="Parameter used to specify the number of processors used in model implementation"
-S or --slice="Number of mesh block to be implemented; A number"
-O or --outfile="File name of out path to store the results; A directroy"
-C or --cartesian='The data coordinate is Cartesian instead of Geodetic coordinate'
''')
   sys.exit(1)

if (len(sys.argv)==1):
   usage()

shortargs='hT:M:B:N:S:O:C'
longargs=['help','topography','model','bin','nprocessor','slice','outfile','cartesian']
opts, args = getopt.getopt( sys.argv[1:], shortargs, longargs )

for opt,val in opts:
   if opt in ('-h','--help'):
      usage()
      continue
   if opt in ('-T','--topography'):
      topo_txt=val
      flag_topo=1
      continue
   if opt in ('-M','--model'):
      model_file=val
      flag_model=1
      continue
   if opt in ('-B','--bin'):
      mesh_dir=val
      continue
   if opt in ('-N','--nprocessor'):
      nproc=int(val)
      continue
   if opt in ('-S','--slice'):
      Nslice=int(val)
      continue
   if opt in ('-O','--outfile'):
      out_dir=val
      continue
   if opt in ('-C','--cartesian'):
      flag_cart=1

print(nproc, Nslice)
#=========End of the subrotinue========#



#-------load mesh parameter file---------#
mesh_par_file = "mesh_par_file.py"
if sys.version_info < (3, ):
   raise Exception("need python3")
elif sys.version_info < (3, 5):
   import importlib
   spec =importlib.machinery.SourceFileLoader("mesh_par", mesh_par_file)
   par = spec.load_module()
else:
   import importlib.util
   spec = importlib.util.spec_from_file_location("mesh_par", mesh_par_file)
   par = importlib.util.module_from_spec(spec)
   spec.loader.exec_module(par)
#--- reference point
ref_lon = par.mesh_ref_lon
ref_lat = par.mesh_ref_lat
ref_alt = par.mesh_ref_alt
ref_ellps = par.mesh_ref_ellps
#===============END of loading parameter file=============#


field_list = [
('nspec','i4'), ('nglob','i4'), ('nspec_irregular','i4'), ('ibool','i4'),
('x','f4'), ('y','f4'), ('z','f4'),
]


#=read the topogaraphy file=#
if(flag_topo==1):
   with open(topo_txt, 'r') as f:
      lines = [ l.split() for l in f.readlines() ]

   topo_lons = np.array([ float(l[0]) for l in lines])
   topo_lats = np.array([ float(l[1]) for l in lines])
   topo_alts = np.array([ float(l[2]) for l in lines])


#====== read in 3-D reference model======#
with open(model_file, 'r') as f:
   lines = [ l.split() for l in f.readlines() if not(l.startswith('#')) ]

model_lat = np.array([float(x[0]) for x in lines])
model_lon = np.array([float(x[1]) for x in lines])
model_dep = np.array([float(x[2]) for x in lines])

npoints = len(model_dep)
model_values = np.zeros((npoints, 5))
model_values[:,0]= np.array([float(x[3]) for x in lines])*1000.0 # Vp, m/s
model_values[:,1]= np.array([float(x[4]) for x in lines])*1000.0 # Vs, m/s
model_values[:,2]= np.array([float(x[5]) for x in lines])*1000.0 # Rho, kg/m^3
model_values[:,3]= np.array([float(x[6]) for x in lines]) # Qp
model_values[:,4]= np.array([float(x[7]) for x in lines]) # Qs


#========interpolate the topography file to each velocity grid==============
if(flag_topo==1):
   altitude=interpolate.griddata((topo_lons, topo_lats), topo_alts, (model_lon, model_lat), method='cubic')
   model_dep = -altitude[:]/1000+model_dep[:]
else:
   model_dep = model_dep[:]


#===========interpolate each SEM slice=====================

#--- convert (lon,lat,alt) to ECEF
ecef = pyproj.Proj(proj='geocent', ellps=ref_ellps)
lla = pyproj.Proj(proj='latlong', ellps=ref_ellps)
x0, y0, z0 = pyproj.transform(lla, ecef, ref_lon, ref_lat, ref_alt)

cosla = np.cos(np.deg2rad(ref_lat))
sinla = np.sin(np.deg2rad(ref_lat))
coslo = np.cos(np.deg2rad(ref_lon))
sinlo = np.sin(np.deg2rad(ref_lon))


def Interpolate(islice):
   input_file = "%s/proc%06d_external_mesh.bin"%(mesh_dir, islice)
   print("Processing:", input_file)
 
   #--- read in SEM mesh
   mesh = {}
 
   with FortranFile(input_file, 'r') as f:
      for field in field_list:
         field_name = field[0]
         data_type = field[1]
         mesh[field_name] = f.read_ints(dtype=data_type)
 
   xgll = mesh['x'][mesh['ibool']-1]
   ygll = mesh['y'][mesh['ibool']-1]
   zgll = mesh['z'][mesh['ibool']-1]

   #  in the case the coordinate of the input model is under geodetic coordinate
   if (flag_cart == 0):
      #====transfering from ENU to ECEF=======#
      xx = x0 + (-sinlo)*xgll + (-sinla*coslo)*ygll + (cosla*coslo)*zgll
      yy = y0 + (coslo)*xgll + (-sinla*sinlo)*ygll + (cosla*sinlo)*zgll
      zz = z0 + (cosla)*ygll + (sinla)*zgll
   
      #====transfering from ECEF to latitude, longtitude, altitude====#
      gll_lon, gll_lat, gll_alt = pyproj.transform(ecef, lla, xx, yy, zz)


      #===Transfering from altitude to depth===#
      gll_depth_km = -1 * gll_alt/1000.0
 
      #======interpolate method=nearest,linear,cubic(1D/2D)======#
      gll_values = interpolate.griddata((model_lon, model_lat, model_dep), model_values, (gll_lon, gll_lat, gll_depth_km), method='nearest')

   else:
   #  if the coordinate of the input model are under the Cartesian coordinate, interpolation can be directly done with coordinate transform 
      gll_lon=xgll
      gll_lat=ygll
      gll_depth_km = -1 * zgll/1000.0
      gll_values = interpolate.griddata((model_lon, model_lat, model_dep), model_values, (gll_lon, gll_lat, gll_depth_km), method='nearest')
 
   #======output model gll file===========#
   out_file = "%s/proc%06d_vp.bin"%(out_dir, islice)
   with FortranFile(out_file, 'w') as f:
      f.write_record(np.array(gll_values[:,0], dtype='f4'))
 
   out_file = "%s/proc%06d_vs.bin"%(out_dir, islice)
   with FortranFile(out_file, 'w') as f:
      f.write_record(np.array(gll_values[:,1], dtype='f4'))
 
   out_file = "%s/proc%06d_rho.bin"%(out_dir, islice)
   with FortranFile(out_file, 'w') as f:
      f.write_record(np.array(gll_values[:,2], dtype='f4'))

   out_file = "%s/proc%06d_qkappa.bin"%(out_dir, islice)
   with FortranFile(out_file, 'w') as f:
      f.write_record(np.array(gll_values[:,3], dtype='f4'))

   out_file = "%s/proc%06d_qmu.bin"%(out_dir, islice)
   with FortranFile(out_file, 'w') as f:
      f.write_record(np.array(gll_values[:,4], dtype='f4'))


start = timeit.default_timer()
if __name__ == '__main__':
   items=[x for x in range(Nslice)]
   p = multiprocessing.Pool(nproc)
   b = p.map(Interpolate,items)
   p.close()
   p.join()
end = timeit.default_timer()
print("Total using time:",str(end-start),'s')
