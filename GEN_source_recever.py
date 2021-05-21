#!/usr/bin/env python
#---import libraries --#
""" 
Prepare the source and receiver file for simulation
"""
import sys
import getopt
import math
import numpy as np
from scipy import interpolate
from scipy.io import FortranFile
import pyproj
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import importlib

#====End of importing libraries====#

#-- default parameters setting ---#
global flag_topo,flag_stat,flag_rece,flag_forc,lag_cmts,ref_lon,ref_lat,ref_alt,ref_ellps
flag_topo=0
flag_stat=0
flag_rece=0
flag_forc=0
flag_cmts=0
topo_txt='topography.txt'
out_file="outfile.dat"
#===End of setting default parameter===#


#-------subroutine 1: FORCESOLUTION parameter------#
def readpara(pathtofile):
    para = {}
    with open(pathtofile, 'r') as f:
        for line in f.readlines():
            line = line.strip('\n')
            content = line.split('!')[0]
            if ':' in content:
                if '0.d0' in content:
                    para[content.split(':')[0]] = 0
                else:
                    para[content.split(':')[0]] = content.split(':')[1].strip()
            else:
                para['flag']=content
    return para
#=============End of subrotine readpara===========#


#--------------subroutine 2: GeodZNE_2_localUYX----------------#
def GeodZNE_2_localUYX(lon, lat,alti):
   ecef = pyproj.Proj(proj='geocent', ellps='WGS84')
   lla = pyproj.Proj(proj='latlong', ellps='WGS84')
   sta_x, sta_y, sta_z = pyproj.transform(lla, ecef, lon, lat, alti)
   x0, y0, z0 = pyproj.transform(lla, ecef, ref_lon, ref_lat, ref_alt)

   cosla = np.cos(np.deg2rad(ref_lat))
   sinla = np.sin(np.deg2rad(ref_lat))
   coslo = np.cos(np.deg2rad(ref_lon))
   sinlo = np.sin(np.deg2rad(ref_lon))

   grd_ee =       -sinlo*(sta_x-x0) +       coslo*(sta_y-y0)
   grd_nn = -sinla*coslo*(sta_x-x0) - sinla*sinlo*(sta_y-y0) + cosla*(sta_z-z0)
   grd_uu =  cosla*coslo*(sta_x-x0) + cosla*sinlo*(sta_y-y0) + sinla*(sta_z-z0)
   return grd_ee,grd_nn,grd_uu
#===============END of subroutine GeodZNE_2_localUYX================#


#-----------subroutine 3: Trans_Matrix--------------#
def Trans_Matrix(lon, lat, alti):
   ratio=1.0/180*math.pi
   sinrLO=math.sin(ref_lon*ratio)
   sinrLA=math.sin(ref_lat*ratio)
   cosrLO=math.cos(ref_lon*ratio)
   cosrLA=math.cos(ref_lat*ratio)
   #averaged radiu of the earth
   R=(6378137+6356752)/2.0

   #
   dX=np.array(np.zeros((3,1)))
   dY=np.array(np.zeros((3,1)))
   dZ=np.array(np.zeros((3,1)))

   sinlo=math.sin(lon*ratio)
   sinla=math.sin(lat*ratio)
   coslo=math.cos(lon*ratio)
   cosla=math.cos(lat*ratio)

   #calculate the normal directions of the geodetic coordinates
   #       Normal of longtitude or East   
   dX[0]= -(R+alti)*cosla*sinlo
   dY[0]=  (R+alti)*cosla*coslo
   dZ[0]=  0

   #       Normal of Latitude or North
   dX[1]= -(R+alti)*sinla*coslo
   dY[1]= -(R+alti)*sinla*sinlo
   dZ[1]=  (R+alti)*cosla

   #       Normal of altitude
   dX[2]= cosla*coslo
   dY[2]= cosla*sinlo
   dZ[2]= sinla

   x =         -sinrLO*dX   +          cosrLO*dY  +      0*dZ
   y =  -sinrLA*cosrLO*dX   +  -sinrLA*sinrLO*dY  + cosrLA*dZ
   z =   cosrLA*cosrLO*dX   +   cosrLA*sinrLO*dY  + sinrLA*dZ
   TransMat=np.hstack((x,y,z))
   TransMat[0][:]=TransMat[0][:]/(x[0]**2+y[0]**2+z[0]**2)**0.5
   TransMat[1][:]=TransMat[1][:]/(x[1]**2+y[1]**2+z[1]**2)**0.5
   TransMat[2][:]=TransMat[2][:]/(x[2]**2+y[2]**2+z[2]**2)**0.5
   return TransMat

#====================End of subroutine Trans_matrix=================#



#----subroutine 4: prepare points force-----#
def GEN_forcs(forcs_fn,out_file):
   para=readpara(forcs_fn)   
   forcs_lat=float(para['latorUTM'])
   forcs_lon=float(para['longorUTM'])
   forcs_dpt=float(para['depth'])
   forcs_E=float(para['component dir vect source E'])
   forcs_N=float(para['component dir vect source N'])
   forcs_Z=float(para['component dir vect source Z_UP'])
   TransMat=Trans_Matrix(forcs_lon, forcs_lat,forcs_dpt)
   print("The geodetic directions of point (%.3f,%.3f,0) under local ENU system origioned at (%.3f,%.3f,0) is:"%(forcs_lon,forcs_lat,ref_lon,ref_lat))
   print(" East:",TransMat[0][:])
   print("North:",TransMat[1][:])
   print("   Up:",TransMat[2][:])
   print(TransMat[0][:]+TransMat[1][:]+TransMat[2][:])
   vector=TransMat[0][:]*forcs_E+TransMat[1][:]*forcs_N+TransMat[2][:]*forcs_Z

   X,Y,Up=GeodZNE_2_localUYX(forcs_lon,forcs_lat,forcs_dpt)

   with open(out_file, "w") as f:
      f.write("%s\n" %(para['flag']))
      f.write("time shift: %s\n" %(para['time shift']))
      f.write("f0: %s\n" %(para['f0']))
      f.write("latorUTM: %f\n" %(Y))
      f.write("longorUTM: %f\n" %(X))
      f.write("depth: %f\n" %(Up))
      f.write("factor force source: %s\n" %(para['factor force source']))
      f.write("component dir vect source E: %f\n" %(vector[0]))
      f.write("component dir vect source N: %f\n" %(vector[1]))
      f.write("component dir vect source Z_UP: %f\n" %(vector[2]))
#=====================END of subroutine GEN_forcs====================#


#------------------------subroutine 5: gen_stat---------------------------#
'''
  If the flag RECEIVERS_CAN_BE_BURIED is set as false in the constant.h.in file,
then the elevation and buried depth should be both set to 0. In that case, the specfem3d 
will place the station on the surface according to (X,Y) value only. If buried depth is not
set as 0, then the vertical coordinate of the station will be calculated as Z = elevation - buried_depth/1000.
'''

def GEN_stat(station_fn,topo_txt,out_file):
   #Read elevation file#
   if (flag_topo == 1 ):
      with open(topo_txt, 'r') as f:
         lines = [ l.split() for l in f.readlines() ]

      topo_lons = np.array([ float(l[0]) for l in lines])
      topo_lats = np.array([ float(l[1]) for l in lines])
      topo_alts = np.array([ float(l[2]) for l in lines])

   #Read station file#
   if ( flag_stat == 1 ):
      with open( station_fn, 'r') as f:
         lines = [ l.split() for l in f.readlines() if not(l.startswith('#')) ]

      sta_name = np.array([str(x[0]) for x in lines])
      sta_net = np.array([str(x[1]) for x in lines])
      sta_lon = np.array([float(x[2]) for x in lines])
      sta_lat = np.array([float(x[3]) for x in lines])
      sta_dep = np.array([float(x[4]) for x in lines])

   if(flag_topo==0):
      model_alti=-sta_dep[:]*1000.0
   else:
      altitudes=interpolate.griddata((topo_lons, topo_lats), topo_alts, (sta_lon, sta_lat), method='cubic')
      model_alti = altitudes[:]-sta_dep[:]*1000.0

   ecef = pyproj.Proj(proj='geocent', ellps='WGS84')
   lla = pyproj.Proj(proj='latlong', ellps='WGS84')
   sta_x, sta_y, sta_z = pyproj.transform(lla, ecef, sta_lon, sta_lat, model_alti)
   x0, y0, z0 = pyproj.transform(lla, ecef, ref_lon, ref_lat, ref_alt)

   cosla = np.cos(np.deg2rad(ref_lat))
   sinla = np.sin(np.deg2rad(ref_lat))
   coslo = np.cos(np.deg2rad(ref_lon))
   sinlo = np.sin(np.deg2rad(ref_lon))

   grd_ee =       -sinlo*(sta_x-x0) +       coslo*(sta_y-y0)
   grd_nn = -sinla*coslo*(sta_x-x0) - sinla*sinlo*(sta_y-y0) + cosla*(sta_z-z0)
   grd_uu =  cosla*coslo*(sta_x-x0) + cosla*sinlo*(sta_y-y0) + sinla*(sta_z-z0)

# notice, if RECEIVERS_CAN_BE_BURIED is set as false in the constant.h.in file of the specfem3d, the last clomun should be zero.
   f=open(out_file,"w")
   for i in range(len(sta_lat)):
      f.write("%10s %5s  %12.3f  %12.3f  0.0 %12.3f\n" %(sta_name[i],sta_net[i],grd_nn[i],grd_ee[i],grd_uu[i]))
   f.close()
#=======================End of subroutine GEN_stat===============================#


#---- subrotinue for I/O of the code -----#
def usage():
   print('''
-h or --help
-T or --topography="Parameter used to specify the name of the grided topography data; If this option is not actived, the earth topography is set as zero"
-R or --receiver="Parameter used to specify the name of the receivers; Data format: name net longitude latitude altitude depth"
-F or --forcesolution="Parameter used to specify the name of point force solusions; For more details please refer to the manual of SPECFEM"
-C or --cmtsolution="Parameter used to specify the mame of CMT solution; For more details please refer to the manual of SPECFEM "
-O or --outfile="arameter used to specify the mame of the outfile"
''')
   sys.exit(1)

if (len(sys.argv)==1):
   usage()

shortargs='hT:R:F:C:O:'
longargs=['help','topography','receiver','forcesolution','cmtsolution','outfile']
opts, args = getopt.getopt( sys.argv[1:], shortargs, longargs )

for opt,val in opts:
   if opt in ('-h','--help'):
      usage()
      continue
   if opt in ('-T','--topography'):
      topo_txt=val
      flag_topo=1
      continue
   if opt in ('-R','--receiver'):
      station_fn=val
      flag_stat=1
      continue
   if opt in ('-F','--forcesolution'):
      forcs_fn=val
      flag_forc=1
      continue
   if opt in ('-C','--cmtsolution'):
      cmts_fn=val
      flag_cmts=1
      continue
   if opt in ('-O','--outfile'):
      out_file=val
      continue

#=========End of the subrotinue========#


#--- load mesh parameter file ---#
mesh_par_file = "mesh_par_file.py"
if sys.version_info < (3, ):
  raise Exception("need python3")
elif sys.version_info < (3, 5):
  spec =importlib.machinery.SourceFileLoader("mesh_par", mesh_par_file)
  par = spec.load_module()
else:
  spec = importlib.util.spec_from_file_location("mesh_par", mesh_par_file)
  par = importlib.util.module_from_spec(spec)
  spec.loader.exec_module(par)

# read mesh parameters #
ref_lon   =  par.mesh_ref_lon
ref_lat   =  par.mesh_ref_lat
ref_alt   =  par.mesh_ref_alt
ref_ellps =  par.mesh_ref_ellps
print('Meshing paramters:',ref_lon , ref_lat, ref_alt,ref_ellps)
#===========End of loading mesh parameters===========#



#--------------generate input station file--------------#
if(flag_stat == 1):
   GEN_stat(station_fn,topo_txt,out_file)  
   print('Notice: If RECEIVERS_CAN_BE_BURIED is set as false in the constant.h.in file of the specfem3d, you should replace the last column with 0!')
   print('        The specfem3d code will place the stations on the surface based on the (X,Y) coordinate.')
#============End of generate input station file==============#


#------------generate FORCESOLUTION file------------#
if ( flag_forc == 1 ) :
   GEN_forcs(forcs_fn,out_file)
   print('Notice: If RECEIVERS_CAN_BE_BURIED is set as false in the constant.h.in file of the specfem3d, you should replace the depth with 0!')
   print('        The specfem3d code will place the source on the surface based on the (X,Y) coordinate.')
#============End of generate input FORCESOLUTION file==============#
