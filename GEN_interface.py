#!/usr/bin/env python
#---import libraries --#
import sys
import getopt
#import warnings
import numpy as np
from scipy import interpolate, ndimage
from netCDF4 import Dataset
import pyproj
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import importlib
#====End of importing libraries====#

#-- default parameters setting ---#
flag_topo=0
flag_inte=0
flag_plot=0
out_file="interface.dat"
#===End of setting default parameter===#


#---- subrotinue for I/O of the code -----#
def usage():
   print('''
-h or --help
-T or --topography="Parameter used to specify the name of the grided topography data; If this option is not actived, the earth topography are set as zero"
-I or --interface="Parameter used to specify the name of grided internal interfaces, such as the sedimentary and Moho discontinuty"
-O or --outfile="Parameter used to specify the name of output file; Dafault file name is interface.dat"
''')
   sys.exit(1)

if (len(sys.argv)==1):
   usage()

shortargs='hT:I:O:'
longargs=['help','topography','interface','outfile']
opts, args = getopt.getopt( sys.argv[1:], shortargs, longargs )

for opt,val in opts:
   if opt in ('-h','--help'):
      usage()
      continue
   if opt in ('-T','--topography'):
      topo_txt=val
      flag_topo=1
      continue
   if opt in ('-I','--interface'):
      disconty=val
      flag_inte=1
      continue
   if opt in ('-O','--outfile'):
      out_file=val
       
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
min_x     =  par.X_MIN
max_x     =  par.X_MAX 
min_y     =  par.Y_MIN
max_y     =  par.Y_MAX 
dx        =  par.dx
dy        =  par.dy
print('Meshing paramters:',ref_lon , ref_lat, ref_alt,ref_ellps, min_x ,min_x,min_y,max_y,dx,dy)
#===========End of loading mesh parameters===========#



#-------interpolate and transfer topography/interface data from geodetic GRD file to the local ENU grid-----#

# generating local ENU grid (x,y)#
x1 = np.arange(min_x, max_x, dx)
y1 = np.arange(min_y, max_y, dy)
x2, y2 = np.meshgrid(x1, y1)

#Read elevation file and depth file#
if (flag_topo == 1 ):
   with open(topo_txt, 'r') as f:
      lines = [ l.split() for l in f.readlines() ]

   topo_lons = np.array([ float(l[0]) for l in lines])
   topo_lats = np.array([ float(l[1]) for l in lines])
   topo_alts = np.array([ float(l[2]) for l in lines])

if ( flag_inte == 1):
   with open(disconty, 'r') as f:
      lines = [ l.split() for l in f.readlines() ]

   disc_lons = np.array([ float(l[0]) for l in lines])
   disc_lats = np.array([ float(l[1]) for l in lines])
   disc_depth = np.array([ float(l[2]) for l in lines])

# in case that only generate free surface #
if ( flag_topo == 1 ) and ( flag_inte == 0 ):
   disc_lons=np.zeros(len(topo_lons))
   disc_lats=np.zeros(len(topo_lats))
   disc_depth=np.zeros(len(topo_alts))
   disc_lons=topo_lons
   disc_lats=topo_lats
   disc_depth=topo_alts/-1000.0

# in case that free surface is at seal level #
if ( flag_topo == 0 ) and ( flag_inte == 1 ):
   topo_lons=np.zeros(len(disc_lons))
   topo_lats=np.zeros(len(disc_lats))
   topo_alts=np.zeros(len(disc_depth))
   topo_lons=disc_lons
   topo_lats=disc_lats
   topo_alts=-1000*disc_depth

# in case that free surface is not at seal level #
if ( flag_topo == 1 ) and ( flag_inte == 1 ):
   topo_alts=topo_alts/1000
   #Interpolate elevation data points to interface data points#
   grd_alts = interpolate.griddata((topo_lons, topo_lats), topo_alts, (disc_lons,disc_lats), method='cubic')
   topo_alts=1000*grd_alts-1000*disc_depth
  
#convert (lon,lat,alt) to ECEF#
ecef = pyproj.Proj(proj='geocent', ellps=ref_ellps)
lla = pyproj.Proj(proj='latlong', ellps=ref_ellps)

print('The results written in %s' %(out_file))
x0, y0, z0 = pyproj.transform(lla, ecef, ref_lon, ref_lat, ref_alt)
xx, yy, zz = pyproj.transform(lla, ecef, disc_lons, disc_lats, topo_alts)

#transform from ECEF to REF_ENU
cosla = np.cos(np.deg2rad(ref_lat))
sinla = np.sin(np.deg2rad(ref_lat))
coslo = np.cos(np.deg2rad(ref_lon))
sinlo = np.sin(np.deg2rad(ref_lon))

print(out_file)
#RotM = np.zeros((3,3))
#RotM[0,:] = [       -sinlo,        coslo,   0.0 ]
#RotM[1,:] = [ -sinla*coslo, -sinla*sinlo, cosla ]
#RotM[2,:] = [  cosla*coslo,  cosla*sinlo, sinla ]
grd_ee =       -sinlo*(xx-x0) +       coslo*(yy-y0)
grd_nn = -sinla*coslo*(xx-x0) - sinla*sinlo*(yy-y0) + cosla*(zz-z0)
grd_uu =  cosla*coslo*(xx-x0) + cosla*sinlo*(yy-y0) + sinla*(zz-z0)

#interpolate uu(ee,nn) to local ENU grid (x,y) and out write out the results
z2 = interpolate.griddata((grd_ee, grd_nn), grd_uu, (x2, y2), method='linear')
with open(out_file, "w") as f:
  for iy in range(len(y1)):
    for ix in range(len(x1)):
        f.write("%+12.5E\n" %(z2[iy,ix]))
#        f.write("%10.3f %10.3f %+12.5E\n" %(x1[ix], y1[iy],z2[iy,ix]))

 #End of interpolating and transfering topography/interface data from geodetic GRD file to the local ENU grid=====#
