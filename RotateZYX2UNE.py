#!/usr/bin/env python
# A python code to rotate the local XYZ coordinates to long,lat,alti #
# From geodetic to ECEF, and here we assume the eath is a perfect sphere 
#  R=(6378137+6356752)/2=6367444.5
#  X  =  [ R + h]*cos(la)cos(lo)            
#  Y  =  [ R + h]*cos(la)sin(lo)      
#  Z  =  [ R + h]*sin(la)
# 
# The partial derivates of la:
#  dX = -[ R + h]*sin(la)cos(lo)dla
#  DY = -[ R + h]*sin(la)sin(lo)dla
#  DZ =  [ R + h]*cos(la)dla
#
# The partial derivates of lo:
#  DX = -[ R + h]*cos(la)sin(lo)dla
#  DY =  [ R + h]*cos(la)cos(lo)dla 
#  DZ = 0 
#
#  The partial derivates of h:
#  DX  =  cos(la)cos(lo)dh           
#  DY  =  cos(la)sin(lo)dh      
#  DZ  =  sin(la)dh
#
#From Ecef to local x-y-z
# 
# dx           -sin(LO)         cos(LO)        0         dX
# dy  =   -sin(LA)cos(LO)  -sin(LA)sin(LO)  cos(LA)  *   DY
# dx	   cos(LA)cos(LO)  cos(LA)sin(LO)   sin(LA)      DZ  
#
# Here, LA, LO are the geodetic latitude and longtidue of the referce point for the loca x-y-z coordinate
#
#  A simply python coda wrriten by Guoliang Li for rotating the simulated waveform from the local x-y-z coordinate to Geodetic coordinate
#  Firstly finished: 28th, December, 2019
#  Polished: 18 th, March, 2021
#
#---import libraries --#
import getopt
from obspy import read
import matplotlib.pyplot as plt
import sys
import numpy as np
import math
import importlib
#====End of importing libraries====#

#-- default parameters setting ---#
global flag_info,stat_name,stat_lat,stat_lon,stat_alt,ref_lon,ref_lat
flag_sac=0
flag_assic=0
flag_info=0
out_dir='./'
#===End of setting default parameter===i


#---subrotine for finding station information---#
def get_info(Z,stat_name,stat_lon,stat_lat,stat_alt):
   lon=-12345
   lat=-12345
   alt=-12345
   a=Z.split('/')[-1].split('.')
   stat=a[0]+'.'+a[1]
   for i in range(len(stat_name)):
      if (stat == stat_name[i]):
         lon=stat_lon[i]
         lat=stat_lat[i]
         alt=stat_alt[i]
   return lon, lat, alt   
#======End of subrotine get_info==============#

#---subrotine for read in assic format waveform---#
def read_assic_data(fn):
   with open(fn, 'r') as f:
      lines = [ l.split() for l in f.readlines() ]

      time_series = np.array([ float(l[0]) for l in lines])
      data = np.array([ float(l[1]) for l in lines])
   b=time_series[0]
   npts=len(time_series)
   delta=(max(time_series)-min(time_series))/(npts-1)
   return data, b ,delta, npts
#=======End of subrotine read_assic_data========#

#---subrotine for writing out the assic format data---#
def write_assic_data(data,b,delta,npts,fn):
   with open(fn, "w") as f:
      for i in range(npts):
         f.write("%10.4f  %+12.5E\n" %(b+delta*i,data[i]))
#========End of subrotine write_assic===========#


#-----------subroutine: Trans_Matrix--------------#
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
   #     partial derivate of longtitude or East   
   dX[0]= -(R+alti)*cosla*sinlo
   dY[0]=  (R+alti)*cosla*coslo
   dZ[0]=  0

   #     partial derivate of Latitude or North
   dX[1]= -(R+alti)*sinla*coslo
   dY[1]= -(R+alti)*sinla*sinlo
   dZ[1]=  (R+alti)*cosla

   #     partial derivate of altitude
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


#---- subrotinue for I/O of the code -----#
def usage():
   print('''
-h or --help
-A or --assic="Parameter used to specify that the input waveforms are in assic format"
-S or --sac="Parameter used to specify that the input waveforms are in assic format"
-I or --outfile="Parameter used to specify the name of stations"
-O or --out_dir="Parameter used to specify the name of output path; A directory"
''')
   sys.exit(1)

if (len(sys.argv)==1):
   usage()

shortargs='hA:S:I:O:'
longargs=['help','assic','sac','information','out_dir']
opts, args = getopt.getopt( sys.argv[1:], shortargs, longargs )

for opt,val in opts:
   if opt in ('-h','--help'):
      usage()
      continue
   if opt in ('-A','--assic'):
      rotate_list=val
      flag_assic=1
      continue
   if opt in ('-S','--sac'):
      rotate_list=val
      flag_sac=1
      continue
   if opt in ('-I','--information'):
      stat_info=val
      flag_info=1
      continue
   if opt in ('-O','out_dir'):
      out_dir=val
      continue


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

#averaged radiu of the earth
R=(6378137+6356752)/2.0

# define the vectors
dX=np.array(np.zeros((3,1)))
dY=np.array(np.zeros((3,1)))
dZ=np.array(np.zeros((3,1)))

#define rotation matrix
TransMat=np.mat(np.zeros((3,3)))

# read in station information
if (flag_info == 1):
   with open(stat_info, 'r') as f:
      lines = [ l.split() for l in f.readlines() ]

   stat_name = np.array([ str(l[0]) for l in lines])
   stat_lon = np.array([ float(l[1]) for l in lines])
   stat_lat = np.array([ float(l[2]) for l in lines])
   stat_alt = np.array([ float(l[3]) for l in lines])
   print(stat_name)
# read in the rotation list#
F=open(rotate_list,"r")
for string in F.readlines():
   fnZ=string.strip('\n')	
   print(fnZ)
   fnY=fnZ.replace('XZ.','XY.')
   fnX=fnZ.replace('XZ.','XX.')

#  the data is in sac format
   if (flag_sac == 1 ) :
      data_Z=read(fnZ,debug_headers=True)
      data_Y=read(fnY,debug_headers=True)
      data_X=read(fnX,debug_headers=True)
      dataZ=data_Z[0].data
      dataY=data_Y[0].data
      dataX=data_X[0].data
      Z_b=data_Z[0].stats.sac.b
      Y_b=data_Y[0].stats.sac.b
      X_b=data_X[0].stats.sac.b
      Z_delta=data_Z[0].stats.sac.delta
      Y_delta=data_Y[0].stats.sac.delta
      X_delta=data_X[0].stats.sac.delta
      Z_npts=data_Z[0].stats.sac.npts
      Y_npts=data_Y[0].stats.sac.npts
      X_npts=data_X[0].stats.sac.npts
      if (flag_info == 0 ):
         stlo=data_Z[0].stats.sac.stlo
         stla=data_Z[0].stats.sac.stla
         stel=data_Z[0].stats.sac.stel
      else:
         stlo,stla,stel=get_info(fnZ,stat_name,stat_lon,stat_lat,stat_alt)

#  the data is in assic format 
   elif (flag_assic == 1):
      dataZ,Z_b,Z_delta,Z_npts=read_assic_data(fnZ)
      dataY,Y_b,Y_delta,Y_npts=read_assic_data(fnY)
      dataX,X_b,X_delta,X_npts=read_assic_data(fnX)
      stlo,stla,stel=get_info(fnZ,stat_name,stat_lon,stat_lat,stat_alt)
   
   #Checking if the begin time, the delta, npts of the input three componets are correct
   if Z_b != Y_b or Y_b != X_b:
      print('The begin time of the input three components are not equal!\n')
      sys.exit() 
   elif Z_delta != Y_delta or Y_delta != X_delta: 
      print('The sampling ratio of the input three components are not equal!\n')
      sys.exit() 
   elif Z_npts != Y_npts or Y_npts != X_npts: 
      print('The npts of the input three components are not equal!\n')
      sys.exit() 
     
   #calculate the vector of the geodetic coordinates under the local x-y-z coordinate
   #=================================================================================

   TransMat=Trans_Matrix(stlo, stla, stel)
      
   dataE=TransMat[0,0]*dataX[:]+TransMat[0,1]*dataY[:]+TransMat[0,2]*dataZ[:]
   dataN=TransMat[1,0]*dataX[:]+TransMat[1,1]*dataY[:]+TransMat[1,2]*dataZ[:]
   dataU=TransMat[2,0]*dataX[:]+TransMat[2,1]*dataY[:]+TransMat[2,2]*dataZ[:]
   
   outU=fnZ.replace('XZ.','XU.').split('/')[-1]
   outN=fnY.replace('XY.','XN.').split('/')[-1]
   outE=fnX.replace('XX.','XE.').split('/')[-1]
   outU=out_dir+'/'+outU
   outN=out_dir+'/'+outN
   outE=out_dir+'/'+outE
   print(outU,outN,outE)
 
   if (flag_sac == 1):
      #rotate the data and wirte out
      data_X[0].data=dataE
      data_Y[0].data=dataN
      data_Z[0].data=dataU

      data_X[0].stats.sac.cmpinc=90
      data_X[0].stats.sac.cmpaz=90
      data_Y[0].stats.sac.cmpinc=90
      data_Y[0].stats.sac.cmpaz=0
      data_Z[0].stats.sac.cmpinc=0
      data_Z[0].stats.sac.cmpaz=0
      
      data_Z.write(outU,format='SAC')
      data_Y.write(outN,format='SAC')
      data_X.write(outE,format='SAC')

   elif (flag_assic == 1):
      write_assic_data(dataU,Z_b,Z_delta,Z_npts,outU)
      write_assic_data(dataN,Y_b,Y_delta,Y_npts,outN)
      write_assic_data(dataE,X_b,X_delta,X_npts,outE)
     
