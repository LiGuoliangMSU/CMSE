#!/usr/bin/env python
"""
Parameters for mesh setup
"""
#--- reference point (WGS84 ellipsoid)
mesh_ref_lon = -119     # in degree
mesh_ref_lat =  75      # in degree
mesh_ref_alt = 0.0      # in meter 
mesh_ref_ellps = 'WGS84' # reference earth model
#---Meshing boudaries
Y_MIN = -650000.0  # South, in meter   
Y_MAX =  650000.0  # North, in meter
X_MIN = -650000.0  # West, in meter
X_MAX =  650000.0  # East, in meter
# grid interval (should be smaller than the actual SEM mesh element size)
dx    =    2000.0  
dy    =    2000.0 
