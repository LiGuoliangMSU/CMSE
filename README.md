## CMSE python based code
# Cartesian Meshing Spherical Earth (CMSE): A code package for generating spherical earth model mesh based on SPECFEM3D_Cartesain 

This new meshing tool CMSE is compatible with the solver of SPECFEM3D_Cartesian code package (https://geodynamics.org/cig/software/specfem3d/). As proved in our SRL manuscript, it combines the flexibility of the Cartesian mesh, which makes the solver more computationally efficient, and the accuracy of global mesh, which is more accurate in wave simulations on regional scale, including in polar regions.

This python based code works along with the SPECFEM3D_Cartesian code package and you can go to link https://geodynamics.org/cig/software/specfem3d/ to download the Fortran based code pacakge. 

# Expaniation of files 
 Four python based codes, that is, GEN_interface.py, GEN_source_recever.py, IMPL_input_model.py and RotateZYX2UNE.py.
 One mesh par file: mesh_par_file.py
 One file: Example

# demonstrations on how to use those Python based code
Go to Example file:
1. 
''' 
    ./GEN_interface.py -T topography.txt -O interface.dat     or   ./GEN_interface.py -T topography.txt -I 24.4km.interface.txt -O interface.dat
'''
