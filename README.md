## CMSE python based code
# Cartesian Meshing Spherical Earth (CMSE): A code package for generating spherical earth model mesh based on SPECFEM3D_Cartesain 

This new meshing tool CMSE is compatible with the solver of SPECFEM3D_Cartesian code package (https://geodynamics.org/cig/software/specfem3d/). As proved in our SRL manuscript, it combines the flexibility of the Cartesian mesh, which makes the solver more computationally efficient, and the accuracy of global mesh, which is more accurate in wave simulations on regional scale, including in polar regions.

This python based code works along with the SPECFEM3D_Cartesian code package and you can go to link https://geodynamics.org/cig/software/specfem3d/ to download the Fortran based code pacakge. 

# Expaniation of files 
 Four python based codes, that is, GEN_interface.py, GEN_source_recever.py, IMPL_input_model.py and RotateZYX2UNE.py. Those python based codes are used to  transfer the geodetic coordinate to the local x-y-z coordinate. For more details explainations, please refer to our SRL manuscripts (a link will be avalibale soon).
 
 One mesh par file: mesh_par_file.py defines the origin of the local x-y-z coordinates and boundaries of the simulated domain. 
 
 One file: Example file contains data used to demonstrate how to use the python based codes

# demonstrations on how to use those Python based code

In each python based code, there is a detailed examplantion about the functions of each parameters. For example, in GEN_interface.py script, you can file the following contents: 

def usage():

   print('''
   
-h or --help

-T or --topography="Parameter used to specify the name of the grided topography data; If this option is not actived, the earth topography are set as zero"

-I or --interface="Parameter used to specify the name of grided internal interfaces, such as the sedimentary and Moho discontinuty"

-O or --outfile="Parameter used to specify the name of output file; Dafault file name is interface.dat"

''')

Next, I will only show how to run those python based codes with the datasets in Example file.

Go to Example file:
1. prepare the interfaces for the simulated region 

```
python ../GEN_interface.py -T topography.txt -O interface.dat     or  
python ../GEN_interface.py -T topography.txt -I 24.4km.interface.txt -O interface.dat
```    
2. The code IMPL_initial_model.py is used to implement the input model. As the used file is too large, I didn't upload it. Here I only post the command how to run it:
```
python ../IMPL_initial_model.py -T topography.txt -M initial_model.txt -B ./DATABASES_MPI -N 144 -K 12 -O ./GLL
```


3. Transfer the receivers and sources from geodetic coordinate to the local x-y-z coornidate
```
python ../GEN_source_recever.py -T topography.txt -R stations.txt -O STATIONS    or
python ../GEN_source_recever.py -T topography.txt -F forcesolution.txt -O FORCESOLUTION
```

4. Rotate the simulated waveforms from local x-y-z coordinate to geodetic coordinate
```
cd ./Roate_example
python ../../RotateZYX2UNE.py -A rotate.lst -I station.lst -O Output
```

# Acknowledgements
We thank the Institute for Cyber-Enabled Research (ICER) at Michigan State University, the Extreme Science and Engineering Discovery Environment (XSEDE supported by NSF grant ACI-1053575) for providing the high-performance computing resources. This research is supported by NSF Grant 1942431 and startup funds of Min Chen at Michigan State University. Additionally, co-author RM was supported by NSF-EAR postdoctoral followship # 1806412.
