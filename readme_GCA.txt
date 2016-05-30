Updated on Sep 14 2012 for the new parameter file type!
Updated on Oct 13, 2012: for both simulated and measured data

Standard input data for inversion: 
1. Forward solver parameter file containing the geometrical and all parameters of the objects (for simulation)
  - for examples forpar_smalldomain.dat
2. Inversion parameter files containing the set of pseudo frequencies, input data file name, regularization parameters, etc.
  - example: invpar_smalldomain.dat
3. Simulated data for the homogeneous medium in the same domain (Laplace transform)
4. Boundary data (simulated or measured) in the frequency domain
5. Simulated data for the object (optional, only for simulated data, used for testing the case of using exact tail only).


Standard folder structure (not required but recommended)
- for each test, create a folder, e.g. "Example1"
- in this folder, create two subfolder "object" and "hom" for object data (for simulated object) and homogeneous medium data. In case of measured data, the "object" folder is not needed.
- run the inversion in the main folder "Example1". 
All results will be stored in this folder.

-----------
COMPILING THE SOURCE CODES:
go to the main GCA folder and run in the terminal: 
make mesh_gen
make forward_solver3D
make laplace
make dat4inv_sim
make gca_scatwave


-----------
HOW TO RUN A TEST WITH SIMULATED DATA:

suppose that the parameter files have been prepared
go to folder Example1 (just an example)

1. Step 1: 
- go to folder object: cd object
- create the mesh for the FEM domain: ..path/mesh_gen ..path/forward_parameter_file.dat density.inp exact_coefficient.m 
(the second and third parameters are optional)
- run the forward solver: ..path/forward_solver3D ..path/forward_parameter_file.dat
- create the Laplace transform of the data of the object (optional, for testing the exact tail only): 
../path/laplace ExactFEM.m NoTimeSteps TimeStep Freq_max lap_ExactFEM_s"Freq_max".m Ngrid 1
Note the output file name  lap_ExactFEM_s"Freq".m must be the same as in the inversion file. Ngrid is not important.

2. Step 2: 
- go to folder ../hom
- copy the exchangeMask.m file from "object" to "hom" in order to speed up the computation: cp ../object/exchangeMask.m exchangeMask.m
- run the forward solver: ..path/forward_solver3D ..path/forward_parameter_file.dat
NOTE: the velocity in the parameter file MUST be 1!!!!
- create the Laplace transform of the data of the homogeneous domain (optional, for testing the exact tail only): 
../path/laplace ExactFEM.m NoTimeSteps TimeStep Freq_max lap_ExactFEM_s"Freq_max".m Ngrid 1

3. Step 3: Create the data for the inversion: 
- go back to the main folder "Example1"
- run ../path/dat4gca_sim ..path/forward_parameter_file.dat ..path/inversion_parameter_file.dat object/SolFEMBoundary.m hom/SolFEMBoundary.m


4. Step 4: Run the Globally convergent algorithm:
- ..path/gca_scatwave ..path/forward_parameter_file.dat ..path/inversion_parameter_file.dat

---------------------
HOW TO RUN A TEST FOR THE SIMULATED DATA:
suppose that the parameter files have been prepared. 

1. Step 1: prepare the data: 
- run the MATLAB function: data_preprocessing_UNCC.m or something similar for data preprocessing
- run the MATLAB function: dat4gca_mea.m for preparing the data for the inversion. 
The parameters should be taken from the inversion parameter files.

Step 2 - Step 4 are the same as the simulated data, see above. NOTE that in step 4, we cannot choose the tail as "exact". We may use a simulated data of some "guess" object for "exact" tail, in this case we have to add step 1 of simulated data. 









