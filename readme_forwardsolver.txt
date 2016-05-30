Updated Sep 14 for new data file type:


how to compile and run the forward solver: 

Before compiling the C routines, make sure that the PETSC_DIR and PETSC_ARCH paths in the Makefile
are correct. 

To compile, run the following jobs:
make mesh_generation  // mesh generation
make forward_solver3D // forward solver
make laplace          // laplace transform

To run the forward solver: 

Step 1: create or modify the input data file .dat, see the structure in the file forpar_new_full.dat:


Step 2: create the grid file: run the following command from the WavES_new:
	./mesh_generation "input parameter file" density.inp coefficient.m 


Step 3: create a folder to store the results:
mkdir "folder name"

and go to the created folder:
cd "folder name"

NOTE: if the folder exists, delete all the old file before running the forward solver!!!
rm *.*


Step 4: from the "folder name", run the forward solver:
../forward_solver3D ../parameters/par2_1obj_h02.dat >/dev/null& 
or if the material is given in a file b.m: 
../forward_solver3D ../parameters/par_1obj_h02.dat 1 b.m >/dev/null& 


Step 5: Do the Laplace transform: 
../laplace inputfile.m NoTimeSteps TimeStep s outputfile.m Nx Ny inputfile_homogeneous.m


For example: 
../laplace Sol_Z04.m 2500 0.002 4 laplapcetr_z04_s4.m 51 61 ../test_hom/Sol_Z04.m 


Extract the boundary scattered waves:

../../extract_scatteredwave SolFEMBoundary.m ../data_hom/SolFEMBoundary.m ScatWaveFEMBoundary.m 1500 4082
argv[1]: input boundary data file with object
argv[2]: input boundary data file with homogeneous medium
argv[3]: output data file of scattered waves on the boundary
argv[4]: number of time steps in the data files
argv[5]: number of boundary grid points 

Update: argv[4] and argv[5] are no longer needed. Only three input parameters. 
../../extract_scatteredwave SolFEMBoundary.m ../data_hom/SolFEMBoundary.m ScatWaveFEMBoundary.m

Visualization of the simulated data in MATLAB:



