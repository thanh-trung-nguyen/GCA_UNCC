/* create the exchangeMask.m file for speeding up the computation in the globally convergent algorithm:

Arguments:argv[1]  name of the input *.dat file for the forward solver

Compile: make create_exchangeMask.

example how to run: ./create_exchangeMask folder_path/parameter_file_name.dat

*/

static char help[] = "";

#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "include/wavesWaveEqOpSeismic.h"
#include "include/waveskraftwerk.hh"
#include "include/wavesSDGeometry.h"
#include "include/wavesScalarEqOpOpt.h"
#include <sys/times.h>
#include <mpi.h>


using namespace std;

typedef double real;
typedef MV_ColMat<real> Mat_real;


//============================================================

int main(int argc, char **argv)
{

	PetscInitialize(&argc, &argv, (char *) 0, help);
	
	int nsd;
	bool USE_FEM, USE_FDM, USE_RHS, USE_DIRICHLET_FEM, 
	     USE_DIRICHLET_FDM, PRINT_FILES, EXCHANGE, USE_ABSORB;

	double maxtime, rhs;
	double guess_velocity = 1; // background velocity.

	int nrSTEPS;

	WavesGridB gg;
	WavesGridB outer_gg;
	WavesSDGeometry sdg;
	WavesSDGeometry sdg_new;
        WavesOutputs out1;

	bool doIncludeCorners = 1;
	int type_of_material1 = 2; 
	double velocity1 = 1.0;

	// load the parameters: call function Configurate from WavesOutput class
	out1.Configurate(argv[1], gg, outer_gg, sdg, sdg_new, nsd, EXCHANGE, USE_FEM, USE_FDM, 
		    USE_DIRICHLET_FEM, USE_DIRICHLET_FDM);

	USE_RHS = false; USE_ABSORB = true; PRINT_FILES = false; 
	nrSTEPS = 100; maxtime = 1; rhs = 1.0; 
	if (argc > 1)
	{
		// initialize the forward solver: 
		WavesScalarEqOpOpt p(gg, outer_gg, sdg, nsd, EXCHANGE, USE_FEM, USE_FDM, 
				     USE_RHS, USE_DIRICHLET_FEM, USE_DIRICHLET_FDM, 
				     USE_ABSORB, PRINT_FILES, nrSTEPS, maxtime, rhs, 
				     type_of_material1, velocity1, guess_velocity, doIncludeCorners);

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		MV_Vector<double> b(gg.getNoNodes()); // vector of material at grid nodes
		
		// calculate the vector of material (coefficient): 
		b = guess_velocity;		

		if (EXCHANGE)
		{
			p.InitFDM();

			if (nsd == 2)
				p.Init2DFEM(b);
			else if (nsd == 3)
				p.Init3DFEM(b);
			p.InitExchangeStructCommon(); // print the exchangeMask.m file
			//p.InitExchangeCommon();
		}
	}
	else
		cout << "Usage:  CONFIGURATION, see arguments argv[1], argv[2], argv[3] argv[4] argv[5]  " << endl;


	PetscFinalize();
	return 0;
}

