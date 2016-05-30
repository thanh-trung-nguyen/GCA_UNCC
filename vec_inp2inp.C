/* 
Function vec_inp2inp.C load the coefficient values from a vector file and the grid file from INP file to create a new INP file with imbeded coefficient values
input: 
argv[1]: the vector file of the coefficient values
argv[2]: the grid file INP consistent with the coefficient file
argv[3]: output file name with extension INP

To compile: make vec_inp2inp

Example to run the program: go to a folder inside the main GCA folder or make a new subfolder, then run:
../vec_inp2inp coef_data.dat gridfile.inp gridfile2.inp

Nguyen Trung Thanh, UNCC, 2013. 
Updated: April 15, 2013: 
*/

static char help[] ="";

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "include/waveskraftwerk.hh"
#include "include/wavesSDGeometry.h"
#include "include/wavesOptional.h" 
#include "include/wavesOutputs.h" 


typedef double real; 
typedef MV_Vector<double> Vec_real;


//=========================================================
//================== THE MAIN PROGRAM =====================
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int main(int argc, char **argv)
{
	PetscInitialize(&argc, &argv, (char *) 0, help);
	int rank;
 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	cout<<" after PETSCInitialize "<<endl;
	// ----------------------------------------

	// input parameters files: 
 	const char* coef_value_file = argv[1];
	const char* inp_grid_file = argv[2];
	const char* out_grid_file = argv[3];


	WavesGridB gg; // in fact outer_gg is the same as gg!!!

	WavesOptional opt; // For saving the coefficient
	WavesOutputs out; 

	gg.scan(inp_grid_file);
 	int Nno = gg.getNoNodes();
	Vec_real Coefficient(Nno);

	out.ReadSolfromFile(coef_value_file,Coefficient);
	
	opt.writeInp3D(out_grid_file, &gg, Coefficient,1); // save the coefficient
	
	

	PetscFinalize();
	return 0;
}
