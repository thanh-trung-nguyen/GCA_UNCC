/* 
function inp2inp.C: 
add the true shape to the reconstruction results. 
Example how to call: 
path/inp2inp inputfile1.inp inputfile1withtrueshape.inp outputfile.inp

Nguyen Trung Thanh, UNCC, 2013. 

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


typedef double real; 
typedef MV_ColMat<real> Mat_real;
typedef MV_ColMat<int> Mat_int;
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;



Vec_real read_inp_file(char *fname)
{
	ifstream inp;
	inp.open(fname);
	int i, j, Nno, Nel;
	double value; 	
	char test[10];
	
	printf("Open the file: %s\n", fname);
	
	// load the number of nodes, number of elements:
	inp >> Nno >> Nel >> i >> i >> i; // the variable i is not used. 

	Vec_real Nodes(Nno);

	//load the nodes' coordinates:
	for (i = 0; i < Nno; i++)
	{
		inp >> j >> value >> value >> value;// NOTE: elements of A is given by (), not []!!!
	}

	//load the element's nodes:
	for (i = 0; i < Nel; i++)
	{
		inp >> j >> j >> test >> value >> value >> value >> value; 
	}
	
	// load the material properties of the nodes:
	inp >> j >> j; //ignore these two data lines
	inp >> test >> j >> test >> test; 
	cout << j << test << endl;

	for (i = 0; i < Nno; i++)
	{
		inp >> j >> Nodes(i);
	}

	inp.close();
	return Nodes;
}




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
 	const char* inpfile1 = argv[1];
	const char* inpfile2 = argv[2];
	const char* outpfile = argv[3];

	ifstream inp;
	inp.open(inpfile1);
	
	int Nno, Nel, i;  // Nno: number of nodes, Nel: number of elements
	// load the number of nodes, number of elements:
	inp >> Nno >> Nel >> i >> i >> i; // the variable i is not used. 
	inp.close();

	Vec_real Coefficient(Nno);
 	WavesGridB gg; // in fact outer_gg is the same as gg!!!
 	WavesOptional opt; 


	Coefficient = read_inp_file(inpfile1);
	gg.scan(inpfile2);
 
        opt.writeInp3D(outpfile, &gg, Coefficient,1); // save the coefficient


	// -----------finish calculation:
 	PetscFinalize();
	return 0;
}
