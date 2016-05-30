/* 
This function extract the scattered wave from the total wave
Given: total wave with objects, total wave without object (incident wave propagation)
Number of time steps
Output: a data file with the scattered wave in the same format as the input files

*/

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <sstream>


#include "include/wavesOutputs.h"
#include "include/waveskraftwerk.hh"
#include "include/wavesSDGeometry.h"
#include "include/wavesOptional.h"


using namespace std;


typedef double real; 
typedef MV_ColMat<real> Mat_real;


//============================================================ 
int main(int argc, char **argv)
{

  int NoTimeSteps,NoBndNodes, i, j;
 char* inputfilename;  char* outputfilename;  char* inputfileHomomedium;
 
 // load the input parameters:
  inputfilename = argv[1];
  inputfileHomomedium = argv[2];
  outputfilename = argv[3];
//  NoTimeSteps = atoi(argv[4]);
//  NoBndNodes = atoi(argv[5]);
  

	// get the number of time steps and number of boundary grid nodes:
	string line; 
	double value; 
	ifstream inpfile;
	inpfile.open(inputfilename);

	// get the number of boundary grid nodes (number of columns in the input file) 
	getline(inpfile, line);
   	istringstream iss(line);
	NoBndNodes = 0; 	
	while (iss >> value)
		NoBndNodes++; 
   	cout << "Number of boundary grid nodes = " << NoBndNodes << endl; 
	// get the number of time steps (number of rows - 1): 
	NoTimeSteps = 0; 
	while ( getline(inpfile, line) )
   		NoTimeSteps++;
   	cout << "Number of time steps = " << NoTimeSteps << endl; 
	inpfile.close();
  
  // load the time domain data from the input file:
  WavesOutputs out_ex; //class containing functions for input-output

  Mat_real data_obj, data_hom;
  data_obj.newsize(NoTimeSteps,NoBndNodes);
  data_hom.newsize(NoTimeSteps,NoBndNodes);
  
  MV_Vector<int> BdNodes(NoBndNodes);	
  MV_Vector<double> data(NoBndNodes); 
  cout << "Load the input file: "<< endl;

  out_ex.load_boundary_data(inputfilename,data_obj,BdNodes);
  out_ex.load_boundary_data(inputfileHomomedium,data_hom,BdNodes);

  out_ex.Write_array_my_part(outputfilename,BdNodes);
 
// subtract the homogeneous medium signal:
	for (i = 0; i < NoTimeSteps; i++)
 	{
 		for (j=0; j < NoBndNodes; j++)
		{
			data(j) = data_obj(i,j) - data_hom(i,j);
		}
		out_ex.Write_array_my_part(outputfilename,data,i);//write the scattered wave to the output file
	} 

  return 0;
}

















