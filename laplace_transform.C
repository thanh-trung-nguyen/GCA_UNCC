/* Do the Laplace transform. Input is a file with structure as a matrix, 
each column represents a time domain signal at a given space point.

 Arguments in this program are:
 argv[1] - name of the data file with extension *.m
 argv[2] - Number of time steps, must be an integer
 argv[3] - time step
 argv[4] - s: value at which the laplace tranform is performed
 argv[5] - name of the output file. 
 argv[6] - Nx: number of columns of the output matrix
 argv[7] - Ny: number of rows.
 argv[8] - (optional) file of the wave for homogeneous medium 
the output file contains only a row vector of Nx elements.
example: in the folder with the file Sol_planeXY_Z004.m, we can run:
../laplace Sol_planeXY_Z04.m 2500 0.002 4 laplapcetr_z04.m 51 61 ../w_hom_h04_fre10/Sol_planeXY_Z04.m 
*/

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
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

  int NoTimeSteps,NoSpacePoints,Nx,Ny, i, j;
  double dt, s;
  char* inputfilename;  char* outputfilename;  char* inputfileHomomedium;
 
 // load the input parameters:
  inputfilename = argv[1];
  NoTimeSteps = atoi(argv[2]);
  dt = atof(argv[3]);
  s = atof(argv[4]);
  outputfilename = argv[5];
  Nx = atoi(argv[6]); // if the resulting vector is re-sized to a matrix Ny * Nx
  Ny = atoi(argv[7]); // if the resulting vector is re-sized to a matrix Ny * Nx
  
  NoSpacePoints = Nx*Ny;

  //check the input:
  cout << inputfilename << endl;  cout << NoTimeSteps << " " << dt << " " << s << " " << endl;
  cout << outputfilename << endl;  cout << Nx << " " << Ny <<  " " << endl;
  
  // load the time domain data from the input file:
  WavesOutputs out_ex; //class containing functions for input-output
  Mat_real time_array;
  time_array.newsize(NoSpacePoints, NoTimeSteps);


  cout << "Load the input file: "<< endl;

  out_ex.ReadSolfromFile(inputfilename,time_array);

  
  //++++++++++add noise in the scattered wave:
  const double RandMax = RAND_MAX + 0.1; // maximum random number but converted to double to avoid the integer division
  srand(time(NULL)); // change the seed of the random function to create different random values for different runs
  double noiselevel = 0.0;
  double   randomvalue; 
 //+++++++++++
 

 // subtract the homogeneous medium signal if provided:
  if (argc > 7)
 	{
		Mat_real array_hommed;
  		array_hommed.newsize(NoSpacePoints,NoTimeSteps);

   		inputfileHomomedium = argv[8];  
		out_ex.ReadSolfromFile(inputfileHomomedium,array_hommed);
		double value;
 		for (i = 0; i < NoSpacePoints; i++)
    			{
      				for (j=0; j < NoTimeSteps; j++)
				{
	  				value = time_array(i,j) - array_hommed(i,j);

					//adding noise:
				        randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
     					time_array(i,j) = value*(1 + randomvalue*noiselevel);

				}
    			}
	}


  MV_Vector<double> Vlaplace(NoSpacePoints);// an intermediate array to get the laplace tranform
  Vlaplace = 0.0;  
  WavesOptional opt;

  cout << "Doing the Laplace tranform:" << endl;
  opt.LaplaceTransform(time_array, dt, s, Vlaplace); //do the Laplace transform, return value to Vlaplace

  // reshape the result to 2D matrix if needed and write the data to the output file:
  cout << "Resizing the array to be a matrix: " << endl;


 if (Ny > 1)
	    {
	      Mat_real LaplaceTr; 
	      LaplaceTr.newsize(Ny,Nx); 
	      LaplaceTr = 0.0;	     
	      for (j = 0; j < Ny; j++)
		{
		  for (i = 0; i < Nx; i++)
		       LaplaceTr(j,i) = Vlaplace(i + Nx*j);
		}
	      out_ex.Write_matrix_to_file(outputfilename,LaplaceTr);
	    }
  else
    {
      out_ex.WriteToFile(outputfilename,Vlaplace);
    }

  return 0;
}

















