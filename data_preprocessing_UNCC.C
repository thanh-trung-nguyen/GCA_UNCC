/* 
function data_preprocessing_UNCC.C is used to preprocess the measured data:
NOT completed yet!!!!



*/


#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "include/wavesMVmtp.h"
#include "include/wavesMVvtp.h"
#include "include/wavesMVblas.h"

using namespace std; 

typedef double real; 
typedef MV_ColMat<real> Mat_real; //this data type is defined in the packages included. 
typedef MV_Vector<double> Vec_real;
const double lightspeed = 0.3; // light speed in meters/nanoseconds

void read_from_file(char *fname, Mat_real & A)
{
	ifstream inp;
	inp.open(fname);
	int i, j;
	printf("Open the file: %s\n", fname);
	cout << "number of rows = " << A.size(0) << "; number of columns = " << A.size(1) << endl;
	
	for (i = 0; i < A.size(0); i++)
	{
		for (j = 0; j < A.size(1); j++)
		{
			inp >> A(i, j);// NOTE: elements of A is given by (), not []!!!

		}
	}

	inp.close();
}

void write_matrix_to_file(char* file, Mat_real& A) // write a matrix to a file
{

	int i, j, n_i, n_j;
	n_i = A.size(0); //number of rows
	n_j = A.size(1); // number of columns

	ofstream outp;
	outp.open(file);

	cout << "opened file  " << file << " No of rows n_i = " << n_i << " No of columns n_j = " << n_j << endl;

	for (i = 0; i < n_i; i++)
	  {
		for (j = 0; j < n_j; j++)
		  {
			outp << A(i,j) << " ";
		  }
		outp << endl;
	  }

	cout << "closed file" << file << endl;
	outp.close();
}

void matrix_shift(Mat_real& A, double b)
{
	int i,j; 
	for (i = 0; i < A.size(0); i++)
	{	for (j=0; j < A.size(1); j++)
			A(i,j) += b;
	}
}

double matrix_mean(Mat_real A)
{
	int i,j;
 	double Mean = 0.0;
	for (i = 0; i < A.size(0); i++)
	{	for (j=0; j < A.size(1); j++)
			Mean += A(i,j);
	}
	return Mean/A.size(0)/A.size(1);
}


double Minim(Vec_real VecValue)
{
  double minim = VecValue(0);
  unsigned int i; 
  for (i=1; i < VecValue.size();i++)
    if( VecValue(i) < minim)
      minim = VecValue(i);
  return minim;
}

double Maxim(Vec_real VecValue)
{
  double maxim = VecValue(0);
  unsigned int i; 
  for (i=1;i < VecValue.size();i++)
    if( VecValue(i) > maxim)
      maxim = VecValue(i);
  return maxim;
}

// do the zero offset for each time-domain signal
void zero_offset_shift(Mat_real& A)
{
	int i,j; 
	for (i = 0; i < A.size(1); i++)
	{	
		double MeanValue = 0.0;
		for (j=0; j < A.size(0); j++)
			MeanValue += A(j,i);
		MeanValue = MeanValue/A.size(0); 
		for (j=0; j < A.size(0); j++)
			A(j,i) -= MeanValue;
		
	}
}

void load_mea_parameters(char* filename, double& MaxTime, int& NoSamples, 
			 int& Nx, int& Ny, double& Xmin, double& Xmax, double& Ymin, double& Ymax, 
			 double& d1, double& d2, double& d3, int& FirstPeak, int& NewNoSamples, double& Threshold)
{
	char text[150];
	ifstream inp;
  	inp.open(filename);
	if (inp.is_open())
	{
		inp >> text >> MaxTime;
		inp >> text >> NoSamples; 
		inp >> text >> Nx;
		inp >> text >> Ny; 
		inp >> text >> Xmin >> Xmax;
		inp >> text >> Ymin >> Ymax;
		inp >> text >> d1; 
		inp >> text >> d2; 
		inp >> text >> d3; 
		inp >> text >> FirstPeak;
		inp >> text >> NewNoSamples;
		inp >> text >> Threshold;

	inp.close();
	}
	else 
		cout << "Error in load_mea_parameters: the file cannot be opened" << endl;
}


//==============MAIN==============================
int main(int argc, char **argv)
{

	char* inputfilename;
  	char* outputfilename;

	char* parameter_file; 
	parameter_file = argv[1];
  	inputfilename = argv[2];
  	outputfilename = argv[3];


	//setup parameters: should be changed if the measurement setup is changed!
 	int NoSamples,NewNoSamples, Nx, Ny, FirstPeak;
	double MaxTime, dist_tr2obj, dist_tr2re, dist_new_source, Xmin, Xmax, Ymin, Ymax, Threshold; 

	load_mea_parameters(parameter_file, MaxTime, NoSamples, 
			    Nx, Ny, Xmin, Xmax, Ymin, Ymax, dist_tr2obj, dist_tr2re, dist_new_source, FirstPeak,NewNoSamples,Threshold);

	int Np = Nx*Ny;
 		
 	// create a matrix to store the data:
 	Mat_real DATA; 
 	DATA.newsize(NoSamples,Np);
 	DATA = 0.0;
 
 	// load the data from the input file to matrix A:
 	cout << "Load the input data from file: " << inputfilename << endl;
 	read_from_file(inputfilename,DATA); //matrix A now contains the data

 	//==== Here is your calculation: 
 	// +++++++++++++++++++++++++++

	//--------Step 1: shift the data to zero-offset: 
	zero_offset_shift(DATA);

	//--------Step 2: shift the data to the correct time zero at a given source distance:
	double dt = MaxTime/(NoSamples+1);
		
 	Mat_real DATA2; 
 	DATA2.newsize(NewNoSamples,Np);
 	DATA2 = 0.0;
	double travel_dist_first_point = sqrt(Xmin*Xmin + Ymin*Ymin + (dist_tr2obj+dist_tr2re)*(dist_tr2obj+dist_tr2re)) + dist_new_source;

	int NewFirstPeak = round(travel_dist_first_point/(dt*lightspeed));

cout << "MaxTime = " << MaxTime << endl;
cout << "Number of samples = " << NoSamples << endl;
cout << "Nx, Ny= " << Nx << Ny << endl;
cout << "Xmin, Xmax, Ymin, Ymax = " << Xmin << " " << Xmax << " " << Ymin << " " << Ymax << endl;
cout << "Distances = " << dist_tr2obj << " " << dist_tr2re << " " << dist_new_source << endl;
cout << "First peak = " << FirstPeak << endl;
cout << "New first peak = " << NewFirstPeak << endl;
cout << "Travel distance from tran. to object to receiver at the top left corner: " << travel_dist_first_point << endl;


	if (NewFirstPeak >= NewNoSamples)
	{
		cout << "Warning: the first peak greater than the number of samples!" << endl;	
		exit(1);
	}
	int i, j;

	if (NewFirstPeak > FirstPeak)
	{
		int Shift = NewFirstPeak - FirstPeak;
		int MaxIdx = NoSamples+Shift; 
		if MaxIdx > NewNoSamples
			MaxIdx = NewNoSamples;

		for (i=0; i < Np; i++)
		{
			for (j=Shift;j < MaxIdx; j++) 
				DATA2(j,i) = DATA(j-Shift,i);
		}
	}
	else if (NewFirstPeak < FirstPeak)
	{
		int Shift = FirstPeak - NewFirstPeak;
		for (i=0; i < Np; i++)
		{
			for (j=1;j < NoSamples-Shift; j++) 
				DATA2(j,i) = DATA(j+Shift,i);
		}
	}


	//------- Step 3: remove the noise before the first peaks:
	// 
	double mindist = dist_tr2obj+dist_tr2re + dist_new_source;
	int minFirstPeak = round(mindist/dt/lightspeed); 
	int minFirstPeakShift = minFirstPeak-10; // 10 is just a test!!!

	//remove all signals before the first peak that come first:
	for (i = 0; i < Np; i++)
	{
		for (j=0; j < minFirstPeakShift; j++) 
			DATA2(j,i) = 0.0;
	}

	// find the first negative peak (use the second negative peak to find the first negative peak):
	MV_Vector<double> tmp(NewNoSamples);
	int fidx; 
	for (i = 0; i < Np; i++)
	{
		tmp = 0.0;
		for (j = minFirstPeakShift; j < NewNoSamples; j++) 
		{
			tmp(j) = DATA2(j,i);//create a vector of large values of data		
		}
		double Min = Minim(tmp); 
		for (j = minFirstPeakShift; j < NewNoSamples; j++) 
		{
			if (fabs(tmp(j)) < -0.7*Min) // 0.7 is just a test!!!!
				tmp(j) = 0.0;
		}

		if (Minim(tmp) < 0)
		{
			// find the first negative peak of vector tmp:
			fidx = minFirstPeakShift;
			while ((tmp(fidx) >= 0) || (tmp(fidx) >= tmp(fidx+1))) 
				fidx += 1;
cout << "test "<< i+1 << " " << fidx+1;

			// next, find the first POSITIVE peak starting from the second negative peak and go backwards: 
			while ((DATA2(fidx,i) <= DATA2(fidx-1,i)) && fidx > 1)		
				fidx -= 1; 
cout << " " << fidx+1; 
			// next, find the first NEGATIVE peak starting from the first POSITIVE peak and go backwards: 
			while ((DATA2(fidx,i) >= DATA2(fidx-1,i)) && fidx > 1)		
				fidx -= 1; 
		
cout << " " << fidx+1; 
			// next, find the first negative value of the first valley! 
			while ((DATA2(fidx,i) < 0) && fidx > 1)		
				fidx -= 1; 
			fidx +=1; //the first negative value of the first valley!

cout << " " << fidx+1 << endl;

			for (j= minFirstPeakShift; j < fidx; j++)
				DATA2(j,i) = 0.0; 	
		}
		else //in case that the signal is smaller than the threshold, we set it to be zero!
		{
			cout << "Data at position " << i << " is too small, set to zero" << endl;
			for (j = minFirstPeakShift; j < NewNoSamples; j++)
				DATA2(j,i) = 0.0; 

		}
	}
	

	//------- Step 4: remove choose only the first peak: 

	for (i = 0; i < Np; i++)
	{
		fidx = 0;
		while ( DATA2(fidx,i) <=0 ) 
			fidx += 1;
		for (j=fidx; j < NewNoSamples; j++)
			DATA2(j,i) = 0.0;

	}


 	// write the data in A to the output file:
 	write_matrix_to_file(outputfilename,DATA);

 	cout << "Data has been written to the output file: " << outputfilename << endl;

 	return 0;	
}



