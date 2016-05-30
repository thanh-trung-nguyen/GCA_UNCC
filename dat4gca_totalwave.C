/* 
This function create the input data for the globally convergent algorithm for simulated data: 





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
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;

//============================================================ 
int main(int argc, char **argv)
{

 	char* forw_par_file = argv[1];
	char* inv_par_file  = argv[2];
	char* obj_dat_file  = argv[3];  
	char* hom_dat_file  = argv[4];  
		
  	WavesOutputs out; 

	// ========= load the boundary grid file name:
 	string bdgridFile = out.Configurate(forw_par_file);
	const char* bdgridfile = bdgridFile.c_str();

	// load the geometrical and temporal parameters:
	int NoTimeSteps;
  	double x_minFEM, y_minFEM, z_minFEM, x_maxFEM, y_maxFEM, z_maxFEM, maxTime;
	
	out.Configure(forw_par_file, x_minFEM, y_minFEM, z_minFEM, x_maxFEM, y_maxFEM, z_maxFEM, maxTime, NoTimeSteps);

	int NoBndNodes = out.load_boundary_grid(bdgridfile);

	// ========= load the inversion parameters:
	double s_min, s_max;
	int Ns; 

	string meadataFile = out.load_inversion_parameters(inv_par_file, s_min, s_max, Ns);
	const char* output_filename = meadataFile.c_str(); // output file name

	// pseudo frequencies:
	Vec_real s(Ns+1); // vector of pseudo-frequencies	
	int i;
	for (i = 0; i <= Ns; i++)
		s(i) = s_max + i*(s_min - s_max)/Ns; 

  
  	// load the time domain data from the input file:

  	Mat_real totalwave, totalwave2, incwave, incwave2;
  	totalwave.newsize(NoTimeSteps,NoBndNodes);
  	totalwave2.newsize(NoTimeSteps,NoBndNodes);
  	incwave.newsize(NoTimeSteps,NoBndNodes);
 	incwave2.newsize(NoTimeSteps,NoBndNodes);

  	MV_Vector<int> BdNodes(NoBndNodes);	
  	cout << "Loading the input file, wait "<< endl;

  	out.load_boundary_data(obj_dat_file,totalwave,BdNodes);
  	out.load_boundary_data(hom_dat_file,incwave,BdNodes);

	
	int j;
	double dt = maxTime/NoTimeSteps;
 	for (i = 0; i < NoTimeSteps; i++)
 	{
 		for (j=0; j < NoBndNodes; j++)
		{
			totalwave2(i,j) = totalwave(i,j)*(i+1)*dt;
			incwave2(i,j) = incwave(i,j)*(i+1)*dt;			
		}
	} 
		
	// load the boundary grid file (for extracting the right measured data): 
	Vec_real X(NoBndNodes), Y(NoBndNodes), Z(NoBndNodes); 
	out.load_boundary_grid(bdgridfile,BdNodes,X,Y,Z); 
	string typeOfBndMeasurement = out.load_inversion_parameters_typebndmea(inv_par_file); // type of boundary measurement


	// compute the boundary condition:
	Vec_real w(NoBndNodes), w2(NoBndNodes), wi(NoBndNodes), wi2(NoBndNodes);
	WavesOptional opt;
	Mat_real BndCond_q, BoundCond, bdc_incwave; 
	BndCond_q.newsize(Ns+1,NoBndNodes);
	BoundCond.newsize(Ns+2,NoBndNodes);
	bdc_incwave.newsize(Ns+2,NoBndNodes); // boundary condition of the incident wave

	for (j = 0; j < NoBndNodes; j++)
		BoundCond(0,j) = BdNodes(j); // boundary indices

	for (i = 0; i <= Ns; i++)
	{
		opt.LaplaceTransform_col(totalwave, dt, s(i), w);
		opt.LaplaceTransform_col(totalwave2, dt, s(i), w2);
		opt.LaplaceTransform_col(incwave, dt, s(i), wi);
		opt.LaplaceTransform_col(incwave2, dt, s(i), wi2);

		for (j = 0; j < NoBndNodes; j++)
		{	BndCond_q(i,j) = -w2(j)/pow(s(i),2)/w(j) - 2*log(w(j))/pow(s(i),3);
			bdc_incwave(i,j) = -wi2(j)/pow(s(i),2)/wi(j) - 2*log(wi(j))/pow(s(i),3);
		}
		if (i==0)
		{
			for (j = 0; j < NoBndNodes; j++)
			{	BoundCond(Ns+1,j) = log(w(j))/pow(s(i),2);//bdc for the first tail
				bdc_incwave(Ns+1,j) = log(wi(j))/pow(s(i),2);
			}
		}
	}

	// compute the average (approximation of the integral:)
	Mat_real bdc_incwave2; 
	bdc_incwave2.newsize(Ns+1,NoBndNodes);
	for (j = 0; j < NoBndNodes; j++)
	{	for (i = 0; i < Ns; i++)
		{	BoundCond(i+1,j) = 0.5*(BndCond_q(i,j) + BndCond_q(i+1,j)); 
			bdc_incwave2(i,j) = 0.5*(bdc_incwave(i,j) + bdc_incwave(i+1,j)); 
		}
		bdc_incwave2(Ns,j) = bdc_incwave(Ns+1,j);
	}	

	// extract the right measured data:
	if (typeOfBndMeasurement.compare("z_max") == 0)//measured at z_max
	{
		cout << "Type of boundary measurement: " << typeOfBndMeasurement << endl;

		for (j=0; j < NoBndNodes; j++)
		{	if (fabs(Z(j) - z_maxFEM) > 0.0000001)
			{
				for (i=1; i <=Ns+1; i++)
					BoundCond(i,j) = bdc_incwave2(i-1,j);
			}	
		}
	}
	out.Write_matrix_to_file(output_filename,BoundCond); // full measurement
	
  	return 0;
}














