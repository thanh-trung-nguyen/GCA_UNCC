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

double Maxim(Vec_real VecValue)
{
  double maxim = VecValue(0);
  unsigned int i; 
  for (i=1;i < VecValue.size();i++)
    if( VecValue(i) > maxim)
      maxim = VecValue(i);
  return maxim;
}

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

  	Mat_real totalwave, incwave, scatwave,scatwave2, incwave2;
  	totalwave.newsize(NoTimeSteps,NoBndNodes);
  	incwave.newsize(NoTimeSteps,NoBndNodes);
    	scatwave.newsize(NoTimeSteps,NoBndNodes);
 	incwave2.newsize(NoTimeSteps,NoBndNodes);
    	scatwave2.newsize(NoTimeSteps,NoBndNodes);

  	MV_Vector<int> BdNodes(NoBndNodes);	
  	cout << "Loading the input file, wait "<< endl;

  	out.load_boundary_data(obj_dat_file,totalwave,BdNodes);
  	out.load_boundary_data(hom_dat_file,incwave,BdNodes);

	// compute the boundary conditions:
	WavesOptional opt;
	Mat_real BndCond_q, BoundCond; 
	BndCond_q.newsize(Ns+1,NoBndNodes);
	BoundCond.newsize(Ns+2,NoBndNodes);
	double dt = maxTime/NoTimeSteps;
	
	int j;
 	for (i = 0; i < NoTimeSteps; i++)
 	{
 		for (j=0; j < NoBndNodes; j++)
		{
			scatwave(i,j) = totalwave(i,j) - incwave(i,j);
			scatwave2(i,j) = scatwave(i,j)*(i+1)*dt;
			incwave2(i,j) = incwave(i,j)*(i+1)*dt;			
		}
	} 
		
	Vec_real ws(NoBndNodes), ws2(NoBndNodes), wi(NoBndNodes), wi2(NoBndNodes);
	
	for (j = 0; j < NoBndNodes; j++)
		BoundCond(0,j) = BdNodes(j);
	for (i = 0; i <= Ns; i++)
	{
		opt.LaplaceTransform_col(scatwave, dt, s(i), ws);
		opt.LaplaceTransform_col(scatwave2, dt, s(i), ws2);
		opt.LaplaceTransform_col(incwave, dt, s(i), wi);
		opt.LaplaceTransform_col(incwave2, dt, s(i), wi2);

		for (j = 0; j < NoBndNodes; j++)
			BndCond_q(i,j) = (ws(j)*wi2(j) - ws2(j)*wi(j))/pow(s(i),2)/wi(j)/(wi(j) + ws(j)) - 2*log(1 + ws(j)/wi(j))/pow(s(i),3);
		if (i==0)
		{
			for (j = 0; j < NoBndNodes; j++)
				BoundCond(Ns+1,j) = log(1 + ws(j)/wi(j))/pow(s(i),2);//bdc for the first tail
		}
	}
	
// 	out.Write_matrix_to_file(output_filename,BndCond_q);

	// compute the average (approximation of the integral:)
	for (i = 0; i < Ns; i++)
	{	for (j = 0; j < NoBndNodes; j++)
			BoundCond(i+1,j) = 0.5*(BndCond_q(i,j) + BndCond_q(i+1,j)); 
	}	



cout << BoundCond(1,0) << " " << BoundCond(2,0) << endl;
// 	out.Write_matrix_to_file(output_filename,BoundCond);


	// load the boundary grid file and extract the right measured data: 
	Vec_real X(NoBndNodes), Y(NoBndNodes), Z(NoBndNodes); 

	out.load_boundary_grid(bdgridfile,BdNodes,X,Y,Z); 
	
	string typeOfBndMeasurement = out.load_inversion_parameters_typebndmea(inv_par_file); // type of boundary measurement


	if (typeOfBndMeasurement.compare("z_max") == 0)//measured at z_max
	{
		cout << "Type of boundary measurement: " << typeOfBndMeasurement << endl;

		Mat_real BoundCond2; 
		BoundCond2.newsize(Ns+2,NoBndNodes);
		BoundCond2 = 0.0;
		Vec_real bdc_tail(NoBndNodes); 
		bdc_tail = -1e10;

		for (j = 0; j < NoBndNodes; j++)
			BoundCond2(0,j) = BoundCond(0,j);

		for (j=0; j < NoBndNodes; j++)
		{	
			if (fabs(Z(j) - z_maxFEM) < 0.0000001)
			{
				for (i=1; i <=Ns+1; i++)
					BoundCond2(i,j) = BoundCond(i,j);
				bdc_tail(j) = BoundCond(Ns+1,j); // extract the bdc for the first tail
			}	
		}
		double MaxTail = Maxim(bdc_tail);
		for (j=0; j < NoBndNodes; j++)
		{	
			if (fabs(Z(j) - z_maxFEM) > 0.0000001)
			{
				BoundCond2(Ns+1,j) = MaxTail;
			}	
		}
		

		out.Write_matrix_to_file(output_filename,BoundCond2);
	}
	else
		out.Write_matrix_to_file(output_filename,BoundCond); // full measurement
	
  	return 0;
}














