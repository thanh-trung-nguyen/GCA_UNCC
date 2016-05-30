/* 
Function glob_con_alg.C is to solve the coefficient identification problem using the globally convergent algorithm
input: 
argv[1]: the parameter file with parameters of the forward solver. For experimental data we also need this parameter file according to the measurement setup. 
argv[2]: the file with parameters for the inversion, e.g. the regularization parameters, pseudo frequencies used,...

To compile: make glob_con_alg

Example to run the program: go to a folder inside the main GCA folder or make a new subfolder, then run:
../glob_con_alg ../parameters/par2_1obj_h02.dat ../parameters/inv_par.dat

Nguyen Trung Thanh, UNCC, 2012. 
*/

static char help[] ="";

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "include/waveskraftwerk.hh"
#include "include/wavesSDGeometry.h"
#include "include/wavesElliptic_GlobConvAlg.h"
#include "include/wavesWaveEqOpElliptic_log.h"
#include "include/wavesPoissonEqOp.h"
#include "include/wavesOptional.h"  
#include "include/wavesSDOperator.h"
#include "include/wavesScalarEqOpOpt.h"
#include <petscksp.h>
#include <sys/times.h>

typedef double real; 
typedef MV_ColMat<real> Mat_real;
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;

const double Pi = 3.141592653589793; 



// +++++ Some vector functions:
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

Vec_real vector_copy(Vec_real Vold)
{
	int N = Vold.size();
	Vec_real Vnew(N);
	for (int i=0; i < N; i++)
		Vnew(i) = Vold(i);
	return Vnew;
}

Vec_real vector_copy(Vec_int Vold)
{
	int N = Vold.size();
	Vec_real Vnew(N);
	for (int i=0; i < N; i++)
		Vnew(i) = Vold(i);
	return Vnew;
}

void vector_addition(Vec_real& V, Vec_real Vadd)// add vector Vadd to V
{
	int N = V.size();
	int N2 = Vadd.size();
	if (N != N2)
		cout << "Error in vector_addition: two vectors must have the same length" << endl;
	else
	{
		for (int i = 0; i < N; i++)
			V(i) += Vadd(i);
	}
}
void vector_multiplication(Vec_real& V, double f)// add vector Vadd to V
{
	int N = V.size();
	for (int i = 0; i < N; i++)
		V(i) *= f;
}

Vec_real vector_subtraction(Vec_real V, Vec_real Vsub)// subtract vector Vadd from V
{	
	int N = V.size();
	int N2 = Vsub.size();
	if (N != N2)
	{
		cout << "Error in vector_subtraction: two vectors must have the same length" << endl;
		return 0.0;
	}
	else
	{
		Vec_real Vnew(N);
		for (int i = 0; i < N; i++)
			Vnew(i) = V(i) - Vsub(i);
		return Vnew;
	}
	
}


// +++++++++++ calculate the weight coefficients in the equation for q (see my note on July 22, 2012): different from the book!!! 
void  ComputeIntegrals(double lambda, double s1, double s2,
                       double& A1, double&  A2, double& A3, double&A4)
{
  double h = s2	- s1; // pseudo-frequency step size
  double CWF = exp(-lambda*h);
  double I0 = (1 - CWF)/lambda;
  double I1 = (s2 - s1*CWF - I0)/lambda;
  double I2 = (s2*s2 - s1*s1*CWF - 2*I1)/lambda;
  double I3 = (pow(s2,3) - pow(s1,3)*CWF - 3*I2)/lambda;	

  A1 = (6*I2 - 4*s2*I1)/I0;
  A2 = (4*I1 - 2*s2*I0)/I0;
  A3 = (4*I3 - 6*s2*I2 + 2*s2*s2*I1)/I0;	
  A4 = 2*I1/I0;

}

// laplace transform of the function f in the boundary condition of the forward problem
// this depends on the function f!!! need to be checked for each problem setup!!!
double lapltr_incident_func(double omega, double s)
{
	double laptr = omega/(s*s + omega*omega)*(1 - exp(-2*Pi*s/omega)); 
	return laptr; 
}


// Compute the boundary values of q and tail from the boundary measured data: data is the scattered wave in the time domain
void ComputeBoundaryCondition_q(char* datafile, char* bdgrid, double Z1, int NoTimeSteps, 
				double dt, Vec_real s, double omega, Mat_real& psi, Vec_int& BdNodes, Vec_real& bdvalue_tail)
{
	
	WavesOutputs out_ex; //class containing functions for input-output
  	int NoBndNodes = psi.size(1); // load the number of boundary grid points
	int Ns = psi.size(0); // number of pseudo-frequencies, exclude the maximum freq!!! 
	
	Mat_real time_data;
  	time_data.newsize(NoTimeSteps,NoBndNodes);

	Mat_real psi2;
	psi2.newsize(Ns+1,NoBndNodes);
  
	Vec_real X(NoBndNodes), Y(NoBndNodes), Z(NoBndNodes); 

	// load the time-domain data on the boundary: 
	out_ex.load_boundary_data(datafile,time_data,BdNodes);

	// load the boundary grid data: 
	out_ex.load_boundary_grid(bdgrid,BdNodes,X,Y,Z); 

	double ds = fabs(s(1) - s(2))/10;
	double s2;
	int i, j;
	WavesOptional opt;

	Vec_real LaplaceTr(NoBndNodes), loga(NoBndNodes), loga2(NoBndNodes);
	Vec_real lapltr_inc_wave(NoBndNodes);

	// calculate the function psi2 at the pseudo-frequencies:
	for (i = 0; i <= Ns; i++)
	{
		// calculate the values of ln() at s(i):
		for (j = 0; j < NoBndNodes; j++)
			lapltr_inc_wave(j) = lapltr_incident_func(omega,s(i))*exp(s(i)*(Z(j)-Z1));
 	
		opt.LaplaceTransform_col(time_data, dt, s(i), LaplaceTr);
		for (j=0; j < NoBndNodes; j++)
			loga(j) = log(1 + LaplaceTr(j)/lapltr_inc_wave(j));

		// calculate the boundary values of the tail:
		if (i==0)
		{			
			for (j = 0; j < NoBndNodes; j++)
				bdvalue_tail(BdNodes(j)-1) = loga(j)/pow(s(0),2);
		}

		// calculate the values of ln() at s(i) + ds2 for approximating the derivative w.r.t. s:
		s2 = s(i) + ds; 
		for (j = 0; j < NoBndNodes; j++)
			lapltr_inc_wave(j) = lapltr_incident_func(omega,s2)*exp(s2*(Z(j)-Z1));
 	
		opt.LaplaceTransform_col(time_data, dt, s2, LaplaceTr);
		for (j=0; j < NoBndNodes; j++)
			loga2(j) = log(1 + LaplaceTr(j)/lapltr_inc_wave(j));

		// calculate the function psi at s(i):		
		for (j=0; j < NoBndNodes; j++)
			psi2(i,j) = (loga2(j) - loga(j))/ds/s(i)/s(i) - 2*loga(j)/pow(s(i),3);
	} // end of for (i)

	// calculate the integrals over the pseudo frequency intervals:
	for (j=0; j < NoBndNodes; j++)
	{	
		for (i = 0; i < Ns; i++)
		psi(i,j) = (psi2(i,j) + psi2(i+1,j))/2; 
	}
}



// compute the right hand side vector of equation for q:
void ComputeRHS_q(double A3n, double A4n, Vec_real Sum_dqdx, Vec_real Sum_dqdy, Vec_real Sum_dqdz, 
		      Vec_real dVdx, Vec_real dVdy, Vec_real dVdz, Vec_real dqdx, Vec_real dqdy, Vec_real dqdz, Vec_real& RHS_vec)
{
	int NrFEMNodes = RHS_vec.size(); 
	int i;
	double v1, v2, v3, v4;
	for (i=0; i < NrFEMNodes; i++)
	{
		v1 = dqdx(i)*dqdx(i) + dqdy(i)*dqdy(i) + dqdz(i)*dqdz(i);
		v2 = Sum_dqdx(i)*Sum_dqdx(i) + Sum_dqdy(i)*Sum_dqdy(i) + Sum_dqdz(i)*Sum_dqdz(i);
		v3 = dVdx(i)*dVdx(i) + dVdy(i)*dVdy(i) + dVdz(i)*dVdz(i);
		v4 = dVdx(i)*Sum_dqdx(i) + dVdy(i)*Sum_dqdy(i) + dVdz(i)*Sum_dqdz(i);
		RHS_vec(i) = -(-A3n*v1 + 2*Sum_dqdz(i) + 2*dVdz(i) - A4n*(v2 + v3 - 2*v4)); // the minus sign is to make it consistent with the setup in the Poisson equation solver in class WavesPoissonEq.
	}

}

// Compute the coefficient matrix for the equation for q:
void ComputeCoeffMatrix_q(double A1n, double A2n, Vec_real dVdx, Vec_real dVdy, Vec_real dVdz, 
			  Vec_real Sum_dqdx, Vec_real Sum_dqdy, Vec_real Sum_dqdz, Vec_real& Ax, Vec_real& Ay, Vec_real& Az)
{
	int NrFEMNodes = Ax.size();
	for (int i = 0; i < NrFEMNodes; i++)
	{
		Ax(i) = -A1n*(dVdx(i) - Sum_dqdx(i)); 
		Ay(i) = -A1n*(dVdy(i) - Sum_dqdy(i)); 
		Az(i) = -(A1n*(dVdz(i) - Sum_dqdz(i)) + A2n); // the minus sign, see comment in ComputeRHS_q 
	}
}


// calculate the function v in the algorithm:
void Compute_v_from_q(double ds, Vec_real q, Vec_real Sum_q_prev, Vec_real V, Vec_real& v)
{
	for (int i=0; i< q.size(); i++)
		v(i) = -ds*q(i) - Sum_q_prev(i) + V(i);

}

// Compute the function w = wi*exp(s^2v):
void Compute_w_from_v(WavesGridB gg, double Z1, double s, double omega, Vec_real v, Vec_real& w)
{
	
	int NrFEMNodes = gg.getNoNodes();
	for (int i = 0; i < NrFEMNodes; i++)
	{
		double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		w(i) = exp(s*s*v(i))*value;
	}
		
}

// calculate the tail from the total wave: 
void Compute_tail_from_w(WavesGridB gg, double Z1, double s, double omega, Vec_real w, Vec_real& V)
{
	
	int NrFEMNodes = gg.getNoNodes();
	for (int i = 0; i < NrFEMNodes; i++)
	{
		double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		V(i) = log(w(i)/value)/s/s;
	}
}

// calculate the tail from the scattered wave:
void Compute_tail_from_ws(WavesGridB gg, double Z1, double s, double omega, Vec_real ws, Vec_real& V)
{
	
	int NrFEMNodes = gg.getNoNodes();
	for (int i = 0; i < NrFEMNodes; i++)
	{	
		double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		V(i) = log(1 + ws(i)/value)/s/s;
	}	
}

void coefficient_truncate(Vec_real& Coefficient, double C_lb, double C_ub)
{
	int NrFEMNodes = Coefficient.size();
	for (int i=0; i< NrFEMNodes; i++)
	{	if (Coefficient(i) < C_lb)
			Coefficient(i) = 1.0; // background value
		else
		{	Coefficient(i) = Coefficient(i);
		if (Coefficient(i) > C_ub)
			Coefficient(i) = C_ub;
		}
	}
}


//update the sum of the previously calculated q_j, j = 0..n-1:
void update_sum_q_prev(double ds, Vec_real q, Vec_real dqdx, Vec_real dqdy, Vec_real dqdz,
			Vec_real& Sum_q_prev, Vec_real& Sum_dqdx_prev, Vec_real& Sum_dqdy_prev,Vec_real& Sum_dqdz_prev)
{
	
	int NrFEMNodes = q.size();
	for (int i = 0; i< NrFEMNodes; i++)
	{
		Sum_q_prev(i) 	 += q(i)*ds;	
		Sum_dqdx_prev(i) += dqdx(i)*ds;	
		Sum_dqdy_prev(i) += dqdy(i)*ds;	
		Sum_dqdz_prev(i) += dqdz(i)*ds;	
	}
}


// compute the true tail function from the simulated solution of the forward problem
void Compute_true_tail(char *fname, char* fname_hom, int NoTimeSteps, double dt, double s, double omega, Vec_real& V)
{
	WavesOptional opt; 
	int NrFEMNodes = V.size();
	Vec_real w(NrFEMNodes), wi(NrFEMNodes); 
cout << "Computing the laplace transform of the object data: " << endl;
	opt.LaplaceTransform(fname, NoTimeSteps, dt, s, w);
cout << "Computing the laplace transform of the object data: " << endl;
	opt.LaplaceTransform(fname_hom, NoTimeSteps, dt, s, wi);
	for (int i = 0; i< NrFEMNodes; i++)
		V(i) = log(1 + (w(i)-wi(i))/wi(i))/s/s;

}

// compute the true tail function from the Laplace transform of the simulated data:
void Compute_true_tail(char *fname, char* fname_hom, double s, Vec_real& V)
{
	WavesOutputs out; 
	int NrFEMNodes = V.size();
	Vec_real ws(NrFEMNodes), wi(NrFEMNodes); 
cout << "Computing the laplace transform of the object data: " << fname << endl;
	out.ReadSolfromFile(fname, ws);
cout << "Computing the laplace transform of the object data: " << fname_hom << endl;
	out.ReadSolfromFile(fname_hom, wi);
	for (int i = 0; i< NrFEMNodes; i++)
		V(i) = log(1 + ws(i)/wi(i))/s/s;

}


// compute the true tail function from the simulated solution of the forward problem, a bit different from the above function
void Compute_true_tail(char *fname, char* fname_hom, WavesGridB gg, double Z1, int NoTimeSteps, double dt, double s, double omega, Vec_real& V)
{
	WavesOptional opt; 
	int NrFEMNodes = V.size();
	Vec_real w(NrFEMNodes), wi(NrFEMNodes); 
cout << "Computing the laplace transform of the object data: " << fname << endl;
	opt.LaplaceTransform(fname, NoTimeSteps, dt, s, w);
cout << "Computing the laplace transform of the incident wave: " << fname_hom << endl;
	opt.LaplaceTransform(fname_hom, NoTimeSteps, dt, s, wi);
	for (int i = 0; i < NrFEMNodes; i++)
	{ 
        	double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
	    	V(i) = log(1 + (w(i)-wi(i))/value)/s/s;
//		wi(i) = value; 
	}
//WavesOutputs out; 
//out.Write_array_my_part("Test_wi2.m",wi,0);

}



// Smoothing the coefficient by averaging: 
void coef_smoothing(WavesGridB gg, double min_value, 
	       double max_value,
               MV_Vector<double>& new_velocity)
{
  
  int i,j,el,n1,n2,n3,n4;
  double value, background;
  background = 1.0;
  Vec_real old_velocity(new_velocity.size());
  old_velocity = vector_copy(new_velocity); 	
 
  int  nsd = gg.getNoSpaceDim();
 
      for (el=0; el < gg.getNoElms();el++)
	{
	  if (nsd == 2)
	    { 
	      n1 = gg.loc2glob(el,0);
	      n2 = gg.loc2glob(el,1);
	      n3 = gg.loc2glob(el,2);

	      if (old_velocity(n1) > min_value &&
		  old_velocity(n2) > min_value &&
		  old_velocity(n3) > min_value)
		{
		  
		  value = (old_velocity(n1) +
			   old_velocity(n2) +
			   old_velocity(n3))/3.0;
		  new_velocity(n1) = value;
		  new_velocity(n2) = value;
		  new_velocity(n3) = value;
		  
		}
	    }
	  else if (nsd == 3)
	    {
	      n1 = gg.loc2glob(el,0);
	      n2 = gg.loc2glob(el,1);
	      n3 = gg.loc2glob(el,2);
	      n4 = gg.loc2glob(el,3);

	      if (old_velocity(n1) > min_value &&
		  old_velocity(n2) > min_value &&
		  old_velocity(n3) > min_value && 
		  old_velocity(n4) > min_value)
		{
		  value = (old_velocity(n1) +
			   old_velocity(n2) +
			   old_velocity(n3) +
			   old_velocity(n4) )/4.0;
		  new_velocity(n1) = value;
		  new_velocity(n2) = value;
		  new_velocity(n3) = value;
		  new_velocity(n4) = value;
	
		}
	    } //nsd == 3
	}
 
 for (int i=0;i <  new_velocity.size();i++)
   {
     if ( new_velocity(i) > max_value)
       new_velocity(i) = max_value;
     
       if ( new_velocity(i) < min_value)
       new_velocity(i) = background;
       
     
   }
}

// Compute the L2-norm of a function given its values at the grid nodes of the FEM mesh
double L2_norm(WavesGridB gg, MV_Vector<double> VecValue)
{
  int i,n,el;
  int ierr;
  int  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();
  int n_1,n_2,n_3,n_4;

  double sol=0.0;

  if(nsd==2){
    FET3n2D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTRI1) {
	F.refill(el);
	real volume = F.area();
	n_1 = F.n1();
	n_2 = F.n2();
	n_3 = F.n3();
	double val1 =  VecValue(n_1)*VecValue(n_1);
	double val2 =  VecValue(n_2)*VecValue(n_2);
	double val3 =  VecValue(n_3)*VecValue(n_3);
      sol+= ((val1 + val2 + val3)/3.0)*volume;
     
      }
    }
  }
  else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());
	n_1 = F.n1();
	n_2 = F.n2();
	n_3 = F.n3();
	n_4 = F.n4();
	
	volume *= 0.25;

       double val1 =  VecValue(n_1)*VecValue(n_1);
       double val2 =  VecValue(n_2)*VecValue(n_2);
       double val3 =  VecValue(n_3)*VecValue(n_3);
       double val4 =  VecValue(n_4)*VecValue(n_4);

      sol+= (val1 + val2 + val3 + val4)*volume;
	
      }
    }
  }

  
  sol = sqrt(sol);
  return(sol);

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



 	const char* forw_par_file = argv[1];
	const char* inv_par_file = argv[2];


// == step 0: load the geometrical and inversion parameters:

 	bool USE_FEM, USE_FDM, USE_DIR_FEM, USE_DIR_FDM, EXCHANGE, USE_ABSORB;
	bool doIncludeCorners=1;

 	double maxTime, omega, dt, t, value;
 	int nsd, NoTimeSteps;
 
 	WavesGridB gg, outer_gg; // in fact outer_gg is the same as gg!!!
 	WavesSDGeometry sdg_FDM, sdg_FEM;	
	WavesOutputs out; 
 
	//========== load the geometrical and forward solver parameters: 
	out.Configure(forw_par_file, gg, outer_gg, sdg_FDM, sdg_FEM, nsd, 
			 EXCHANGE, USE_FEM, USE_FDM, USE_DIR_FEM, USE_DIR_FDM, USE_ABSORB,
			 NoTimeSteps, maxTime, omega); 

	// ========= initialize the forward solver: 
	bool USE_RHS = false; bool PRINT_FILES = false; 
	int type_of_mat = 3; double velocity = 3.0; double guess_velocity = 1.0; // these parameters can be arbitrary!!!
	double rhs = 1.0;

	WavesScalarEqOpOpt p_wave_eq(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
				     USE_RHS, USE_DIR_FEM, USE_DIR_FDM, USE_ABSORB, PRINT_FILES, 
				     NoTimeSteps, maxTime, rhs, type_of_mat, velocity, guess_velocity, doIncludeCorners);

	// ========= initialize the poisson equation solver:
	bool USE_FEM_P = true; bool USE_RHS_P = true; bool USE_DIR_FEM_P = true;	bool PRINT_FILES_P = false; 
	double EVALPT_X = 0.0;	double EVALPT_Y = 0.0;	double EVALPT_Z = 0.0;	
	int noQpts = 4; // number of quadrature points for each element, see file poisson_UNCC_h04.dat of Larisa

	WavesPoissonEqOp  p_Poisson(gg, USE_FEM_P, USE_RHS_P, USE_DIR_FEM_P, PRINT_FILES_P, EVALPT_X, EVALPT_Y, EVALPT_Z, rhs, noQpts);
	WavesElliptic_gradient  q_grad(&sdg_FEM); // for calculating the spatial derivatives
		
	// ========= load the inversion parameters:
	double s_min, s_max, C_lb, C_ub;
	int Ns, i, k; // number of pseudo-frequencies

	string meadataFile = out.load_inversion_parameters(inv_par_file, s_min, s_max, Ns);
	const char* meadatafile = meadataFile.c_str();

	string incwaveFile = out.load_inversion_parameters(inv_par_file); // file with incident wave
	const char* incwavefile = incwaveFile.c_str();

	Vec_int  MaxIter(Ns);
	Vec_real Lambda_CWF(Ns), RegPar(Ns);

	out.load_inversion_parameters(inv_par_file, Lambda_CWF, RegPar, MaxIter, C_lb, C_ub);


// == step 0.1: load the measured data, calculate the boundary values: 

	string bdgridFile = out.Configurate(forw_par_file); //the boundary grid file name
	const char* bdgridfile = bdgridFile.c_str();

	int NoBoundaryNodes = out.load_boundary_grid(bdgridfile);

	double Z_max = sdg_FDM.getZ(sdg_FDM.getN_k()-1);  
	
	dt = maxTime/NoTimeSteps; 
	
	Vec_real pseudoFreq(Ns+1); // vector of pseudo-frequencies	
	for (i = 0; i <= Ns; i++)
		pseudoFreq(i) = s_max + i*(s_min - s_max)/Ns; 

	Mat_real psi; //!!!!!!!!!!!!!!! boundary values for q_n for n = 1, 2, ... Ns.
	psi.newsize(Ns,NoBoundaryNodes);

	Vec_int BoundaryNodes(NoBoundaryNodes);
	
	int NrFEMNodes = gg.getNoNodes(); // number of FEM nodes
	Vec_real bdvalue_tail(NrFEMNodes); //boundary values for tail: if we use model 2 to approximate the first tail
	bdvalue_tail = -1.0;	

//	ComputeBoundaryCondition_q(meadatafile, bdgridfile, Z_max, NoTimeSteps, dt, pseudoFreq, omega, psi,BoundaryNodes,bdvalue_tail);
	
	// new boundary data structure: 
	out.load_boundary_data_inversion(meadatafile,psi,BoundaryNodes, bdvalue_tail);
		

// ==  INITIALIZATION of the functions and solvers:
	
	Vec_real V(NrFEMNodes), dVdx(NrFEMNodes), dVdy(NrFEMNodes), dVdz(NrFEMNodes), V_old(NrFEMNodes);  // tail function and its derivatives
	Vec_real v(NrFEMNodes), w(NrFEMNodes), ws(NrFEMNodes), Wi(NrFEMNodes); // the functions v, w and Wi of homogeneous medium
	Vec_real q(NrFEMNodes), dqdx(NrFEMNodes), dqdy(NrFEMNodes), dqdz(NrFEMNodes);
	Vec_real Sum_q_prev(NrFEMNodes), Sum_dqdx_prev(NrFEMNodes), Sum_dqdy_prev(NrFEMNodes), Sum_dqdz_prev(NrFEMNodes);
	Vec_real Ax(NrFEMNodes), Ay(NrFEMNodes), Az(NrFEMNodes); // coefficient vector in the convection term of eq. for q
	Vec_real RHS_vec(NrFEMNodes); // right hand side vector in equation for q
	Vec_real boundary_values(NoBoundaryNodes); 
	Vec_real Coefficient(NrFEMNodes), Coeff_old(NrFEMNodes); //coefficient to be reconstructed




/*
// TEST++++++++++++++++++++

NrFEMNodes = gg.getNoNodes();
Vec_real initialguess(NrFEMNodes);
Vec_real exactsol(NrFEMNodes);
initialguess = 1.0;
MV_Vector<int> markNodes;
  opt.findBoundaryNodes(gg,markNodes); 

Ax = 0.0; Ay = 0.0;  Az = 0.0; RHS_vec = 0.0;
 
// boundary values:
for (i=0; i < NrFEMNodes; i++)
{	
	double x = gg.getCoor(i,0); 
	double y = gg.getCoor(i,1);
	double z = gg.getCoor(i,2);
	//exactsol(i) = x*x + y*y - 2*z*z;
	exactsol(i) = sin(x+y-sqrt(2)*z);
	if (markNodes(i) == 1)
		initialguess(i) = exactsol(i);
} 
//p_Poisson.InitGlobConvFEM(Ax, Ay, Az, 0.0);  // create the stiffness matrix
//p_Poisson.PoissonEqSolverDirichletFEMnew(initialguess, sdg_FEM, RHS_vec); // solve the Poisson eq.
	p_Poisson.InitPoissonFEM();
 	p_Poisson.PoissonEqSolverFEM(initialguess);
	p_Poisson.extract_solution(V); // extract the solution to vector q, q now is the UPDATED one!
	p_Poisson.DirichletFEM(initialguess,V);  

	
out.Write_array_my_part("Test1.m",initialguess,0);	
out.Write_array_my_part("Test2.m",V,0);	
out.Write_array_my_part("Test3.m",exactsol,0);	

// end TEST+++++++++++++++++++++
*/

	// calculate the Laplace transform of the simulated incident wave (used in calculating the tail function)
//	opt.LaplaceTransform(incwavefile, NoTimeSteps, dt, s_max, Wi); // input is in the time domain
	out.ReadSolfromFile(incwavefile, Wi); // input data is in the Laplace transform domain


	//--------!!! Initialization of the tail--------------------: 

	// way 1: zero = homogeneous medium
//	V = 0.0; dVdx = 0.0; dVdy = 0.0; dVdz = 0.0; 

	
	// way 2: true tail 
	string truesolFile = out.load_inversion_parameters_fp(inv_par_file); 
	const char* truesolfile = truesolFile.c_str();

	out.ReadSolfromFile(truesolfile, w);
	for (i = 0; i< NrFEMNodes; i++)
		V(i) = log(w(i)/Wi(i))/s_max/s_max;

//	Compute_true_tail(truesolfile, incwavefile, s_max, V);
	//Compute_true_tail(truesolfile, incwavefile, NoTimeSteps, dt, s_max, omega, V);//time domain input files
	//Compute_true_tail(truesolfile, incwavefile, gg, Z_max,  NoTimeSteps, dt, s_max, omega, V); // time domain input files
	q_grad.compute_gradient(V, dVdx,  dVdy,  dVdz); // calculate the gradient of V



/*	// way 3: using model 2: tail is computed by solving the Laplace equation
	// -------------------
	Ax = 0.0; Ay = 0.0; Az = 0.0; RHS_vec = 0.0;	
	p_Poisson.InitGlobConvFEM(Ax, Ay, Az, 0.0);  // create the stiffness matrix
	p_Poisson.PoissonEqSolverDirichletFEMnew(bdvalue_tail, sdg_FEM, RHS_vec); // solve the Poisson eq.
//	p_Poisson.InitPoissonFEM();
// 	p_Poisson.PoissonEqSolverFEM(bdvalue_tail);
	p_Poisson.extract_solution(V); // extract the solution to vector q, q now is the UPDATED one!
	p_Poisson.DirichletFEM(bdvalue_tail,V);  
	q_grad.compute_gradient(V, dVdx,  dVdy,  dVdz); 
	//--------
*/



//out.Write_array_my_part("bdvtail.m",bdvalue_tail,0);	

cout << endl << "Test tail: Max = " << Maxim(V) << " Min = " << Minim(V) << endl << endl;
out.Write_array_my_part("Test_tail.m",V,0);	



	// initialization of q at the highest frequency:
	q = 0.0; dqdx = 0.0; dqdy =0.0; dqdz = 0.0; 
	Sum_dqdx_prev = 0.0; Sum_dqdy_prev = 0.0; Sum_dqdz_prev = 0.0; Sum_q_prev = 0.0;

	bool STOP;
	int Iter; 
	double A1n, A2n, A3n, A4n; // integral coefficients
	double ds = fabs(pseudoFreq(0)-pseudoFreq(1));
	char C_filename[30];
	Vec_real   Vec_Error(NrFEMNodes);
	Vec_real   E_FEM(NrFEMNodes);
	Vec_real bd_vector(NrFEMNodes);  //vector of boundary value and initial guess for eq. for q
	bd_vector = 0.0; 
	Coefficient = 0.1;
	
	// check the presence of the exchangeMask.m file: 
	const char* exchangeMaskfile = "exchangeMask.m"; //test if the exchangeMask file is available
	ifstream ifile(exchangeMaskfile);
	if (!ifile)
	{	cout << "ERROR: the exchangeMask.m file must be copied to the same folder" << endl;
		return 0;
	}


	WavesOptional opt(gg,USE_FEM,EXCHANGE,USE_FDM); // For saving the coefficient

// ===========================================================================================================
// ++++++++++ MAIN LOOP:
//============================================================================================================

	for (int n=0; n < Ns; n++) //loop w.r.t. pseudo frequencies
	{	STOP = false;
		Iter = 0; 

		ComputeIntegrals(Lambda_CWF(n), pseudoFreq(n+1), pseudoFreq(n), A1n, A2n, A3n, A4n);

		// update the average of the previous q_j, j=0,...,n-1:
		update_sum_q_prev(ds, q, dqdx, dqdy, dqdz, Sum_q_prev, Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev);
		
		// extract the boundary values from psi: 
		for (i = 0; i < NoBoundaryNodes; i++)
			bd_vector(BoundaryNodes(i)-1) = psi(n,i);

cout << "max of bdc: " << Maxim(bd_vector) << " Min of bdc: " << Minim(bd_vector) << endl;

	
		while (STOP == false) // internal loop at each pseudo frequency
		{

			Iter++; 
			V_old = vector_copy(V); // for stopping only!
			Coeff_old = vector_copy(Coefficient); // for stopping only!

//==== STEP 1: SOLVE the Poisson equation for q_n: NOTE that the equation to be solved using the p_Poisson functions are 	
// ( grad q, grad v)  +  ((Ax,Ay,Az)*grad q, v)  +  (\epsilon q, v) =   (RHS, v) + (dq_dn, v). 
// Or equivalent form: \Delta q - (Ax,Ay,Az)*grad q - epsilon q = RHS with the Dirichlet b.c.

			ComputeRHS_q(A3n,A4n,Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev, dVdx, dVdy, dVdz, dqdx, dqdy, dqdz, RHS_vec); 
			ComputeCoeffMatrix_q(A1n, A2n, dVdx, dVdy, dVdz, Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev, Ax, Ay, Az);

			p_Poisson.InitGlobConvFEM(Ax, Ay, Az, RegPar(n));  // create the stiffness matrix
			p_Poisson.PoissonEqSolverDirichletFEMnew(bd_vector, sdg_FEM, RHS_vec); // solve the Poisson eq.
			p_Poisson.extract_solution(q); // extract the solution to vector q, q now is the UPDATED one!
			p_Poisson.DirichletFEM(bd_vector,q);  
			q_grad.compute_gradient(q, dqdx,  dqdy,  dqdz); // calculate the gradient of q

// === end of step 1.

//**************

// STEP 2: update the coefficient:
 
			Compute_v_from_q(ds, q, Sum_q_prev, V, v);
			Compute_w_from_v(gg, Z_max, pseudoFreq(n+1), omega, v, w);

//out.ReadSolfromFile(truesolfile, w);

			// calculate the coefficient using the FEM method:
			p_wave_eq.Init3DParameter();
     			p_wave_eq.ParameterFEM2(w, pseudoFreq(n+1), Coefficient);
		
cout << "Test Coefficient: max = " << Maxim(Coefficient) << ", min = " << Minim(Coefficient) << endl;
			coefficient_truncate(Coefficient, C_lb, C_ub); // apply the lower and upper bounds to the coefficient
//			coef_smoothing(gg,C_lb,C_ub,Coefficient); 

    			// save the new coefficient in the file of the form C_n_i.inp:
			i=sprintf(C_filename,"%s%d%s%d%s","Coef_",n+1,"_",Iter,".inp"); 
			cout << "File name of the coefficient: "; printf ("%s\n",C_filename);

			opt.writeInp3D(C_filename, &gg, Coefficient,1); // save the coefficient

// ===== end of step 2

// ****************

// === STEP 3: update the tail function: 

//test the computation of tail with exact parameters:
//out.ReadSolfromFile("object/exact_coefficient.m",Coefficient);

/*			// solve the forward problem, calculate the Laplace transform of the scattered wave:
			WavesScalarEqOpOpt p_wave_eq2(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
						 USE_RHS, USE_DIR_FEM, USE_DIR_FDM, USE_ABSORB, PRINT_FILES, 
						 NoTimeSteps, maxTime, rhs, type_of_mat, velocity, guess_velocity, doIncludeCorners);	
			p_wave_eq2.InitFDM();
			p_wave_eq2.Init3DFEM(Coefficient);
			p_wave_eq2.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
		
cout << " Test coefficient: Max = " << Maxim(Coefficient) << endl;
			w = 0.0; // laplace transform of the scattered wave
			t = 0.0;
			//main loop of the forward solver:
			for (k = 0; k < NoTimeSteps; k++)
			{   	      
				p_wave_eq2.PlaneWaveBackFDM(k,omega);
				p_wave_eq2.WaveEqSolverFEMforDifMat(t, k, type_of_mat, guess_velocity); 

				p_wave_eq2.ApplyExchangeCommon();
				p_wave_eq2.ApplySwap();

 				//Extract the solution at the current time step
				E_FEM = 0.0;
				p_wave_eq2.SaveSolutionsFEM(k,E_FEM);

cout << "max of E_FEM: " << Maxim(E_FEM) << endl;

				t += dt;  
				value = exp(-s_max*t)*dt;
				for (i = 0; i < NrFEMNodes; i++)		
					w(i) += E_FEM(i)*value;

			}	

*/




		// initialize the forward solver: 
		WavesScalarEqOpOpt p_wave_eq2(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
				     USE_RHS, USE_DIR_FEM, USE_DIR_FDM, 
				     USE_ABSORB, PRINT_FILES, NoTimeSteps, maxTime, rhs, 
				     type_of_mat, velocity, guess_velocity, doIncludeCorners);

		// time paramters:
		dt = p_wave_eq2.InitTime();
		t = 0.0;

		p_wave_eq2.InitFDM();
		p_wave_eq2.Init3DFEM(Coefficient);
		p_wave_eq2.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
	
		MV_Vector<real>   E_FEM(NrFEMNodes);
		w = 0.0;

		//============= ++++++++++ MAIN LOOP FORWARD SOLVER ++++++++++ ++=====================

		for (k = 0; k < NoTimeSteps; k++)
		{   	p_wave_eq2.PlaneWaveBackFDM(k,omega);
			p_wave_eq2.WaveEqSolverFEMforDifMat(t, k, type_of_mat, velocity);
			p_wave_eq2.ApplyExchangeCommon();
			p_wave_eq2.ApplySwap();

			//save global FEM solution
			p_wave_eq2.SaveSolutionsFEM(k,E_FEM);

			t += dt;  
			double value = exp(-s_max*t)*dt;
			for (int i = 0; i < NrFEMNodes; i++)		
				w(i) += E_FEM(i)*value;

		}


cout << "Max coeff = " << Maxim(Coefficient) << " min = " << Minim(Coefficient) << endl;
cout << "Freq = " << s_max << endl;
cout << "dt = " << dt << " maxTime = " << maxTime << " NoTimeSteps: " << NoTimeSteps << endl;


out.Write_array_my_part("Test_w.m",w,n); 
out.ReadSolfromFile(truesolfile, w); // test the accuracy of the calculation of tail: 

			ws = vector_subtraction(w,Wi); // subtract the Laplace transform of the incidient wave from the total wave

			Compute_tail_from_ws(gg, Z_max, s_max, omega, ws, V);
			q_grad.compute_gradient(V, dVdx,  dVdy,  dVdz);// calculate the spatial derivatives of V


//cout << endl << "Test tail: Max = " << Maxim(V) << " Min = " << Minim(V) << endl << endl;			
out.Write_array_my_part("Test_tail2.m",V,0);	

// check the stopping rules:
			Vec_Error = vector_subtraction(Coefficient,Coeff_old);
			double Error = L2_norm(gg,Vec_Error);
cout << endl << "Test the error: " << "==== " << Error << endl << endl; 
			if (( Iter >= MaxIter(n) ) || (Error < 0.01))
				STOP = true;


		} // end of while loop
	} // end of for n = ... loop

// ============================================================================== 
// ------- end of MAIN LOOP
//===============================================================================

cout << endl << " +++++++++++++++++++++++++++++++++" << endl <<" Finish the test " << endl << " +++++++++++++++++++++++++++++++++++ " <<endl; 



	// -----------finish calculation:
 	PetscFinalize();
	return 0;
}
