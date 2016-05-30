/* 
Function gca2_scatwave.C solves the coefficient identification problem using the NEW globally convergent algorithm
input: 
argv[1]: the parameter file with parameters of the forward solver. For experimental data we also need this parameter file according to the measurement setup. 
argv[2]: the file with parameters for the inversion, e.g. the regularization parameters, pseudo frequencies used,...

To compile: make gca_scatwave

Example to run the program: go to a folder inside the main GCA folder or make a new subfolder, then run:
../gca_scatwave ../parameters/par2_1obj_h02.dat ../parameters/inv_par.dat

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
int Maxim(Vec_int VecValue)
{
  int maxim = VecValue(0);
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
	{	double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		w(i) = exp(s*s*v(i))*value;
	}
		
}

// calculate the tail from the total wave: 
void Compute_tail_from_w(WavesGridB gg, double Z1, double s, double omega, Vec_real w, Vec_real& V)
{
	int NrFEMNodes = gg.getNoNodes();
	for (int i = 0; i < NrFEMNodes; i++)
	{	double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		V(i) = log(w(i)/value)/s/s;
	}
}

// calculate the tail from the scattered wave:
void Compute_tail_from_ws(WavesGridB gg, double Z1, double s, double omega, Vec_real ws, Vec_real& V)
{	
	int NrFEMNodes = gg.getNoNodes();
	for (int i = 0; i < NrFEMNodes; i++)
	{	double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
		V(i) = log(1 + ws(i)/value)/s/s;
	}	
}


// truncate the coefficient and remove the "artifacts" outside a given domain
void coefficient_truncate(Vec_real& Coefficient, double C_lb, double C_ub, WavesGridB& gg, Vec_real Subdomain)
{
	int NrFEMNodes = Coefficient.size();
	double Max = Maxim(Coefficient); cout << "Max of coefficient = " << Max << endl;
	for (int i=0; i< NrFEMNodes; i++)
	{	
		double x = gg.getCoor(i,0); 
		double y = gg.getCoor(i,1);
		double z = gg.getCoor(i,2); 
		if (x < Subdomain(0) || x > Subdomain(1) ||
		    y < Subdomain(2) || y > Subdomain(3) ||
		    z < Subdomain(4) || z > Subdomain(5))
			Coefficient(i) = 1.0;
		else
		{	if (Coefficient(i) < C_lb)
				Coefficient(i) = 1.0; // background value
			if (Coefficient(i) > C_ub)
				Coefficient(i) = C_ub;			
		}
	}
}

// compute the coefficient given the function v using the FDM 
void compute_coefficient_FDM(Vec_real& C, Vec_real v, WavesSDGeometry& sdg, double s)
{
	int Nx = sdg.getN_i();
	int Ny = sdg.getN_j();
	int Nz = sdg.getN_k();
	double dx = sdg.getDx();
	double dy = sdg.getDy();
	double dz = sdg.getDz();
	int i, j, k, node, nxl, nxu, nyl, nyu, nzl, nzu;
	double secder, grad;

	for (k=1; k< Nz-1; k++)
	{	for (j = 1; j < Ny-1; j++)
		{	for (i=1;i < Nx-1; i++)
			{	node = i + j*Nx + k*Nx*Ny;
				nxl = (i-1) + j*Nx + k*Nx*Ny; // the neighboring nodes
				nxu = (i+1) + j*Nx + k*Nx*Ny;
				nyl = i + (j-1)*Nx + k*Nx*Ny;
				nyu = i + (j+1)*Nx + k*Nx*Ny;
				nzl = i + j*Nx + (k-1)*Nx*Ny;
				nzu = i + j*Nx + (k+1)*Nx*Ny;
				
				secder = (v(nxl) - 2*v(node) + v(nxu))/dx/dx 
				       + (v(nyl) - 2*v(node) + v(nyu))/dy/dy
				       + (v(nzl) - 2*v(node) + v(nzu))/dz/dz;
				grad = 	pow((v(nxu) - v(nxl))/2/dx,2) + pow((v(nyu) - v(nyl))/2/dy,2) + pow((v(nzu) - v(nzl))/2/dz,2);
					
				C(node) = 1.0 + secder + grad*s*s + 2*s*(v(nzu) - v(nzl))/2/dz;
				if (C(node) < 1.0)
					C(node) = 1.0; 	

			}
		}

	}
}

void smooth_coefficient_FDM(Vec_real& C, WavesSDGeometry& sdg)
{
	int Nx = sdg.getN_i();
	int Ny = sdg.getN_j();
	int Nz = sdg.getN_k();
	double dx = sdg.getDx();
	double dy = sdg.getDy();
	double dz = sdg.getDz();
	int i, j, k, node, nxl, nxu, nyl, nyu, nzl, nzu;
	double secder, grad;
	
	int NrFEMNodes = C.size();
	Vec_real Cnew(NrFEMNodes);
	Cnew = 1.0;

	for (k=1; k< Nz-1; k++)
	{	for (j = 1; j < Ny-1; j++)
		{	for (i=1;i < Nx-1; i++)
			{	node = i + j*Nx + k*Nx*Ny;
				nxl = (i-1) + j*Nx + k*Nx*Ny; // the neighboring nodes
				nxu = (i+1) + j*Nx + k*Nx*Ny;
				nyl = i + (j-1)*Nx + k*Nx*Ny;
				nyu = i + (j+1)*Nx + k*Nx*Ny;
				nzl = i + j*Nx + (k-1)*Nx*Ny;
				nzu = i + j*Nx + (k+1)*Nx*Ny;								
				Cnew(node) = (C(nxl) + C(nxu) + C(nyl) + C(nyu) + C(nzl) + C(nzu) + 3*C(node))/9;

			}
		}

	}
	C = vector_copy(Cnew);
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
	opt.LaplaceTransform(fname, NoTimeSteps, dt, s, w);
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
	out.ReadSolfromFile(fname, ws);
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
	opt.LaplaceTransform(fname, NoTimeSteps, dt, s, w);
	opt.LaplaceTransform(fname_hom, NoTimeSteps, dt, s, wi);
	for (int i = 0; i < NrFEMNodes; i++)
	{      	double value = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1));
	    	V(i) = log(1 + (w(i)-wi(i))/value)/s/s;
	}

}


// Coefficient averaging: 
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
   {     if ( new_velocity(i) > max_value)
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
	// ----------------------------------------

	// input parameters files: 
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
	Vec_real Subdomain(6);

	out.load_inversion_parameters(inv_par_file, Lambda_CWF, RegPar, MaxIter, C_lb, C_ub, Subdomain);
	string firsttail = out.load_inversion_parameters_firsttail(inv_par_file);

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
	bdvalue_tail = 0.0;	

	// load boundary values for q and first Tail: 
	out.load_boundary_data_inversion(meadatafile,psi,BoundaryNodes, bdvalue_tail); 	// new boundary data structure in Laplace domain

	// ==  INITIALIZATION of the functions and solvers:
	
	Vec_real V(NrFEMNodes), dVdx(NrFEMNodes), dVdy(NrFEMNodes), dVdz(NrFEMNodes), V_old(NrFEMNodes), V2(NrFEMNodes);  // tail function and its derivatives
	Vec_real v(NrFEMNodes), w(NrFEMNodes), ws(NrFEMNodes), Wi(NrFEMNodes); // the functions v, w and Wi of homogeneous medium
	Vec_real q(NrFEMNodes), dqdx(NrFEMNodes), dqdy(NrFEMNodes), dqdz(NrFEMNodes);
	Vec_real Sum_q_prev(NrFEMNodes), Sum_dqdx_prev(NrFEMNodes), Sum_dqdy_prev(NrFEMNodes), Sum_dqdz_prev(NrFEMNodes);
	Vec_real Ax(NrFEMNodes), Ay(NrFEMNodes), Az(NrFEMNodes); // coefficient vector in the convection term of eq. for q
	Vec_real RHS_vec(NrFEMNodes); // right hand side vector in equation for q
	Vec_real boundary_values(NoBoundaryNodes); 
	Vec_real Coefficient(NrFEMNodes), Coeff_old(NrFEMNodes); //coefficient to be reconstructed
	Coefficient = 1.0; 

	WavesOptional opt(gg,USE_FEM,EXCHANGE,USE_FDM); // For saving the coefficient

	// Laplace transform of the simulated incident wave (used in calculating the tail function):
	out.ReadSolfromFile(incwavefile, Wi); // input data is in the Laplace transform domain

	//--------Compute the first tail--------------------: 
	if (firsttail.compare("exact") == 0) // true tail 
	{	string truesolFile = out.load_inversion_parameters_fp(inv_par_file); 
		const char* truesolfile = truesolFile.c_str();

		out.ReadSolfromFile(truesolfile, w);
		for (i = 0; i< NrFEMNodes; i++)
			V(i) = log(w(i)/Wi(i))/s_max/s_max;
		q_grad.compute_gradient(V, dVdx,  dVdy,  dVdz); // calculate the gradient of V

		remove("ExactTail.m");	out.Write_array_my_part("ExactTail.m", V, 0); // for Matlab visualization
		opt.writeInp3D("ExactTail.inp", &gg, V,1); // for AVS-viewer

	}
	else if (firsttail.compare("model2") == 0) // using model 2: tail is computed by solving the Laplace equation
	{	//for (i = 0; i< NrFEMNodes; i++)
		//	bdvalue_tail(i) = bdvalue_tail(i) + log(Wi(i))/s_max/s_max; // bdc for p a.w. the total wave		

		Ax = 0.0; Ay = 0.0; Az = 0.0; RHS_vec = 0.0;	
		p_Poisson.InitGlobConvFEM(Ax, Ay, Az, 0.0);  // create the stiffness matrix
		p_Poisson.PoissonEqSolverDirichletFEMnew(bdvalue_tail, sdg_FEM, RHS_vec); // solve the Poisson eq.
		p_Poisson.extract_solution(V); // extract the solution!
		p_Poisson.DirichletFEM(bdvalue_tail,V);  

		remove("FirstTail_model2.m"); out.Write_array_my_part("FirstTail_model2.m",V,0);
		opt.writeInp3D("FirstTail_model2.inp", &gg, V,1); // for AVS-viewer

		//for (i = 0; i< NrFEMNodes; i++)
		//	V(i) = V(i) - log(Wi(i))/s_max/s_max; // convert back to the tail V of the scattered wave
		
		V2 = vector_copy(V);
		//vector_multiplication(V2,-1.0/s_max);
		q_grad.compute_gradient(V2, dVdx,  dVdy,  dVdz); 		
	}
	else if (firsttail.compare("upperbound") == 0)
	{	Coefficient = -C_ub*Minim(bdvalue_tail)*3.0e3; cout << "Test coeff: " << Maxim(Coefficient) << endl; 
		coefficient_truncate(Coefficient, C_lb, C_ub, gg, Subdomain);
		WavesScalarEqOpOpt p_wave_eq2(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
		     				USE_RHS, USE_DIR_FEM, USE_DIR_FDM, 
		     				USE_ABSORB, PRINT_FILES, NoTimeSteps, maxTime, rhs, 
		     				type_of_mat, velocity, guess_velocity, doIncludeCorners);
		
		dt = p_wave_eq2.InitTime(); // important!!!
		p_wave_eq2.InitFDM();
		p_wave_eq2.Init3DFEM(Coefficient);
		p_wave_eq2.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
	
		MV_Vector<real>   E_FEM(NrFEMNodes);
		w = 0.0; t = 0.0;

		//=====+++ MAIN LOOP FORWARD SOLVER +++=====
		for (k = 0; k < NoTimeSteps; k++)
		{   	p_wave_eq2.PlaneWaveBackFDM(k,omega);
			p_wave_eq2.WaveEqSolverFEMforDifMat(t, k, type_of_mat, velocity);
			p_wave_eq2.ApplyExchangeCommon();
			p_wave_eq2.ApplySwap();
			p_wave_eq2.SaveSolutionsFEM(k,E_FEM); // extract the solution at time t
			t += dt;  
			double value = exp(-s_max*t)*dt;
			for (i = 0; i < NrFEMNodes; i++)		
				w(i) += E_FEM(i)*value;
		}
		ws = vector_subtraction(w,Wi); // subtract the Laplace transform of the incidient wave from the total wave
		Compute_tail_from_ws(gg, Z_max, s_max, omega, ws, V);
		V2 = vector_copy(V);
		vector_multiplication(V2,-1.0/s_max);	// test the tail
		q_grad.compute_gradient(V2, dVdx,  dVdy,  dVdz);// calculate the spatial derivatives of V	

		remove("FirstTail_upperbound.m"); out.Write_array_my_part("FirstTail_upperbound.m",V,0);
		opt.writeInp3D("FirstTail_upperbound.inp", &gg, V,1); // for AVS-viewer

	}
	else
	{		V = 0.0; dVdx = 0.0; dVdy = 0.0; dVdz = 0.0; // default choice of the first tail
	}

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
	Coefficient = 1.0;
	
	// check the presence of the exchangeMask.m file: 
	const char* exchangeMaskfile = "exchangeMask.m"; //test if the exchangeMask file is available
	ifstream ifile(exchangeMaskfile);
	if (!ifile)
	{	cout << "ERROR: the exchangeMask.m file must be copied to the same folder" << endl;
		return 0;
	}


// some test output files:

remove("Test_q2.m"); remove("Test_v2.m"); remove("Test_w2.m"); remove("Test_c2.m"); remove("Test_error2.m"); remove("Test_tail2.m"); 

// error:

int Niter = 20;
Vec_real Error(Niter); Error = 0.0;  

// ===========================================================================================================
// ++++++++++ MAIN LOOP:
//============================================================================================================

for (int ii =0; ii < Niter; ii++) // repeating the whole process twice:
{	V_old = vector_copy(V); // for stopping only!
	Coeff_old = vector_copy(Coefficient); // for stopping only!

	// **** Step 1: solve the equation of function q:
	for (int n=0; n < Ns; n++) //loop w.r.t. pseudo frequencies
	{	//cout << "Frequency s= " << pseudoFreq(n+1) << endl;
		ComputeIntegrals(Lambda_CWF(n), pseudoFreq(n+1), pseudoFreq(n), A1n, A2n, A3n, A4n); // compute coefficient of eq. for q:

		// update the average of the previous q_j, j=0,...,n-1:
		update_sum_q_prev(ds, q, dqdx, dqdy, dqdz, Sum_q_prev, Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev);
		
		// extract the boundary values from psi: 
		for (i = 0; i < NoBoundaryNodes; i++)
			bd_vector(BoundaryNodes(i)-1) = psi(n,i);
		
		//****: SOLVE the Poisson equation for q_n: NOTE that the equation to be solved using the p_Poisson functions are 	
		// ( grad q, grad v)  +  ((Ax,Ay,Az)*grad q, v)  +  (\epsilon q, v) =   (RHS, v) + (dq_dn, v). 
		// Or equivalent form: \Delta q - (Ax,Ay,Az)*grad q - epsilon q = RHS with the Dirichlet b.c.

		ComputeRHS_q(A3n,A4n,Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev, dVdx, dVdy, dVdz, dqdx, dqdy, dqdz, RHS_vec); 
		ComputeCoeffMatrix_q(A1n, A2n, dVdx, dVdy, dVdz, Sum_dqdx_prev, Sum_dqdy_prev, Sum_dqdz_prev, Ax, Ay, Az);

		p_Poisson.InitGlobConvFEM(Ax, Ay, Az, RegPar(n));  // create the stiffness matrix
		p_Poisson.PoissonEqSolverDirichletFEMnew(bd_vector, sdg_FEM, RHS_vec); // solve the Poisson eq.
		p_Poisson.extract_solution(q); // extract the solution to vector q, q now is the UPDATED one!
		p_Poisson.DirichletFEM(bd_vector,q);  
		q_grad.compute_gradient(q, dqdx,  dqdy,  dqdz); // calculate the gradient of q
	} // end of FOR loop
	
	// **** Step 2: compute the coefficient:
	Compute_v_from_q(ds, q, Sum_q_prev, V, v);
	out.Write_array_my_part("Test_v2.m",v,ii);	

	//Compute the coefficient using FEM equation:
	Compute_w_from_v(gg, Z_max, s_min, omega, v, w);

	p_wave_eq.Init3DParameter();
	p_wave_eq.ParameterFEM2(w, s_min, Coefficient);

	// using the FDM approximation: 
//	compute_coefficient_FDM(Coefficient, v, sdg_FEM, s_min);
	smooth_coefficient_FDM(Coefficient,sdg_FEM);

//	coef_smoothing(gg, C_lb, C_ub, Coefficient);
	coefficient_truncate(Coefficient, C_lb, C_ub, gg, Subdomain); // apply the lower and upper bounds to the coefficient
	out.Write_array_my_part("Test_c2.m",Coefficient,ii);				

 	// save the new coefficient in the file of the form C_n_i.inp:
	i=sprintf(C_filename,"%s%d%s","Coef_",ii,".inp"); 
	opt.writeInp3D(C_filename, &gg, Coefficient,1); // save the coefficient
	
	// **** STEP 3: update the tail function: 
	if (firsttail.compare("exact") != 0)  // if exact first tail is used, no update of tail is performed!
	{	// initialize the forward solver: 
		WavesScalarEqOpOpt p_wave_eq2(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
		     				USE_RHS, USE_DIR_FEM, USE_DIR_FDM, 
		     				USE_ABSORB, PRINT_FILES, NoTimeSteps, maxTime, rhs, 
		     				type_of_mat, velocity, guess_velocity, doIncludeCorners);

		
		dt = p_wave_eq2.InitTime(); // important!!!
		p_wave_eq2.InitFDM();
		p_wave_eq2.Init3DFEM(Coefficient);
		p_wave_eq2.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
	
		MV_Vector<real>   E_FEM(NrFEMNodes);
		w = 0.0;
		t = 0.0;

		//=====+++ MAIN LOOP FORWARD SOLVER +++=====

		for (k = 0; k < NoTimeSteps; k++)
		{   	p_wave_eq2.PlaneWaveBackFDM(k,omega);
			p_wave_eq2.WaveEqSolverFEMforDifMat(t, k, type_of_mat, velocity);
			p_wave_eq2.ApplyExchangeCommon();
			p_wave_eq2.ApplySwap();

			//save global FEM solution
			p_wave_eq2.SaveSolutionsFEM(k,E_FEM);
			t += dt;  
			double value = exp(-s_max*t)*dt;
			for (i = 0; i < NrFEMNodes; i++)		
				w(i) += E_FEM(i)*value;

		}

		out.Write_array_my_part("Test_tail2.m",V,ii);	
		ws = vector_subtraction(w,Wi); // subtract the Laplace transform of the incidient wave from the total wave
		Compute_tail_from_ws(gg, Z_max, s_max, omega, ws, V);
		V2 = vector_copy(V);
		//vector_multiplication(V2,1.0/pseudoFreq(n+1));	// test the tail
		q_grad.compute_gradient(V2, dVdx,  dVdy,  dVdz);// calculate the spatial derivatives of V
	} //end of if firsttail
	// ----end of step 3---


	// **** stopping rules:
	double L2NormC = L2_norm(gg,Coefficient); cout << "L2-norm of coefficient: " << L2NormC << endl;
	Vec_Error = vector_subtraction(Coefficient,Coeff_old);
	Error(ii) = L2_norm(gg,Vec_Error)/(L2_norm(gg,Coeff_old) + 1e-15);	
	cout << " Iteration " << ii << " ******* Relative error =  " << Error(ii) << endl;
/*	if ((ii > 0) && (Maxim(Coefficient) > Maxim(Coeff_old)))
	{	cout << "Stopped due to increasing value of coefficient" << endl;
		break;
	}
/*	if ((ii > 0) && (Error(ii) > 1.2*Error(ii-1)))
	{	cout << "Stopped due to increasing of error" << end;
		break;
	}	
*/
} // end of repeating the whole process
out.Write_array_my_part("Test_error2.m",Error,0);	

// ============================================================================== 
// ------- end of MAIN LOOP
//===============================================================================

cout << endl << " +++++++++++++++++++++++++++++++++" << endl 
     << " Finish the test " << endl 
     << " +++++++++++++++++++++++++++++++++++ " <<endl; 

	// -----------finish calculation:
 	PetscFinalize();
	return 0;
}
