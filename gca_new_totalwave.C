/*
Function gca_new_totalwave.C implements the new globally convergent algorithm 


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

Vec_int vector_copy(Vec_int Vold)
{
	int N = Vold.size();
	Vec_int Vnew(N);
	for (int i=0; i < N; i++)
		Vnew(i) = Vold(i);
	return Vnew;
}

Vec_real vector_addition(Vec_real V, Vec_real Vadd)// add vector Vadd to V
{
	int N = V.size();
	Vec_real Vnew(N); 
	Vnew = 0.0;
	int N2 = Vadd.size();
	if (N != N2)
		cout << "Error in vector_addition: two vectors must have the same length" << endl;
	else
	{	
		for (int i = 0; i < N; i++)
			Vnew(i) = V(i) + Vadd(i);
	}
	return Vnew;
}

Vec_real vector_subtraction(Vec_real V, Vec_real Vsub)// subtract vector Vadd from V
{	
	int N = V.size();
	Vec_real Vnew(N);
	Vnew = 0.0;
	int N2 = Vsub.size();
	if (N != N2)
		cout << "Error in vector_subtraction: two vectors must have the same length" << endl;
	else
	{
		for (int i = 0; i < N; i++)
			Vnew(i) = V(i) - Vsub(i);
	}
	return Vnew;
}

Vec_real vector_multiplication(Vec_real V, double f)// V*f
{
	int N = V.size();
	Vec_real Vnew(N);
	for (int i = 0; i < N; i++)
		Vnew(i) = V(i)*f;
	return Vnew;
}

// the weight coefficients in the equation for q
void  ComputeIntegrals(double lambda, double s1, double s2, double& A1, double&A4)
{
  double h = s2	- s1; // pseudo-frequency step size
  double CWF = exp(-lambda*h);
  double I0 = (1 - CWF)/lambda; 
  double I1 = (s2 - s1*CWF - I0)/lambda; 
  double I2 = (s2*s2 - s1*s1*CWF - 2*I1)/lambda; 
 
  A1 = (6*I2 - 4*s2*I1)/I0;
  A4 = 2*I1/I0;

}

// laplace transform of the function f(t) = sin(t*omega), t\in(0,2pi/omega), = 0 for t > 2pi/omega.
double lapltr_incident_func(double omega, double s)
{
	double laptr = omega/(s*s + omega*omega)*(1 - exp(-2*Pi*s/omega)); 
	return laptr; 
}

// compute q of the incident wave
Vec_real Compute_q_from_incidentwave(double s, double omega, WavesGridB& gg, WavesSDGeometry& sdg_FDM)
{
	double Z1 = sdg_FDM.getZ(sdg_FDM.getN_k()-1);
	int NrFEMNodes = gg.getNoNodes();
	Vec_real q0(NrFEMNodes);
	double ds = 0.01;  
	double s2 = s + ds; 
	double w, w2;
	for (int i = 0; i < NrFEMNodes; i++)
	{	w = lapltr_incident_func(omega,s)*exp(s*(gg.getCoor(i,2)-Z1))/s/s;
		w2 = lapltr_incident_func(omega,s2)*exp(s2*(gg.getCoor(i,2)-Z1))/s2/s2;
		q0(i) = (w2 - w)/ds;
	}
	return q0;
}

// compute the right hand side vector of equation for q using the total wave:
Vec_real ComputeRHS_q_totalwave(double A4n, Vec_real dVdx, Vec_real dVdy, Vec_real dVdz)
{
	int NrFEMNodes = dVdx.size(); 
	Vec_real RHS_vec(NrFEMNodes); 
	for (int i=0; i < NrFEMNodes; i++)
		RHS_vec(i) = A4n*(dVdx(i)*dVdx(i) + dVdy(i)*dVdy(i) + dVdz(i)*dVdz(i));// the + sign, see comment in step 1.
	return RHS_vec;
}


// Compute the coefficient matrix for the equation for q using the total wave:
void ComputeCoeffMatrix_q_totalwave(double A1n, Vec_real dVdx, Vec_real dVdy, Vec_real dVdz, 
			            Vec_real& Ax, Vec_real& Ay, Vec_real& Az)
{
	int NrFEMNodes = Ax.size();
	for (int i = 0; i < NrFEMNodes; i++)
	{	Ax(i) = -A1n*dVdx(i); 
		Ay(i) = -A1n*dVdy(i); 
		Az(i) = -A1n*dVdz(i); // the minus sign, see comment in ComputeRHS_q 
	}
}

// Compute the function w = exp(s^2v) (using the total wave):
Vec_real Compute_w_from_v_totalwave(double s, Vec_real v)
{	
	int NrFEMNodes = v.size();
	Vec_real w(NrFEMNodes);
	for (int i = 0; i < NrFEMNodes; i++)
		w(i) = exp(s*s*v(i));
	return w;
}

// calculate the tail from the total wave (using the total wave): 
Vec_real Compute_tail_from_w_totalwave(double s, Vec_real w)
{
	int NrFEMNodes = w.size();
	Vec_real V(NrFEMNodes);
	double s2 = s*s; 
	for (int i = 0; i < NrFEMNodes; i++)
	{	V(i) = log(w(i))/s2;
	}
	return V;
}


// compute the coefficient given the function v using the FDM 
void compute_coefficient_FDM_totalwave(Vec_real& C, Vec_real v, WavesSDGeometry& sdg, double s)
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
					
				C(node) = secder + grad*s*s;
				if (C(node) < 1.0)
					C(node) = 1.0; 	

			}
		}

	}
}


// truncate the coefficient and remove the "artifacts" outside a given domain
void coefficient_truncate(Vec_real& Coefficient, double C_lb, double C_ub, WavesGridB& gg, Vec_real Subdomain)
{
	int NrFEMNodes = Coefficient.size();
	//double Max = Maxim(Coefficient); cout << "Max of coefficient = " << Max << endl;
	for (int i=0; i< NrFEMNodes; i++)
	{	double x = gg.getCoor(i,0); 
		double y = gg.getCoor(i,1);
		double z = gg.getCoor(i,2); 
		if (x < Subdomain(0) || x > Subdomain(1) ||
		    y < Subdomain(2) || y > Subdomain(3) ||
		    z < Subdomain(4) || z > Subdomain(5))
			Coefficient(i) = 1.0;
		else
		{	if (Coefficient(i) < C_lb)
				Coefficient(i) = 1.0; // background value
			else if (Coefficient(i) > C_ub)
				Coefficient(i) = C_ub;			
		}
	}
}

void average_coefficient_FDM(Vec_real& C, WavesSDGeometry& sdg)
{
	int Nx = sdg.getN_i();
	int Ny = sdg.getN_j();
	int Nz = sdg.getN_k();
	int i, j, k, node, nxl, nxu, nyl, nyu, nzl, nzu;
	
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
				Cnew(node) = (C(nxl) + C(nxu) + C(nyl) + C(nyu) + C(nzl) + C(nzu) + C(node))/7.0;
			}
		}
	}
	C = vector_copy(Cnew);
}

// Compute the L2-norm of a function given its values at the grid nodes of the FEM mesh
double L2_norm(WavesGridB gg, MV_Vector<double> VecValue)
{
  int el;
  int  nsd = gg.getNoSpaceDim();
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


        // Load the geometrical and inversion parameters:

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
	dt = maxTime/NoTimeSteps; // time step

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

	// ===== load the measured data, calculate the boundary values: 

	string bdgridFile = out.Configurate(forw_par_file); //the boundary grid file name
	const char* bdgridfile = bdgridFile.c_str();

	int NoBoundaryNodes = out.load_boundary_grid(bdgridfile);

	Vec_real pseudoFreq(Ns+1); // vector of pseudo-frequencies	
	for (i = 0; i <= Ns; i++)
		pseudoFreq(i) = s_max + i*(s_min - s_max)/Ns; 

	Mat_real psi; //!!!!!!!!!!!!!!! boundary values for q_n for n = 1, 2, ... Ns.
	psi.newsize(Ns,NoBoundaryNodes);

	Vec_int BoundaryNodes(NoBoundaryNodes);
	
	int NrFEMNodes = gg.getNoNodes(); // number of FEM nodes


	// ==  INITIALIZATION of the functions and solvers:
	
	Vec_real V(NrFEMNodes), dVdx(NrFEMNodes), dVdy(NrFEMNodes), dVdz(NrFEMNodes), V2(NrFEMNodes);  // tail function and its derivatives
	Vec_real v(NrFEMNodes), w(NrFEMNodes), Wi(NrFEMNodes), w2(NrFEMNodes), wmax(NrFEMNodes); // the functions v, w and Wi of homogeneous medium
	Vec_real q(NrFEMNodes), E_FEM(NrFEMNodes);
	Vec_real Ax(NrFEMNodes), Ay(NrFEMNodes), Az(NrFEMNodes); // coefficient vector in the convection term of eq. for q
	Vec_real RHS_vec(NrFEMNodes); // right hand side vector in equation for q
	Vec_real Coefficient(NrFEMNodes), Coeff_old(NrFEMNodes); //coefficient to be reconstructed
	Coefficient = 1.0; 

	WavesOptional opt(gg,USE_FEM,EXCHANGE,USE_FDM); // For saving the coefficient
	// Laplace transform of incident wave (used in calculating the first tail):
	out.ReadSolfromFile(incwavefile, Wi); // input data is in the Laplace transform domain


	Vec_real bdc_tail(NrFEMNodes); //boundary values for tail: if we use model 2 to approximate the first tail
	bdc_tail = Compute_tail_from_w_totalwave(s_max,Wi);// initial guess for the first tail!	

	// load boundary values for q and first Tail: 
	out.load_boundary_data_inversion(meadatafile,psi,BoundaryNodes, bdc_tail); 	// new boundary data structure in Laplace domain


	//--------Compute the first tail--------------------: 
	if (firsttail.compare("exact") == 0) // true tail 
	{	string truesolFile = out.load_inversion_parameters_fp(inv_par_file); 
		const char* truesolfile = truesolFile.c_str();

		out.ReadSolfromFile(truesolfile, w);
		V = Compute_tail_from_w_totalwave(s_max,w);
		
		remove("ExactTail.m");	out.Write_array_my_part("ExactTail.m", V, 0); // for Matlab visualization
		opt.writeInp3D("ExactTail.inp", &gg, V,1); // for AVS-viewer
	}
	else if (firsttail.compare("model2") == 0) // using model 2: tail is computed by solving the Laplace equation
	    {	Ax = 0.0; Ay = 0.0; Az = 0.0; RHS_vec = 0.0;	
		p_Poisson.InitGlobConvFEM(Ax, Ay, Az, 0.0);  // create the stiffness matrix
		p_Poisson.PoissonEqSolverDirichletFEMnew(bdc_tail, sdg_FEM, RHS_vec); // solve the Poisson eq.
		p_Poisson.extract_solution(V); // extract the solution!
		p_Poisson.DirichletFEM(bdc_tail,V); // apply the bdc to the computed solution  

		remove("FirstTail_model2.m"); out.Write_array_my_part("FirstTail_model2.m",V,0);
		opt.writeInp3D("FirstTail_model2.inp", &gg, V,1); // for AVS-viewer	
	    }
	    else if (firsttail.compare("upperbound") == 0)
	    {	Coefficient = C_ub; cout << "Test coeff: " << Maxim(Coefficient) << endl; 
		coefficient_truncate(Coefficient, C_lb, C_ub, gg, Subdomain);
		
		dt = p_wave_eq.InitTime(); // important!!!
		p_wave_eq.InitFDM();
		p_wave_eq.Init3DFEM(Coefficient);
		p_wave_eq.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
	
		w = 0.0; t = 0.0;

		//=====+++ MAIN LOOP FORWARD SOLVER +++=====
		for (k = 0; k < NoTimeSteps; k++)
		{   	p_wave_eq.PlaneWaveBackFDM(k,omega);
			p_wave_eq.WaveEqSolverFEMforDifMat(t, k, type_of_mat, velocity);
			p_wave_eq.ApplyExchangeCommon();
			p_wave_eq.ApplySwap();
			p_wave_eq.SaveSolutionsFEM(k,E_FEM); // extract the solution at time t
			t += dt;  
			double value = exp(-s_max*t)*dt;
			for (i = 0; i < NrFEMNodes; i++)		
				w(i) += E_FEM(i)*value;
		}
		V = Compute_tail_from_w_totalwave(s_max, w);

		remove("FirstTail_upperbound.m"); out.Write_array_my_part("FirstTail_upperbound.m",V,0);
		opt.writeInp3D("FirstTail_upperbound.inp", &gg, V,1); // for AVS-viewer

	    }
	    else
	    {	//V = vector_copy(bdc_tail); 
		V = Compute_tail_from_w_totalwave(s_max,Wi);
		remove("FirstTail_hom.m"); out.Write_array_my_part("FirstTail_hom.m",V,0);
	    }
	


	bool STOP;
	int Iter; 
	double A1n, A4n; // integral coefficients
	double ds = fabs(pseudoFreq(0)-pseudoFreq(1));
	char C_filename[30];
	Vec_real   Vec_Error(NrFEMNodes);
	Vec_real bdc_q(NrFEMNodes);  //vector of boundary value and initial guess for q

	bdc_q = Compute_q_from_incidentwave(s_max, omega, gg, sdg_FDM); // first initial guess for tail from the homogeneous medium
	Coefficient = 1.0;
	
	// check the presence of the exchangeMask.m file: 
	const char* exchangeMaskfile = "exchangeMask.m"; //test if the exchangeMask file is available
	ifstream ifile(exchangeMaskfile);
	if (!ifile)
	{	cout << "ERROR: the exchangeMask.m file must be copied to the same folder" << endl;
		return 0;
	}


	// some test output files:

	remove("Test_q_new.m"); remove("Test_v_new.m"); remove("Test_w_new.m"); remove("Test_c_new.m"); remove("Test_tai_new.m"); 
	remove("Test_L2Coeff_new.m"); remove("Test_MaxCoeff_new.m"); remove("Test_RelativeError_new.m"); remove("Test_L2Error_new.m");
        remove("Test_error_new.m");

	Vec_real L2Error(Maxim(MaxIter)); Vec_real RelativeError(Maxim(MaxIter));
	Vec_real MaxCoeff(Maxim(MaxIter)); Vec_real L2Coeff(Maxim(MaxIter));
	Vec_real Coeff_n(NrFEMNodes), V_n(NrFEMNodes), Error(Ns); 



// ============================================
// ++++++++++ MAIN LOOP:
//=============================================

for (int ii =0; ii < 1; ii++) // repeating the whole process:
{	for (int n=0; n < Ns; n++) //loop w.r.t. pseudo frequencies
	{	STOP = false;
		Iter = 0; 
		L2Error = 0.0; RelativeError = 0.0; MaxCoeff = 0.0; L2Coeff = 0.0;
		double s_low = pseudoFreq(n+1); 
		double s_up = pseudoFreq(n);
		cout << endl << "FREQUENCY INTERVAL [s_n, s_{n-1}] = (" << s_low << ", " << s_up << ")" << ", n = " << n << endl;
		ComputeIntegrals(Lambda_CWF(n), s_low, s_up, A1n, A4n); // compute coefficient of eq. for q:
		
		// extract the boundary values from psi: 
		for (i = 0; i < NoBoundaryNodes; i++)
			bdc_q(BoundaryNodes(i)-1) = psi(n,i);
		
		V_n = vector_copy(V);
		Coeff_n = vector_copy(Coefficient); 

		while (STOP == false) // internal loop at each pseudo frequency
		{	Coeff_old = vector_copy(Coefficient); // for stopping only!

			V2 = vector_copy(V);
			//if (firsttail.compare("exact") != 0)	
			vector_multiplication(V2,-1.0/s_up);	// test using the derivative of tail
			q_grad.compute_gradient(V2, dVdx,  dVdy,  dVdz);// calculate the spatial derivatives of tail
		
			//**** STEP 1: find q_n: NOTE that the equation to be solved using the p_Poisson functions is: 	
			// ( grad q, grad v)  +  ((Ax,Ay,Az)*grad q, v)  +  (\epsilon*q, v) =   (RHS, v) + (dq_dn, v). 
			// Or equivalent form: \Delta q - (Ax,Ay,Az)*grad q - epsilon*q = -RHS with the Dirichlet b.c.

			RHS_vec = ComputeRHS_q_totalwave(A4n, dVdx, dVdy, dVdz); 
			ComputeCoeffMatrix_q_totalwave(A1n, dVdx, dVdy, dVdz, Ax, Ay, Az);

			p_Poisson.InitGlobConvFEM(Ax, Ay, Az, -RegPar(n));  // create the stiffness matrix
			p_Poisson.PoissonEqSolverDirichletFEMnew(bdc_q, sdg_FEM, RHS_vec); // solve the Poisson eq.
			p_Poisson.extract_solution(q); // extract the solution to vector q, q now is the UPDATED one!
			p_Poisson.DirichletFEM(bdc_q,q);  //assign the bdc to the computed solution
			out.Write_array_my_part("Test_q_new.m",q,Iter);	
			// === end of step 1.

			//*** STEP 2: update the coefficient
			v = vector_addition(vector_multiplication(q,-ds),V);
			out.Write_array_my_part("Test_v_new.m",v,Iter);	

			//Compute the coefficient using FEM:
			w = Compute_w_from_v_totalwave(s_low, v);
			p_wave_eq.Init3DParameter();
     			p_wave_eq.ParameterFEM2(w, s_low, Coefficient);
			
			//compute_coefficient_FDM_totalwave(Coefficient, v, sdg_FEM, s_low); // using the FDM
			
			coefficient_truncate(Coefficient, C_lb, C_ub, gg, Subdomain); // apply the lower and upper bounds to the coefficient
			//average_coefficient_FDM(Coefficient,sdg_FEM); // coefficient averaging
			out.Write_array_my_part("Test_c_new.m",Coefficient,Iter);				

    			// save the new coefficient in the file of the form C_n_i.inp:
			i=sprintf(C_filename,"%s%d%s%d%s","Coef_new_",n+1,"_",Iter,".inp"); 
			opt.writeInp3D(C_filename, &gg, Coefficient,1); // save the coefficient
			// ===== end of step 2


			// **** STEP 3: update the tail function: 
			// initialize the forward solver: 
			WavesScalarEqOpOpt p_wave_eq2(gg, outer_gg, sdg_FDM, nsd, EXCHANGE, USE_FEM, USE_FDM, 
			     				USE_RHS, USE_DIR_FEM, USE_DIR_FDM, 
			     				USE_ABSORB, PRINT_FILES, NoTimeSteps, maxTime, rhs, 
			     				type_of_mat, velocity, guess_velocity, doIncludeCorners);
		
			dt = p_wave_eq2.InitTime(); // important!!!
			p_wave_eq2.InitFDM();
			p_wave_eq2.Init3DFEM(Coefficient);
			p_wave_eq2.InitExchangeReadMask(); // read exchangeMask.m file.
	
			w = 0.0; t = 0.0; w2 = 0.0; wmax = 0.0;

			//=====+++ MAIN LOOP FORWARD SOLVER +++=====
			for (k = 0; k < NoTimeSteps; k++)
			{   	p_wave_eq2.PlaneWaveBackFDM(k,omega);
				p_wave_eq2.WaveEqSolverFEMforDifMat(t, k, type_of_mat, velocity);
				p_wave_eq2.ApplyExchangeCommon();
				p_wave_eq2.ApplySwap();
				p_wave_eq2.SaveSolutionsFEM(k,E_FEM); // extract the solution

				t += dt;  
				double value1 = exp(-s_up*t)*dt;
				double value2 = exp(-s_low*t)*dt;
				double value3 = exp(-s_max*t)*dt; 
		
				for (i = 0; i < NrFEMNodes; i++)		
				{	w(i) += E_FEM(i)*value1;
					w2(i)+= E_FEM(i)*value2; // this is for the next frequency interval
					wmax(i) +=E_FEM(i)*value3;
				}
			}

			out.Write_array_my_part("Test_tai_new.m",V,Iter);	// save the previous tail
			V = Compute_tail_from_w_totalwave(s_up, w);
			// ----end of step 3---

			// --- stopping rules:::::::::::::::::
			Vec_Error = vector_subtraction(Coefficient,Coeff_old);
			L2Error(Iter) = L2_norm(gg,Vec_Error); 
			RelativeError(Iter) = L2Error(Iter)/(L2_norm(gg,Coeff_old) + 1e-15);	
			MaxCoeff(Iter) = Maxim(Coefficient); 
			L2Coeff(Iter) = L2_norm(gg,Coefficient);

			cout << " Iteration " << Iter << ", Max Coeff = " << MaxCoeff(Iter) << ", L2-error =  " << L2Error(Iter) << endl;
/*			cout <<  endl << " Relative L2 error =  " << RelativeError(Iter) 
			     << ", L2-error =  " << L2Error(Iter) << endl;
			cout << "L2 norm error tail = " << L2_norm(gg, vector_subtraction(V,V_old)) << endl; 
			cout << " Relative L2 norm error tail = " << L2_norm(gg, vector_subtraction(V,V_old))/L2_norm(gg,V_old)<< endl; 
*/
			if ((Iter > 0) && (RelativeError(Iter) >= RelativeError(Iter-1)))
				STOP = true; 
			if (Iter >= MaxIter(n)-1 || L2Error(Iter) <=1e-5)
				STOP = true;
			Iter++; 

		} // end of while loop

		// compute the new tail for the new frequency interval:
		V = Compute_tail_from_w_totalwave(s_low, w2);

		out.Write_array_my_part("Test_L2Error_new.m",L2Error,n);	
		out.Write_array_my_part("Test_RelativeError_new.m",RelativeError,n);	
		out.Write_array_my_part("Test_MaxCoeff_new.m",MaxCoeff,n);	
		out.Write_array_my_part("Test_L2Coeff_new.m",L2Coeff,n);	

		cout << endl << "after iteration: " << n << endl;
		cout << "L2 norm error of tail: " << 	L2_norm(gg, vector_subtraction(V,V_n)) << endl; 	
		cout << "L2 norm error of coeff: " << 	L2_norm(gg, vector_subtraction(Coefficient,Coeff_n)) << endl; 	
		Error(n) = fabs(Maxim(Coefficient) - Maxim(Coeff_n));
		cout << " max error of coeff: " << Error(n) << endl; 	
		cout << "Relative max error of coeff: " << fabs(Maxim(Coefficient) - Maxim(Coeff_n))/Maxim(Coeff_n) << endl; 	


	} // end of for n = ... loop
	V = Compute_tail_from_w_totalwave(s_max, wmax); //if we repeat the whole process, recompute the tail at the maximum frequency
	out.Write_array_my_part("Test_error_new.m",Error,ii);	
} // end of repeating the whole process

// ================================================== 
// ------- end of MAIN LOOP
//===================================================

	cout << endl << " +++++++++++++++++++++++++++++++++" << endl 
	<< " Finish the test " << endl 
     	<< " +++++++++++++++++++++++++++++++++++ " <<endl; 

	// -----------finish calculation:
 	PetscFinalize();
	return 0;
}
