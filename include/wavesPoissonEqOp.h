// Copyright (C) 2000 Larisa Beilina
//
// This file is part of WavES project.
//
// WavES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with WavES. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2000-01-01 by Larisa Beilina
// Last changed: 2012-05-29 by Larisa Beilina



#ifndef _WAVESPOISSONEQOP_
#define _WAVESPOISSONEQOP_

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include "wavesSDDefs.h"
#include "waveskraftwerk.hh"
#include "wavesfindBoundaryNodes.h"
#include "wavesutil2.h"
#include "wavesFuncDefs.h"
#include "wavesBCOperators.h"
#include "wavesOptional.h"
#include "wavesOutputs.h"
#include <petscksp.h>

typedef double real; 
typedef MV_ColMat<real> Mat_real;
//enum overlap_code{over1,over2,over3,over4,over5};

class WavesPoissonEqOp
{

protected:


   WavesOptional opt;
  
  Grid gg;
 
 bool USE_FEM;

 bool USE_RHS;

 bool USE_DIRICHLET_FEM;

 bool PRINT_FILES;
 int  ierr;
 double EVALPT_X; 
 double EVALPT_Y;
 double EVALPT_Z;



 double rhs;
 int noQpts;

 int nbn;
 int* boundindex;
 // Scalar* bvalues;
 double* bound_values;
 
 int nbn1;
 int* boundindex1;
 // Scalar* bvalues1;
 double* bound_values1;
 
int nbn2;
 int* boundindex2;
 // Scalar* bvalues2;
 double* bound_values2;
 
 double* bvalues;
 


  ofstream ofs; 
  int  nsd,n_i,n_j,n_k;
  bool doIncludeCorners;
  
 public: 

  WavesPoissonEqOp(Grid &gg, 
			 bool& USE_FEM,
			 bool &USE_RHS,
			 bool &USE_DIRICHLET_FEM,
			 bool &PRINT_FILES,
			 double &EVALPT_X, 
			 double &EVALPT_Y,
			 double &EVALPT_Z,
			 double &rhs, int& noQpts_) ;
  

 enum overlap_code{over1,over2,over3,over4,over5};
 typedef real (*Fcn3d) (real x,real y,real z);
 MV_Vector<int> u1nodes;
 MV_Vector<int> u2nodes;
 


// initialization for FEM solution 
  Mat A, B, Adual;
  Vec grad_;

  Vec lumpedAmass;
  Vec Fn;
 
  Vec u01;
  Vec u02;
  Vec u1; 
  Vec u2;
  Vec puas_sol;
  Vec exsol;
  int FEM_initialised;
  
  int Poisson_FEM_initialized;
  int bvalues_initialised;

  /// destructor 
  ~WavesPoissonEqOp();

  //void  Configurate(char *fname);

  void constructGrad_x(Mat&  Gradient_x);
  void constructGrad_y(Mat&  Gradient_y);
  void constructGrad_z(Mat&  Gradient_z);

  void GradFEM(Mat& Gradient, MV_Vector<double>& q, Vec& u_2 );
  void Convert( Vec& u_2, double* solfdm);

  void constructStiffnessMatrixOpt( Mat& A);
  
  void constructStiffnessMatrixOpt(Mat& A, 
				   MV_Vector<double>& grad_q_x,  
				   MV_Vector<double>& grad_q_y, 
				   MV_Vector<double>& grad_q_z, double eps_);

  void getTriangleQuadrature(const int noQpts, 
			     MV_ColMat<real>& points, 
			     MV_Vector<real>& weights);
  
  void getTetraQuadrature2(const int noQpts,
			   MV_ColMat < real > &points,
			   MV_Vector < real > &weights);
  
  void InitPoissonFEM();
  
  void InitGlobConvFEM(MV_Vector<double>& grad_q_x,
		       MV_Vector<double>& grad_q_y, 
		       MV_Vector<double>& grad_q_z, double eps);
  
  
  void ApplyDirichletFEM(Vec& exact);
  void ApplyDirichletFEM(MV_Vector<double>& exact);  
  void DirichletFEM(MV_Vector<double>& exsol, real* function);

  int dirichletBCxxx(Mat&  Matrice );

  void PoissonEqInitBoundaryValues(MV_Vector<double>& exsol, 
				    double* initbvalues);

  void FEM_Boundary(MV_Vector<double>& exsol,
		    MV_Vector<double>&  func_obs);

  void evalRHS(Vec& Fn_, MV_Vector<double>& rhs_);
  
  void PoissonEqSolverFEM(MV_Vector<double>& exsol);

  void PoissonEqSolverDirichletFEMnew(MV_Vector<double>& exsol, 
				      WavesSDGeometry& sdg, 
				      MV_Vector<double>& grad_tail_square);
  
  void writeSol(char* filename);
  void writeSol(char* filename, Vec& sol);
  void extract_solution(MV_Vector<double>& Sol); // added by Thanh
  void DirichletFEM(MV_Vector<double> exsol, MV_Vector<double>& function);
  
  void  Configurate(char *fname);
  //  void ComputeL2NormPoisson();
  

};



#endif

