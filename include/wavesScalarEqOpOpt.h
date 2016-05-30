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
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESSCALAREQOPOPT_H
#define __WAVESSCALAREQOPOPT_H

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include "wavesSDDefs.h"
#include "waveskraftwerk.hh"
#include "wavesfindBoundaryNodes.h"
#include "wavesutil2.h"
#include "wavesBCOperators.h"
#include "wavesOptional.h"
#include "wavesOutputs.h"
#include "wavesSDIndexes.h"
#include "wavesWaveEquationOp.h"
#include "wavesWaveEqOpSeismic.h"
#include "wavesPlaneWaveOpexpl.h"
#include "wavesFuncDefs.h"
#include "wavesReflectedWave.h"
#include "wavesWaveEqOpPar.h"

typedef double real;
typedef MV_ColMat<real> Mat_real;

class WavesScalarEqOpOpt
{

	protected:

		WavesOptional opt;
		//==========================================00

		Grid gg;

		Grid outer_gg;
		WavesSDGeometry sdg;
		bool EXCHANGE;

		bool USE_FEM;
		bool USE_FDM;
		bool USE_RHS;
		bool USE_DIRICHLET_FEM;
		bool USE_DIRICHLET_FDM;
		bool PRINT_FILES;
		bool USE_ABSORB;

		WavesSDBoundary leftBoundary;
		WavesSDBoundary rightBoundary;
		WavesSDBoundary lowBoundary;
		WavesSDBoundary topBoundary;
		WavesSDBoundary perBoundary;
		WavesSDBoundary backBoundary;
		WavesSDBoundary corners;
		DirichletOp dirichletBC;
		DirichletOp cornerBC;

		WavesWaveEqBoundary abs1;
		WavesWaveEqBoundary abs2;
		WavesWaveEqBoundary abs3;
		WavesWaveEqBoundary abs4;
		WavesWaveEqBoundary abs5;
		WavesWaveEqBoundary abs6;

		WavesSDIndexes outerWithHole;
		WavesSDIndexes outerBoundary1;
		WavesSDIndexes outerBoundary2;

		SDBoundaryType sdb;

		//for output 2D only for FDM
		AVSOutputOp avs2dFDM;

		AssignTimeFunctionOp initialize2DTimeFunction;
		AssignTimeFunctionOp initialize3DTimeFunction;
		AssignTimeFunctionOp initialize3D_1;
		AssignTimeFunctionOp initrhs1;
		AssignTimeFunctionOp initrhs2;

		// from WavesSDOperator.h

		ApplyFunctionOp add_sol_and_function;

		int nrSTEPS, type_of_material;
		int ierr;

		double maxtime, rhs, velocity, guess_velocity;

		int nbn;
		int* boundindex;
		PetscScalar* bvalues;
		double* bound_values;

		int* innerIndexArray1;
		int* innerIndexArray2;
		int *outerIndexArray1;
		int *outerIndexArray2;

		ofstream ofs;
		int nsd;

		int n_i, n_j, n_k;


	public:

		WavesScalarEqOpOpt(Grid &gg_, Grid &outer_gg_, WavesSDGeometry &sdg_, int &nsd, bool& EXCHANGE_, bool& USE_FEM_, bool& USE_FDM_, bool &USE_RHS_, bool &USE_DIRICHLET_FEM_, bool &USE_DIRICHLET_FDM_, bool &USE_ABSORB_, bool &PRINT_FILES_, int &nrSTEPS_,
				double &maxtime_, double &rhs_, int &type_of_material_, double &velocity_, double &guess_velocity_, bool &doIncludeCorners_);

		enum overlap_code
		{
				over1,
				over2,
				over3,
				over4,
				over5
		};
		typedef real (*Fcn3d)(real x, real y, real z);
		MV_Vector<int> u1nodes;
		MV_Vector<int> u2nodes;

// initialization for FEM solution 
		Mat A, Adual;
		Vec C;
		Mat AA;

		Vec lumpedAmass;
		Vec Fn;

		Vec uhelp;
		Vec uhelpC;
		Vec u01;
		Vec u02;
		Vec u1;
		Vec u2;
                Vec k2;

		Vec u1_x;
		Vec u1_y;
		Vec u1_z;
		Vec u2_x;
		Vec u2_y;
		Vec u2_z;

		double dt;
		int FEM_initialised;

		int bvalues_initialised;

		bool doIncludeCorners;

		//==========================  
		// for FDM initialization
		//==========================

		real* u_Outer;
		real* v_Outer;

		real* b_Outer;

		real* u_x_Outer;
		real* v_x_Outer;

		real* u_y_Outer;
		real* v_y_Outer;

		real* u_z_Outer;
		real* v_z_Outer;

		real* tmpOuter;
		real* fn_Outer;

		real* u_exact;

		double InitTime();

		void print_GID_MESH(char* filename);
		void print_GID_MESH_FEM(WavesGridB& grid, char* filename);
		void print_GID_MESH_FDM(char* filename);

		void print_mesh_common(char* filename, WavesGridB& grid);
		void print_mesh_common(char* filename);

		void ApplySwap();
		void ApplySwap_LENS();

		void InitFDM();

		int print_GID_VEC_FDM(char* filename, int sch);
		int print_GID_VEC_FEM(char* filename, int sch);
		int print_GID_VEC_COMMON(char* filename, int sch);

		int print_GID_VEC_FEM_LENS(char* filename, int sch);

		int printCommon(char* filename, int sch);

		void ApplyRHS(double& t);
		void ApplyAdjRHS();

		void AdjScalarFDM(int k, real* adj_f);
		void AdjScalarFEM(double t, int k, int type_of_material, MV_Vector<double>& Dif_sol, MV_Vector<double>& velocity);

		void AdjFEM_LENS(double t, int k, MV_Vector<double>& Dif_sol, MV_Vector<double>& velocity);

		int nno;

/// destructor 
		~WavesScalarEqOpOpt();

		void PlaneWaveFEMforDifMat(double t, int k, int type_of_material, double velocity, double y_fix);

		void findNodesLENS(MV_Vector<int>& markNodes);
		void findNodesEarth(MV_Vector<int>& markNodes);
		void findNodesEarthFDM(MV_Vector<int>& markNodes, int code);
		void findNodesPhotonFDM(MV_Vector<int>& markNodes, int code);

		void findNodesPhotonFDMY(MV_Vector<int>& markNodes, double fix, double fix_min, double fix_max, int code);

		void findNodesPhotonFDMX(MV_Vector<int>& markNodes, double fix, double fix_min, double fix_max, int code);

		void PlaneWaveFEM_LENS(double t, int k, MV_Vector<int>& markNodes);

		void PlaneWaveFEM_Earth(double t, int k, MV_Vector<int>& markNodes);

		void InitExchangeCommon();
		void InitExchangeCommon(WavesScalarEqOpOpt& p);
		void InitExchangeStructCommon();

		void InitExchangeReadMask();

		void ExFDMtoFEM();

		void ApplyExchangeCommon();
		void ApplyExchange();
		void ApplyExchange1(WavesScalarEqOpOpt& p);

		void ApplyExchange_LENS();

		int myVecSwap(Vec u1_, Vec u2_);

		int dropZeroWavesElements(Mat& Ain, Mat& Aout, double droptol);

		void writeSol();

		int MatMult0(Mat AA, Vec xx, Vec yy);

		void constructIterationMatrixOpt(Mat& A_, Vec& lumpedAmass_);

		void getTetraQuadrature(const int noQpts, MV_ColMat<real>& points, MV_Vector<real>& weights);

		void getTriangleQuadrature(const int noQpts, MV_ColMat<real>& points, MV_Vector<real>& weights);

		void evalRHS(Vec& Fn_, Fcn3d rhs_);

		void evalRHSforDifMat(Vec& Fn_, double t, Fcn3dtime rhs_, Vec& lumpedAmass_, int type_of_material, double velocity);

		void evalRHSforDifMat(Vec& Fn_, double t, Fcn3dtime rhs_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity);

		void evalAdjRHSforDifMat(Vec& Fn_, double t, MV_Vector<double>& Dif_sol, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity);

		void applyDirichletCondition(Vec& u, int nbn, int* boundindex, double* bvalues);

		void copyVal(const real *x, real *y, const int *IX, const int *IY, const int &n);

		void InitDirichletFEM();

		void InitDirichletFEM(Grid& grid, Grid& grid_outer);

		void InitDirichletFEM(Mat_real& BoundaryCoord);

		void InitFEMDifMat(int type_of_material, MV_Vector<double>& velocity);

		void Init3DFEM();
		void Init3DFEM(MV_Vector<double>& b);
		void Init2DFEM(MV_Vector<double>& b);
		void Init2DGLK(MV_Vector<double>& b);
		void Init3DFEM_LENS();
		void Init3DFEMDifMat_LENS(MV_Vector<double>& new_velocity);

		void Init2DParameter();
		void Init3DParameter();

	

		void InitPlaneWave();
		void InitPlaneWave1();
		void InitPlaneWave_GLK();
		void InitPlaneWave_LENS();

		void PlaneWaveLeftFDM(int k);
		void PlaneWaveLeftFDM(int k, double omega);
		void PlaneWaveBackFDM(int k);
		void PlaneWaveBackFDM(int k, double omega);
		void PlaneWaveBackFDM(int k, double omega, double& param);

		void PlaneWavePerFDM(int k);
		void PlaneWaveTopFDM(int k, double omega);
		void PlaneWaveTopFDM(int k, double omega, double& param);
		void PlaneWaveTop1FDM(int k, double omega, double& param);

		void PlaneWaveBotFDM(int k, double omega);
		void PlaneWaveBotFDM(int k, double omega, double& param);

		void apply_reflwave_low(int& k, double omega, MV_Vector<double>& array);
		void PlaneWaveTopFDMSeismic(int k);

		void AdjPlaneWaveTopFDM(int k, real *Difsol);

		void AdjPlaneWaveFDM(int k, real *Difsol, WavesSDBoundary& inner);

		void AdjPlaneWaveLeftFDM(int k, real *Difsol);
		void AdjPlaneWaveFDMSeismic(int k, real *Difsol);
		void ReflAdjPlaneWaveFDM(int k, real *Difsol);

		void PlaneWaveTopFDM_(int k);

		void applyPlaneWaveFEM(MV_Vector<int>& markNodes, double t, double y_fix);

		void WaveEqSolverFEMforDifMat(double t, int k, int type_of_material, double velocity);

		void WaveEqGLK(double t, int k, MV_Vector<double>& bb);

		void WaveEqSolverFEMforDifMat(double t, int k, int type_of_material, MV_Vector<double>& velocity);

		void IterationMatrix2D(Mat& A_, Vec& lumpedAmass_, int type_of_material);

		void IterationMatrix2DDifMat(Mat& A_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& bb);

		void IterationMatrix2DGLK(Mat& A_, Vec& C_, Vec& lumpedAmass_, MV_Vector<double>& bb);

		void IterationMatrix2DGLKnew(Mat& A_, Vec& C_, Vec& lumpedAmass_, MV_Vector<double>& bb);

		void IterationMatrix3D(Mat& A_, Vec& lumpedAmass_, int type_of_material);

		void IterationMatrix3Dnew(Mat& A_, Vec& lumpedAmass_, int type_of_material);

		void IterationMatrix3DDifMatnew(Mat& A_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity);

		void IterationMatrix3DDifMat(Mat& A_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity);

		void  ParameterMatrix2D( Mat& A_, 
					 Vec& lumpedAmass_,
					 int type_of_material);
		
		void ParameterMatrix3D( Mat& A_, 
					Vec& lumpedAmass_,
					int type_of_material);

		void ParameterFEM(MV_Vector<double>& V_tilde, double& omega,  
				  MV_Vector<double>& param);
		
		void ParameterFEM2(MV_Vector<double> w, double omega,  
				      MV_Vector<double>& param); // added by Thanh, calculate the coefficient given w


		void PlaneFDMtoFEM(int k);
		void PlaneFEM(double t, int k, int type_of_material, double velocity);
		
		void getFileNameElastic1(char* buf);
		void getFileNameElastic2(char* buf);

		void InitPuassonFEM();
		void PuassonEqSolverFEM();

		void Save_Sol_FEM(int nr_of_it, MV_Vector<double>& Save_u_fem, MV_Vector<int>& Vec_Glob_Num);

		void Save_Sol_FEM(int nr_of_it, Mat_real& Save_u_fem);

		void Save_Sol_FEM(MV_Vector<double>& Save_u_fem);

		void Save_Sol_FDM(int nr_of_it, MV_Vector<double>& Save_u_fdm, MV_Vector<int>& Vec_Glob_Num);

		void Save_Sol_FDM(int nr_of_it, MV_Vector<double>& Save_u_fdm);

		void SaveSolutionsFEM(int nr_of_it, MV_Vector<double>& Save_fem);

		void ScalarFDM(int k);
		void ScalarFDM_LENS(int k);

		void ApplyLaplasOp(real* y);

		void ApplyAvsOut(char* filename);
		void ApplyAvsOut(int k);

		int nrFEMnod()
		{
			return gg.getNoNodes();
		}
		void printAVS(int& k);

		void print_GID(char* filename, int sch);
		void ComputeL2NormPoisson();
		int print_AVS_FEM(int sch);
		int print_AVS_FDM(int sch);

		int print_AVS_ElType(int sch, int EL_TYPE);

};

#endif

