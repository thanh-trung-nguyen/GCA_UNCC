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
// Changed by Vladimir Timonov: 2012-01-09 
// Changed by Thanh Nguyen: 2014-02-15


#ifndef __WAVESOPTIONAL_H
#define __WAVESOPTIONAL_H

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include "wavesSDDefs.h"
#include "waveskraftwerk.hh"
#include "wavesSDGeometry.h"
#include "wavesfindBoundaryNodes.h"
#include "wavesSDOperator.h"
#include "wavesSDIndexes.h"
#include "wavesSDMaskIndexes.h"
#include "wavesFuncDefs.h"
#include "wavesBCOperators.h"

typedef double real;
typedef MV_ColMat<real> Mat_real;
typedef MV_ColMat<int> Mat_int;

enum overlap_code
{
		over1,
		over2,
		over3,
		over4,
		over5
};

class WavesOptional
{
	protected:

		Grid outer_gg;
		Grid gg;

		bool EXCHANGE;
		bool USE_FEM;
		bool USE_FDM;
		bool USE_RHS;
		bool USE_SYSTEM;
		bool USE_DIRICHLET;
		bool USE_DIRICHLET_FEM;
		bool USE_ABSORB;
		bool PRINT_FILES;
		int ierr, nrSTEPS;
		double EVALPT_X;
		double EVALPT_Y;
		double EVALPT_Z;
		double maxtime;
		double dt, rhs;

		int nbn;
		int* boundindex;
		PetscScalar* bvalues;
		double* bound_values;
		int* maska;
		int *outerIndexArray0;
		int *outerIndexArray1;
		int *outerIndexArray2;
		int* innerIndexArray1;
		int* innerIndexArray2;
		int* innerIndexArray0;

		MV_Vector<int> u0nodes;
		MV_Vector<int> u1nodes;
		MV_Vector<int> u2nodes;

		MV_Vector<int> exchangeMask;

		ofstream ofs;
		int nsd, n_i, n_j, n_k;

	public:

		WavesSDGeometry sdg;

		WavesSDIndexes outerWithHole;
		WavesSDIndexes outerBoundary1;
		WavesSDIndexes outerBoundary2;

		WavesOptional(WavesSDGeometry &sdg, Grid &outer_gg, Grid &gg, bool& EXCHANGE, bool& USE_FEM, bool &USE_FDM, double &EVALPT_X, double &EVALPT_Y, double &EVALPT_Z, double &maxtime);

		WavesOptional(WavesSDGeometry &sdg, Grid &outer_gg, Grid &gg, bool& EXCHANGE, bool& USE_FEM, bool &USE_FDM, double &maxtime);

		WavesOptional(Grid &gg, bool& USE_FEM, double &EVALPT_X, double &EVALPT_Y, double &EVALPT_Z);

		WavesOptional(Grid &gg, bool& USE_FEM, bool& EXCHANGE, bool &USE_FDM);

		WavesOptional(Grid &gg);
		WavesOptional();
		enum overlap_code
		{
				over1,
				over2,
				over3,
				over4,
				over5
		};

		void extractNodesFDM(WavesSDGeometry& ref, WavesSDGeometry& coarse, MV_Vector<int>& NodesFDM);

		void copyValues(double *x, real *y, const int *IX, const int *IY, const int &n);

		void sort(int *ia, const int &n);

		void Write_To_File(char *file_name, MV_Vector<double>& Ar1);

		void InitExchangeFEM(Grid& grid, Grid& grid_outer);

		void InitExchangeConvergFEM(Grid& grid);
		void InitExchangeFEMfromFile(Grid& grid, Grid& grid_outer);

		void initExchangeFEMtoFDM(WavesSDIndexes& innerHole);
		void initExchangeFDM();

		void initExchangeStructFDM();
		void initExchangeMaskFDM();
		void initExchangeFDM(WavesSDIndexes& innerHole);
		void initExchangeFDMtoFEM();

		void ApplyExchangeFDMtoFEM(Vec& u2, real* uOuter);
		void ApplyExchangeFEMtoFDM(Vec& u2, real* uOuter);

		void ApplyExchange(Vec& u2, real* uOuter);
		void ApplyExchangeConverg(Vec& u2, real* uOuter, WavesSDIndexes& innerHole);

		int bvalues_initialised;

		/// destructor 
		~WavesOptional();

		void LaplaceTransform(WavesGridB& gg, Mat_real& array, double omega, double dt, MV_Vector<real>& LaplaceTr);

		void LaplaceTransform_(WavesGridB& gg, Mat_real& array, double omega, double dt, MV_Vector<real>& LaplaceTr);

		void LaplaceTransform(Mat_real time_array, double dt, double s, MV_Vector<real>& LaplaceTr); //added by Thanh.
		void LaplaceTransform_col(Mat_real time_array, double dt, double s, MV_Vector<real>& LaplaceTr); //added by Thanh.
		void LaplaceTransform(char* datafile, int NoTimeSteps, double dt, double s, MV_Vector<real>& LaplaceTr); //added by Thanh
		void LaplaceTransform(char* datafile, int NoTimeSteps, double dt, MV_Vector<double> s, Mat_real& LaplaceTr);// added by Thanh

		double Compute_space_int(WavesGridB& gg, MV_Vector<real>& difference, int rank);

		double Compute_space_int(MV_Vector<real>& Values_at_Elms);

		real computeMinElmSize(Grid& grid);
		real computeMaxElmSize(Grid& grid);

		void BoundaryNodes(MV_Vector<int>& bndNodes);

		int *maskaFDM(MV_Vector<int>& bndNodes);
		int *maskaStructFDM(MV_Vector<int>& bndNodes);

		real *fromFEMtoFDM(Vec& u2);

		real sortMakeIndex(MV_Vector<real>& values, MV_Vector<int>& index, int sort_range);

		void findElmsCommonNodes(WavesGridB& grid, Mat_int& elmNeighbor, Mat_int& El2Nodes, Mat_int& El2NotCommonNode);

		void findElmsCommonNodes3D(WavesGridB& grid, Mat_int& El2Nodes, Mat_int& El2NotCommonNode);

		void findNeighborElms(WavesGridB& grid, Mat_int& elmNeighbor);
		void findNeighborElms3D(WavesGridB& grid, Mat_int& elmNeighbor);

		void extractNodeNumbers(MV_Vector<int>& exchangeMask, MV_Vector<int>& unodes, int num);

		void computeSortVectorfromGrid(MV_Vector<int>& unodes, MV_Vector<real>& sortkey);

		void makeSortedNodesIndex(MV_Vector<int>& exchangeMask, MV_Vector<int>& unodes, int num);

		void getElmsAtFaces(const WavesNeighborFE& n2e, const int n1, const int n2, const int n3, const int n4, const int e, int& e3, int& e4);

		void computeExNodesOuterbord(Grid& grid, MV_Vector<int>& bndNodes, MV_Vector<int>& markNodes, int &outer_bord_code);

		void computeExNodesCommon(Grid& grid, MV_Vector<int>& bndNodes, MV_Vector<int>& markNodes, const int &bord_code1, const int &bord_code2);

		void extractBoundNodes(Grid& grid_outer, Grid& grid, MV_Vector<int>& bndNodes);

		void makeExchangeNodes(Grid& grid, const overlap_code &ovc, MV_Vector<int>& bndNodes);

		void myOpen(ofstream &ofs_);
		void MTVOpen(ofstream &ofs_);
		void writeMatlabDataInPnt(MV_Vector<real>& solFEM, MV_Vector<real>& solFDM, MV_Vector<real>& timedata);

		void writeSolInMatlabFile(Vec& u1, Vec& u2);

		int getElmAtFace(const WavesNeighborFE& n2e, const int n1, const int n2, const int n3, const int e);

		int getFaceNo(const WavesGridB& grid, const int n1, const int n2, const int n3, const int e, bool safe);

		void findBoundaryNodes(WavesGridB& grid, MV_Vector<int>& boundaryMarkers);

		void findB_Elms_Nodes(WavesGridB& grid, MV_Vector<int>& boundaryElms, MV_Vector<int>& boundaryMarkers);

		void copyVal(const real *x, real *y, const int *IX, const int *IY, const int &n);

		void sortFDM(int *ia, const int &n);
		void getFileNameFDM(char* buf);
		void getFileName(char* buf);
		void getFileNameElastic1(char* buf);
		void getFileNameElastic2(char* buf);

		int findInitDisturbPoint(MV_Vector<real>& evalpt);

		void print_gid_mesh(char *file);
		void print_gid_mesh_FEM(WavesGridB& gg, char *file);
		void print_gid_mesh_Amira(WavesGridB& gg, char *file);

		void print_gid_jump(WavesGridB& gg, char *file, MV_Vector<int>& jump);

		void print_gid_elements_FEM(WavesGridB& gg, char *file);
		void print_gid_nodes_FEM(WavesGridB& gg, char *file);

		int print_common_2DFEM(char *file, real* u_array, real* u_1_, real* u_2_, int k);

		int print_common_3DFEM(char* file, real* u_array, real* u_1_, real* u_2_, real* u_3_, int k);

		void getMTVName(char* buf);
		int writeInp(char *file, Grid *grid, Vec *u, int nsys);
		int writeInpElType(char *file, Grid *grid, Vec *u, int nsys, int type_of_mat, int NODE_DATA);

		int writeInpElType2D(char *file, Grid *grid, Vec *u, int nsys, int type_of_mat, int NODE_DATA);

		int writeInp2D(char *file, Grid *grid, Vec *u, int nsys, int p3d);

		int writeInp2D(char *file, Grid *grid, MV_Vector<double>& b, int nsys, int p3d);
		int writeInp3D(char *file, Grid *grid, MV_Vector<double>& b, int nsys);
		int writeInpFDM2D(char *file, WavesSDGeometry *grid, double *u, int nsys, int p3d);

		void ApplyAvsOut(char* filename, Vec& sol);
		void ApplyAvsOutElType(int sch, Vec& sol, int EL_TYPE);

		void ApplyAvsOut(int k, Vec& u2);
		void ApplyAvsOutElastic(int k, Vec& u21, Vec& u22, Vec& ucommon, double* uOuter);

		int writeMTV2D(char *file, Grid *grid, Vec *u, int nsys);
		int writeMTV4D(char *file, Grid *grid, Vec *u, int nsys);
		void PlotMTVFDMGrid();
		void ApplyPlotMTVOut(int k, real* vOuter, Vec& u2);
		bool Inp_to_RES_Gid(Vec& u2, char *inp_file, int results, int k);

		double InitTime();

		int nrOuternod()
		{
			return sdg.getNoNodes();
		}
		int nrFEMnod()
		{
			return gg.getNoNodes();
		}

		void printAVS(char* filename);
		void print_GID_MESH(char* filename);

		int print_2D_onepoint(char *file, real& u_array, real& u_1_, real& u_2_, int k);

		int print_onepoint(char *file, real& u_array, real& u_1_, real& u_2_, int k);

		int print_3D_onepoint(char *file, real& u_array, real& u_1_, real& u_2_, real& u_3_, int k);

		int print_1Dgid_result(char *file, MV_Vector<double>& results, int k);

		int print_params(char *file, MV_Vector<double>& results, int k);

		int print_params_int(char *file, MV_Vector<int>& results, int k);

		int print_2Dgid_result(char *file, real* u_1_, real* u_2_, real* u_array, int k);

		int print_2Dgid_result(char *file, MV_Vector<real>& u1, MV_Vector<real>& u2, MV_Vector<real>& u3, int k);

		int print_2Dgid_result_common(char *file, real* u_1_, real* u_2_, real* u_array, int k);

		int print_3Dgid_result(char *file, real* u_array, real* u_1_, real* u_2_, real* u_3_, int k);

		int print_3Dgid_result(char *file, MV_Vector<double>& u_array, MV_Vector<double>& u_1_, MV_Vector<double>& u_2_, MV_Vector<double>& u_3_, int k);

		int print_3Dgid_result_common( char *file, 
							      real*  u_array,
							      real* u_1_,
							      real* u_2_,
							      real* u_3_,
							      int k);

		int print_3Dgid_LENS(char *file, real* u_1_, real* u_2_, real* u_3_, int k);

		void print_GID(char* filename);

		void ComputeL2Norm(double* func_to_norm, double& l2norm);

		double L2Norma(Grid& gg, MV_Vector<double>& Gradient);
		double Norma(Grid& gg, MV_Vector<double>& Gradient);

		double GetLength2D(int n1, int n2);
		double GetArea3D(int n1, int n2, int n3);

		void GetNormal2D(int n1, int n2, MV_Vector<double>& normal);
		void GetNormal3D(int A, int B, int C, MV_Vector<double>& normal);

		void NormalDir2D(int n1, int n2, int n3, MV_Vector<double>& normal, int& normal_minus, int& normal_plus);
		void NormalDir3D(int n1, int n2, int n3, int n4, MV_Vector<double>& normal, int& normal_minus, int& normal_plus);

		void makeNormalArray2D(Mat_real& NormalArray2D);
		void makeNormalArray3D(Mat_real& NormalArray3D);
};

#endif

