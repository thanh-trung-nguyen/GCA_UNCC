//KRAFTWERK

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
// Second changed: 2012-01-09 by Vladimir Timonov
// Last modified by Thanh Nguyen, 2014-01-15.

#ifndef __WAVESOUTPUTS_H
#define __WAVESOUTPUTS_H

#include <stdio.h>
#include <errno.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <sys/times.h>
#include <iostream>
#include <string.h>
#include <fstream>

#include <sys/types.h>
#include <stdio.h>

#include "iomanip"
#include "assert.h"

#include "waveskraftwerk.hh"
#include "wavesSDGeometry.h"
#include "wavesSDOperator.h"
#include "wavesFEcomp.h"
#define WIDTH  16

typedef double real;
typedef MV_ColMat<real> Mat_real;
typedef MV_ColMat<int> Mat_int;

class WavesOutputs
{

	public:

		Grid gg;
		WavesSDGeometry sdg;
		MV_Vector<double> bb;
		MV_Vector<int> b;
		double norma;
		int it;

		WavesOutputs();
		~WavesOutputs();
		WavesOutputs(Grid& gg_, MV_Vector<double>& bb_);
		WavesOutputs(Grid& gg_, MV_Vector<int>& b_);
		WavesOutputs(Grid& gg_);
		WavesOutputs(WavesSDGeometry& sdg_);
		WavesOutputs(double norma_, int it_);

		WavesOutputs(double x_min, double x_max, double y_min, double y_max, double z_min, double z_max);

		void CreateVecGlobNum(char *fname, MV_Vector<int>& Vec_Glob_Num);

		void Write_const(char *fname, double& alfa);
		void Write_const(char *fname, int alfa);

		void WriteMyPart(char* file, double func, int iter);
		void WriteMyPart(char* file, int nr, double x, double y, double z, int iter, double norma);

		int writeInp(char *file, Grid *grid, double* u, int nsys);

		void ReadConst(char *fname, double& alfa);
		void ReadSolfromFile(char *fname, Mat_real& Ar1);
		void ReadSolfromFile(char *fname, MV_Vector<double>& Ar1);
		void ReadSolfromFile(char *fname, MV_Vector<int>& Ar1);

		void WriteToFile(char *fname, real* Ar1, int nno);
		void WriteToFile(char *fname, MV_Vector<double>& Ar1);
		void WriteToFile(char *fname, MV_Vector<bool>& Ar1);
		void WriteToFile(char *fname, MV_Vector<int>& Ar1);
		void WriteSolutionsToFile(char *fname, Mat_real& Ar1);

		int print_GID_FDM(char* filename, int sch, MV_Vector<double>& E1_new, MV_Vector<double>& E2_new);

		int print_GID_FDM3D(char* filename, int sch, MV_Vector<double>& E1_new, MV_Vector<double>& E2_new, MV_Vector<double>& E3_new);

		int print_GID_FDM(char* filename, int sch, Mat_real& E);

		void WriteToGidFile1(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num);
		void WriteToGidFile(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num, int code = 0);

		void WriteWavesElementArrays(WavesGridB& gg, Mat_int& Gid_Elms, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int& sch);

		void WriteToGidFile2(char *fname, WavesGridB& gg, MV_Vector<int>& Vec_glob_Num);

		void WriteToGidFile(char *fname, WavesSDGeometry& sdg, MV_Vector<int>& Vec_glob_Num);

		void GidWrite(char *fname, WavesSDGeometry& sdg, MV_Vector<int>& Vec_glob_Num);

		bool Inp_to_RES_Gid(char *inp_file, int results, int k);

		void WriteToFile(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code);
		void WriteToFile2(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code);
		void WriteToFile002(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code);
		void WriteToFile0025(char *fname, WavesSDGeometry& sdg, WavesGridB& gg, int code);

		double L2Norma();
		double L2Norma(MV_Vector<int>& glob_nodes, MV_Vector<double>& values_in_glob_nodes);
		int offset(int j);
		int offset(int n_i, int n_j, int iter);
		void Write_my_part(char* file);
		void Write_my_part(char* file, double& value, int iter);
		void Write_array_my_part(char* file, Mat_real& array_part, int iter);
		void Write_vector_my_part(char* file, Mat_real& array_part, int iter);
		void Write_array_my_part(char* file, MV_Vector<double>& array_part, int iter);

		// +++++new functions added by Thanh	
		void Write_array_my_part(char* file, MV_Vector<int> array_part); // added by Thanh
		void Write_matrix_to_file(char* file, Mat_real& array_part); // write a matrix to a file. Thanh added
		void load_boundary_data(char *fname, Mat_real& Ar1, MV_Vector<int>& BoundNodes);
		void Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 	WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
			 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
			 bool &USE_DIR_FDM, bool& USE_ABSORB, bool &USE_RHS,
			 bool &PRINT_FILES, int &NoTimeSteps, double &maxTime, double &rhs_code, 
			 int& NoObjects, int &TypeOfMat1, int &TypeOfMat2, int &TypeOfMat3, 
			 double& velocity1, double& velocity2, double& velocity3,
			 double& noiselevel, double& freq, MV_Vector<double>& Zposition); //added by Thanh

		void Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 		WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
			 	bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
		 		bool &USE_DIR_FDM);


		void Configurate(char *fname, Grid& grid, Grid& outer_gg, 
		 	WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
			 bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
			 bool &USE_DIR_FDM, bool& USE_ABSORB,
			 int &NoTimeSteps, double &maxTime, double& freq); //added by Thanh

		string Configurate(char *fname); // added by Thanh
		string Configurate(char *fname, int idx);

		//load the type of material to the grid:
		void VecGetTypeofMat(WavesGridB& gg, int nsd, int type_of_material, double velocity, MV_Vector<double>& b);
		void VecGetTypeofMat(WavesGridB& gg, int nsd,
			     		int type_of_material1, int type_of_material2,
		     			double velocity1, double velocity2, 
		 	    		MV_Vector<double>& b);
		void VecGetTypeofMat(WavesGridB& gg, int nsd,
		   			  int type_of_material1, int type_of_material2, int type_of_material3,
		   			  double velocity1, double velocity2, double velocity3,
		   			  MV_Vector<double>& b);

		int  load_boundary_grid(char *bdgrid);
		void load_boundary_grid(char *bdgrid, MV_Vector<int>& BoundNodes);
		void load_boundary_grid(char *bdgrid, MV_Vector<int>& BoundNodes, MV_Vector<double>& X, MV_Vector<double>& Y, MV_Vector<double>& Z);

		string load_inversion_parameters(char* fname, double& s_min, double& s_max, int& Ns); // load inversion parameters
		void load_inversion_parameters(char* fname, MV_Vector<double>& Lambda_CWF, MV_Vector<double>& RegPar, 
			     			MV_Vector<int>& MaxIter, double& C_lb, double& C_ub);
		void load_inversion_parameters(char* fname, MV_Vector<double>& Lambda_CWF, MV_Vector<double>& RegPar, 
			     			MV_Vector<int>& MaxIter, double& C_lb, double& C_ub,MV_Vector<double>& Subdomain);
		string load_inversion_parameters(char* fname); //load the file name with incident wave
		string load_inversion_parameters_fp(char* fname);
		string load_inversion_parameters_typebndmea(char* fname);
		string load_inversion_parameters_firsttail(char* fname);

		// new input file type:
		void Configure(char *fname, Grid& grid, Grid& outer_gg, 
		 		WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 		bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, 
		 		bool &USE_DIR_FDM, bool& USE_ABSORB, bool &USE_RHS,
		 		bool &PRINT_FILES, int &NoTimeSteps, double &maxTime, double &rhs_code, 
		 		double& noiselevel, double& freq, int& NoObjects, MV_Vector<double>& Zposition);
		void Configure(char *fname, Grid& grid, Grid& outer_gg, WavesSDGeometry& sdg, WavesSDGeometry& sdg_new, int& nsd, 
		 		bool& USE_EXCHANGE, bool& USE_FEM, bool& USE_FDM, bool &USE_DIR_FEM, bool &USE_DIR_FDM, bool& USE_ABSORB, 
				int &NoTimeSteps, double &maxTime, double& freq);

		void Configure(char *fname, double& x_minFEM, double& y_minFEM, double& z_minFEM, 
		  		    double& x_maxFEM, double& y_maxFEM, double& z_maxFEM, double& maxTime, int& NoTimeSteps);

		string boundary_grid_file_name(char *fname); // load the boundary grid file name
		string FEM_grid_file_name(char *fname); //load the FEM grid file name
		string Configure(char *fname, int idx); // load file names to save the solution at some positions
		string type_of_boundary_condition(char *fname);
		void load_FEM_parameters(char *fname, int& nsd, int& NxFEM, int& NyFEM, int& NzFEM, 
					double& dxFEM, double& dyFEM, double& dzFEM, 
					double& x_minFEM, double& y_minFEM, double& z_minFEM, int& NoObjects);
		string load_object_type(char *fname, int objidx);
		void load_object_parameters(char *fname, int objidx, int& TypeOfMat, double& velocity, MV_Vector<double>& Coord);
		void load_object_parameters(char *fname, int objidx, int& TypeOfMat, double& velocity);

		void load_boundary_data_inversion(char* fname, MV_ColMat<double>& psi, MV_Vector<int>& BndNodes, MV_Vector<double>& bdvalue_tail); 
		void load_grid_file(char *fname, MV_ColMat<double>& Nodes, MV_ColMat<int>& Elements);
		void write2tecplot(char* file, MV_ColMat<double> Nodes, MV_ColMat<int> Elements);
		// +++++ ----------


	
		void ReadConst(char *fname, int& alfa);

		double x_max_;
		double x_min_;
		double y_max_;
		double y_min_;
		double z_min_;
		double z_max_;


};

#endif

