/* function forward_solver3D.C: solve the forward wave equation problem

Compile: make forward_solver3D
How to call: 
./path/forward_solver3D parameter_file.dat 
./path/forward_solver3D parameter_file.dat arg2 coefficient_file
./path/forward_solver3D parameter_file.dat arg2 coefficient_file arg4

where: 
- parameter_file.dat: file contains geometrical, mesh, material information
- arg2: 1 or 0. If arg2 = 1, the coefficient value is read from coefficient_file.
- coefficient_file: file containing the values of the coefficient at the corresponding FEM mesh. 
- arg4: 1 or 0. If arg4 = 0, the solution WILL NOT be written to the file ExactFEM.m

NOTE: The first call is equivalent to arg2 = 0, arg4 = 1.
Last updated: Dec 17, 2012.

*/

static char help[] = "";

#include <iostream>
#include <string>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include "include/wavesWaveEqOpSeismic.h"
#include "include/waveskraftwerk.hh"
#include "include/wavesSDGeometry.h"
#include "include/wavesScalarEqOpOpt.h"
#include <sys/times.h>
#include <mpi.h>


using namespace std;

typedef double real;
typedef MV_ColMat<real> Mat_real;
typedef MV_Vector<real> Vec_real; 

double Maxim(MV_Vector<double> Gradient)
{
	double maxim = Gradient(0);
	for (int i = 0; i < Gradient.size(); i++)
		if (Gradient(i) > maxim)
			maxim = Gradient(i);
	return maxim;
}
double Minim(MV_Vector<double> Gradient)
{
	double minim = Gradient(0);
	for (int i = 0; i < Gradient.size(); i++)
		if (Gradient(i) < minim)
			minim = Gradient(i);
	return minim;
}


//============================================================

int main(int argc, char **argv)
{

	PetscInitialize(&argc, &argv, (char *) 0, help);
	
	const double RandMax = RAND_MAX + 0.1; // maximum random number but converted to double to avoid the integer division
	int nsd, k;
	bool USE_FEM, USE_FDM, USE_RHS, USE_DIRICHLET_FEM, 
	     USE_DIRICHLET_FDM, PRINT_FILES, EXCHANGE, USE_ABSORB;

	double maxTime, rhs, t, dt, noiselevel, omega, randomvalue;
	double guess_velocity = 1.0; // background velocity.
	double y_fix = 1.0; //test the FEM method, check!!!
	MV_Vector<double> Z(4); // location to save data

	int NoTimeSteps, NoObjects, i;

	WavesGridB gg;
	WavesGridB outer_gg;
	WavesSDGeometry sdg;
	WavesSDGeometry sdg_new;
        WavesOutputs out1;

	bool doIncludeCorners = 1;
	char* paramfile = argv[1]; 

	// load the parameters: call function Configure from WavesOutput class
	out1.Configure(paramfile, gg, outer_gg, sdg, sdg_new, nsd, EXCHANGE, USE_FEM, USE_FDM, 
		    USE_DIRICHLET_FEM, USE_DIRICHLET_FDM, USE_ABSORB, USE_RHS, PRINT_FILES, 
		    NoTimeSteps, maxTime, rhs, noiselevel, omega, NoObjects, Z);


	// ---load the boundary nodes:
 	string bdgridFile = out1.boundary_grid_file_name(paramfile);
	const char* bdgrid = bdgridFile.c_str();

	int Nbnodes = out1.load_boundary_grid(bdgrid); // number of boundary node
 	MV_Vector<int> BoundNodes(Nbnodes);
	out1.load_boundary_grid(bdgrid,BoundNodes); // load the boundary grid nodes from the grid file

 	MV_Vector<double> BoundaryData(Nbnodes);

	//-----------------------
	string File1 = out1.Configure(paramfile,1); //file name of the extract data
	const char* file1 = File1.c_str();

	string File2 = out1.Configure(paramfile,2); //file name of the extract data
	const char* file2 = File2.c_str();

	string File3 = out1.Configure(paramfile,3); //file name of the extract data
	const char* file3 = File3.c_str();

	string File4 = out1.Configure(paramfile,4); //file name of the extract data
	const char* file4 = File4.c_str();

	// load the properties of objects:
	int TypeOfMat = 2; 
	double velocity = 1.0; 
	

	if (argc > 1)
	{
		// initialize the forward solver: 
		WavesScalarEqOpOpt p(gg, outer_gg, sdg, nsd, EXCHANGE, USE_FEM, USE_FDM, 
				     USE_RHS, USE_DIRICHLET_FEM, USE_DIRICHLET_FDM, 
				     USE_ABSORB, PRINT_FILES, NoTimeSteps, maxTime, rhs, 
				     TypeOfMat, velocity, guess_velocity, doIncludeCorners);

		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);


		MV_Vector<double> b(gg.getNoNodes()); // vector of material at grid nodes


		WavesOptional opt(gg);	
		WavesOutputs out(gg,b);	
		WavesOutputs out_ex(sdg_new);
		WavesOutputs out_planeXY; // extract the solution on a X-Y plane (Z= const). added by Thanh
		
		// calculate the vector of material (coefficient): 
		b = guess_velocity;		
		if (argc > 3 && atoi(argv[2]) == 1)
		{
		    out_ex.ReadSolfromFile(argv[3], b);
		}
		else
		  for (i=0; i < NoObjects; i++)
		    {
			out1.load_object_parameters(paramfile, i, TypeOfMat, velocity); // the first object
			out1.VecGetTypeofMat(gg, nsd, TypeOfMat, velocity, b);
		    }
		

		// save the vector of material:
		if (rank == 0)
			out.WriteToFile((char *) "Exact_param.m", b);

		if (rank == 0)
		{
			if (nsd == 2)
				int ierr = opt.writeInp2D((char *) "density.inp", &gg, b, 1, -1);
			else if (nsd == 3)
				int ierr = opt.writeInp3D((char *) "density.inp", &gg, b, 1);
		}

		// time paramters:
		dt = p.InitTime();
		t = 0.0;
		MV_Vector<real> thetimesteps(NoTimeSteps);
		tms utime0, tottime;
		times(&tottime);

		const char* exchangeMaskfile = "exchangeMask.m"; //test if the exchangeMask file is available
		ifstream ifile(exchangeMaskfile);

		if (EXCHANGE)
		{
			p.InitFDM();

			if (nsd == 2)
				p.Init2DFEM(b);
			else if (nsd == 3)
				p.Init3DFEM(b);

			if (ifile) 
			{ // The file exists, 
				
				p.InitExchangeReadMask(); // read exchangeMask.m file. The file must be generated by the 
			}
			else
				//p.InitExchangeCommon();
				p.InitExchangeStructCommon(); // print the exchangeMask.m file
		}

		if (USE_FDM)
			p.InitPlaneWave_GLK();
		if (USE_FEM)
			p.Init2DFEM(b);
		else if  (nsd ==3)
		  p.Init3DFEM(b);


		// find nodes where save data in FDM domain

		int n_i = sdg.getN_i();
		int n_j = sdg.getN_j();
		
		MV_Vector<real>   E_FEM(gg.getNoNodes());
		E_FEM = 0.0;
		
		
		// to save the solution on some XY planes:Thanh
		int XYsize = n_i * n_j; 
		MV_Vector<int> Nodes_XY1(XYsize); // nodes on a X-Y plane (z = const)
		MV_Vector<int> Nodes_XY2(XYsize); // nodes on a X-Y plane (z = const)
		MV_Vector<int> Nodes_XY3(XYsize); // nodes on a X-Y plane (z = const)
		MV_Vector<int> Nodes_XY4(XYsize); // nodes on a X-Y plane (z = const)

		MV_Vector<double> Sol_XY1(XYsize); // solution values on a X-Y plane (z = const)
		MV_Vector<double> Sol_XY2(XYsize); // solution values on a X-Y plane (z = const)
		MV_Vector<double> Sol_XY3(XYsize); // solution values on a X-Y plane (z = const)
		MV_Vector<double> Sol_XY4(XYsize); // solution values on a X-Y plane (z = const)
		
   		sdg.coord2node3Dz(Z(0), Nodes_XY1); // calculate the node indices on the plane {z=z1}
		sdg.coord2node3Dz(Z(1), Nodes_XY2);
		sdg.coord2node3Dz(Z(2), Nodes_XY3);
		sdg.coord2node3Dz(Z(3), Nodes_XY4);

		const char* FEMfile = "ExactFEM.m"; 
		const char* SolFEMboundary = "SolFEMBoundary.m"; 
		
		srand(time(NULL)); // change the seed of the random function to create different random values for different runs
		//load the type of boundary condition:
		string TypeOfBC = out1.type_of_boundary_condition(paramfile);

		// remove old files:
		ifstream ifile0(FEMfile); 
		if (ifile0) {  	remove(FEMfile);}
		ifstream ifile1(file1); 
		if (ifile1){ remove(file1);}
		ifstream ifile2(file2); 
		if (ifile2){ remove(file2);}
		ifstream ifile3(file3); 
		if (ifile3){ remove(file3);}
		ifstream ifile4(file4); 
		if (ifile4){ remove(file4);}
		ifstream ifile5(SolFEMboundary); 
		if (ifile5){ remove(SolFEMboundary);}
	
		out_ex.Write_array_my_part(SolFEMboundary,BoundNodes);

		//======================================================================
		// Step 1. Compute the exact solution of the Wave equation and
		//         keep results in the arrays Forward_u_fem and Forward_u_fdm
		//======================================================================
		//main loop

		for (k = 0; k < NoTimeSteps; k++)
		{   

			if (EXCHANGE)
			{

				// =====  worked before with plane wave initialized at the left of FDM ===========
				//cout << "*****************works fdm***********************" << endl;
	      			if (TypeOfBC.compare("planewave_z_max") == 0)
					p.PlaneWaveBackFDM(k,omega);


				// in the case of normal wave equation
				p.WaveEqSolverFEMforDifMat(t, k, TypeOfMat, velocity);

				p.ApplyExchangeCommon();

				p.ApplySwap();

			}

			if (USE_FEM)
			{

				p.PlaneWaveFEMforDifMat(t, k, TypeOfMat, velocity, y_fix);

				p.ApplySwap();
			}

			if (USE_FDM)
			{

				p.PlaneWaveTopFDMSeismic(k);

				p.ApplySwap();
			}

			// write solution at the obs.points

			p.Save_Sol_FDM(k,Sol_XY1,Nodes_XY1); // extract the solution at the nodes given by Nodes_XY1
			p.Save_Sol_FDM(k,Sol_XY2,Nodes_XY2);
			p.Save_Sol_FDM(k,Sol_XY3,Nodes_XY3);
			p.Save_Sol_FDM(k,Sol_XY4,Nodes_XY4);

			//save global FEM solution
			E_FEM = 0.0;
			p.SaveSolutionsFEM(k,E_FEM);


			// add noise to the exact solution
			
			for (i =0; i < gg.getNoNodes(); i++)
			  {
			    randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
			    E_FEM(i) = E_FEM(i)*(1 + randomvalue*noiselevel);
			    
			  }

			for (i =0; i < XYsize; i++)
			  {
			    randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
			    Sol_XY1(i) = Sol_XY1(i)*(1 + randomvalue*noiselevel);
			    
			    randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
			    Sol_XY2(i) = Sol_XY2(i)*(1 + randomvalue*noiselevel);
			    
			    randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
			    Sol_XY3(i) = Sol_XY3(i)*(1 + randomvalue*noiselevel);
			    
			    randomvalue = 2*(rand()/RandMax) - 1.0;   // randomvalue is a random number at the interval [-1;1)
			    Sol_XY4(i) = Sol_XY4(i)*(1 + randomvalue*noiselevel);			
			    
			  }
 
			// extract the solution at the boundary of the FEM domain (for inversion)
			for (i = 0; i < Nbnodes; i++)
				BoundaryData(i) = E_FEM(BoundNodes(i)-1);


			// write global solution to compute Laplace transform
			if (rank == 0)
			  {
			    if (argc <=  4 || (argc > 4 && atoi(argv[4]) != 0))
				out_ex.Write_array_my_part(FEMfile,E_FEM,k);

			    out_planeXY.Write_array_my_part(file1,Sol_XY1,k);
			    out_planeXY.Write_array_my_part(file2,Sol_XY2,k);
			    out_planeXY.Write_array_my_part(file3,Sol_XY3,k);
			    out_planeXY.Write_array_my_part(file4,Sol_XY4,k);
			    out_ex.Write_array_my_part(SolFEMboundary,BoundaryData,k+1);
			  }


			thetimesteps(k) = t;
			t += dt;  

		}

		times(&utime0);

		cout << "##### total time = " << (utime0.tms_utime - tottime.tms_utime) / 100.0 << " #####\n";
		cout << "For " << int(maxTime / dt) << " steps\n" << endl;
		cout << "Min computational total_time/nrSteps*nrNodes is " 
		     << ((utime0.tms_utime - tottime.tms_utime) / 100.0) / (NoTimeSteps * sdg.getNoNodes()) << endl;

	}
	else
	{
		cout << "Usage:  CONFIGURATION, see arguments argv[1], argv[2], argv[3]  " << endl;
	}



	PetscFinalize();
	return 0;
}

