/* This program generates 3D mesh in the FEM domain. 
The geometrical parameters is given in the input data file with extension .dat
Input: 	argv[1]: parameter file name
	argv[2]: (optional) the *.inp file name to save the coefficient with grid
	argv[3]: (optional) the file name to save the coefficient at nodes 
Output: the grid file with the name provided in the parameter file
Compile: make mesh_gen

Example of use: 

./mesh_gen forpar_new_full.dat 
./mesh_gen forpar_new_full.dat density.inp coefficient.m

@Nguyen Trung Thanh UNCC 2012. 

*/

static char help[] ="";

#include "include/wavesSDOperator.h"
#include "include/wavesSDGeometry.h"
#include "include/wavesSDIndexes.h"
#include "include/wavesGridB.h"
#include "include/wavesutil2.h"
#include "include/wavesFEcomp.h"
#include "include/wavesfindBoundaryNodes.h"
#include "include/wavesAddGrids.h"
#include "include/wavesOptional.h"
#include "include/wavesOutputs.h"

#include <iostream>
#include <sys/times.h>
#include <mpi.h>

using namespace std;


typedef double real; 
typedef MV_ColMat<real> Mat_real;
typedef MV_ColMat<int> Mat_int;
typedef MV_Vector<real> Vec_real;

void makeGrid(WavesGridB& square_grid,
              WavesGridB& triangle_grid,
	      WavesGridB& new_grid)
{
  // find the size of the grid.
  int nel = square_grid.getNoElms();
  int nno = square_grid.getNoNodes();
  int maxnne = square_grid.getMaxNoNodesInElm ();
  int nsd = square_grid.getNoSpaceDim();
  int nbind = 0;
  double eps = 1e-6;

  MV_Vector<real> glob_pt1(2);
  MV_Vector<real> glob_pt2(2);
  MV_Vector<real> glob_pt3(2);
  MV_Vector<real> glob_pt4(2);
  glob_pt1 = 0.0;
  glob_pt2 = 0.0;
  glob_pt3 = 0.0;
  glob_pt4 = 0.0;

  int i,d,e,n;
  int n0 = 0;
  int n1 = 0;
  int n2 = 0;
  int n3 = 0;
  int ecount = 0;
  ElementType etype = ELMTRI1;
   
  maxnne=3;
  new_grid.redim( nsd, triangle_grid.getNoNodes(),
		   triangle_grid.getNoElms(),
		   maxnne, nbind, etype);

  //loop for square_grid
  
    for(e=0;e<nel;e++)
      {
	
	int n0e = square_grid.loc2glob( e,0 );
	int n1e = square_grid.loc2glob( e,1 );
	int n2e = square_grid.loc2glob( e,2 );
	int n3e = square_grid.loc2glob( e,3 );
	
	glob_pt1(0) = square_grid.getCoor(n0e,0);
	glob_pt1(1) = square_grid.getCoor(n0e,1); 
	glob_pt2(0) = square_grid.getCoor(n1e,0);
	glob_pt2(1) = square_grid.getCoor(n1e,1);
	
	glob_pt3(0) = square_grid.getCoor(n2e,0);
	glob_pt3(1) = square_grid.getCoor(n2e,1); 
	glob_pt4(0) = square_grid.getCoor(n3e,0);
	glob_pt4(1) = square_grid.getCoor(n3e,1); 

	    
        int marker = 0;
	for (i = 0; i < triangle_grid.getNoNodes();i++)
	  { 
              
	    if ( fabs(glob_pt1(0) - triangle_grid.getCoor(i,0)) < eps &&
		 fabs(glob_pt1(1) - triangle_grid.getCoor(i,1)) < eps)
	      {  n0 = i; marker++; break;}
	  }
    
	for (i = 0; i < triangle_grid.getNoNodes();i++)
	  { 
	    if ( fabs(glob_pt2(0) - triangle_grid.getCoor(i,0)) < eps &&
		 fabs(glob_pt2(1) - triangle_grid.getCoor(i,1)) < eps)
	      {     marker++;  n1 = i; break;}
	  }
	
	for (i = 0; i < triangle_grid.getNoNodes();i++)
	  { 
	    if ( fabs(glob_pt3(0) - triangle_grid.getCoor(i,0)) < eps &&
		 fabs(glob_pt3(1) - triangle_grid.getCoor(i,1)) < eps)
	      {marker++;  n2 = i;break;}
	  }
	
	for (i = 0; i < triangle_grid.getNoNodes();i++)
	  { 
	    if ( fabs(glob_pt4(0) - triangle_grid.getCoor(i,0)) < eps &&
		 fabs(glob_pt4(1) - triangle_grid.getCoor(i,1)) < eps)
	      { marker++;  n3 = i; break;}
	  }
	 
	//if all nodes belongs to square from square_grid
	if (marker == 4 )
	  { 
	    // cout<<"n1 "<<n1<<"n0 "<<n0<<"n2 "<<n2<<endl;
	    int mattype = 0;
	    new_grid.putLoc2glob(ecount,0,n1);
	    new_grid.putLoc2glob(ecount,1,n0);
	    new_grid.putLoc2glob(ecount,2,n2);
	    new_grid.setMaterialType(ecount,mattype);
	    ecount++;
	    cout<<ecount<<" n1 "<<n1<<"n0 "<<n0<<"n2 "<<n2<<endl;
	    new_grid.putLoc2glob(ecount,0,n3);
	    new_grid.putLoc2glob(ecount,1,n2);
	    new_grid.putLoc2glob(ecount,2,n0);
	    new_grid.setMaterialType(ecount,mattype); 
	    ecount++;
	    cout<<ecount<<" n3 "<<n3<<"n2 "<<n2<<"n0 "<<n0<<endl;
	  }
      }  //for all elements in square_grid


  for( n=0; n < triangle_grid.getNoNodes();n++)
    for(d=0;d<nsd;d++)
      new_grid.putCoor(n,d,triangle_grid.getCoor(n,d));
}


     

void makeGrid(WavesSDGeometry& sgrid,
	      WavesGridB& ugrid,ElementType etype,
	      int mattype)
{
  // find the size of the grid.
  int nel = sgrid.getNoElms();
  int nno = sgrid.getNoNodes();
  int maxnne = sgrid.getMaxNoNodesInElm ();
  int nsd = sgrid.getNoSpaceDim();
  int nbind = 0;


  int d,e,n;
  int ecount = 0;
  if(etype==ELMQUAD1) {
    maxnne=4;
    ugrid.redim ( nsd, nno, nel, maxnne, nbind, etype);
    for(e=0;e<nel;e++){
      int n0 = sgrid.loc2glob(e,0);
      int n1 = sgrid.loc2glob(e,1);
      int n2 = sgrid.loc2glob(e,2);
      int n3 = sgrid.loc2glob(e,3);
      //      int mattype = sgrid.getMaterialType(e);
      //      int mattype = 0;
      ugrid.putLoc2glob(ecount,0,n0);
      ugrid.putLoc2glob(ecount,1,n1);
      ugrid.putLoc2glob(ecount,2,n3);
      ugrid.putLoc2glob(ecount,3,n2);
      ugrid.setMaterialType(ecount,mattype);
      ecount++;
    }
  }
  else  if(etype==ELMTRI1) {
    maxnne=3;
    ugrid.redim ( nsd, nno, 2*nel, maxnne, nbind, etype);
    for(e=0;e<nel;e++){
      int n0 = sgrid.loc2glob(e,0);
      int n1 = sgrid.loc2glob(e,1);
      int n2 = sgrid.loc2glob(e,2);
      int n3 = sgrid.loc2glob(e,3);
      //      int mattype = sgrid.getMaterialType(e);
      // int mattype = 0;
      ugrid.putLoc2glob(ecount,0,n0);
      ugrid.putLoc2glob(ecount,1,n1);
      ugrid.putLoc2glob(ecount,2,n2);
      ugrid.setMaterialType(ecount,mattype);
      ecount++;
      ugrid.putLoc2glob(ecount,0,n3);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n1);
      ugrid.setMaterialType(ecount,mattype); 
      ecount++;
    }
  }
  else if ( etype==ELMTET1) {
    maxnne=4; // tetrahedron
    ugrid.redim ( nsd, nno, 6*nel, maxnne, nbind, etype); 
    for(e=0;e<nel;e++){
      int n0 = sgrid.loc2glob(e,0);
      int n1 = sgrid.loc2glob(e,1);
      int n2 = sgrid.loc2glob(e,2);
      int n3 = sgrid.loc2glob(e,3);
      int n4 = sgrid.loc2glob(e,4);
      int n5 = sgrid.loc2glob(e,5);
      int n6 = sgrid.loc2glob(e,6);
      int n7 = sgrid.loc2glob(e,7);
      
      //int mattype;// = sgrid.getMaterialType(e);
      //  mattype = 1;

      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n1);
      ugrid.putLoc2glob(ecount,3,n0);

      ecount++;
      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n0);
      ugrid.putLoc2glob(ecount,3,n4);
      ecount++;

      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n4);
      ugrid.putLoc2glob(ecount,3,n6);
      ecount++;

      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n3);
      ugrid.putLoc2glob(ecount,3,n1);
      ecount++;
    
      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n7);
      ugrid.putLoc2glob(ecount,3,n3);
      ecount++;
    
      ugrid.setMaterialType(ecount,mattype);
      ugrid.putLoc2glob(ecount,0,n5);
      ugrid.putLoc2glob(ecount,1,n2);
      ugrid.putLoc2glob(ecount,2,n6);
      ugrid.putLoc2glob(ecount,3,n7);      
      ecount++;
    }
  }
  else if ( etype==ELMHEX1) {
    maxnne=8; // cube
    ugrid.redim ( nsd, nno, nel, maxnne, nbind, etype); 
    for(e=0; e<nel ;e++){
      
      int n0 = sgrid.loc2glob(e,0);
      int n1 = sgrid.loc2glob(e,1);
      int n2 = sgrid.loc2glob(e,2);
      int n3 = sgrid.loc2glob(e,3);
      int n4 = sgrid.loc2glob(e,4);
      int n5 = sgrid.loc2glob(e,5);
      int n6 = sgrid.loc2glob(e,6);
      int n7 = sgrid.loc2glob(e,7);

      //      int mattype;// = sgrid.getMaterialType(e);
      // mattype = 1;

      ugrid.setMaterialType(e,mattype);
      //      int i;
      //      for(i=0;i<maxnne;i++)
      ugrid.putLoc2glob(e,0,n0);
      ugrid.putLoc2glob(e,1,n1);
      ugrid.putLoc2glob(e,3,n2);
      ugrid.putLoc2glob(e,2,n3);
      ugrid.putLoc2glob(e,4,n4);
      ugrid.putLoc2glob(e,5,n5);
      ugrid.putLoc2glob(e,7,n6);
      ugrid.putLoc2glob(e,6,n7);
    }
  }

  for(n=0;n<nno;n++)
    for(d=0;d<nsd;d++)
      ugrid.putCoor(n,d,sgrid.getCoor(n,d));
}

int checkInsideSphere(double x, double y, double z, double center_x, double center_y, 
	    double center_z, double radius)
{
  int check;

       double diff_x = fabs(x - center_x);
       double diff_y = fabs(y - center_y);
       double diff_z = fabs(z - center_z);
       
       double d = sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);

      if (d <= radius)  
	  check= 1;
      else 
	  check = 0;

  return check;
}

int checkInsideSphere(WavesGridB& grid, int n, double center_x, double center_y, 
	    double center_z, double radius)
{
  double x,y,z;
  
       x =  grid.getCoor(n,0);
       y =  grid.getCoor(n,1);
       z =  grid.getCoor(n,2);

  return checkInsideSphere(x,y,z,center_x,center_y,center_z,radius);
}


int checkInsideCylinderX(double x, double y, double z, double center_y, double center_z, 
	    double radius, double x_min, double x_max)
{
	int check;
 
       double diff_1 = fabs(y - center_y);
       double diff_2 = fabs(z - center_z);
       
       double d = sqrt(diff_1*diff_1 + diff_2*diff_2);

      if ((d <= radius) && (x >= x_min) && (x <= x_max))  
	  check= 1;
      else 
	  check = 0;

  return check;
}
int checkInsideCylinderX(WavesGridB& grid, int n, double center_y, double center_z, 
	    double radius, double x_min, double x_max)
{
 	double x,y,z;
  
       x =  grid.getCoor(n,0);
       y =  grid.getCoor(n,1);
       z =  grid.getCoor(n,2);

  return checkInsideCylinderX(x,y,z,center_y,center_z,radius,x_min,x_max);
}


int checkInsideCylinderY(double x, double y, double z, double center_x, double center_z, 
	    double radius, double y_min, double y_max)
{
	int check;

       double diff_1 = fabs(x - center_x);
       double diff_2 = fabs(z - center_z);
       
       double d = sqrt(diff_1*diff_1 + diff_2*diff_2);

      if ((d <= radius) && (y >= y_min) && (y <= y_max))  
	  check= 1;
      else 
	  check = 0;


  return check;
}

int checkInsideCylinderY(WavesGridB& grid, int n, double center_x, double center_z, 
	    double radius, double y_min, double y_max)
{
  	double x,y,z;
  
       x =  grid.getCoor(n,0);
       y =  grid.getCoor(n,1);
       z =  grid.getCoor(n,2);

	return checkInsideCylinderY(x,y,z,center_x,center_z,radius,y_min,y_max);
}

int checkInsideCylinderZ(double x, double y, double z, double center_x, double center_y, 
	    double radius, double z_min, double z_max)
{
	int check;

       double diff_1 = fabs(y - center_y);
       double diff_2 = fabs(x - center_x);
       
       double d = sqrt(diff_1*diff_1 + diff_2*diff_2);

      if ((d <= radius) && (z >= z_min) && (z <= z_max))  
	  check= 1;
      else 
	  check = 0;


  return check;
}
int checkInsideCylinderZ(WavesGridB& grid, int n, double center_x, double center_y, 
	    double radius, double z_min, double z_max)
{
  	double x,y,z;
  
       x =  grid.getCoor(n,0);
       y =  grid.getCoor(n,1);
       z =  grid.getCoor(n,2);

  return checkInsideCylinderZ(x,y,z,center_x,center_y,radius,z_min,z_max);
}

int checkInsidePrism(double x, double y, double z, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
       int check;

      if ((x >= x_min) && (x <= x_max) && (y >= y_min) && (y <= y_max) && (z >= z_min) && (z <= z_max))  
	  check = 1;
      else 
	  check = 0;


  return check;
}


int checkInsidePrism(WavesGridB& grid, int n, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max)
{
       double x,y,z;
  
       x =  grid.getCoor(n,0);
       y =  grid.getCoor(n,1);
       z =  grid.getCoor(n,2);

	return checkInsidePrism(x,y,z,x_min,x_max,y_min,y_max,z_min,z_max);
}

// check a point with coordinate (x,y,z) is inside an object given by its parameters in Coord:
int checkInsideObject(string objecttype, MV_Vector<double> Coord, double x, double y, double z)
{
	double center_x, center_y, center_z, radius, x_min, x_max, y_min, y_max, z_min, z_max; 		
	if (objecttype.compare("sphere") == 0)
	{
	  	center_x = Coord(0); center_y = Coord(1); center_z = Coord(2); radius = Coord(3); 
		return checkInsideSphere(x,y,z,center_x,center_y,center_z,radius);
	}
	else if ((objecttype.compare("cylinder_x") == 0) || (objecttype.compare("cylinder-x") == 0))
	{
		center_y = Coord(0); center_z = Coord(1); radius = Coord(2); 
		x_min = Coord(3); x_max = Coord(4);
		return checkInsideCylinderX(x,y,z,center_y,center_z,radius,x_min,x_max);
	}
	else if ((objecttype.compare("cylinder_y") == 0) || (objecttype.compare("cylinder-y") == 0))
	{
		center_x = Coord(0); center_z = Coord(1); radius = Coord(2); 
		y_min = Coord(3); y_max = Coord(4);
		return checkInsideCylinderY(x,y,z,center_x,center_z,radius,y_min,y_max);
	}
	else if ((objecttype.compare("cylinder_z") == 0) || (objecttype.compare("cylinder-z") == 0))
	{
		center_x = Coord(0); center_y = Coord(1); radius = Coord(2); 
		z_min = Coord(3); z_max = Coord(4);
		return checkInsideCylinderZ(x,y,z,center_x,center_y,radius,z_min,z_max);
	}
	else if (objecttype.compare("prism") == 0)
	{
		x_min = Coord(0); x_max = Coord(1);
		y_min = Coord(2); y_max = Coord(3);
		z_min = Coord(4); z_max = Coord(5); 
	return checkInsidePrism(x,y,z,x_min,x_max,y_min,y_max,z_min,z_max);
	}
	else
		return 0;
}

void SetMatType(WavesGridB& grid, string objecttype, int TypeOfMat, MV_Vector<double> Coord)
{
  	int nel = grid.getNoElms();
 	int nsd = grid.getNoSpaceDim();
  	int e;
  	int code1, code2, code3,code4;
  	int n0, n1, n2, n3;
  	double center_x, center_y, center_z, radius, x_min, x_max, y_min, y_max, z_min, z_max; 	

  	if (nsd==3)
  	{
		if (objecttype.compare("sphere") == 0)
		{
		  	center_x = Coord(0); center_y = Coord(1); center_z = Coord(2); radius = Coord(3); 
		        for(e=0;e<nel;e++)
      			{
	  			n0 = grid.loc2glob(e,0);
	  			n1 = grid.loc2glob(e,1);
	  			n2 = grid.loc2glob(e,2);
	  			n3 = grid.loc2glob(e,3);
	          		code1 = checkInsideSphere(grid, n0, center_x, center_y, center_z, radius);
	          		code2 = checkInsideSphere(grid, n1, center_x, center_y, center_z, radius);
	          		code3 = checkInsideSphere(grid, n2, center_x, center_y, center_z, radius);
		  		code4 = checkInsideSphere(grid, n3, center_x, center_y, center_z, radius);
				if (code1 == 1 && code2 == 1 && code3 == 1 && code4 == 1)
	  	 	      		grid.setMaterialType(e, TypeOfMat);

	  		}
 		}
		else if ((objecttype.compare("cylinder_x") == 0) || (objecttype.compare("cylinder-x") == 0))
	  	{
			center_y = Coord(0); center_z = Coord(1); radius = Coord(2); 
			x_min = Coord(3); x_max = Coord(4);
		        for(e=0;e<nel;e++)
      			{
	  			n0 = grid.loc2glob(e,0);
	  			n1 = grid.loc2glob(e,1);
	  			n2 = grid.loc2glob(e,2);
	  			n3 = grid.loc2glob(e,3);
				code1 = checkInsideCylinderX(grid, n0, center_y, center_z, radius, x_min, x_max);
				code2 = checkInsideCylinderX(grid, n1, center_y, center_z, radius, x_min, x_max);
				code3 = checkInsideCylinderX(grid, n2, center_y, center_z, radius, x_min, x_max);
				code4 = checkInsideCylinderX(grid, n3, center_y, center_z, radius, x_min, x_max);
				if (code1 == 1 && code2 == 1 && code3 == 1 && code4 == 1)
	  	 	      		grid.setMaterialType(e, TypeOfMat);

	  		}

		}		    
		else if ((objecttype.compare("cylinder_y") == 0) || (objecttype.compare("cylinder-y") == 0))
		{
			center_x = Coord(0); center_z = Coord(1); radius = Coord(2); 
			y_min = Coord(3); y_max = Coord(4);
		        for(e=0;e<nel;e++)
      			{
	  			n0 = grid.loc2glob(e,0);
	  			n1 = grid.loc2glob(e,1);
	  			n2 = grid.loc2glob(e,2);
	  			n3 = grid.loc2glob(e,3);
				code1 = checkInsideCylinderY(grid, n0, center_x, center_z, radius, y_min, y_max);
				code2 = checkInsideCylinderY(grid, n1, center_x, center_z, radius, y_min, y_max);
				code3 = checkInsideCylinderY(grid, n2, center_x, center_z, radius, y_min, y_max);
				code4 = checkInsideCylinderY(grid, n3, center_x, center_z, radius, y_min, y_max);
				if (code1 == 1 && code2 == 1 && code3 == 1 && code4 == 1)
	  	 	      		grid.setMaterialType(e, TypeOfMat);

	  		}
		}
		else if ((objecttype.compare("cylinder_z") == 0) || (objecttype.compare("cylinder-z") == 0))
		{
			center_x = Coord(0); center_y = Coord(1); radius = Coord(2); 
			z_min = Coord(3); z_max = Coord(4);
		        for(e=0;e<nel;e++)
      			{
	  			n0 = grid.loc2glob(e,0);
	  			n1 = grid.loc2glob(e,1);
	  			n2 = grid.loc2glob(e,2);
	  			n3 = grid.loc2glob(e,3);
				code1 = checkInsideCylinderZ(grid, n0, center_x, center_y, radius, z_min, z_max);
				code2 = checkInsideCylinderZ(grid, n1, center_x, center_y, radius, z_min, z_max);
				code3 = checkInsideCylinderZ(grid, n2, center_x, center_y, radius, z_min, z_max);
				code4 = checkInsideCylinderZ(grid, n3, center_x, center_y, radius, z_min, z_max);
				if (code1 == 1 && code2 == 1 && code3 == 1 && code4 == 1)
	  	 	      		grid.setMaterialType(e, TypeOfMat);

	  		}
		}
		else if (objecttype.compare("prism") == 0)
		{

			x_min = Coord(0); x_max = Coord(1);
			y_min = Coord(2); y_max = Coord(3);
			z_min = Coord(4); z_max = Coord(5); 
		        for(e=0;e<nel;e++)
      			{
	  			n0 = grid.loc2glob(e,0);
	  			n1 = grid.loc2glob(e,1);
	  			n2 = grid.loc2glob(e,2);
	  			n3 = grid.loc2glob(e,3);
				code1 = checkInsidePrism(grid, n0, x_min, x_max, y_min, y_max, z_min, z_max); 
				code2 = checkInsidePrism(grid, n1, x_min, x_max, y_min, y_max, z_min, z_max);
				code3 = checkInsidePrism(grid, n2, x_min, x_max, y_min, y_max, z_min, z_max);
				code4 = checkInsidePrism(grid, n3, x_min, x_max, y_min, y_max, z_min, z_max);
				if (code1 == 1 && code2 == 1 && code3 == 1 && code4 == 1)
	  	 	      		grid.setMaterialType(e, TypeOfMat);

	  		}
		}
		else
			cout << "ERROR: The object type " << objecttype << " is not recognized" << endl;
      	
  	}

}



void SetMat(WavesGridB& grid, string objecttype, double velocity, MV_Vector<double> Coord, MV_Vector<double>& MatVec)
{
  	int nno = grid.getNoNodes();
  	int nsd = grid.getNoSpaceDim();
	int code, i; 
	double center_x, center_y, center_z, radius, x_min, x_max, y_min, y_max, z_min, z_max; 		
	
	if (nsd==3)
	{
		if (objecttype.compare("sphere") == 0)
		{
		  	center_x = Coord(0); center_y = Coord(1); center_z = Coord(2); radius = Coord(3); 
			for(i = 0; i< nno; i++)
		      	{	  
				code = checkInsideSphere(grid, i, center_x, center_y, center_z, radius);
				if (code == 1)
				{
				      MatVec(i)=velocity;
				}
			      
			}
 		}
		else if ((objecttype.compare("cylinder_x") == 0) || (objecttype.compare("cylinder-x") == 0))
	  	{
			center_y = Coord(0); center_z = Coord(1); radius = Coord(2); 
			x_min = Coord(3); x_max = Coord(4);
			for(i = 0; i< nno; i++)
		      	{	  
				code = checkInsideCylinderX(grid, i, center_y, center_z, radius, x_min, x_max);
				if (code == 1)
				{
				      MatVec(i)=velocity;
				}
			      
			}

		}		    
		else if ((objecttype.compare("cylinder_y") == 0) || (objecttype.compare("cylinder-y") == 0))
		{
			center_x = Coord(0); center_z = Coord(1); radius = Coord(2); 
			y_min = Coord(3); y_max = Coord(4);
			for(i = 0; i< nno; i++)
		      	{	  
				code = checkInsideCylinderY(grid, i, center_x, center_z, radius, y_min, y_max);
				if (code == 1)
				{
				      MatVec(i)=velocity;
				}
			      
			}
		}
		else if ((objecttype.compare("cylinder_z") == 0) || (objecttype.compare("cylinder-z") == 0))
		{
			center_x = Coord(0); center_y = Coord(1); radius = Coord(2); 
			z_min = Coord(3); z_max = Coord(4);
			for(i = 0; i< nno; i++)
		      	{	  
				code = checkInsideCylinderZ(grid, i, center_x, center_y, radius, z_min, z_max);
				if (code == 1)
				{
				      MatVec(i)=velocity;
				}
			      
			}
		}
		else if (objecttype.compare("prism") == 0)
		{
			x_min = Coord(0); x_max = Coord(1); 
			y_min = Coord(2); y_max = Coord(3);
			z_min = Coord(4); z_max = Coord(5);
			for(i = 0; i< nno; i++)
		      	{	  
				code = checkInsidePrism(grid, i, x_min, x_max, y_min, y_max, z_min, z_max); 
				if (code == 1)
				{
				      MatVec(i)=velocity;
				}
			      
			}
		}
		else
			cout << "ERROR:  The object type " << objecttype << " is not recognized" << endl;
	}
}

Vec_real linspace(double x1, double x2, int NoPoints)
{
	Vec_real points(NoPoints);
	double dx = (x2-x1)/NoPoints; 

	for (int i = 0; i < NoPoints; i++)
		points(i) = x1+ i*dx;
	return points; 
}

double mean_value(Vec_real x)
{
	double mean = 0.0; 
	for (int i = 0; i < x.size(); i++)
		mean += x(i);

	return mean/x.size();
}
	
void SetMat2(WavesGridB& grid, string objecttype, double velocity, MV_Vector<double> Coord, double dx, double dy, double dz, MV_Vector<double>& MatVec)
{
  int nno = grid.getNoNodes();
  int nsd = grid.getNoSpaceDim();
  int code, code1, code2, code3, code4, code5, code6, code7, code8, n, ix, iy, iz; 
  double x, y, z, xl, yl, zl, dxl, dyl, dzl; 
  int nohalfNodes = 2; 
  int nlocalNodes = 2*nohalfNodes+1; // at each grid node	
	
  dxl = dx/nlocalNodes; 
  dyl = dy/nlocalNodes; 
  dzl = dz/nlocalNodes; 

  if (nsd==3)
  {
	Vec_real localcoeff(nlocalNodes*nlocalNodes*nlocalNodes); 
	for (n=0; n < nno; n++)
	{   localcoeff = 1.0;  	
	    x =  grid.getCoor(n,0);
       	    y =  grid.getCoor(n,1);
       	    z =  grid.getCoor(n,2);
	    code1 = checkInsideObject(objecttype, Coord, x-dx/2, y-dy/2, z-dz/2);
	    code2 = checkInsideObject(objecttype, Coord, x+dx/2, y-dy/2, z-dz/2);
            code3 = checkInsideObject(objecttype, Coord, x-dx/2, y+dy/2, z-dz/2);
	    code4 = checkInsideObject(objecttype, Coord, x+dx/2, y+dy/2, z-dz/2);
	    code5 = checkInsideObject(objecttype, Coord, x-dx/2, y-dy/2, z+dz/2);
	    code6 = checkInsideObject(objecttype, Coord, x+dx/2, y-dy/2, z+dz/2);
	    code7 = checkInsideObject(objecttype, Coord, x-dx/2, y+dy/2, z+dz/2);
            code8 = checkInsideObject(objecttype, Coord, x+dx/2, y+dy/2, z+dz/2);
	    if ((code1==1)||(code2==1) ||(code3==1)||(code4==1)||(code5==1)||(code6==1)||(code7==1)||(code8==1)) 
		{    for (iz = 0; iz < nlocalNodes; iz++)
		    {	zl = z + dzl*(iz-nohalfNodes); 
			for (iy = 0; iy < nlocalNodes; iy++)
			{   yl = y + dyl*(iy-nohalfNodes);
			    for (ix = 0; ix < nlocalNodes; ix++)
			    {   xl = x + dxl*(ix-nohalfNodes);				  
		   		code = checkInsideObject(objecttype, Coord, xl, yl, zl);
				if (code == 1)
				localcoeff(ix + nlocalNodes*iy + nlocalNodes*nlocalNodes*iz) = velocity;	
			    }
			}
		    }
		    if (MatVec(n) == 1)      
			MatVec(n) = mean_value(localcoeff);	
		}							
	}
  }
}



double Maxim(MV_Vector<double> Gradient)
{
	double maxim = Gradient(0);
	for (int i = 0; i < Gradient.size(); i++)
		if (Gradient(i) > maxim)
			maxim = Gradient(i);
	return maxim;
}
//===========================================================================
  

int main(int argc, char **argv)
{
   if (argc > 1)
   {
 	// get the input parameters from the .dat file:
   	double dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM;
  
  	int nsd, NxFEM, NyFEM, NzFEM;
 	int NoObjects, i;
	char* fname = argv[1];   	

	WavesOutputs out1; 
	//load the geometrical parameters:
	out1.load_FEM_parameters(fname, nsd, NxFEM, NyFEM, NzFEM, dxFEM, dyFEM, dzFEM, 
				x_minFEM, y_minFEM, z_minFEM, NoObjects);
 	
  	// create the mesh for the homogeneous medium: 
 	ElementType etype = ELMTET1;
 
 	// Create the FDM grid in the FEM domain: 
 	WavesSDGeometry outer1(NxFEM, NyFEM, NzFEM, nsd, dxFEM, dyFEM, dzFEM, x_minFEM, y_minFEM, z_minFEM);
 
 	// initialize FEM mesh using structure of FDM mesh (outer1) 
 	WavesGridB tet1;
 	makeGrid(outer1,tet1,etype,1);
 	WavesOptional opt(tet1);

 	// add the objects into the mesh:
 	MV_Vector<double> Markers; 
 	Markers.newsize(tet1.getNoNodes());

	string objecttype; 
	MV_Vector<double> Coord(6); 
	Coord = 0.0; 
	double velocity; 
	int TypeOfMat; 
 	
	MV_Vector<double> MatVec(tet1.getNoNodes()); //vector of coefficient values at nodes
	MatVec = 1.0;
	cout << "Number of objects: " << NoObjects << endl;	
	for (i = 0; i < NoObjects; i++)
	{	
		objecttype = out1.load_object_type(fname, i); 
		cout << "Object no. " << i+1 << ": object type: "  <<  objecttype << endl;

		out1.load_object_parameters(fname, i, TypeOfMat, velocity, Coord);
 		cout << "Type of Material and velocity: " << TypeOfMat << " " << velocity << endl; 

		SetMatType(tet1, objecttype, TypeOfMat, Coord);
		SetMat(tet1, objecttype, velocity, Coord, MatVec);
		//SetMat2(tet1, objecttype, velocity, Coord, dxFEM, dyFEM, dzFEM, MatVec);
		for (int j=0;j < 6; j++)
			cout << Coord(j) << " "; 
		cout << endl;
	}
	cout << "Max of Material " << Maxim(MatVec) << endl;
 	WavesOutputs out(tet1, Markers);
 
	string FEMGridFileName = out1.FEM_grid_file_name(fname); 
	const char* FEMGridFile = FEMGridFileName.c_str();
	 
 	tet1.print(FEMGridFile);

	if (argc >2)
	{
		char* fname1 = argv[2]; 
		opt.writeInp3D(fname1, &tet1, MatVec,1);
	}
	if (argc > 3) // save the coefficient:
	{
		char* fname2 = argv[3];		
		ofstream output;
		output.open(fname2);
		if (output.is_open())
		{
			for (int i = 0; i < MatVec.size(); i++)
				output << MatVec(i) << " "; 
			output.close();
		}
	}
 
	// extract the boundary nodes:
	double x, y, z, x_maxFEM, y_maxFEM, z_maxFEM; 
	int ii, bd_idx;
	bd_idx = 0;  

	x_maxFEM = x_minFEM + (NxFEM-1)*dxFEM; 
	y_maxFEM = y_minFEM + (NyFEM-1)*dyFEM;
	z_maxFEM = z_minFEM + (NzFEM-1)*dzFEM;

	for (ii = 0; ii < tet1.getNoNodes(); ii++)
	{
		x = tet1.getCoor(ii,0);
		y = tet1.getCoor(ii,1);
		z = tet1.getCoor(ii,2);
		if ((x - x_minFEM)*(x - x_maxFEM)*(y - y_minFEM)*(y - y_maxFEM)*(z-z_minFEM)*(z - z_maxFEM) == 0)
		{
			bd_idx +=1;
		}	
	}

	cout << "Number of boundary nodes: " << bd_idx << endl;

	// save the boundary grid file:
	string bdgridFile = out.boundary_grid_file_name(argv[1]);
	const char* BoundaryGridFile = bdgridFile.c_str();

	FILE* fp = fopen(BoundaryGridFile,"w");

	fprintf(fp, "%i \n", bd_idx);

	for (ii = 0; ii < tet1.getNoNodes(); ii++)
	{
		x = tet1.getCoor(ii,0);
		y = tet1.getCoor(ii,1);
		z = tet1.getCoor(ii,2);
		if ( (x - x_minFEM)*(x - x_maxFEM)*(y - y_minFEM)*(y - y_maxFEM)*(z - z_minFEM)*(z - z_maxFEM) == 0 )
			fprintf(fp, "%i %f %f %f\n", ii+1, x, y, z);
	}
	fclose(fp);
	cout << "Boundary grid written to the file: " << BoundaryGridFile << endl;

   }
   else 
	cout << "Not enough input " << endl;

   return 0;
}

