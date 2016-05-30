                                                                 // -*- C++ -*-
#ifndef WavesAddGrids_h_IS_INCLUDED
#define WavesAddGrids_h_IS_INCLUDED

#include <include/wavesGridB.h>

/*<AddGrids:*/
class WavesAddGrids
{
public:
  WavesAddGrids (real tol=0.05) {tolerance=tol;};
  ~WavesAddGrids () {};				
  bool addGrids( WavesGridB& grid1, WavesGridB& grid2, WavesGridB& gridadd );
  bool extractGrid( WavesGridB& g1, WavesGridB& g2, MV_Vector<int>& partition);
  bool createMirrorGrid( WavesGridB& grid, WavesGridB& newgrid, 
			 MV_Vector<real>& plane_);
  void symmetrizeAndRotateGrid (WavesGridB& grid, WavesGridB& outgrid, int rotations);
  void rotateGrid( WavesGridB& grid, real alpha);
  void setTolerance(const real tol){tolerance=tol;}
private:
  real tolerance;
};

/*
NAME:  WavesAddGrids - various grid modification algorithms

KEYWORDS:

  grid generation, symmetry, reflection, rotation


DESCRIPTION:

  The class enables construction of grid by adding two grids
  Construction of symmetric and rotational symmetric grids.


CONSTRUCTORS AND INITIALIZATION:

  There is one constructor without arguments.
  No initialization is needed.

MEMBER FUNCTIONS:

  "addGrids" - adding two grids together. It is assumed that the two grids 
  are disjoint except for common nodes at the boundary.
  The numbering of nodes from grid 1 are kept. The nodes in grid 2 
  which are common with grid 1 are given the grid 1 nodal numbers. 
  The remaining nodes are given increasing node numbers starting where 
  the grid 1 numbering stopped.
  The element numbering of grid 1 will be kept in the added grid,
  and the elements which was in grid 2 gets the following element
  numbers in the same order.
  The indicator names from grid 1 will overwrite the names of grid 2.
  If an indicator is set in one of grid 1 and grid 2 then it will also
  be set for the nodes in the added grid. The existing version
  works only for homogeneous grids, i.e. all elements must
  be of the same type.

  "extractGrid" - given a grid, create a new grid which only 
  contains a portion of the elements of the input grid. The elements
  are selected by a vector containing the element numbers of the input grid
  selected elements. The element numbering of the new grid will be
  after the order the elements are given in the vector. The node
  numbering are given by the order which the required nodes are
  part of elements in the new grid.

  "createMirrorGrid" - construct a mirrored grid. The new grid will have
  twice as many elements as the input grid.
  Detect the nodes which lies on the plane, these nodes will be common with
  the original and the new mirrored part which is added to the grid.
  Also the indicators and the subdomain types of the reflected nodes and 
  elements will be copied.

  "symmetrizeAndRotateGrid" - Assumptions: A 2D or 3D grid is in a 
  sector with center at origin. One part of the boundary of the grid 
  has its nodes at the y=0 plane. The boundary at the sector must be 
  straight. The size of the angle of the sector is M_PI/rotations.
  If these assumptions are fulfilled then the outgrid will be constructed
  by first mirroring the initial grid in the y=0 line/plane. The resulting
  symmetric grid is then rotated and added to itself (rotations-1) times. 
  The output grid takes the form of an annulus or torus or 
  cylinder shape, depending on the shape of the input grid.

  "rotateGrid" - rotate a grid around origin,
  used by "symmetrizeAndRotateGrid".

EXAMPLES:

  Construction a grid mirrored in a plane.

  Suppose that a finite element grid has been constructed (of type "WavesGridB").
  Given a hyper plane, i.e., a point in 1D, a line in 2D and a plane in 3D,
  of the form 

  A*x + D = 0 in 1D.
  A*x + B*y + D = 0 in 2D.
  A*x + B*y + C*z + D = 0 in 3D.

  It is necessary that  A*A+B*B+C*C is non-zero.

 The coefficients in the representation of the hyper plane is put in a 
 "Ptv(real)" object defining the plane. 

| 
| #include <WavesGridB.h>
| #include <readOrMakeGrid.h>
| #include <AddGrids.h>
| 
| int main (int , const char**  )
| { 
|   Handle(WavesGridB) grid,mirrorgrid;; 
|   real xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1;
|   int xdiv = 4,ydiv = 4,zdiv = 4;
|   String preproinfo;
| 
|   AddGrids cmg;
| 
|   grid.rebind( new WavesGridB() );
|   mirrorgrid.rebind( new WavesGridB() );
|   preproinfo = aform("PREPROCESSOR=PreproBox");
|   preproinfo += aform("/d=2 [%g,%g]x[%g,%g]",xmin,xmax,ymin,ymax);
|   preproinfo+=aform("/d=2 elm_tp=ElmT3n2D div=[%d,%d] grad=[1,1]",xdiv,ydiv);
|   s_o<<"input:"<<preproinfo.chars()<<endl;
|   readOrMakeGrid (grid(),preproinfo.chars());
|   Ptv(real) plane2(3);
|   plane2(1) = 0.0;  plane2(2) = 1.0;  plane2(3) = 0.0;
|   // will mirror the grid on the line x*0 + y*1 + 0 = 0, i.e. y=0
|   cmg.AddGrids( grid(), mirrorgrid(), plane2);
|   grid().print("FILE=T3.grid");
|   mirrorgrid().print("FILE=T3.mirror.grid");
| }

*/


#endif


