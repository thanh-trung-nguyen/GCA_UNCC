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
// Last changed: 2012-05-15 by Larisa Beilina

 
/*
  Class for hierarchical representation of a finite element grid
*/

#ifndef WavesGridFEHier_h_IS_INCLUDED
#define WavesGridFEHier_h_IS_INCLUDED

#include "wavesGridB.h"
#include "wavesFEcomp.h" // computations for triangles and tetrahedrals
#include "wavesMVvtp.h"
#include "wavesMVvind.h"
#include "wavesMVmtp.h"

#include <stdio.h>
//----------------------------------------------------------------------
/*<WavesGridFEHier:*/

class WavesGridFEHier : public WavesGridB
{
private:
  // The number in row e, column j, refers to the element number on the other
  // side of node(1D)/edge(2D)/face(3D) number j of element e.
  MV_ColMat<int> elm_neigh;// contains neighbor elements of elements
  // new_nodes is in 2D and 3D used as temporary storage for new node numbers. 
  // In 2D row and column positions agree with the ordering in
  // elm_neigh. In 1D new_nodes will contain the number of new nodes in the
  // element. In 3D the column position refers to the local edge number 
  // where a new node will be inserted.
  MV_ColMat<int> new_nodes; // number for new nodes at marked edges
  MV_Vector<int> pos_child; // position for first child of element in next_grid
  MV_Vector<int> parent;    // element no of parent element in prev_grid

  WavesGridFEHier*  prev_grid; // parent grid,NULL if this is the coarsest grid
  WavesGridFEHier*  next_grid; // child grid, NULL if this is the finest grid
  
// some parameters for refinement alternatives
  
  bool regular_uniform_ref; // true for use of regular refinement method
  void    setRegularUniformRef ( const bool flag ){ regular_uniform_ref = flag; }
  bool getRegularUniformRef () const        { return regular_uniform_ref; }
  int     refinement_method; // refinement method, # of edges to be divided.
  void    setRefinementMethod(int method)  { refinement_method = method; }
  int     getRefinementMethod() const  { return refinement_method;       }
  int     reg_bisect_short;  // Used in 2D refinement if all edges   
                             // of a triangle is marked for refinement
  // 0 regular: refine triangle into 4 triangles of same shape
  // 1 bisection: 4 triangles by bisection at longest edge
  // 2 shortest: choose the alternative of 0) and 1) which gives shortest edge
//--------------------------------------------------
// set, get and increase element and node counters used in the node numbering
  
  int current_elm_number;      // counter for new element
  int current_node_number;     // counter for new nodes

  void setCurrElmNo(const int e) {     current_elm_number  = e;  }
  void setCurrNodeNo(const int n){     current_node_number = n;  }
  int  getCurrNodeNo()const { return current_node_number;   }
  int  getCurrElmNo()const  { return current_elm_number;    }
  int  incCurrElmNo()      { return ++current_elm_number;  }
  int  incCurrNodeNo()     { return ++current_node_number; }

  bool new_nodes_are_set;   // true if node numbers of the new nodes are set
  void   setNewNodesSet(const bool flag)   { new_nodes_are_set = flag;}

//--------------------------------------------------
// Information on neighbor elements and grid hierarchy, 
// accessed by public functions below.

  bool parent_info;         // information about parent grid is computed
  void   setParentInfo  ( bool flag )     { parent_info   = flag;}
  bool children_info;       // information about children grid is computed
  void   setChildrenInfo( bool flag )     { children_info = flag;}
  bool neighbor_computed;   // true if initNeighbor is initiated 
  void   setNeighborsComputed(bool flag)  { neighbor_computed = flag;}
  bool finest_grid;  // true if this grid is the finest grid in hierarchy
  void   setFinestGrid(bool flag)         { finest_grid = flag;  }
  bool coarsest_grid;// true if this grid is the coarsest grid in hierarchy
  void   setCoarsestGrid(bool flag)       { coarsest_grid = flag;}
  void   setParentGrid( WavesGridFEHier* grid )   { prev_grid= grid;     }
  void   setChildGrid ( WavesGridFEHier* grid )   { next_grid = grid;    }
  
  // set that child_elm in this grid has par_elm in parent grid as parent elm
  void putParent(int child_elm, int par_elm ){ parent(child_elm) = par_elm; }
  void putParent1(int child_elm, int par_elm ){ parent(child_elm-1)=par_elm-1;}

//-------------------------------------------------------------
  // functions used for numbering of new nodes
  
  // put a node on edge edge_no in element e and give it a global node number
  void numberingNewNodes( int e, int edge_no ); // used in 2D and 3D
  void extraNumberingForT10 (); // ElmT10n3D refinement,. for future use.

  // mark (and number) nodes on the number of edges given in mark_elms.
  // allows varying number of added nodes for elements.
  bool markElements( const MV_Vector<int>& mark_elms); 
  // mark (and number) for uniform refinement
  bool markAllElements ();
  // version where refinement of element is ON or OFF, number of edges to
  // be refined is given by the parameter refinement_method.
  //  bool markElementsBool( const VecSimple(bool)& mark_elms); 
  bool markElementsBool( const MV_Vector<bool>& mark_elms); 
  //  bool markElementsBool( int* mark_elms ); 
  // mark edge recursively in element e
  void    markEdgeRec( int e, int edge_no ); 
  // mark edges recursively in element e, use parameter refinement_method
  void    markAllEdgesRec( int e ); 
  // mark (number) the no_mark_edges longest edges recursively in element e
  void    markAllEdgesRec( int e, int no_mark_edges ); 

// True if edge ed1 is ordered before edge ed2, used if edge lengthes are equal
  bool test_equal_with_geometry;
  bool lessThan ( const int ed1, const int ed2, const int element);

//-------------------------------------------------------------
// After that position and numbers of new nodes are done construct the
// new grid. Return true if the refinement was successful.
  bool refine ( WavesGridFEHier* newgrid ); // chooses between refine1D, refine2D and refine3D
  bool refine1D ( WavesGridFEHier* newgrid ); // 1D refinement
  bool refine2D ( WavesGridFEHier* newgrid ); // 2D refinement
  bool refine3D ( WavesGridFEHier* newgrid ); // 3D refinement

// initialization functions where some commonalities are put:

  void initRefine ();   // initialization in public refine versions
  // initialization and redim of new grid in refine1D, refine2D and refine3D
  void initRefinedGrid ( WavesGridFEHier* newgrid, int new_nno, int new_nel );
  // common initialization and redim of elm_neigh in initNeighbor1D-3D
  void initNumbering (); 
  void initGridFEHier ();   // initialization in constructor
  void checkReadyForRefine ();//check that numbering and neighborcomp. are done
  // -------------------------------------------------------------

// functions used in numbering, refinement and neighbor computations:
  
  // return 1, 2 or 3 depending on length of edges between nodes 5-10, 6-8, 7-9
  int  getShortestDiagonal( int e, MV_ColMat<real>& coornew); // 3D
  //  int  getShortestDiagonal( int e, float* coornew); // 3D
  // get the local face number of element e which have nodes n1,n2 and n3.
  int  getElmAtFace(const int n1, const int n2, const int n3, const int e);//3D
  // get the local face number of element e which have midpoint nodes n1 and n2
  int  getElmAtFace(const int n1, const int n2, const int e);// for T10
  // get the two elm numbers of neighbors which share face and edge n1 and  n2
  // with element e, where n3 and n4 are the two other nodes in the element.
  void getElmsAtFaces(const int n1, const int n2, const int n3, const int n4,
		      const int e, int& e3, int& e4); // 3D
  // put the list of all the elements and local edge numbers (including e and
  // edge_no) on stacks, elmstack and ednostack, increase stackcounter.
  void getElmsAtEdge(int e, int edge_no); // 3D, put result in elmstack

// functions for setting grid data,
  // set that element nb is neighbor element on the other side of 
  // node/edge/face (in 1D,2D and 3D) ed_no of element e.
  void putElmNeighbor(int e, int ed_no, int nb){ elm_neigh(e,ed_no) = nb; }
  void putElmNeighbor1(int e, int ed_no, int nb){ elm_neigh(e-1,ed_no-1) = nb;}
  int getNewNodeNo(int e, int ed_no) { return new_nodes(e,ed_no);}
  int getNewNodeNo1(int e, int ed_no) { return new_nodes(e-1,ed_no-1);}
  /*
  // some "putLoc2Glob-macros" to save space in refinement code
  void putLoc2Glob12 (int glob_node1, int glob_node2, int e){
    nodpek[0][e-1] = glob_node1; nodpek[1][e-1] = glob_node2;}
  void putLoc2Glob123 (int glob_node1, int glob_node2, int glob_node3, int e){
    nodpek[0][e-1] = glob_node1; nodpek[1][e-1] = glob_node2; 
    nodpek[2][e-1] = glob_node3;}
  void putLoc2Glob456 (int glob_node4, int glob_node5, int glob_node6, int e){
    nodpek[3][e-1] = glob_node4; nodpek[4][e-1] = glob_node5; nodpek[5][e-1] = glob_node6;}
  void putLoc2Glob1234 (const int glob_node1, const int glob_node2,
			const int glob_node3, const int glob_node4,const int e)
  { nodpek[0][e-1] = glob_node1;    nodpek[1][e-1] = glob_node2;
    nodpek[2][e-1] = glob_node3;    nodpek[3][e-1] = glob_node4;  }
  void putLoc2Glob5678910 (const int glob_node5, const int glob_node6,
			   const int glob_node7, const int glob_node8,
			   const int glob_node9, const int glob_node10,
			   const int e)
  {nodpek[4][e-1] = glob_node5; nodpek[5][e-1] = glob_node6; 
  nodpek[6][e-1] = glob_node7;  nodpek[7][e-1] = glob_node8; 
  nodpek[8][e-1] = glob_node9; nodpek[9][e-1]= glob_node10;}
  */
  
  // some "putLoc2Glob-macros" to save space in refinement code
  void putLoc2Glob12 (int glob_node1, int glob_node2, int e){
    nodpek(e,0) = glob_node1; nodpek(e,1) = glob_node2;}
  void putLoc2Glob123 (int glob_node1, int glob_node2, int glob_node3, int e){
    nodpek(e,0) = glob_node1; nodpek(e,1) = glob_node2; nodpek(e,2) = glob_node3;}
  void putLoc2Glob456 (int glob_node4, int glob_node5, int glob_node6, int e){
    nodpek(e,3) = glob_node4; nodpek(e,4) = glob_node5; nodpek(e,5) = glob_node6;}

  void putLoc2Glob1234 (const int glob_node1, const int glob_node2,
			const int glob_node3, const int glob_node4,const int e);
  //  { nodpek(e,0) = glob_node1;    nodpek(e,1) = glob_node2;
  //    nodpek(e,2) = glob_node3;    nodpek(e,3) = glob_node4;  }

  void putLoc2Glob5678910 (const int glob_node5, const int glob_node6,
			   const int glob_node7, const int glob_node8,
			   const int glob_node9, const int glob_node10,
			   const int e)
  {nodpek(e,4) = glob_node5; nodpek(e,5) = glob_node6; nodpek(e,6) = glob_node7;
   nodpek(e,7) = glob_node8; nodpek(e,8) = glob_node9; nodpek(e,9)= glob_node10;}

// functions for checking which are longest edge on face 2 and 3 of tetrahedron
// (order has been computed in advance)

  int firstOf256 ( MV_Vector<int>& order );//find which no comes first in order
  int firstOf346 ( MV_Vector<int>& order );// when 1 is first.
  //  int firstOf256 ( int* order );//find which no comes first in order
  //  int firstOf346 ( int* order );// when 1 is first.
  // get first of ed2, ed5, and ed6 which after reordering will be edges 2,5,6
  int firstOf256 ( MV_Vector<int>& order,
		   const int ed2, const int ed5, const int ed6); 
  // get first of ed3, ed4, and ed6 which after reordering will be edges 3,4,6
  int firstOf346 ( MV_Vector<int>& order,
		   const int ed3, const int ed4, const int ed6); 
  void initNeighbor1D(); // computation of neighbor data in 1D
  // computation of neighbor data + longest edge first
  void initNeighbor2D( bool alwaysSortEdgesIn2D = false );
  void initNeighbor3D(); // computation of neighbor data + longest edge data

// common initialization in the changeGrid routines
  void initL2QorQ2L( WavesGridFEHier* newgrid, const int new_nno );
// called from changeGridL2Q or changeGridQ2L
//  WavesGridFEHier* changeGridB2toB3(); // Change from ElmB2n1D  to ElmB3n1D  elms
//  WavesGridFEHier* changeGridB3toB2(); // Change from ElmB3n1D  to ElmB2n1D  elms
//  WavesGridFEHier* changeGridT3toT6(); // Change from ElmT3n2D  to ElmT6n2D  elms
//  WavesGridFEHier* changeGridT6toT3(); // Change from ElmT6n2D  to ElmT3n2D  elms
//  WavesGridFEHier* changeGridT4toT10();// Change from ElmT4n3D  to ElmT10n3D elms
//  WavesGridFEHier* changeGridT10toT4();// Change from ElmT10n3D to ElmT4n3D  elms
  
  WavesGridFEHier* changeGridB2toB3( WavesGridFEHier* newgrid ); // Change from ElmB2n1D  to ElmB3n1D  elms
  WavesGridFEHier* changeGridB3toB2( WavesGridFEHier* newgrid ); // Change from ElmB3n1D  to ElmB2n1D  elms
  WavesGridFEHier* changeGridT3toT6( WavesGridFEHier* newgrid ); // Change from ElmT3n2D  to ElmT6n2D  elms
  

  WavesGridFEHier* changeGridT6toT3( WavesGridFEHier* newgrid ); // Change from ElmT6n2D  to ElmT3n2D  elms
  
  WavesGridFEHier* changeGridT4toT10( WavesGridFEHier* newgrid );// Change from ElmT4n3D  to ElmT10n3D elms
  
  WavesGridFEHier* changeGridT10toT4( WavesGridFEHier* newgrid );// Change from ElmT10n3D to ElmT4n3D  elms
  
  // some scratch variables used in the refinement process
  //  VecSimple(int) edge_face2, edge_face3; // after initialization
  MV_Vector<int> edge_face2;
  MV_Vector<int> edge_face3;
  // edgeface2 contains the longest edge number in face 2 which of edges 2,5,6
  // edgeface3 contains the longest edge number in face 3 which of edges 3,4,6
  // (edges are sorted so that edge number 1 is longest in face 1 and 4)
  //  VecSimple(real) lengthsscratch; // used in getEdgeLengthOrder
  MV_Vector<real> lengthsscratch;    // used in getEdgeLengthOrder
  MV_Vector<int> orderscratch;// used in calling getEdgeLengthOrder
  int stackcounter;         // size of stacks elmstack and ednostack
  // where neighbors are put in the recursive refinement of T4 elements
  //  VecSimple(int) elmstack;  // used in the node numbering of T4 elements
  //  VecSimple(int) ednostack; // used in the node numbering of T4 elements 
  MV_Vector<int> elmstack;
  MV_Vector<int> ednostack;

public:
  
  WavesGridFEHier ();
  WavesGridFEHier ( const WavesGridB& grid );     // make a WavesGridFEHier from a GridFE
  WavesGridFEHier ( const WavesGridFEHier& grid ); // construct a copy
  
  ~WavesGridFEHier();
  
  bool ok () const;
  void operator = ( const WavesGridFEHier& grid );
  void operator = ( const WavesGridB& grid );
  void print_gid_mesh_FEM( WavesGridFEHier& gg,char *file);
// the main refinement methods

  // refinement with fixed number ref_method of subdivided edges in 
  // refined element
  //  WavesGridFEHier* refine(const VecSimple(bool)& mark_elms,const int ref_method);
  //  WavesGridFEHier* refine(int* mark_elms,const int ref_method);
  WavesGridFEHier* refine( WavesGridFEHier* newgrid, const MV_Vector<bool>& mark_elms,const int ref_method);
  // refinement where the number of subdivided edges in refined element
  // may vary (mixed method)
  //  WavesGridFEHier* refine(const VecSimple(int)& mark_elms);
  WavesGridFEHier* refine( WavesGridFEHier* newgrid, const MV_Vector<int>& mark_elms);
  WavesGridFEHier* refine0( WavesGridFEHier* newgrid, MV_Vector<int>& mark_elms);
  // WavesGridFEHier* refine(int* mark_elms);

  // Refine all elements. Regular refinement is used if regref is true 
  // and ref_method is 3 in 2D or 6 in 3D, respectively. 
  // To refine all edges set ref_method=3 in 2D, and ref_method=3D.
  WavesGridFEHier* refineUniformly( WavesGridFEHier* newgrid, const int ref_method=3, bool regref=false );

// refine all edges with length longer than break_ratio*|longest edge in grid|.
  WavesGridFEHier* refineLongEdges( WavesGridFEHier* newgrid, const real break_ratio );

// make a new grid from with linear or quadratic elements and vice versa,
// the new grid will be connected to the position in the hierarchy which 
// the old grid had. The old grid will still exist (if addressed by a handle)
// but it will be disconnected from the hierarchy.

  WavesGridFEHier* changeGridL2Q( WavesGridFEHier* newgrid ); // linear to quadratic grid
  WavesGridFEHier* changeGridQ2L( WavesGridFEHier* newgrid ); // quadratic to linear grid

  //k void makeHigherOrderGrid( WavesGridFEHier* newgrid, const int order );
  //k void makeHigherOrderGrid( WavesGridB* newgrid, 
  //k		    const MV_Vector<int>& elm_order,
  //k		    const bool KEEP_LINEAR_NODES_CONTINUOUS=true);
  //k  void checkCurvedElementBoundaries(MV_Vector<bool>& elm_is_curved);
  void    setRegBisectOrShort(int method)  { reg_bisect_short = method; }
  int     getRegBisectOrShort()     { return reg_bisect_short;          }

// members of utility class for geometry and finite element computations
// for triangles (2D) and tetrahedras (3D), normally at most one of 
// FET3 and FET4 is used depending on dimension of the grid.
  //????
  FET3n2D FET3; // triangle (2D) computations
  FET4n3D FET4; // tetrahedral (3D) computations 

// make a copy of whole or part of the grid hierarchy,
// return a pointer to the finest grid of the copied hierarchy

  WavesGridFEHier* makeCopyOfHierarchy( MV_Vector<WavesGridFEHier>& copiedgrids, 
				   const int firstlevel, 
				   const int lastlevel);

// In 2D and 3D initNeighbor reorder the nodes in the element so that the 
// longest edge will be first. It also computes the numbers of the element
// neighbors. To avoid unnecessary computation check with "isNeighborsComputed"
// to see if neighbor information already is computed. If the parameter
// "alwaysSortEdgesIn2D" is true then in 2D sorting of triangle edges are
// always forced, if not sorting is only done for the coarsest grid.
// (normally the triangle edges are in correct order after 2D refinement).
// For other dimensions alwaysSortEdgesIn2D plays no role.
  void initNeighbor( bool alwaysSortEdgesIn2D = false ); 
  
// methods for removal of data and disconnection of grids from hierarchy
// note that the disconnected grids will be deleted unless no handle 
// addresses the grid
  void removeNeighborInfo (); // rem. neighbor data: elm_neigh and new_nodes
  void removeNewNodeInfo ();  // removal of the new_nodes data
  void removeHierInfo ();     // rem. of children, parent, prev_ and next_grid
  void removeChildrenInfo (); // removal of children and next_grid
  void removeChildGrid ();    // removal of child grid and its info
  void removeParentInfo ();   // removal of parent, prev_gri
  void removeParentGrid ();   // removal of parent grid and its info
  
// check value of parameters
  bool isParentInfo() const       { return   parent_info;     }
  bool isChildrenInfo() const     { return children_info;     }
  bool isNeighborsComputed() const{ return neighbor_computed; }
  bool isNewNodesSet() const      { return new_nodes_are_set; }
  bool isThisFinestGrid() const   { return finest_grid;       }
  bool isThisCoarsestGrid() const { return coarsest_grid;     }

  // there are two different methods in determining which of 
  // two equally long edges should be ordered first:
  // Either by geometry or by the global node number.
  void setEqualTest2Geometry ( const bool flag){ test_equal_with_geometry = flag;}
  bool isEqualTestGeometry () const { return test_equal_with_geometry;}
  
// check with isThisFinestGrid and isThisCoarsestGrid if grid exists
  WavesGridFEHier* getParentGrid() const { return prev_grid; }// get coarser grid
  WavesGridFEHier* getChildGrid() const { return next_grid; }// get finer grid
  
  int getGridLevelNo()const;     // get level number of this grid in the hierarchy, 
                                // where the coarsest level has level number 1.
  int getFinestGridLevelNo()const;   // level number of the finest grid level
  const WavesGridFEHier* getFinestGrid()const;//get pointer to finest grid in hierarchy
  const WavesGridFEHier* getCoarsestGrid()const;//get pointer to coarsest grid in hierarchy
// Initially the grid itself is both the finest and coarsest grid.

  // Return 0 if there is no parent grid or if parent info is not computed
  int  getParent(const int e)const; 
  // get parent element on a rel_level coarser grid, where rel_level=1
  // will give identical result as getParent(int).
  int  getParent(const int rel_level, const int e) const;
  // getParent_eff: same as the getParent function but more efficient, no check
  // is performed if parent grid exists and if parent info has been computed
  int  getParent_eff(const int e)const{ return parent( e ); } //efficient, no check
  int  getParent_eff( const int rel_level, const int e)const; // efficient, no check

  int  getParent_eff1( const int e ) const { return parent( e )+1; } 
  // efficient, no check

  // Get the number of element (children) an element is divided into
  // in the child grid (next_grid). Minimum no is 1 (not refined element)
  // normally max is 4 children in 2D and 8 children in 3D.
  int  getNoChildren(const int e)const { return pos_child(e+1)-pos_child(e); }
  int  getNoChildren1(const int e)const { return pos_child(e)-pos_child(e-1); }
//c  int  getNoChildren(int e){ return number_of_children(e); }
  int  getChild(const int e, const int child_no)const ; 
// return 0 if no child info is computed
  // Same as getChild but no checking whether a finer grid exists, children
  // info is computed, and child_no is wrong.
  int  getChild_eff(const int e, const int child_no) const 
  {return pos_child( e ) + child_no;}
  // get the element number of the first child element
  int  getChild_eff1( const int e, const int child_no) const 
  {return pos_child( e-1 ) + child_no;}
  // get the element number of the first child element
  int  getFirstChild(const int e) const { return pos_child( e ); }
  int  getFirstChild1(const int e) const { return pos_child( e-1 )+1; }

  bool checkParentInfo() const;   // check that  parent  info is consistent
  bool checkChildrenInfo() const; // check that children info is consistent

  // computes the order of lengths of the local edges of a tetrahedron, 
  // where after call, order(1) will contain the local edge number of longest
  // edge, order(2) next longest, etc.
  //  void getEdgeLengthOrder( VecSimple(int)& order, int e ); // 3D
  void getEdgeLengthOrder( MV_Vector<int>& order, int e ); // 3D
  // get the local edge number on which two global node numbers n1 and n2.
  // The algorithm assumes that n1 and n2 are two of the four nodes
  // vertex nodes of tetrahedron e, but only the nonoptimized version
  // will detect the error that n1 and n2 not are nodes.
  //  In the optimized version a more efficient version is used which
  // do not detect the mistake.
  int  getLocEdgeNo( int n1, int n2, int e ); // 3D
  // get the global node numbers of the two vertex nodes lying on edge edge_no
  void getNodesOnEdge  ( int e, int edge_no, int& n1, int& n2 ); // 3D
  // get the local face number of the face which contain nodes n1,n2,n3.
  // If safe=true then a (slower) algorithm is used which detects whether any
  // of n1,n2,n3 in fact is not nodes of element e.
  int  getFaceNo( int n1, int n2, int n3, int e, bool safe=false ); // 3D
  // get the three global node numbers which are on face f of element e.
  void getGlobNodesOnFace( int& n1, int& n2, int& n3, int f, int e ); // 3D

// Computations of length of element edges:
  int getLongestEdgeNo(int e);  // edge number of longest edge in element e
  int getSecondLongestEdgeNo1( int e ); // middle longest edge, 1 longest
  int getShortestEdgeNo(int e);     // number of shortest edge in element e
  int getShortestEdgeNo1(int e); // number of shortest edge in elm e, 1 longest

// Get the element number on the neighbor element on the other side of 
// node (1D), edge (2D) or face (3D). Make sure that the neighbor info has
// been computed by using "isNeighborsComputed". Use "initNeighbor" to
// compute the neighbors.
  int  getElmNeighbor( const int e, const int ef_no )const 
  { return elm_neigh(e,ef_no);}
  int  getElmNeighbor1( const int e, const int ef_no )const 
  { return elm_neigh(e-1,ef_no-1);}
// Is element elnei a neighbor of element e? return true if so.
  bool isElmNeighbor( const int e, const int elnei )const;

// clear (set to OFF) indicators for nodes in the interior of domain.
  void removeInteriorIndicators();   
  bool checkElementOrientation() const; // check that areas/volumes are positive
  bool isElmNeighborConsistent() const; // true if neighbor info are consistent

// Mesh improvement (smoothing) based on criterion qualcrit.
// The quality criteria which determines the quantity which the
// smoothing algorithm tries to optimize are given in FEcomp.
// Only nodes without any set indicator can be moved.
// By calling getMinEta with the neighbor list of elements of a node 
// with different tentative new positions of the node, the position of the
// the node which gives best (largest value) of getMinEta is chosen.
  //k  void makeSmoothing (int nrepetitions=1, int qualcrit=1);
//k  void smoothNode (int node, int qualcrit=1); // smooth a single node
// smooth neighbor nodes of a node
//k  void smoothNeighbors (int node, int qualcrit=1); 
//k  void smoothNeighborsExt (int node, int qualcrit=1); 

// compute the smallest (worst) value of the quality criterion of the 
// elements in element_list.
  //  real getMinEta( VecSimplest(int)& element_list, int qualcrit=1);
  real getMinEta( MV_Vector<int>& element_list, int qualcrit=1);
  //  real getMinEta( int* element_list, int qualcrit=1);

  //  void makeStatistics (Os os); // compute some statistics of grid elements
  void makeStatistics (); // compute some statistics of grid elements
  
// Findin which element in level "level" the point x is contained in,
// return false if the point was not found inside any element.
// The search algorithm uses the tree structure of  parent and child grids
// and that the refinements are nested. For efficiency in a extensive use 
// of findElm call "computeBoundingBoxes" prior to use (but only once).

  /*
  bool findElm( Ptv(real)& x, int level, int& element );
  bool findElm( Ptv(real)& x, int& element) {
    return findElm (x, getGridLevelNo(), element );}

// The next version of findElm uses the following extra piece of information
// to faster find the correct element, (much faster if coarsest grid has
// many elements): The given point is inside element elm0 on level level0. 
// Return the element in level1 which the point is inside.
// If level1 < level0 it is more efficient to use getParent (or getParent_eff)
  int findElm( Ptv(real)& x, int level0, int elm0, int level1 );

// Warning: findElm assumes that nodes has not been moved from their placement
// in the refinement.

// Search through all elements in grid until the element which x is inside
// is found. This is the algorithm used on the coarsest grid.
// Return true if element was found.
  bool simpleSearch( Ptv(real)& x, int& element);
// Is the point x inside a given element?
  bool isInsideElm( Ptv(real)& x, int element);
  
// Precomputation of bounding boxes for elements makes the test whether a
// a point _could_ be inside an element faster.
// To obtain this efficiency call computeBoundingBoxes before findElm are used.
// If it is found that a point can be inside an element, barycentric
// coordinates are computed to find if the point is inside or outside.

private:
  MatSimple(real) bounding_box;
  bool bound_box_computed;
  void setBoundingBoxesComputed( bool flag ) { bound_box_computed = flag; }
public:
  void computeBoundingBoxes (); // compute bounding boxes for elements
  void removeBoundingBoxes ();  // remove the computed bounding boxes
  bool isBoundingBoxesComputed() { return bound_box_computed; }
  
// It is possible to compute statistics from the use of findElm.
// Call initSearchStat or initSearchStatAllGrids to set counters to zero.
// The number of used boxtests and local element coordinates are counted.

private:
  int box_counter, in_elm_counter;
public:
  void initSearchStat () { box_counter = in_elm_counter = 0; }
  void initSearchStatAllGrids (); // set counters to zero
  void printSearchStatAllGrids (Os os); // print statistics
  int getNoBoxTests () { return box_counter; } // box tests
  int getNoElmTests () { return in_elm_counter; } // comp. loc. coordinates
  */

  bool findElm( MV_Vector<real>& x, int level, int& element );
  bool findElm( MV_Vector<real>& x, int& element) {
    return findElm (x, getGridLevelNo(), element );}

// The next version of findElm uses the following extra piece of information
// to faster find the correct element, (much faster if coarsest grid has
// many elements): The given point is inside element elm0 on level level0. 
// Return the element in level1 which the point is inside.
// If level1 < level0 it is more efficient to use getParent (or getParent_eff)
  int findElm( MV_Vector<real>& x, int level0, int elm0, int level1 );

// Warning: findElm assumes that nodes has not been moved from their placement
// in the refinement.

// Search through all elements in grid until the element which x is inside
// is found. This is the algorithm used on the coarsest grid.
// Return true if element was found.
  bool simpleSearch( MV_Vector<real>& x, int& element);
// Is the point x inside a given element?
  bool isInsideElm( MV_Vector<real>& x, int element);
  
// Precomputation of bounding boxes for elements makes the test whether a
// a point _could_ be inside an element faster.
// To obtain this efficiency call computeBoundingBoxes before findElm are used.
// If it is found that a point can be inside an element, barycentric
// coordinates are computed to find if the point is inside or outside.

private:
MV_Vector<real> node_indicators;
real max_aniso;

void setNewNodeCoordinates(const int m1, 
			   const int m2,
			   const int mid, 
			   WavesGridFEHier* newgrid);

public:
void initNodeIndicatorsFromElmIndicators(MV_Vector<real>& elm_indicators,
						 real max_anisotropy_=0.5);
void removeNodeIndicators();


private:
  MV_ColMat<real> bounding_box;
  bool bound_box_computed;
  void setBoundingBoxesComputed( bool flag ) { bound_box_computed = flag; }
public:
  void computeBoundingBoxes (); // compute bounding boxes for elements
  void removeBoundingBoxes ();  // remove the computed bounding boxes
  bool isBoundingBoxesComputed() { return bound_box_computed; }
  
// It is possible to compute statistics from the use of findElm.
// Call initSearchStat or initSearchStatAllGrids to set counters to zero.
// The number of used boxtests and local element coordinates are counted.

private:
  int box_counter, in_elm_counter;
public:
  void initSearchStat () { box_counter = in_elm_counter = 0; }
  void initSearchStatAllGrids (); // set counters to zero
  //  void printSearchStatAllGrids (Os os); // print statistics
  void printSearchStatAllGrids (); // print statistics
  int getNoBoxTests () { return box_counter; } // box tests
  int getNoElmTests () { return in_elm_counter; } // comp. loc. coordinates

// print grid in plotmtv format
  //  void printGrid(Os os); //print the grid
  void printGrid(char* file); //print the grid
// print in plotmtv format the boundary and the interfaces between different 
// subdomaintypes
  //  void printBoundary(Os os);
  void printBoundary(char* file);
  //k void printBoundary(char* file, int crve_srf);
  //k void printParamgrid(char* file, int curvsurf);
  //k void printOneSurfaceGrid(char* file, int curvsurf);
  // crve_srf0 both crv and srf, 1 only surf, 2 only curve
// print only elements of subdomain type sbdtype
  //  void printGrid(Os os, const int sbdtype) const; 
  void printGrid(char* file, const int sbdtype); 

  // print number of nodes and elements of the grids in the hierarchy
  //  void presentGridHier (Os os);
  void presentGridHier ();

};

/*Class:WavesGridFEHier

NAME:  WavesGridFEHier - finite element grid

KEYWORDS:

  finite element, grid, mesh, GridFE, WavesGridFEHier, adaptive refinement,
  triangle, tetrahedron

DESCRIPTION:

  The class represents a collection of functions for refinement and 
  representation of a hierarchy of finite element grids.

  Class "WavesGridFEHier" is derived from GridFE.

CONSTRUCTORS AND INITIALIZATION:

 "WavesGridFEHier" can be constructed from a GridFE, if all the elements of the 
 grid are of the same type, where the allowed element types are 
 "ElmB2n2D", "ElmB3n2D", "ElmT3n2D", "ElmT6n2D", "ElmT4n3D" and "ElmT10n3D".
 Since "WavesGridFEHier" is derived from "GridFE" can the "GridFE" preprocessors 
 also be used directly to create an initial "WavesGridFEHier" grid.

BASIC FUNCTIONALITY:
 
 WavesGridFEHier is an extension of the GridFE class, for description
 of the members which are inherited see the documentation in GridFE.

 The refinement algorithms described below will construct refined grids which 
 will be nested in the original "parent" grid, i.e., each of the elements 
 in the refined "child" grid will lie inside an element of the 
 parent grid.

 When a grid is refined it will be connected to the old grid from
 which it was refined. This connection consists of pointers and information
 on which "parent" element in the new grid an element in the refined grid
 lies in. The old grid will for each element contain the element numbers 
 of the "children" elements in the refined grid into which the element 
 was refined. This parent-child information describes a tree structure
 which is useful e.g., for searching in the grid (see "findElm"). The methods
 "getParent", "getChild" and other methods are used to extract the tree 
 structure information.
 
 The sequence of refined grids will define a grid hierarchy where
 each of its components is a "WavesGridFEHier". In addition, by inheritance 
 from "GridFE" each of the grid can also be used as a "GridFE" grid in
 all the normal uses of a "GridFE" grid, such as assembly of systems, 
 defining fields, etc.  
 Each of the grids in the hierarchy is given a grid level number, 
 where the coarsest grid (without parent grid) has level 1, the level 
 is then increased for each finer grid. The grid level number is obtained by 
 the method "getGridLevelNo".

 In the refinement algorithms and for many other applications such as 
 finding outer or internal boundaries, computing error estimates, etc., 
 it is important to have information on which neighbor elements an element has.
 The neighbor element information is initialized by "initNeighbor",
 and contain the element number of the neighbor elements.
 The neighbor element is found by using "getElmNeighbor".
 In 1D the element on the left side is in column 1, and the element on the
 right side of the element is in column 2 of the array which "getElmNeighbor"
 access. In 2D, triangle edges have a local ordering where edge 1 goes between
 local nodes 1 and 2, edge 2 goes between local node 2 and 3, and edge 3 is
 between nodes 3 and 1. In 3D the following ordering of faces are used:
 Face 1 contains nodes 1, 2 and 4; face 2 contains nodes 2, 3 and 4;
 face 3 contains nodes 1, 3 and 4; face 4 contains nodes 1, 2 and 3.

 The local numbering of edges in a tetrahedron is as follows: 
 edge 1 is between node 1 and 2, edge 2 is between node 2 and 3,
 edge 3 is between node 3 and 1, edge 4 is between node 1 and 4,
 edge 5 is between node 2 and 4, edge 6 is between node 3 and 4.

 If an element does not have a neighbor at a node (in 1D), edge (2D) or
 face (3D), the corresponding column position in elm_neigh accessed by
 "getElmNeighbor" will contain the value 0.

----------------------------------------
 Refinement algorithms:

 The algorithm used to refine a grid consists of two phases.
 The first phase is to determine where new nodes should be inserted
 and which global node number these new nodes are going to get in the new grid.
 The second phase consists of using the node numbers introduced in phase one
 and refine the elements individually into their children elements.
 It is crucial that the first phase has succeeded in choosing the
 positions and numbers of the new nodes (which are shared by many elements) 
 so that the refined grid will be consistent (regular) and that the shape 
 of the elements in the refined grids will not degenerate after
 repeated refinement.

 For the linear elements "ElmT3n2D" (triangles) and "ElmT4n3D" (tetrahedra)
 the new nodes will be inserted on midpoints of edges in the grid. 
 The number of edges in each element to be refined is determined
 by the argument of the "refine" functions. The n longest edges
 in the element will be refined, where 0<=n<=3 in 2D and 0<=n<=6 in 3D.
 In addition an edge can be refined because the neighbor element wants
 to refine that edge. In addition, to ensure that the shape of the refined 
 elements will not degenerate by repeated refinement, the longest edge (in 2D)
 in a triangle will be refined, (if n was larger than 1 it was already
 refined). Correspondingly for 3D the longest edge of the element and 
 the longest edges of faces with any refined edges of tetrahedra are refined. 

 The refinement for "ElmT6n2D" elements are similar to the case 
 of "ElmT3n2D", but new nodes are inserted (and numbered), also between vertex 
 nodes and  midpoints nodes (which already was present from the old grid).
 New nodes will also be inserted on new interior edges of the elements.
 Refinement of quadratic element in 3D ("ElmT10n3D") is not yet implemented,
 but it is possible to use the function "changeL2Q" to change an "ElmT4n3D" 
 grid to an "ElmT10n3D" grid.

 Setting of indicators and subdomain types for the refined grid:

 All the nodes in the old grid will also be present in the new grid 
 (with the same numbering) and their old indicator values are copied.
 For the new inserted nodes in the refined grid an indicator 
 value will be set ON if both the nodes on the edge where the new node
 was inserted had that indicator value ON. Note that for a refined grid
 this means that nodes inserted on an edge, which goes between
 two nodes on the boundary, but not belonging to the boundary could get
 an indicator value set ON (if both the nodes on the boundary had the same
 indicator value set to ON). If this is not what the user wanted this can 
 be solved by either using different indicators on different part of the 
 boundary, or by calling "removeInteriorIndicators" after refinement, 
 (but this will also remove other, perhaps wanted, indicators in the interior).

 The element subdomain type of refined elements are given
 the same (integer) value as the parent element.
 
----

 Here are the details of the subdivision of elements (the second phase)
  of the refinement of  "ElmT3n2D" triangles. In phase one the nodes of 
 the element has been
 permuted so that edge 1 is the longest edge. The following cases can 
 occur: 0, 1, 2, or 3 of the edges has been subdivided by new nodes. 
 Let Ti(k,l,m) be a created refined element with local node numbers k,l and m.

 0 marked edges: Copy the old element which gives element T1=(1,2,3), see fig.
 _________________________________________________________
 1 marked edge: Divide the element into 2 elements, there is only one way to do
 this. Create elements T1=(3,1,m1) and T2=(2,3,m1), see figure below. 
 By the permutation of edges and the method of marking of new nodes at 
 the longest edges m1 will lie on edge 1.

 3______2
  |\ 2 /     m1 is the midpoint of edge 1 (the edge between node 1 and 2).
  |1\/m1
  |/ 
 1
 _________________________________________________________
 2 marked edges: (here the nodes are on edges 1 and 2, could also be 1 and 3).
 Divide the element into 3 new elements, there are 
 2 alternative ways to do this. We choose the alternative which will introduce
 the shortest edge of the two alternatives, 
 (rule of thumb: divide long edges and add short edges):
 Choose alternative 1 if edge 1-2 is longer than edge 2-3. Alternative 1
 is hence always chosen since edge 1 is the longest edge of the element.
 alt 1: T1=(3,1,m1), T2=(3,m1,m2), T3=(m1,2,m2),
 alt 2: T1=(3,1,m2), T2=(m2,1,m1), T3=(2,m2,m1),
 with the number convention below

 3__m2__2
  |    /
  |  /m1
  |/
 1
 _________________________________________________________
 3 marked edges: Divide into 4 new elements. There are 4 ways to do this, 
 (one regular and three bisections). Only 2 ways if bisection with
 the two shorter sides as base are excluded. 
 a) Regular refinement: T1=(1,m1,m3), T2=(m1,2,m2), T3=(m3,m2,3), T4=(m2,m3,m1)
 b) Bisection refinement with edge 1 as "base"
 T1=(1,m1,m3), T2=(m1,2,m2), T3=(m1,3,m3), T4=(3,m1,m2)
 The choice between regular or bisection refinement is determined by setting
 "setRegBisectOrShort". The choices for this parameters is 0) always regular 
 refinement, 1) always bisection refinement, or 2) the alternative which
 adds the shortest alternative edge to the new grid.

  3__m2__2
   |    /
 m3|  /m1
   |/
  1

 Observe that all new and old elements use counter clockwise orientation.

 ----
 Refinement in 3D:

 In 3D, bisection of elements are always used, except for uniform refinement
 of all edges, in which case regular refinement can be used. 
 By phase 1 of the refinement there is for every tetrahedron to be refined
 always a node dividing the longest edge of the tetrahedron.
 The tetrahedron is bisected into two tetrahedra by connecting the
 node on the longest edge to the opposing edge to the longest edge,
 (i.e. an edge between the midpoints of edge 1 and 6 with the local edge 
 numbering). The two new tetrahedra are then further subdivided if any of 
 the 5 other edges also have midpoints nodes to be inserted. The order 
 in which these further bisections are done will depend on the length of the 
 edges. For the refinement to result in a regular grid it is necessary that
 the edges subdividing a face shared by two tetrahedras are drawn consistently
 for both the involved tetrahedras.
   To use regular refinement in 3D it is necessary that all elements and all
 edges are refined. In regular refinement the 8 new tetrahedra are constructed 
 by first cutting of 4 tetrahedras similar to the parent tetrahedron
 from the parent tetrahedron, and making a regular
 refinement of the triangular face. The four remaining tetrahedra are
 constructed by adding an edge between the nodes of two opposing edges.
 There are three such possible tetrahedral interior edges. The 
 alternative which adds the shortest edge length is chosen.

-------
 Refinement in 1D:

 Phase 1 of the refinement is trivial since all the new nodes are
 interior to the elements and can hence be determined locally in the
 construction of the refined elements. Up to 4 nodes can be inserted
 into an element (a segment).

MEMBER FUNCTIONS:

  refine - Given a array with the length of the number of elements, 
    the vector contains information on which of the elements should be 
    refined. There are two versions of "refine". In the first version
    the vector contains bools, which are true if the corresponding 
    element should be refined.
    In the second version of "refine" the array contains integers, which
    values determines the (minimum) number of edges of the element to be 
    refined. Some additional edges can be subdivided by the rules to ensure
    stability of refinement, (the longest edge is always refined if any, etc).
    In the bool vector version, the number of subdivided edges in refined 
    elements is determined by the integer parameter "ref_method".
    In 1D "ref_method" means how many times the interval element is 
    sub divided in this refinement, the maximum allowed number is 4.
    For 2D "ref_method" is between 1 and 3, and for 3D "ref_method" is between
    1 and 6. The return of "refine" is a pointer to the refined grid.

  refineUniformly - refine all elements, the parameter "ref_method" is used as
    in "refine", if the parameter "regref" is true, then regular refinement
    will be used if in addition ref_method==6 for 3D grids or 
    ref_method==3 for 2D grids.

  refineLongEdges - refine all edges with length longer than 
    break_ratio*|longest edge|, where |longest edge| is the length of the
    longest edge in the grid. This refinement function will construct
    grids where all elements will be of approximately the same size.

  changeGridL2Q - from a grid with piece-wise linear elements construct
    a grid with piece-wise quadratic elements. Works in 1D, 2D and 3D.
  changeGridQ2L - The opposite of "changeGridL2Q".

  setRegBisectOrShort - Only 2D grids
  getRegBisectOrShort - Only 2D grids
      Set and get the parameter determining which refinement alternative 
    to be used in the case that all 3 edges are refined. The value 0 
    of the parameter means, "always regular refinement" (4 triangles of same 
    shape as original), value 1 means "always bisection", 
    (longest edge as base), and value 2 means that the alternative of regular 
    and bisection refinement which introduces the shortest new edge is used.

  ok - Simply calls GridFE::ok 

  makeCopyOfHierarchy - make a copy of a part of the gridhierarchy, the two
    parameters determines the levels in the old grid which will be
    the first and last levels in the copied grid hierarchy. A pointer
    to the finest grid of the copied hierarchy is returned.

  initNeighbor - In 2D and 3D "initNeighbor" reorder the nodes in the element 
    so that the longest edge is first. It also computes the numbers of the 
    element neighbors. In the refinement process "initNeighbor" is called 
    for the old grid, but "initNeighbor" not automatically called for the 
    new refined grid (until the new grid is refined), so if e.g. the
    method "getElmNeighbor" needs to be used before the grid
    has been refined then "initNeighbor" must be called.

-- some methods to remove information on neighbors, new_nodes structure,
   grid hierarchy and parent-child information.

  removeNeighborInfo - remove neighbor data, i.e., elm_neigh and new_nodes
  removeNewNodeInfo -  remove the new_nodes data
  removeHierInfo -     remove the coupling children, parent,prev_ and next_grid
  removeChildrenInfo - removal of pointer to next_grid and children info
  removeChildGrid -    removal of pointer to child grid and its information
  removeParentInfo -   removal of pointer to parent, prev_grid
  removeParentGrid -   removal of parent grid and its information
  
   The methods removeHierInfo, removeChildGrid and removeParentGrid
   do remove handle/pointer to the other grids in the hierarchy, 
   this may lead to that those grids are deleted depending whether 
   other handles are referring to that grid.

-- methods for getting status of grid hierarchy information

  isParentInfo - is information about the parent grid computed so that 
    getParent can be used?
  isChildrenInfo - is information about the child grid computed so that 
    getChild can be used?
  isNeighborsComputed - is element-element info computed?
  isNewNodesSet - is the array new_nodes filled?
  isThisFinestGrid - return true if the grid has no child grid.
  isThisCoarsestGrid - return true if the grid has no parent grid.
  
  getParentGrid - get pointer to the parent grid (a coarser grid)
  getChildGrid -  get pointer to the child grid (a finer grid)
    return NULL if no grid or safer: check with "isThisCoarsestGrid"/
    "isThisFinestGrid" before calling "getParentGrid" or "getChildGrid".
  
  getGridLevelNo - return the level number in the hierarchy of this grid, 
    where the coarsest level (without parent grid) has level number 1 and
    the level number for each added child grid to the hierarchy is increased 
    by 1.
  getFinestGridLevelNo - return the level number of the finest grid level
  getFinestGrid - get a pointer to the finest grid in the hierarchy
  getCoarsestGrid - get a pointer to the coarsest grid in the hierarchy
 
  getParent - parent(e) returns the element number of the parent element of e,
    0 is returned if no parent grid.
  getParent - with two arguments "rel_level" returns the parent element 
    rel_level levels up, (rel_level==0 returns the input element number)
  getParent_eff - more efficient version, can be used if it has been checked
    that the input element number is legal, that the grid has parent grid 
    and that the parent information is computed.
  getNoChildren - returns the number of elements this element has been 
    divided into in the child grid.
  getChild -  return the element number of child number "child_no" 
    of element "e".
  getChild_eff - more efficient version of getChild (skips some checking).
  getFirstChild - return the element number of the first child of an element.

  checkParentInfo - check that the parent info is consistent, i.e., that
    an element which has parent element p is also a child of p.
  checkChildrenInfo - check that children info is consistent, i.e. that
    the child of an element has this element as its parent.

  getElmNeighbor - get the element number of the neighbor on local 
    node/edge/face in 1D/2D/3D, respectively. The value 0 is returned if no 
    neighbor, i.e. the element is on the outer boundary.
  isElmNeighbor - does element "e" have element "elnei" as one of its
     neighbors?

-- some functions used in the refinement, but which might also be of more 
     general usage.

    -- functions only applicable in 3D:
  getEdgeLengthOrder - compute the order of the lengths of the edges.
    If the integer sequence (5,2,1,6,4,3) is returned this means that
    edge 5 was longest, then edge 2 and so on. If the lengths are equal
    the priority is uniquely decided either by the midpoint of the edge
    or by the node numbers . See the function "lessThan".
  getLocEdgeNo - return the local edge number of the edge with global node 
    numbers n1 and n2
  getNodesOnEdge - get the two global node numbers of the nodes which are
    on the local edge "edge_no" in element "e"
  getFaceNo - return the (local) face number (1-4) which contains global
    nodes n1, n2, and n3. The bool parameter "safe" chooses between two 
    different  algorithms. For the faster algorithm with safe=false then 
    wrong face is returned IF n1, n2, n3, are not nodes of the element.
  getGlobNodesOnFace - get the global node numbers of the nodes which are on 
    the local face f in element e.

-- functions which can be used both in 2D and 3D
  getLongestEdgeNo - return the local edge number of the longest edge
  getSecondLongestEdgeNo1 - return the second longest edge number, where in 2D 
    it is assumed that edge number 1 is longest, (no assumption in 3D).
  getShortestEdgeNo - get the edge number of the shortest edge in element e
  getShortestEdgeNo1 - get the edge number of the shortest edge in element e,
     where in 2D it is assumed that edge number 1 is longest, 
     (no assumption in 3D).

  removeInteriorIndicators - remove all the indicators which
    are not on the boundary of the domain, useful if refinement has caused
    unwanted indicators. (A node inserted at the midpoint of an edge
    will get an indicator set to ON if both the endpoints nodes 
    have that indicator set to ON.

  checkElementOrientation - computes area/volume of elements and checks that 
    area/volume are positive as they should. (Otherwise the element can be
    degenerate or the nodes are in wrong order.
  isElmNeighborConsistent - return true if the element neighbor information
    is self consistent. (a neighbor to an element must have the element as
    its neighbor)

  makeStatistics - compute some statistics about the elements of the grid

-- In smoothing of a grid the parameter "qualcrit" refers to the quality 
-- criteria used in FEcomp. 

  makeSmoothing: Mesh improvement according to the criterion "qualcrit".
    Move nodes without indicators so that the form of the elements are improved
  smoothNode - move one node
  smoothNeighbors - smooth the neighbor nodes of a node
  getMinEta - return the smallest value of all the elements in the list of
    elements 
 
-- methods for finding in which element a point is inside
-- general: safety region, box test, or compute the barycentric coordinates

  findElm - given a point "x", find in which element the point is in level 
    "level" without a priori info on where the node is.
    If no level is given the search is performed on this grid.

  findElm - given a point "x" which is in element "elm0" on level "level0"
    return in which element in level "level1" the point is inside.

  simpleSearch - return the element number a given point is inside
  isInsideElm - return true if the node is inside a given element
  
-- statistics for search

  box_counter, in_elm_counter - counters
  initSearchStat - initialization for search statistics set counters to zero
  initSearchStatAllGrids - initialization
  printSearchStatAllGrids - print the result of the search
  getNoBoxTests - return number of boxtests
  getNoElmTests - return number of barycentric tests

-- bounding boxes for triangles to make the search faster
 
  bounding_box - array containing coordinates of the element bounding boxes
  bound_box_computed - boolean telling if bounding boxes are computed
  removeBoundingBoxes - free the space occupied by the bounding boxes
  setBoundingBoxesComputed - set flag
  isBoundingBoxesComputed - check flag
  computeBoundingBoxes - computation of bounding boxes 

  printGrid - print a picture of the grid in .mtv-format, 1D, 2D, 3D
  printGrid(sbdtype) - mtv plot of elements with subdomain sbdtype, 2D and 3D
  printBoundary - mtv plot of the boundary of the grid, 2D and 3D
  presentGridHier - write info on the grids in the hierarchy

EXAMPLES:

  #include <WavesGridFEHier.h>
  #include <readOrMakeGrid.h>
  int main (int nargs, const char** args)
  { 
    initDIFFPACK (nargs, args);  

    Handle(WavesGridFEHier) grid; 
    grid.rebind( new WavesGridFEHier() );
    // constructing a grid with preprocessor
    real xmin = 0.0,  xmax = 1.0,  ymin = 0.0,  ymax = 1.0;
    int  xdiv = 4,  ydiv = 4;

    String preproinfo = aform("PREPROCESSOR=PreproBox");
    preproinfo+=aform("/d=2 [%g,%g]x[%g,%g]",xmin,xmax,ymin,ymax);
    preproinfo+=aform("/d=2 elm_tp=ElmT3n2D div=[%d,%d] grad=[1,1]",xdiv,ydiv);
  
    readOrMakeGrid (grid(),preproinfo.chars());// create an initial grid
    grid().presentGridHier(s_o);             // present the grid
    grid.rebind(grid().refineUniformly());   // uniform refinement
    grid().presentGridHier(s_o);             // present the grid hierarchy
    grid.rebind(grid().refineLongEdges(.5)); // refine the 50% longest edges
    grid().presentGridHier(s_o);             // present the grid hierarchy
  }

*/

#endif
