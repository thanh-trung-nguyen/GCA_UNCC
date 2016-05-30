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
      

//  Class for representation of a hierarchical finite element grid.
//  Code implementation of member functions.



#include "include/wavesGridFEHier.h"
//-----------------------------------------------------------------------------
WavesGridFEHier:: WavesGridFEHier ()
//-----------------------------------------------------------------------------
{
// initialize the WavesGridB data
  //  nsd = nno = nel = nne = maxnne = nbind = 0;
  nsd = nno = nel = maxnne = nbind = 0;
  //  elm_tp = "Dummy element";
  onemat = false;
  //  uniform_mesh = false;
  //  same_elm_tp = true;
  current_node = 0;
  //  associated_bfg = NULL;

  initGridFEHier();  // initialize the GridFEHier specific part
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initGridFEHier ()
//-----------------------------------------------------------------------------
{
  next_grid = NULL;
  max_aniso=0.5;
  neighbor_computed = new_nodes_are_set = false;
  finest_grid = coarsest_grid =           true; 
  parent_info = children_info =           false;
  regular_uniform_ref =                   false;
  refinement_method = 1;
  reg_bisect_short = 2; // choose shortest edge alternative, used in 2D
  test_equal_with_geometry = true;
  //  test_equal_with_geometry = false;
  current_elm_number = current_node_number = 0;
  bound_box_computed = false;
  FET3.attach(this); // at most one of FET3 and FET4 are going to be used
  FET4.attach(this); // for the same grid
  //k  getSurfaceTriangulation().attachGrid(this);

  lengthsscratch.newsize(6);
  orderscratch.newsize(6);
  initSearchStat ();
  //  lengthsscratch.redim(6); //global arrays used in getEdgeLengthOrder in 3D
  //  orderscratch.redim(6);
}

//-----------------------------------------------------------------------------
WavesGridFEHier:: WavesGridFEHier (const WavesGridB& grid)
//-----------------------------------------------------------------------------
{
  if (!grid.ok())
    //    fatalerrorFP("WavesGridFEHier:: WavesGridFEHier", "Empty grid is given.");
    cout<<" fatalerrorFP WavesGridFEHier:: WavesGridFEHier Empty grid is given\n"<<flush;

  if (this == &grid) 
    return;

  // purify (operator= calls redim which tests this->nsd==nsd_ etc., hence
  // these variables must be initialized!)
  //  nsd = nno = nel = nne = maxnne = nbind = 0;
  nsd = nno = nel = maxnne = nbind = 0;

  *this = grid;
}

//-----------------------------------------------------------------------------
WavesGridFEHier:: WavesGridFEHier (const WavesGridFEHier& grid)
//-----------------------------------------------------------------------------
{
  if (!grid.ok())
    //    fatalerrorFP("WavesGridFEHier:: WavesGridFEHier", "Empty grid is given.");
    cout<<" fatalerrorFP WavesGridFEHier:: WavesGridFEHier Empty grid is given\n"<<flush;

  if (this == &grid)
    return;

  //  nsd = nno = nel = nne = maxnne = nbind = 0;
  nsd = nno = nel = maxnne = nbind = 0;

  *this = grid; 

}

//-----------------------------------------------------------------------------
WavesGridFEHier:: ~WavesGridFEHier()
//-----------------------------------------------------------------------------
{ 
//  if (verbose_mode >= 0)
   printf("WavesGridFEHier:: ~WavesGridFEHier Deleted an object with %d nodes and %d elements\n",nno,nel);
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: ok () const
//-----------------------------------------------------------------------------
{
// Check the GridFE parts ...

  bool test = WavesGridB::ok();

// continue to check the WavesGridFEHier specific parts...

// consistency of the bool flags, hierarchy etc..

  return test;  // if no errors are found, this statement will be reached
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: operator = (const WavesGridB& grid)
//-----------------------------------------------------------------------------
{
  if (!grid.ok())
    //    fatalerrorFP("WavesGridFEHier:: operator =","Empty grid is given as argument.");
    printf(" fatalerrorFP WavesGridFEHier:: operator = Empty grid is given\n");

  /*
  String etyp = grid.getElmType(1);
  if( grid.allElmsOfSameType() == false &&
      "ElmT4n3D" != etyp && "ElmT10n3D" != etyp &&
      "ElmT3n2D" != etyp &&  "ElmT6n2D" != etyp && 
      "ElmB2n1D" != etyp &&  "ElmB3n1D" != etyp )
      fatalerrorFP("WavesGridFEHier:: operator =(const GridFE&)\n",
		   "Only elements ElmB2n1D, ElmB3n1D, ElmT3n2D, ElmT6n2D, "
		   "ElmT4n3D and ElmT10n3D are implemented\n"
		   "so far, this is element type %s.\n"
		   "It might be possible to use Prepro::box2triangle "
		   "to change from ElmB4n2D and ElmB8n3D elements?",
		   etyp.chars()); 
		   */

  WavesGridB::operator = (grid); // copy WavesGridB part

  initGridFEHier ();  // initiate GridFEHier specific data
}


//-----------------------------------------------------------------------------
void WavesGridFEHier:: operator = (const WavesGridFEHier& grid)
//-----------------------------------------------------------------------------
{ 
  if (!grid.ok())
    //    fatalerrorFP("GridFEHier:: operator =","Empty grid is given as argument.");
    printf(" fatalerrorFP WavesGridFEHier:: operator = Empty grid is given\n");

  WavesGridB::operator = (grid);

  initGridFEHier ();

// copy WavesGridFEHier specific data

  neighbor_computed = grid.neighbor_computed;
  new_nodes_are_set = grid.new_nodes_are_set;
  finest_grid = grid.finest_grid;
  coarsest_grid = grid.coarsest_grid;
  parent_info = grid.parent_info;
  children_info = grid.children_info;
  regular_uniform_ref = grid.regular_uniform_ref;
  refinement_method = grid.refinement_method;
  reg_bisect_short = grid.reg_bisect_short;

  current_elm_number = grid.current_elm_number;
  current_node_number = grid.current_node_number;
  bound_box_computed = grid.bound_box_computed;

  if( bound_box_computed ){
    bounding_box.newsize(nel,2*nsd);
    bounding_box = grid.bounding_box;
  }
  if( neighbor_computed ){
    //if already exist, delete it (if not correct size)
    elm_neigh.newsize( nel, nsd+1 );
    elm_neigh = grid.elm_neigh;
    if( grid.edge_face2.size() || grid.edge_face3.size() ){
      edge_face2.newsize(grid.edge_face2.size());
      edge_face3.newsize(grid.edge_face3.size());
      edge_face2 = grid.edge_face2;
      edge_face3 = grid.edge_face3;
    }
  }

  if( new_nodes_are_set ){
    new_nodes.newsize( nel, grid.new_nodes.size(2) );
    new_nodes = grid.new_nodes;
  }
  if( children_info ){
    pos_child.newsize( nel+1 );
    pos_child = grid.pos_child;
  }
  if( parent_info ){
    parent.newsize( nel );
    parent = grid.parent;
  }

  //  if( prev_grid.ok())
  //    prev_grid.rebind(grid.prev_grid());
  prev_grid = grid.prev_grid;

  next_grid = grid.next_grid;

}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: makeCopyOfHierarchy(  MV_Vector<WavesGridFEHier>& copiedgrids,
					       const int firstlevel, 
					       const int lastlevel )
//-----------------------------------------------------------------------------
{
// Create a copy of a part of the grid hierarchy.
// firstlevel is the coarsest level in "this" which is copied
// and which will be the coarsest grid in the copy.
// lastlevel is the finest level which is copied
// and which will be the finest grid in the copy.
//
// One application of this member function is to create a sequence
// of hierarchical grids which have some common coarse grids.

  //  Handle(WavesGridFEHier) oldgrid;
  //  oldgrid.rebind(this);
  WavesGridFEHier* oldgrid = this;

  //  if( firstlevel < 1 ){
  if( firstlevel < 0 ){
    //    errorFP("WavesGridFEHier:: makeCopyOfHierarchy\n",
    //	    "Requested first level %d is to coarse, must be > 0.", firstlevel);
    printf(" WavesGridFEHier:: makeCopyOfHierarchy Requested first level %d is too coarse, must be > 0.\n",firstlevel);
    return this;
    //    return oldgrid.getPtr();
  }

  if( lastlevel > getFinestGridLevelNo () ){
    //    errorFP("WavesGridFEHier:: makeCopyOfHierarchy\n",
    //	    "Requested last level %d is to fine, only %d levels exists.",
    //	    lastlevel,getFinestGridLevelNo());
    printf("WavesGridFEHier:: makeCopyOfHierarchy, Requested last level %d is to fine, only %d levels exists.\n",lastlevel,getFinestGridLevelNo ());
    //    return oldgrid.getPtr();
    return this;
  }

  if( firstlevel > lastlevel ){
    //    errorFP("WavesGridFEHier:: makeCopyOfHierarchy\n",
    //	    "First level is larger than last level.");
    printf("WavesGridFEHier:: makeCopyOfHierarchy First level is larger than last level\n");
    return this;
    //    return oldgrid.getPtr();
  }

  int i;
  
  //  oldgrid.rebind( getCoarsestGrid() );
  oldgrid = (WavesGridFEHier*) this->getCoarsestGrid();
// before the first level
  for( i=0; i<firstlevel; i++ )
    //    oldgrid.rebind(oldgrid().getChildGrid());
    oldgrid = oldgrid->getChildGrid();

  WavesGridFEHier* g;
  WavesGridFEHier* gg;
  int j=0;
  for( i=firstlevel; i<=lastlevel; i++ ){
    //    g = new WavesGridFEHier(oldgrid()); // copy one of the grids 
    //    g = new WavesGridFEHier(*oldgrid); // copy one of the grids 
    copiedgrids(j) = *oldgrid;
    g = &copiedgrids(j);
    j++;
    if ( i == firstlevel ){
      gg = g; // connect to new grid
      gg->removeParentGrid(); // make sure this is the coarsest grid 
    }
    else {
      gg->setChildGrid( g ); // set child pointer
      g->setParentGrid( gg ); // set parent pointer
      gg = g; // connect to new grid
    } 

    if ( i == lastlevel )
      gg->removeChildGrid(); // now the finest grid in the new hierarchy 
    else
      //      oldgrid.rebind( oldgrid().getChildGrid() );// go one step up in hierarchy
      oldgrid=oldgrid->getChildGrid();// go one step up in hierarchy
  }
  return gg;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getGridLevelNo() const
//-----------------------------------------------------------------------------
{ //recursive call
  return ( isThisCoarsestGrid() ) ? 1 : getParentGrid()->getGridLevelNo() + 1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getFinestGridLevelNo() const
//-----------------------------------------------------------------------------
{
  return getFinestGrid()->getGridLevelNo();
}

//-----------------------------------------------------------------------------
const WavesGridFEHier* WavesGridFEHier:: getFinestGrid () const
//-----------------------------------------------------------------------------
{ //recursive call
  return ( isThisFinestGrid() ) ? this : getChildGrid()->getFinestGrid();
}

//-----------------------------------------------------------------------------
const WavesGridFEHier* WavesGridFEHier:: getCoarsestGrid () const
//-----------------------------------------------------------------------------
{ //recursive call
  return ( isThisCoarsestGrid() ) ? this : getParentGrid()->getCoarsestGrid();
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeNeighborInfo ()
//-----------------------------------------------------------------------------
// Removal of neighbor data, i.e., elm_neigh and new_nodes
{
  elm_neigh.newsize(0,0);
  new_nodes.newsize(0,0);
  edge_face2.newsize(0);
  edge_face3.newsize(0);
  setNeighborsComputed( false );
  setNewNodesSet( false );
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeNewNodeInfo()
//-----------------------------------------------------------------------------
// Removal of new_nodes data showing the node numbering  or marking of edges
{
  new_nodes.newsize(0,0);
  setNewNodesSet( false );
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeChildrenInfo()
//-----------------------------------------------------------------------------
// Removal of the information about where the children element are.
{
  pos_child.newsize(0);
  setChildrenInfo( false );
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeChildGrid()
//-----------------------------------------------------------------------------
// Removal of pointer to the child grid and 
// info about where the children elements are.
{
  next_grid = NULL;
  setFinestGrid( true ); // This grid is now finest
  removeChildrenInfo();    // without child grid no need for children info
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeParentInfo()
//-----------------------------------------------------------------------------
// removal of info on where the parent elements are
{
  parent.newsize(0);
  setParentInfo( false );
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeParentGrid()
//-----------------------------------------------------------------------------
// removal of parent info and prev_grid
{
  //  prev_grid.detach();
  prev_grid=NULL; // should maybe delete the grid, any references
  setCoarsestGrid( true ); // if no parent then this grid is coarsest
  removeParentInfo(); // No need for info about parent elms if no parent grid
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeHierInfo()
//-----------------------------------------------------------------------------
// removal of children, parent, prev_ and next_grid
{
  removeChildGrid();
  removeParentGrid();
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initRefine()
//-----------------------------------------------------------------------------
{
  /*
  if( "ElmT4n3D" != elm_tp && "ElmT10n3D" != elm_tp &&
      "ElmT3n2D" != elm_tp && "ElmT6n2D" != elm_tp && 
      "ElmB2n1D" != elm_tp && "ElmB3n1D" != elm_tp )
      fatalerrorFP("WavesGridFEHier:: initRefine\n",
		   "Only elements ElmB2n1D, ElmB3n1D, ElmT3n2D, ElmT6n2D, "
//		   " ElmT4n3D and ElmT10n3D are implemented so far,\n"
		   " and ElmT4n3D are implemented so far,\n"
		   " this is element type %s",elm_tp.chars());
		   */


  //?????????????
  if( ! isNeighborsComputed() ) // is the element-element info computed?
    initNeighbor();            // required for numbering of new nodes
  else if (nsd==3) {//in 3D it is necessary to make sure that n2e info is comp.
    //    additional_info.initNeighbor(*this,true,false,false);
    //    WavesNeighborFE& neighbor = additional_info.neighbor;
    if( neighbor.nodeSize() == 0 )
      neighbor.init(*this,true,false,false);
    //    if( neighbor.nodeSize() == 0 )
    //      neighbor.init(*this,true,false,false);
  }
  if( ! isThisFinestGrid() ) {
    //    warningFP("WavesGridFEHier:: initRefine", 
    //	      "\nA refined grid already exists, will replace it.");
    printf("WavesGridFEHier:: initRefine, refined grid already exists, will replace it.\n");
    removeChildGrid();
  }
}

//-----------------------------------------------------------------------------
//WavesGridFEHier* WavesGridFEHier:: refine(const VecSimple(int)& mark_elms)
WavesGridFEHier* WavesGridFEHier:: refine ( WavesGridFEHier* newgrid, const MV_Vector<int>& mark_elms )
//-----------------------------------------------------------------------------
{
  initRefine();
  /*  printf(" edgeface2 has %d elements",edge_face2.size());
  printf(" edgeface3 has %d elements",edge_face3.size());
  for(int i=0;i< edge_face2.size();i++)
    printf(" i=%d, edge_face2(i)=%d, edge_face3(i)=%d\n",i,edge_face2(i),edge_face3(i));
    */
  markElements( mark_elms );
  if( refine( newgrid ) ) 
    return getChildGrid();
  else
    //    errorFP("WavesGridFEHier:: refine", "Refinement failed.");
    printf("WavesGridFEHier:: refine Refinement failed\n.");
  return NULL;
}

//-----------------------------------------------------------------------------
//WavesGridFEHier* WavesGridFEHier:: refine(const VecSimple(bool)& mark_elms,
WavesGridFEHier* WavesGridFEHier:: refine ( WavesGridFEHier* newgrid, const MV_Vector<bool>& mark_elms,
				const int ref_method )
//-----------------------------------------------------------------------------
{
  initRefine();
  setRefinementMethod( ref_method );
  markElementsBool( mark_elms );
  if( refine( newgrid ) ) 
    return getChildGrid();
  else
    //    errorFP("WavesGridFEHier:: refine", "Refinement failed.");
    printf("WavesGridFEHier:: refine Refinement failed\n.");
  return NULL;
}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: refineUniformly( WavesGridFEHier* newgrid, const int ref_method, bool regref)
//-----------------------------------------------------------------------------
// refine uniformly, if regref is true and all edges in 2D and 3D will
// be refined, then regular subdivision is chosen
{ //CPUclock cpu; cpu.initTime();
  //  cout<<"refineUniformly begin\n"<<flush;
  initRefine(); //cout<<"\n after initrefine = " << cpu.getInterval()<<"\n";
  //  cout<<"refineUniformly after initRefine\n"<<flush;
  setRefinementMethod( ref_method );
  markAllElements(); //cout<<"\n after markall = " << cpu.getInterval()<<"\n";
  //  cout<<"refineUniformly after markAllElements"<<flush;
  if( (nsd == 3 && ref_method == 6) || (nsd == 2 && ref_method == 3) )
    setRegularUniformRef(regref); // bisection or regular refinement?
  if( refine( newgrid ) ){ // this is where the actual refinement takes place
    setRegularUniformRef(false);// reset value for next refinement
//    cout<<"\n after refine = " << cpu.getInterval()<<"\n";
   return getChildGrid();
  }
  else
    //    errorFP("WavesGridFEHier:: refineUniformly", "refinement failed.");
    printf("WavesGridFEHier:: refineUniformly Refinement failed\n.");
  return NULL;
}
//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: refineLongEdges( WavesGridFEHier* newgrid, const real break_ratio )
//-----------------------------------------------------------------------------
// Let maxlen be the longest edge in the whole grid
// refine all edges which are longer or equal than break_ratio of maxlen
// if ElmT6n2D or ElmT10n3D only the vertex nodes will be used.
{
  int e;
  //  VecSimple(int) markings(nel);
  //  printf("DEBUG begin refineLongEdges\n");
  MV_Vector<int> markings(nel);
  if ( nsd==3 ) {
    FET4.refill(0);
    real maxlen = FET4.maxLength(); // find longest edge in grid
    for( e=0; e<nel; e++ ) {
      FET4.refill(e);
      //      maxlen = max(maxlen,FET4.maxLength());
      maxlen = (maxlen>FET4.maxLength())?maxlen:FET4.maxLength();
    } 
    real b2 = sqr(break_ratio*maxlen); // use square for efficiency
    for( e=0; e<nel; e++ ) {
      FET4.refill(e);
      int p = 0;
      if( FET4.l1sq() >= b2 ) p++;
      if( FET4.l2sq() >= b2 ) p++;
      if( FET4.l3sq() >= b2 ) p++;
      if( FET4.l4sq() >= b2 ) p++;
      if( FET4.l5sq() >= b2 ) p++;
      if( FET4.l6sq() >= b2 ) p++;
      markings(e) = p; // at least p edges in element e will be subdivided
      //      printf("DEBUG element %d markings %d refineLongEdges\n",e,p);
    }
  }
  else if ( nsd==2 ) {
    FET3.refill(0);
    real maxlen = FET3.maxLength(); 
    for( e=1; e<nel; e++ ) {
      FET3.refill(e);
      //      maxlen = max(maxlen,FET3.maxLength());
      maxlen = (maxlen>FET3.maxLength())?maxlen:FET3.maxLength();
    } 
    real b2 = sqr(break_ratio*maxlen);
    for( e=0; e<nel; e++ ){
      FET3.refill(e);
      int p = 0;
      if( FET3.l1sq() >= b2 ) p++;
      if( FET3.l2sq() >= b2 ) p++;
      if( FET3.l3sq() >= b2 ) p++;
      markings(e) = p; // at least p edges in element e will be subdivided
    }
  }
  else
    //    errorFP("WavesGridFEHier:: refineLongEdges","dimension %d not implemented",nsd);
    printf("WavesGridFEHier:: refineLongEdges Dimension %d not implemented\n",nsd);
  return refine( newgrid, markings ); // refine with "mixed method" 
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: refine(WavesGridFEHier* newgrid)
//-----------------------------------------------------------------------------
// return true if refinement successful
{
  switch ( nsd ) {
  case 1: return refine1D (newgrid);
  case 2: return refine2D (newgrid);
  case 3: return refine3D (newgrid);
//  default: fatalerrorFP("WavesGridFEHier:: refine"," nsd = %d not implemented",nsd);
  default: printf("WavesGridFEHier:: refine Dimension %d not implemented\n",nsd);
  } 
  return false;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initRefinedGrid ( WavesGridFEHier* newgrid,int new_nno, int new_nel )
//-----------------------------------------------------------------------------
// dimensioning the new grid: new WavesGridFEHier, children, parent
// common for 1D,2D and 3D linear and quadratic
{
//  newgrid->redim( nsd, new_nno, new_nel, maxnne, nbind, nne, onesbd );
  ElementType etype = getElementType(0);
  newgrid->redim( nsd, new_nno, new_nel, maxnne, nbind,etype );
//  newgrid->setElmType( elm_tp );
//  newgrid->setNonUniformMesh();
  newgrid->setRefinementMethod( getRefinementMethod() );
  newgrid->setRegBisectOrShort( getRegBisectOrShort() );// only used in 2D
// set hierarchy relations of old grid and refined grid
  newgrid->setParentGrid( this );  

  //k  newgrid->getSurfaceTriangulation().attachGeometry( this->getSurfaceTriangulation().getGeometryPtr() );

  setChildGrid ( newgrid );
  setFinestGrid( false );
  newgrid->setFinestGrid( true );
  newgrid->setCoarsestGrid( false );
  setChildrenInfo( true );
  newgrid->setParentInfo( true );
  pos_child.newsize( nel+1 ); // pos_child and parent will be filled later
  newgrid->parent.newsize( new_nel );
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: checkReadyForRefine()
//-----------------------------------------------------------------------------
// Check that the grid is ready for refinement
{
  if ( ! isNeighborsComputed() )
    //    fatalerrorFP("WavesGridFEHier:: checkReadyForRefine", 
    //		 "Neighbor data is not computed.");
    printf("WavesGridFEHier:: checkReadyForRefine Neighbor data is not computed/n");
  if ( ! isThisFinestGrid() ) 
    //    fatalerrorFP("WavesGridFEHier:: checkReadyForRefine", 
    //		 "Refined grid already exists.");
    printf("WavesGridFEHier:: checkReadyForRefine Refined grid already exists/n");
  if( ! isNewNodesSet() ) 
    //    fatalerrorFP("WavesGridFEHier:: checkReadyForRefine",
    //		 "Edges must be marked before refinement.");
    printf("WavesGridFEHier:: checkReadyForRefine Edges must be marked before refinemen/n");
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: refine1D( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// return true if refinement successful
{
  checkReadyForRefine (); // check if the grid is ok
  /*
  if( "ElmB2n1D" != elm_tp && "ElmB3n1D" != elm_tp )
      fatalerrorFP("WavesGridFEHier:: refine",
		   "Not implemented yet, element type = %s",elm_tp.chars() );
		   */
  int b,e,i;

// Find the number of nodes and elements in the refined grid.
// The new_nodes gives the number of new nodes in interval

  int new_nel = nel;
  //  for( e=1; e<=nel; e++ ) // should work for both linear and quadratic elms
  //    new_nel += new_nodes( e, 1 );
  for( e=0; e<nel; e++ ) // should work for both linear and quadratic elms
    new_nel += new_nodes( e, 0 );

  setCurrNodeNo(nno);
  //why +1 ?
    int new_nno = getCurrNodeNo()+1 + new_nel - nel;
 

  //  int new_nno = getCurrNodeNo() + new_nel - nel; 
  //  if( "ElmB3n1D" == elm_tp ) 
  if( maxnne==3 ) 
    new_nno = getCurrNodeNo()+1 + 2*(new_nel - nel); 

// dimensioning the new grid: new WavesGridFEHier, children, parent
  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  initRefinedGrid( newgrid, new_nno, new_nel );

// copy the name of the boundary indicators
  /*
  for (i = 1; i <= nbind; i++)
    newgrid->putBoIndName( getBoIndName(i), i );

// copy the old nodes which will be present in the refined grid
// and copy boundary indicators

  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */

  // for( i=1; i<=nno; i++ )
  //  newgrid->coord(i,1) = coord(i,1); 
 
  for( i=0; i < nno; i++ )
    newgrid->coord(i,0) = coord(i,0); 
 
   
 /*  if ( nbind ) 
    for( i=1; i<=nno; i++ )
      for( b=1; b<=nbind; b++ )
	if( ind(i,b) == '1' ) 
	  nind(i,b) = '1';
	  */
// Start the numbering of the elements in the new grid so that the first 
// element is given element number 1.
  setCurrElmNo(-1); 

  int e1,e2,e3,e4,e5; //numbers for new elements
  int m1,m2,m3,m4,m5,m6,m7,m8;    //numbers for new nodes
  int n1,n2,n12;      //node numbers in old grid
  real x1,x2,x12;     //coordinates in old grid

  //  if( "ElmB2n1D" == elm_tp ) {
  if( maxnne==2 ) {
    for( e=0; e < nel; e++ ) {//Refine the elements independently and set data
       pos_child( e ) = getCurrElmNo()+1; //???????
    

      n1 = loc2glob( e, 0 );  
      n2 = loc2glob( e, 1 );
      x1 = coord( n1, 0 );  
      x2 = coord( n2, 0 ); 
      
      int number_of_new_elements = new_nodes( e, 0 ) + 1;
      switch ( number_of_new_elements ) {
//----------------------------------------------------------------
      case 1: // no refinement: just copy the info to the new element
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
	newgrid->putLoc2Glob12(n1, n2, e1);
	break;
//----------------------------------------------------------------
      case 2: // the element is divided into two elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();
	newgrid->putLoc2Glob12( n1, m1, e1 );
	newgrid->putLoc2Glob12( m1, n2, e2 );
	newgrid->coord(m1,0) = (x1 + x2)/2.0;

	cout<<" e1 = incCurrElmNo() "<<e1<<"e2 "<<e2<<"  m1  "<<m1<<" n1 "<<n1<<" n2  "<<n2<<"  coord(m1,0) "<<(x1+x2)/2.0<<endl;


	/*	for( b=1; b<=nbind; b++ )
	  if( ind(n1,b) == '1' && ind(n2,b) == '1') // if both ends are set...
	    nind(m1,b) = '1';  // then set midpoint
	    */
	break;
//----------------------------------------------------------------
      case 3: // the element is divided into three elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();	m2 = incCurrNodeNo();
	newgrid->putLoc2Glob12( n1, m1, e1 );
	newgrid->putLoc2Glob12( m1, m2, e2 );
	newgrid->putLoc2Glob12( m2, n2, e3 );
	newgrid->coord(m1,0) = x1*0.6666666666666666 + x2*0.3333333333333333;
	newgrid->coord(m2,0) = x1*0.3333333333333333 + x2*0.6666666666666666;
	/*	for( b=1; b<=nbind; b++ )
	  if( ind(n1,b) == '1' && ind(n2,b) == '1')
	    nind(m1,b) = nind(m2,b) = '1';
	    */
	break;
//----------------------------------------------------------------
      case 4: // the element is divided into four elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();  m3 = incCurrNodeNo();
	newgrid->putLoc2Glob12( n1, m1, e1 );
	newgrid->putLoc2Glob12( m1, m2, e2 );
	newgrid->putLoc2Glob12( m2, m3, e3 );
	newgrid->putLoc2Glob12( m3, n2, e4 );
	newgrid->coord(m1,0) = x1*0.75 + x2*0.25;
	newgrid->coord(m2,0) = (x1 + x2)*0.5;
	newgrid->coord(m3,0) = x1*0.25 + x2*0.75;
	/*	for( b=1; b<=nbind; b++ )
	  if( ind(n1,b) == '1' && ind(n2,b) == '1')
	    nind(m1,b) = nind(m2,b) = nind(m3,b) = '1';
	    */
	break;
//----------------------------------------------------------------
      case 5: // the element is divided into five elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
	newgrid->putParent( e5 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();
	m3 = incCurrNodeNo();  m4 = incCurrNodeNo();
	newgrid->putLoc2Glob12( n1, m1, e1 );
	newgrid->putLoc2Glob12( m1, m2, e2 );
	newgrid->putLoc2Glob12( m2, m3, e3 );
	newgrid->putLoc2Glob12( m3, m4, e4 );
	newgrid->putLoc2Glob12( m4, n2, e5 );
	newgrid->coord(m1,0) = x1*0.8 + x2*0.2;
	newgrid->coord(m2,0) = x1*0.6 + x2*0.4;
	newgrid->coord(m3,0) = x1*0.4 + x2*0.6;
	newgrid->coord(m4,0) = x1*0.2 + x2*0.8;
	/*	for( b=1; b<=nbind; b++ )
	  if( ind(n1,b) == '1' && ind(n2,b) == '1')
	    nind(m1,b) = nind(m2,b) = nind(m3,b) = nind(m4,b) = '1';*/
	break;
//------------------------------------------------------
      default: 
//------------------------------------------------------
	//	fatalerrorFP("WavesGridFEHier:: refine1D",
	//		     "There must be 1-5 new elements.");
	printf("WavesGridFEHier:: refine1D There must be 1-5 new elements.\n");
      } //end of switch
    }
  }
  //  else if( "ElmB3n1D" == elm_tp ) {
  else if( maxnne==3 ) {
    for( e=0; e < nel; e++ ) {//Refine the elements independently and set data
      pos_child( e ) = getCurrElmNo()+1; 
      
      n1  = loc2glob( e, 0 ); 
      n12 = loc2glob( e, 1 ); 
      n2  = loc2glob( e, 2 );
      x1  = coord( n1, 0 );  
      x12 = coord( n12,0 );
      x2  = coord( n2, 0 );

      int number_of_new_elements = new_nodes( e, 0 ) + 1;
      switch ( number_of_new_elements ) {
//----------------------------------------------------------------
      case 1: // no refinement: just copy the info to the new element
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
	newgrid->putLoc2Glob123( n1, n12, n2, e1 );
	break;
//----------------------------------------------------------------
      case 2: // the element is divided into two elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n1,  m1, n12, e1 );
	newgrid->putLoc2Glob123( n12, m2, n2,  e2 );
	newgrid->coord(m1,0) = (x1 + x12)/2;
	newgrid->coord(m2,0) = (x2 + x12)/2;
	/*	for( b=1; b<=nbind; b++ ){
	  if( ind(n1,b) == '1' && ind(n12,b) == '1') //if both ends are set...
	    nind(m1,b) = '1';
	  if( ind(n12,b) == '1' && ind(n2,b) == '1') //if both ends are set...
	    nind(m2,b) = '1';
           }
	    */
      
	break;
//----------------------------------------------------------------
      case 3: // the element is divided into three elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();
	m3 = incCurrNodeNo();  m4 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n1, m1,  m2, e1 );
	newgrid->putLoc2Glob123( m2, n12, m3, e2 );
	newgrid->putLoc2Glob123( m3, m4,  n2, e3 );
	x1 = coord(n1,1);  x2 = coord(n2,1); 
	newgrid->coord(m1,0) = x1*0.6666666666666666 + x12*0.3333333333333333;
	newgrid->coord(m2,0) = x1*0.3333333333333333 + x12*0.6666666666666666;
	newgrid->coord(m3,0) = x2*0.3333333333333333 + x12*0.6666666666666666;
	newgrid->coord(m4,0) = x2*0.6666666666666666 + x12*0.3333333333333333;
	/*	for( b=1; b<=nbind; b++ ) {
	  if( ind(n1,b) == '1' && ind(n12,b) == '1')
	    nind(m1,b) = nind(m2,b) = '1';
	  if( ind(n12,b) == '1' && ind(n2,b) == '1')
	    nind(m3,b) = nind(m4,b) = '1';
	}
	    */
	break;
//----------------------------------------------------------------
      case 4: // the element is divided into four elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();  m3 = incCurrNodeNo();
	m4 = incCurrNodeNo();  m5 = incCurrNodeNo();  m6 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n1,  m1, m2,  e1 );
	newgrid->putLoc2Glob123( m2,  m3, n12, e2 );
	newgrid->putLoc2Glob123( n12, m4, m5,  e3 );
	newgrid->putLoc2Glob123( m5,  m6, n2,  e4 );
	newgrid->coord(m1,0) = x1*0.75 + x12*0.25;
	newgrid->coord(m2,0) = x1*0.5  + x12*0.5;
	newgrid->coord(m3,0) = x1*0.25 + x12*0.75;
	newgrid->coord(m4,0) = x2*0.25 + x12*0.75;
	newgrid->coord(m5,0) = x2*0.5  + x12*0.5;
	newgrid->coord(m6,0) = x2*0.75 + x12*0.25;
	/*	for( b=1; b<=nbind; b++ ){
	  if( ind(n1,b) == '1' && ind(n12,b) == '1')
	    nind(m1,b) = nind(m2,b) = nind(m3,b) = '1';
	  if( ind(n12,b) == '1' && ind(n2,b) == '1')
	    nind(m4,b) = nind(m5,b) = nind(m6,b) = '1';
	}
	    */
	break;
//----------------------------------------------------------------
      case 5: // the element is divided into five elements
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
	newgrid->putParent( e5 = incCurrElmNo(), e );
	m1 = incCurrNodeNo();  m2 = incCurrNodeNo();
	m3 = incCurrNodeNo();  m4 = incCurrNodeNo();
	m5 = incCurrNodeNo();  m6 = incCurrNodeNo();
	m7 = incCurrNodeNo();  m8 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n1, m1,  m2, e1 );
	newgrid->putLoc2Glob123( m2, m3,  m4, e2 );
	newgrid->putLoc2Glob123( m4, n12, m5, e3 );
	newgrid->putLoc2Glob123( m5, m6,  m7, e4 );
	newgrid->putLoc2Glob123( m7, m8,  n2, e5 );
	newgrid->coord(m1,0) = x1*0.8 + x12*0.2;
	newgrid->coord(m2,0) = x1*0.6 + x12*0.4;
	newgrid->coord(m3,0) = x1*0.4 + x12*0.6;
	newgrid->coord(m4,0) = x1*0.2 + x12*0.8;
	newgrid->coord(m5,0) = x2*0.2 + x12*0.8;
	newgrid->coord(m6,0) = x2*0.4 + x12*0.6;
	newgrid->coord(m7,0) = x2*0.6 + x12*0.4;
	newgrid->coord(m8,0) = x2*0.8 + x12*0.2;
	/*	for( b=1; b<=nbind; b++ ){
	  if( ind(n1,b) == '1' && ind(n12,b) == '1')
	    nind(m1,b) = nind(m2,b) = nind(m3,b) = nind(m4,b) = '1';
	  if( ind(n12,b) == '1' && ind(n2,b) == '1')
	    nind(m5,b) = nind(m6,b) = nind(m7,b) = nind(m8,b) = '1';
         }
	    */
	break;
//------------------------------------------------------
      default: 
//------------------------------------------------------
	//	fatalerrorFP("WavesGridFEHier:: refine1D",
	//		     "There must be 1-5 new elements.");
	printf("WavesGridFEHier:: refine1D There must be 1-5 new elements.\n");
      } //end of switch
    }
  }
  else
//    fatalerrorFP("WavesGridFEHier:: refine1D","Not implemented element");
    printf("WavesGridFEHier:: refine1D  Not implemented element\n");
  if( new_nno != getCurrNodeNo()+1 )
    //    fatalerrorFP("WavesGridFEHier:: refine1D","new_nno != getCurrNodeNo()");
    printf("WavesGridFEHier:: refine1D new_nno != getCurrNodeNo()\n");

  //  pos_child( nel+1 ) = getCurrElmNo()+1; 
  pos_child( nel ) = getCurrElmNo()+1; 

// set the subdomain types for the new grid:
// copy the info from the parent elements

  if( ! onemat )
    for(e=0;e<new_nel;e++) 
      newgrid->setMaterialType(e,getMaterialType(newgrid->getParent_eff(e)));

  removeNewNodeInfo(); // remove the data showing the marked edges
  newgrid->setParentInfo( true );
  setChildrenInfo( true );

  //#ifdef ARRAY_RANGECHECK
  if( ! newgrid->checkElementOrientation() )
    printf("Wrong orientation.\n");
  if( ! checkChildrenInfo() ) 
    printf("children info failed.\n");
  if( ! newgrid->checkParentInfo() ) 
    printf("parent info failed.\n");
  //#endif
  return true; //everything went fine
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: refine2D( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// return true if refinement successful
{
  checkReadyForRefine (); // check if the grid is ok
  /*  if( "ElmT3n2D" != elm_tp && "ElmT6n2D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: refine",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
  int b,e,i;

// Find the number of nodes and elements in the refined grid.
// The new nodes have been numbered before.
// A triangle which have n marked/numbered edges will be refined
// into n+1 elements.

// The new_nodes structure is 0 where no edge refinement is taking place.

  int new_nel = nel;
  for( e=0; e<nel; e++ )
    if( new_nodes( e, 0 ) ) { 
      new_nel++;
      if( new_nodes( e, 1 ) ) new_nel++;
      if( new_nodes( e, 2 ) ) new_nel++;
    }

// current node number has been incremented during the numbering process
  int new_nno = getCurrNodeNo()+1;
  //  if( "ElmT6n2D" == elm_tp ) // interior nodes are added below
  if( maxnne==6 ) // interior nodes are added below
    new_nno += new_nel - nel;

// dimensioning the new grid: new WavesGridFEHier, children, parent

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  initRefinedGrid( newgrid, new_nno, new_nel );
  /*
// copy the name of the boundary indicators

  for (i = 1; i <= nbind; i++) newgrid->putBoIndName( getBoIndName(i), i );

// copy the old nodes which will be present in the refined grid
// and copy boundary indicators

  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
  for( i=0; i<nno; i++ ){ // copy all the nodes from old grid
    newgrid->coord(i,0) = coord(i,0); 
    newgrid->coord(i,1) = coord(i,1); 
  }
  /*
  if ( nbind ) 
    for( i=1; i<=nno; i++ )
      for( b=1; b<=nbind; b++ )
	nind(i,b) = ind(i,b);
	*/
// Go through all elements and compute position of the new nodes
// at the midpoints of edges which have been assigned a node number.
// Also set indicators for these nodes

  int n1,n2,n3,n4,n5,n6;
  //  if( "ElmT3n2D" == elm_tp )
  if( maxnne=3 )
    for( e=0; e<nel; e++ )
      if ( (n4 = new_nodes( e, 0 )) ){
	n1 =  loc2glob( e, 0 );
	n2 =  loc2glob( e, 1 );
	n3 =  loc2glob( e, 2 );
	n5 = new_nodes( e, 1 );
	n6 = new_nodes( e, 2 );
	if( ! newgrid->coord(n4,0) ){ // to avoid computing twice
	  newgrid->coord(n4,0) = (coord(n1,0)+coord(n2,0))/2;
	  newgrid->coord(n4,1) = (coord(n1,1)+coord(n2,1))/2;
	  /*
	  for( b=1; b<=nbind; b++ )
	    if( ind(n1,b) == '1' && ind(n2,b) == '1') //if both ends are set...
	      nind(n4,b) = '1';
	      */
	}
	if( ((n5) && ! newgrid->coord(n5,0)) ){ 
	  newgrid->coord(n5,0) = (coord(n2,0)+coord(n3,0))/2;
	  newgrid->coord(n5,1) = (coord(n2,1)+coord(n3,1))/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( ind(n2,b) == '1' && ind(n3,b) == '1' )
	      nind(n5,b) = '1';
	      */
	}
	if( ((n6) && ! newgrid->coord(n6,0)) ){
	  newgrid->coord(n6,0) = (coord(n3,0)+coord(n1,0))/2;
	  newgrid->coord(n6,1) = (coord(n3,1)+coord(n1,1))/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( ind(n3,b) == '1' && ind(n1,b) == '1' )
	      nind(n6,b) = '1';
	      */
	}
      }

// Start the numbering of the elements in the new grid so that the first 
// element is given element number 0.
  setCurrElmNo(-1); 
  //  bool T6 = ( "ElmT6n2D" == elm_tp ) ? true : false;
  bool T6 = ( maxnne==6 ) ? true : false;

  int e1,e2,e3,e4; //numbers for new elements
  int m1,m2,m3;
  int mm1,mm2,mm3, nn1,nn2,nn3; // used for quadratic elements

// what should be done if 3 edges are marked:
  bool REGULAR, TEST_SHORTEST;
  if( getRegularUniformRef () || reg_bisect_short==0 ){
    REGULAR = true; // always regular refinement
    TEST_SHORTEST = false; 
  }
  else {
    if( reg_bisect_short == 1 ){
      REGULAR = false; // always bisection refinement
      TEST_SHORTEST = false;
    }
    else // case reg_bisect_short == 2
      TEST_SHORTEST = true;// REGULAR will be determined by shortest edge
  }

  real l1,l2;
  int number_of_new_elements;
  if( ! T6 ) {  
    for( e=0; e < nel; e++ ) { // Refine elements independently and set data
      pos_child( e ) = getCurrElmNo()+1; 
      number_of_new_elements = 1;
      if( (m1 = new_nodes( e, 0 )) ) {
	number_of_new_elements++;
	if( (m2 = new_nodes( e, 1 )) ) // m2 and m3 only if m1 is set
	  number_of_new_elements++;
	if( (m3 = new_nodes( e, 2 )) ) 
	  number_of_new_elements++;
      }
      switch ( number_of_new_elements ) {
//----------------------------------------------------------------
      case 1: // no refinement: just copy the info to the new element
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
	newgrid->putLoc2Glob123(loc2glob(e,0),loc2glob(e,1),loc2glob(e,2), e1);
	break;
//----------------------------------------------------------------
      case 2:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	n3 = loc2glob( e, 2 );
	newgrid->putLoc2Glob123( n3, loc2glob(e,0), m1, e1 );
	newgrid->putLoc2Glob123( loc2glob(e,1), n3, m1, e2 );
	break;
//----------------------------------------------------------------
      case 3:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	n3 =  loc2glob( e, 2 );
	if( m2 ) { // edge 1 and edge 2 marked
	  newgrid->putLoc2Glob123( n3, loc2glob(e,0), m1, e1 );
	  newgrid->putLoc2Glob123( n3, m1, m2, e2 );
	  newgrid->putLoc2Glob123( m1, loc2glob(e,1), m2, e3 );
	}
	else { // edge 1 and edge 3 marked
	  newgrid->putLoc2Glob123( loc2glob(e,1), n3, m1, e1 );
	  newgrid->putLoc2Glob123( m1, n3, m3, e2 );
	  newgrid->putLoc2Glob123( loc2glob(e,0), m1, m3, e3 );
        }
	break;
//----------------------------------------------------------------
      case 4:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
// it is assumed that edge 1 is the longest edge, these triangles are common
	n3 = loc2glob( e, 2 );
	newgrid->putLoc2Glob123( loc2glob( e, 0 ), m1, m3, e1 );
	newgrid->putLoc2Glob123( m1, loc2glob( e, 1 ), m2, e2 );
	if ( TEST_SHORTEST ) { // choose the edge which is shortest
          // an alternative is to check if largest angle is obtuse
	  l1 = sqr(newgrid->coord(m2,0) - newgrid->coord(m3,0)) +
	    sqr(newgrid->coord(m2,1) - newgrid->coord(m3,1));
	  l2 = sqr(newgrid->coord(m1,0) - newgrid->coord(n3,0)) +
	    sqr(newgrid->coord(m1,1) - newgrid->coord(n3,1));
	  REGULAR = ( l1*0.999999 <= l2 ) ? true : false;
	  // the factor 0.9999999 is to avoid arbitrariness in the equal case
	}
	if( REGULAR ) { // regular refinement
	  newgrid->putLoc2Glob123( m3, m2, n3, e3 );
	  newgrid->putLoc2Glob123( m2, m3, m1, e4 );
	} else { // three edges drawn to midpoint of longest edge
	  newgrid->putLoc2Glob123( m1, n3, m3, e3 );
	  newgrid->putLoc2Glob123( n3, m1, m2, e4 );
	}
	break;
//------------------------------------------------------
      default: 
//------------------------------------------------------
	//	fatalerrorFP("WavesGridFEHier:: refine","There must be 1-4 new elements.");
	printf("WavesGridFEHier:: refine There must be 1-4 new elements.\n");
      } //end of switch
    }
  }
  else { // T6 is true,  quadratic elements
    for( e=0; e < nel; e++ ) { // Refine elements independently and set data
      pos_child( e ) = getCurrElmNo()+1; 
      number_of_new_elements = 1;
      if( (nn1 = new_nodes( e, 0 )) ) {
	number_of_new_elements++;
	if( (nn2 = new_nodes( e, 1 )) ) 
	  number_of_new_elements++;
	if( (nn3 = new_nodes( e, 2 )) ) 
	  number_of_new_elements++;
      }
      switch ( number_of_new_elements ) {
//----------------------------------------------------------------
      case 1: // no refinement: just copy the info to the new element
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
	newgrid->putLoc2Glob123(loc2glob(e,0),loc2glob(e,1),loc2glob(e,2), e1);
	newgrid->putLoc2Glob456(loc2glob(e,3),loc2glob(e,4),loc2glob(e,5), e1);
	break;
//----------------------------------------------------------------
      case 2:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	m1 = loc2glob( e, 3 );
	n3 = loc2glob( e, 2 );
	mm1 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n3, loc2glob(e,0), m1, e1 );
	newgrid->putLoc2Glob456( loc2glob(e,5), nn1, mm1, e1 );
	newgrid->putLoc2Glob123( loc2glob(e,1), n3, m1, e2 );
	newgrid->putLoc2Glob456( loc2glob(e,4), mm1, new_nodes(e,3), e2 );
	break;
//----------------------------------------------------------------
      case 3:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	n3 =  loc2glob( e, 2 );
	mm1 = incCurrNodeNo();
	m1 = loc2glob( e, 3 );
	if( nn2 ) { // edge 1 and edge 2 are marked
	  m2 = loc2glob( e, 4 );
	  mm2 = incCurrNodeNo();
	  newgrid->putLoc2Glob123( n3, loc2glob(e,0), m1, e1 );
	  newgrid->putLoc2Glob456( loc2glob(e,5) , nn1, mm2, e1 );
	  newgrid->putLoc2Glob123( n3, m1, m2, e2 );
	  newgrid->putLoc2Glob456( mm2, mm1, new_nodes(e,4), e2 );
	  newgrid->putLoc2Glob123( m1, loc2glob(e,1), m2, e3 );
	  newgrid->putLoc2Glob456( new_nodes(e,3), nn2, mm1, e3 );
	}
	else { // edge 1 and edge 3 are marked
	  m3 = loc2glob( e, 5 );
	  mm3 = incCurrNodeNo();
	  newgrid->putLoc2Glob123( loc2glob(e,1), n3, m1, e1 );
	  newgrid->putLoc2Glob456( loc2glob(e,4) , mm1, new_nodes(e,3), e1 );
	  newgrid->putLoc2Glob123( m1, n3, m3, e2 );
	  newgrid->putLoc2Glob456( mm1, nn3, mm3, e2 );
	  newgrid->putLoc2Glob123( loc2glob(e,0), m1, m3, e3 );
	  newgrid->putLoc2Glob456( nn1, mm3, new_nodes(e,5), e3 );
        }
	break;
//----------------------------------------------------------------
      case 4:
//----------------------------------------------------------------
	newgrid->putParent( e1 = incCurrElmNo(), e );
	newgrid->putParent( e2 = incCurrElmNo(), e );
	newgrid->putParent( e3 = incCurrElmNo(), e );
	newgrid->putParent( e4 = incCurrElmNo(), e );
// assume that 1 is the longest edge
	n1 = loc2glob( e, 0 );	n2 = loc2glob( e, 1 );	n3 = loc2glob( e, 2 );
	m1 = loc2glob( e, 3 );	m2 = loc2glob( e, 4 );	m3 = loc2glob( e, 5 );
	mm1 = incCurrNodeNo();	mm2 = incCurrNodeNo();	mm3 = incCurrNodeNo();
	newgrid->putLoc2Glob123( n1, m1, m3, e1 );
	newgrid->putLoc2Glob456( nn1, mm1, new_nodes(e,5), e1 );
	newgrid->putLoc2Glob123( m1, n2, m2, e2 );
	newgrid->putLoc2Glob456( new_nodes(e,3), nn2, mm2, e2 );
	if ( TEST_SHORTEST ) {// choose shortest edge alternative 
	  l1 = sqr(newgrid->coord(m2,0) - newgrid->coord(m3,0)) +
	    sqr(newgrid->coord(m2,1) - newgrid->coord(m3,1));
	  l2 = sqr(newgrid->coord(m1,0) - newgrid->coord(n3,0)) +
	    sqr(newgrid->coord(m1,1) - newgrid->coord(n3,1));
	  REGULAR = ( l1*0.999999 <= l2 ) ? true : false;
	}
	if( REGULAR ) { //regular if shorter distance
	  newgrid->putLoc2Glob123( m3, m2, n3, e3 );
	  newgrid->putLoc2Glob456( mm3, new_nodes( e, 4 ), nn3, e3 );
	  newgrid->putLoc2Glob123( m2, m3, m1, e4 );
	  newgrid->putLoc2Glob456( mm3, mm1, mm2, e4 );
	} else { // three edges drawn to midpoint of longest edge
	  newgrid->putLoc2Glob123( m1, n3, m3, e3 );
	  newgrid->putLoc2Glob456( mm3, nn3, mm1, e3 );
	  newgrid->putLoc2Glob123( n3, m1, m2, e4 );
	  newgrid->putLoc2Glob456( mm3, mm2, new_nodes( e, 4 ), e4 );
	}
	break;
//------------------------------------------------------
      default: 
//------------------------------------------------------
	//	fatalerrorFP("WavesGridFEHier:: refine","There must be 1-4 new elements.");
	printf("WavesGridFEHier:: refine There must be 1-4 new elements.\n");
      } //end of switch
    }
  } // end of quadratic elements (T6)

  //  pos_child( nel+1 ) = getCurrElmNo()+1; 
  pos_child( nel ) = getCurrElmNo()+1; 

// set subdomain types for the new grid copy the info from the parent elements
  if( ! onemat )
    //    for(e=1;e<=new_nel;e++)
    //      newgrid->subdomain_type(e) = subdomain_type(newgrid->getParent_eff(e));
    for(e=0;e<new_nel;e++)
      newgrid->setMaterialType(e,getMaterialType(newgrid->getParent_eff(e)));

  real x1,y1,x2,y2,x3,y3;

// newgrid->coord is initialized to 0

// compute positions and boundary indicators for the new midpoint nodes
  if( T6 ) // already done for T3 case
    //    for( e=1; e<=new_nel; e++) {
    for( e=0; e<new_nel; e++) {
      if( getNoChildren( newgrid->getParent_eff( e ) ) > 1 ) {
	// the element has at least one new node
	n1 = newgrid->loc2glob( e, 0 );
	n2 = newgrid->loc2glob( e, 1 );
	n3 = newgrid->loc2glob( e, 2 );
	n4 = newgrid->loc2glob( e, 3 );
	n5 = newgrid->loc2glob( e, 4 );
	n6 = newgrid->loc2glob( e, 5 );
	x1 = newgrid->coord(n1,0); y1 = newgrid->coord(n1,1);
	x2 = newgrid->coord(n2,0); y2 = newgrid->coord(n2,1);
	x3 = newgrid->coord(n3,0); y3 = newgrid->coord(n3,1);
	if( ! newgrid->coord(n4,0) ){ // to avoid computing things twice
	  newgrid->coord(n4,0) = (x1 + x2)/2;
	  newgrid->coord(n4,1) = (y1 + y2)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n1,b) == '1'&& nind(n2,b) == '1' )
	      nind(n4,b) = '1';
	      */
	}
	if( ! newgrid->coord(n5,0) ){
	  newgrid->coord(n5,0) = (x2 + x3)/2;
	  newgrid->coord(n5,1) = (y2 + y3)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n2,b) == '1' && nind(n3,b) == '1' )
	      nind(n5,b) = '1';
	      */
	}
	if( ! newgrid->coord(n6,0) ){
	  newgrid->coord(n6,0) = (x3 + x1)/2;
	  newgrid->coord(n6,1) = (y3 + y1)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n3,b) == '1' && nind(n1,b) == '1' )
	      nind(n6,b) = '1';
	      */
	}
      }
    }

  removeNewNodeInfo(); // remove the data showing the marked edges
  newgrid->setParentInfo( true );
  setChildrenInfo( true );
  //#ifdef ARRAY_RANGECHECK
  if( ! newgrid->checkElementOrientation() )
    printf("Wrong orientation.\n");
  if( ! checkChildrenInfo() ) 
    printf("children info failed.\n");
  if( ! newgrid->checkParentInfo() ) 
    printf("parent info failed.\n");
  //#endif
  return true; //everything went fine
}

void WavesGridFEHier:: putLoc2Glob1234 (const int glob_node1, const int glob_node2,
				   const int glob_node3, const int glob_node4,const int e)
{ nodpek(e,0) = glob_node1;    nodpek(e,1) = glob_node2;
nodpek(e,2) = glob_node3;    nodpek(e,3) = glob_node4;  }


//-----------------------------------------------------------------------------
bool WavesGridFEHier:: refine3D( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
{
// return true if refinement successful
  checkReadyForRefine ();  // check if the grid is ok

  /*  if( "ElmT4n3D" != elm_tp && "ElmT10n3D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: refine3D",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
  int b,e,i;

// Find the number of nodes and elements in the refined grid.
// The new nodes have been numbered before.
// A triangle which have n marked/numbered edges will be refined
// into n+2 elements if edge 6 is refined otherwise into n+1 elements

// The new_nodes structure is 0 where no edge refinement is taking place.

  int new_nel = nel;
  for( e=0; e<nel; e++ )
    if( new_nodes( e, 0 ) ) { // if edge 2-6 are refined then so are edge 1
      new_nel++;
      if( new_nodes( e, 1 ) )  new_nel++;
      if( new_nodes( e, 2 ) )  new_nel++;
      if( new_nodes( e, 3 ) )  new_nel++;
      if( new_nodes( e, 4 ) )  new_nel++;
      if( new_nodes( e, 5 ) ){ new_nel++; new_nel++; }
    }

// current node number was incremented during the node numbering process
  int new_nno = getCurrNodeNo()+1; //????

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  initRefinedGrid( newgrid, new_nno, new_nel ); // newsize and initialization

// copy the name of the boundary indicators
  /*
  for (i = 1; i <= nbind; i++) newgrid->putBoIndName( getBoIndName(i), i );

// copy the old nodes and indicators which will be present in refined grid

  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
  for( i=0; i<nno; i++ ){
    newgrid->coord(i,0) = coord(i,0); 
    newgrid->coord(i,1) = coord(i,1); 
    newgrid->coord(i,2) = coord(i,2); 
  }
  /*
  if ( nbind )
    for( i=1; i<=nno; i++ )
      for( b=1; b<=nbind; b++ )
	nind(i,b) = ind(i,b);
	*/
  int e1,e2,e3,e4,e5,e6,e7,e8; // numbers for new elements
  int n1,n2,n3,n4,n5,n6,n7,n8,n9,n10; // numbers for nodes in element
  int mid,m1,m2,node2,node3;

  setCurrElmNo(-1); // first element is given element number 0.

  //  bool T10 = ( "ElmT10n3D" == elm_tp ) ? true : false;
  bool T10 = ( maxnne == 10 ) ? true : false;

  //  if (T10) errorFP("WavesGridFEHier::refine3D","ElmT10n3D not yet implemented");
  if (T10) 
    printf("WavesGridFEHier::refine3D ElmT10n3D not yet implemented\n");

  //#ifdef ARRAY_RANGECHECK
  //  VecSimple(int) edgeorder;
  //  edgeorder.newsize(6);
  MV_Vector<int> edgeorder(6);
  //#endif

  int refcase;
  for( e=0; e < nel; e++ ) { // Refine the elements independently and set data
    pos_child( e ) = getCurrElmNo()+1; 
    n1 = loc2glob(e,0);  // these are needed for all cases
    n2 = loc2glob(e,1);
    n3 = loc2glob(e,2);
    n4 = loc2glob(e,3);
    n5 = new_nodes( e, 0 ); // always set if any refinement

    refcase = 1;
    if( n5 ){ // if not n5 then none of the others
      n6  = new_nodes( e, 1 );
      n7  = new_nodes( e, 2 );
      n8  = new_nodes( e, 3 );
      n9  = new_nodes( e, 4 );
      n10 = new_nodes( e, 5 );
      for(i=1;i<=6;i++){
	switch ( i ) { // set nodes on edge 1-6
	case 1: mid=n5;  m1=n1;  m2=n2; break;
	case 2: mid=n6;  m1=n2;  m2=n3; break;
	case 3: mid=n7;  m1=n3;  m2=n1; break;
	case 4: mid=n8;  m1=n1;  m2=n4; break;
	case 5: mid=n9;  m1=n2;  m2=n4; break;
	case 6: mid=n10; m1=n3;  m2=n4; break;
	}
	
	if( mid ) refcase++;  
	if( mid && ! newgrid->coord(mid,0) ){ // to avoid recomputing things
// Compute position and indicators of the new nodes at midpoints of edges 
	  //k	  if(max_aniso != 0.5)
	 //k   setNewNodeCoordinates(m1,m2,mid,newgrid);
	 //k else 
	  {
	    newgrid->coord(mid,0) = (coord(m1,0) + coord(m2,0))/2;
	    newgrid->coord(mid,1) = (coord(m1,1) + coord(m2,1))/2;
	    newgrid->coord(mid,2) = (coord(m1,2) + coord(m2,2))/2;
	  }
	  /*	  for( b=1; b<=nbind; b++ )
	    if(ind(m1,b) == '1' && ind(m2,b) == '1') // if both ends are set...
	      nind(mid,b) = '1'; // then set midpoint
	      */
	}
      }
    }
    //#ifdef ARRAY_RANGECHECK  // some tests in nonoptimized version
    getEdgeLengthOrder( edgeorder, e );
    if( edgeorder(0) != 0 ) 
      //      errorFP("WavesGridFEHier::refine3D","1 must be longest");
      printf("WavesGridFEHier::refine3D  1 must be longest\n");
    if( refcase > 1 && !n5 ) 
      //      errorFP("WavesGridFEHier::refine3D","n5 must be set");
      printf("WavesGridFEHier::refine3D n5 must be set\n");
    //#endif

    switch ( refcase ) {
//----------------------------------------------------------------
    case 1 : // no refinement: just copy the info to the new element
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
      newgrid->putLoc2Glob1234( n1, n2, n3, n4  , e1 );
      break;
//----------------------------------------------------------------
    case 2 :
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
      newgrid->putParent( e2 = incCurrElmNo(), e ); // e is parent to e2
//	cout<<e<<" 5"<<"\n";
      switch ( edge_face2(e)+1 ) {
      case 2:  newgrid->putLoc2Glob1234( n3, n2, n4, n5  , e1 );  break;
      case 5:  newgrid->putLoc2Glob1234( n2, n4, n3, n5  , e1 );  break;
      case 6:  newgrid->putLoc2Glob1234( n4, n3, n2, n5  , e1 );  break;
	//      default:  errorFP("WavesGridFEHier::refine3D","wrong number in edge_face2");
      default:  printf("WavesGridFEHier::refine3D wrong number in edge_face2\n");
      }
      switch ( edge_face3(e)+1 ) {
      case 3:  newgrid->putLoc2Glob1234( n1, n3, n4, n5  , e2 );  break;
      case 4:  newgrid->putLoc2Glob1234( n4, n1, n3, n5  , e2 );  break;
      case 6:  newgrid->putLoc2Glob1234( n3, n4, n1, n5  , e2 );  break;
	//      default:  errorFP("WavesGridFEHier::refine3D","wrong number in edge_face3"); 
      default:  cout<<"WavesGridFEHier::refine3D wrong number in edge_face3\n";
      }
      break;
//----------------------------------------------------------------
    case 3 :
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
      newgrid->putParent( e2 = incCurrElmNo(), e ); // e is parent to e2
      newgrid->putParent( e3 = incCurrElmNo(), e ); // e is parent to e3
      for( node2=1; node2<6; node2++ ) 
	if( new_nodes( e, node2) )  
	  break;
      switch ( node2+1 ) {
      case 2:
//	  cout<<e<<" 56"<<"\n";
	switch ( edge_face3(e)+1 ) {
	case 3:  newgrid->putLoc2Glob1234( n1, n3, n4, n5  , e1 );  break;
	case 4:  newgrid->putLoc2Glob1234( n4, n1, n3, n5  , e1 );  break;
	case 6:  newgrid->putLoc2Glob1234( n3, n4, n1, n5  , e1 );  break;
	  //	default:  errorFP("WavesGridFEHier::refine3D","wrong number in edge_face3");
      default:  cout<<"WavesGridFEHier::refine3D wrong number in edge_face3\n";
	}
	newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e2 );
	newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e3 );
	break;
      case 3:
//	  cout<<e<<" 57"<<"\n";
	switch ( edge_face2(e)+1 ) {
	case 2:  newgrid->putLoc2Glob1234( n3, n2, n4, n5  , e1 );  break;
	case 5:  newgrid->putLoc2Glob1234( n2, n4, n3, n5  , e1 );  break;
	case 6:  newgrid->putLoc2Glob1234( n4, n3, n2, n5  , e1 );  break;
	  //	default:  errorFP("WavesGridFEHier::refine3D","wrong number in edge_face2");
      default:  cout<<"WavesGridFEHier::refine3D wrong number in edge_face2\n";
	}
	newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e2 );
	newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e3 );
	break;
      case 4:
//	  cout<<e<<" 58"<<"\n";
	switch ( edge_face2(e)+1 ) {
	case 2:  newgrid->putLoc2Glob1234( n3, n2, n4, n5  , e1 );  break;
	case 5:  newgrid->putLoc2Glob1234( n2, n4, n3, n5  , e1 );  break;
	case 6:  newgrid->putLoc2Glob1234( n4, n3, n2, n5  , e1 );  break;
	  //	default:  errorFP("WavesGridFEHier::refine3D","wrong number in edge_face2");
      default:  cout<<"WavesGridFEHier::refine3D wrong number in edge_face2\n";
	}
	newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e2 );
	newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e3 );
	break;
      case 5:
//	  cout<<e<<" 59"<<"\n";
	switch ( edge_face3(e)+1 ) {
	case 3:  newgrid->putLoc2Glob1234( n1, n3, n4, n5  , e1 );  break;
	case 4:  newgrid->putLoc2Glob1234( n4, n1, n3, n5  , e1 );  break;
	case 6:  newgrid->putLoc2Glob1234( n3, n4, n1, n5  , e1 );  break;
	  //	default: errorFP("WavesGridFEHier::refine3D","wrong number in edge_face3");
      default:  cout<<"WavesGridFEHier::refine3D wrong number in edge_face2\n";
	}
	newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e2 );
	newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e3 );
	break;
      case 6:
//	  cout<<e<<" 510"<<"\n";
	newgrid->putParent( e4 = incCurrElmNo(), e ); // e is parent to e4
	newgrid->putLoc2Glob1234( n1, n3, n10, n5  , e1 );
	newgrid->putLoc2Glob1234( n2, n4, n10, n5  , e2 );
	newgrid->putLoc2Glob1234( n3, n2, n10, n5  , e3 );
	newgrid->putLoc2Glob1234( n4, n1, n10, n5  , e4 );
	break;
	//      default: errorFP("WavesGridFEHier::refine3D",
	//		       "case 3 did not find second node");
      default: cout<<"WavesGridFEHier::refine3D case 3 did not find second node\n";
      }
      break;
//----------------------------------------------------------------
    case 4 : // find the two edges except edge 1 to be refined
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e ); // e is parent to e1
      newgrid->putParent( e2 = incCurrElmNo(), e ); // e is parent to e2
      newgrid->putParent( e3 = incCurrElmNo(), e ); // e is parent to e3
      newgrid->putParent( e4 = incCurrElmNo(), e ); // e is parent to e4
      if( n10 )
	newgrid->putParent( e5 = incCurrElmNo(), e ); // e is parent to e5
// first find if second node is n6 n7 n8 or n9
// there after find third node which comes after second node
      for( node2=1; node2<6; node2++ ) 
	if( new_nodes( e, node2 ) )  
	  break;
      switch ( node2+1 ) {
      case 2:
	for( node3=2; node3<6; node3++ ) 
	  if( new_nodes(e,node3) )  
	    break;
	switch ( node3+1 ) {
	case 3: // case 5 6 7 
//	    cout<<e<<" 567"<<"\n";
	  newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e1 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e2 );
	  newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e3 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e4 );
	  break;
	case 4:// case 5 6 8
//	    cout<<e<<" 568"<<"\n";
	  newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e1 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e2 );
	  newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e3 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e4 );
	  break;
	case 5:// case 5 6 9
//	    cout<<e<<" 569"<<"\n";
	  switch ( edge_face3(e)+1 ) {
	  case 3:  newgrid->putLoc2Glob1234( n1, n3, n4, n5  , e1 );  break;
	  case 4:  newgrid->putLoc2Glob1234( n4, n1, n3, n5  , e1 );  break;
	  case 6:  newgrid->putLoc2Glob1234( n3, n4, n1, n5  , e1 );  break;
	    //	  default:  errorFP("WavesGridFEHier::refine3D","wrong num in edge_face3");
	  default:  cout<<"WavesGridFEHier::refine3D wrong num in edge_face3\n";
	  }
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e2 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e3 );
	    newgrid->putLoc2Glob1234( n6, n4, n5, n9  , e4 );
	  } else {
	    newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e3 );
	    newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e4 );
	  }
	  break;
	case 6: // case 5 6 10
//	    cout<<e<<" 5610"<<"\n";
	  newgrid->putLoc2Glob1234( n1, n3, n10, n5  , e1 ); // common elms
	  newgrid->putLoc2Glob1234( n4, n1, n10, n5  , e2 );
	  newgrid->putLoc2Glob1234( n5, n3, n10, n6  , e3 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e4 );
	    newgrid->putLoc2Glob1234( n4, n6, n5, n10 , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n10  , e4 );
	    newgrid->putLoc2Glob1234( n10, n2, n5, n6  , e5 );
	  }
	  break;
	default:
	  //	  errorFP("WavesGridFEHier:: refine3D","node 5 6,did not find third node");
	  cout<<"WavesGridFEHier:: refine3D node 5 6,did not find third node\n";
	}
	break;
      case 3:
	for( node3=3; node3<6; node3++ ) 
	  if( new_nodes( e, node3 ) )  
	    break;
	switch(node3+1){
	case 4: // case 5 7 8
//	    cout<<e<<" 578"<<"\n";
	  switch ( edge_face2(e)+1 ) {
	  case 2:  newgrid->putLoc2Glob1234( n3, n2, n4, n5  , e1 );  break;
	  case 5:  newgrid->putLoc2Glob1234( n2, n4, n3, n5  , e1 );  break;
	  case 6:  newgrid->putLoc2Glob1234( n4, n3, n2, n5  , e1 );  break;
	  default:  
	    //	    errorFP("WavesGridFEHier::refine3D","wrong number in edge_face2");
	    cout<<"WavesGridFEHier::refine3D wrong number in edge_face2\n";
	  }
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e2 ); // common elms
	  if ( edge_face3(e) == 2 ) {
	    newgrid->putLoc2Glob1234( n4, n7, n5, n8  , e3 );
	    newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e4 );
	  } else {
	    newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e3 );
	    newgrid->putLoc2Glob1234( n8, n3, n5, n7  , e4 );
	  }
	  break;
	case 5: // case 5 7 9
//	    cout<<e<<" 579"<<"\n";
	  newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e1 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e2 );
	  newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e3 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e4 );
	  break;
	case 6: // case 5 7 10
//	    cout<<e<<" 5710"<<"\n";
	  newgrid->putLoc2Glob1234( n2, n4, n10, n5  , e1 );
	  newgrid->putLoc2Glob1234( n3, n2, n10, n5  , e2 );
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10  , e3 );
	  if ( edge_face3(e) == 2 ) {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e4 );
	    newgrid->putLoc2Glob1234( n7, n4, n5, n10 , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n10  , e4 );
	    newgrid->putLoc2Glob1234( n1, n10, n5, n7  , e5 );
	  }
	  break;
	default:
	  //	  errorFP("WavesGridFEHier:: refine3D","node 5 7,did not find third node");
	  cout<<"WavesGridFEHier:: refine3D node 5 7,did not find third node\n";
	}
	break;
      case 4:
	if( n9 ) { // case 5 8 9 
//	    cout<<e<<" 589"<<"\n";
	  newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e1 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e2 );
	  newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e3 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e4 );
	} else {	  // case 5 8 10
//	    cout<<e<<" 5810"<<"\n";
	  newgrid->putLoc2Glob1234( n2, n4, n10, n5  , e1 );
	  newgrid->putLoc2Glob1234( n3, n2, n10, n5  , e2 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8  , e3 );
	  if ( edge_face3(e) == 3 ) {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e4 );
	    newgrid->putLoc2Glob1234( n3, n8, n5, n10 , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n10 , e4 );
	    newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e5 );
	  }
	}
	break;
      case 5: // case 5 9 10
//	  cout<<e<<" 5910"<<"\n";
	newgrid->putLoc2Glob1234( n1, n3, n10, n5  , e1 );
	newgrid->putLoc2Glob1234( n4, n1, n10, n5  , e2 );
	newgrid->putLoc2Glob1234( n4, n5, n10, n9  , e3 );
	if( edge_face2(e) == 4 ) {
	  newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e4 );
	  newgrid->putLoc2Glob1234( n9, n3, n5, n10 , e5 );
	} else {
	  newgrid->putLoc2Glob1234( n2, n3, n5, n10 , e4 );
	  newgrid->putLoc2Glob1234( n10, n2, n9, n5 , e5 );
	}
	break;
      case 6:
	//      default: errorFP("WavesGridFEHier::refine3D",
	//		       "case 4 did not find second node");
      default: cout<<"WavesGridFEHier::refine3D case 4 did not find second node\n";
      }
      break;
//----------------------------------------------------------------
    case 5: // find the two edges which are not refined
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e );
      newgrid->putParent( e2 = incCurrElmNo(), e );
      newgrid->putParent( e3 = incCurrElmNo(), e );
      newgrid->putParent( e4 = incCurrElmNo(), e );
      newgrid->putParent( e5 = incCurrElmNo(), e );
      if ( n10 )
	newgrid->putParent( e6 = incCurrElmNo(), e ); // e is parent to e6
      for( node2=1; node2<6; node2++ ) 
	if( ! new_nodes(e,node2) )  
	  break;
      switch ( node2+1 ) {
      case 2:
	for( node3=2; node3<6; node3++ ) 
	  if( ! new_nodes( e, node3 ) ) 
	    break;
	switch ( node3+1 ) {
	case 3: // not node 6 and 7 => case 5 8 9 10
//	    cout<<e<<" 58910"<<"\n";
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8  , e1 );
	  if ( edge_face3(e) == 3 ) {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e2 );
	    newgrid->putLoc2Glob1234( n3, n8, n5, n10 , e3 );
	  } else {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n10 , e2 );
	    newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e3 );
	  }
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9  , e4 );
	  if( edge_face2(e) == 4 ) {
	    newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e5 );
	    newgrid->putLoc2Glob1234( n9, n3, n5, n10 , e6 );
	  } else {
	    newgrid->putLoc2Glob1234( n2, n3, n5, n10 , e5 );
	    newgrid->putLoc2Glob1234( n10, n2, n9, n5 , e6 );
	  }
	  break;
	case 4: // not node 6 and 8 => case 5 7 9 10
//	    cout<<e<<" 57910"<<"\n";
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10  , e1 );
	  if ( edge_face3(e) == 2 ) {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e2 );
	    newgrid->putLoc2Glob1234( n4, n5, n7, n10 , e3 );
	  } else {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n10  , e2 );
	    newgrid->putLoc2Glob1234( n1, n10, n5, n7  , e3 );
	  }
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9  , e4 );
	  if( edge_face2(e) == 4 ) {
	    newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e5 );
	    newgrid->putLoc2Glob1234( n3, n5, n9, n10 , e6 );
	  } else {
	    newgrid->putLoc2Glob1234( n2, n3, n5, n10 , e5 );
	    newgrid->putLoc2Glob1234( n10, n2, n9, n5 , e6 );
	  }
	  break;
	case 5: // not node 6 and 9 => case 5 7 8 10
//	    cout<<e<<" 57810"<<"\n";
	  newgrid->putLoc2Glob1234( n2, n4, n10, n5  , e1 );
	  newgrid->putLoc2Glob1234( n3, n2, n10, n5  , e2 );
	  switch ( edge_face3(e)+1 ) {
	  case 3:
	    newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e3 );
	    newgrid->putLoc2Glob1234( n7, n4, n5, n10 , e4 );
	    newgrid->putLoc2Glob1234( n4, n7, n5, n8  , e5 );
	    newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e6 );
	    break;
	  case 4:
	    newgrid->putLoc2Glob1234( n8, n3, n5, n7  , e3 );
	    newgrid->putLoc2Glob1234( n3, n8, n5, n10 , e4 );
	    newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e5 );
	    newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e6 );
	    break;
	  case 6:
	    newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e3 );
	    newgrid->putLoc2Glob1234( n1, n10, n5, n7 , e4 );
	    newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e5 );
	    newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e6 );
	    break;
	    //	  default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face3");
	    cout<<"WavesGridFEHier:: refine3D  illegal no in edge_face3\n";
	  }
	  break;
	case 6: // not node 6 and 10 => case 5 7 8 9
//	    cout<<e<<" 5789"<<"\n";
	  newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e1 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e2 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e3 );// common elms
	  if ( edge_face3(e) == 2 ) {
	    newgrid->putLoc2Glob1234( n5, n4, n7, n8  , e4 );
	    newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e4 );
	    newgrid->putLoc2Glob1234( n8, n3, n5, n7  , e5 );
	  }
	  break;
	default:
	  //	  errorFP("WavesGridFEHier:: refine3D","node 6,did not find second node");
	  printf("WavesGridFEHier:: refine3D node 6,did not find second node e=%d a\n",e);
	}
	break;
      case 3:
	for( node3=3; node3<=5; node3++ ) //!!!!!!!!!!!!
	  if( ! new_nodes( e, node3 ) ) 
	    break;
	switch ( node3+1 ) {
	case 4: // not node 7 and 8 => case 5 6 9 10
//	    cout<<e<<" 56910"<<"\n";
	  newgrid->putLoc2Glob1234( n1, n3, n10, n5  , e1 );
	  newgrid->putLoc2Glob1234( n4, n1, n10, n5  , e2 );
	  switch ( edge_face2(e)+1 ) {
	  case 2:
	    newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e3 );
	    newgrid->putLoc2Glob1234( n6, n4, n5, n9  , e4 );
	    newgrid->putLoc2Glob1234( n4, n6, n5, n10 , e5 );
	    newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e6 );
	    break;
	  case 5:
	    newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e3 );
	    newgrid->putLoc2Glob1234( n9, n3, n5, n10 , e4 );
	    newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e5 );
	    newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e6 );
	    break;
	  case 6:
	    newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e3 );
	    newgrid->putLoc2Glob1234( n10, n2, n5, n6 , e4 );
	    newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e5 );
	    newgrid->putLoc2Glob1234( n2, n10, n5, n9 , e6 );
	    break;
	    //	  default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face2");
	  default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face2\n";
	  }
	  break;
	case 5: // not node 7 and 9 => case 5 6 8 10
//	    cout<<e<<" 56810"<<"\n";
	  newgrid->putLoc2Glob1234( n5, n3, n10, n6  , e1 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e2 );
	    newgrid->putLoc2Glob1234( n4, n6, n5, n10 , e3 );
	  } else {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n10  , e2 );
	    newgrid->putLoc2Glob1234( n10, n2, n5, n6  , e3 );
	  }
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8  , e4 );
	  if ( edge_face3(e) == 3 ) {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e5 );
	    newgrid->putLoc2Glob1234( n3, n8, n5, n10 , e6 );
	  } else {
	    newgrid->putLoc2Glob1234( n3, n1, n5, n10 , e5 );
	    newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e6 );
	  }
	  break;
	case 6: // not node 7 and 10 => case 5 6 8 9
//	    cout<<e<<" 5689"<<"\n";
	  newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e1 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e2 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e3 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e4 );
	    newgrid->putLoc2Glob1234( n4, n5, n6, n9  , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e4 );
	    newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e5 );
	  }
	  break;
	default:
	  //	  errorFP("WavesGridFEHier:: refine3D","node 6,did not find second node");
	  printf("WavesGridFEHier:: refine3D node 6,did not find second node e=%d b\n",e);
	  printf("newnodes %d %d %d %d %d %d\n",new_nodes(e,0),new_nodes(e,1),new_nodes(e,2),new_nodes(e,3),new_nodes(e,4),new_nodes(e,5));
	  printf("node3+1=%d\n\n",node3+1);
	}
	break;
      case 4:
	if( ! n9 ) { // not node 8 and 9 => case 5 6 7 10
//	    cout<<e<<" 56710"<<"\n";
	  newgrid->putLoc2Glob1234( n5, n3, n10, n6  , e1 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e2 );
	    newgrid->putLoc2Glob1234( n4, n6, n5, n10 , e3 );
	  } else {
	    newgrid->putLoc2Glob1234( n4, n2, n5, n10  , e2 );
	    newgrid->putLoc2Glob1234( n10, n2, n5, n6  , e3 );
	  }
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10  , e4 );
	  if ( edge_face3(e) == 2 ) {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e5 );
	    newgrid->putLoc2Glob1234( n7, n4, n5, n10 , e6 );
	  } else {
	    newgrid->putLoc2Glob1234( n1, n4, n5, n10  , e5 );
	    newgrid->putLoc2Glob1234( n1, n10, n5, n7  , e6 );
	  }
	} else { // not node 8 and 10 => case 5 6 7 9
//	  cout<<e<<" 5679"<<"\n";
	  newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e1 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e2 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e3 );
	  if ( edge_face2(e) == 1 ) {
	    newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e4 );
	    newgrid->putLoc2Glob1234( n6, n4, n5, n9  , e5 );
	  } else {
	    newgrid->putLoc2Glob1234( n5, n3, n9, n6  , e4 );
	    newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e5 );
	  }
	}
	break;
      case 5: // not node 9 and 10 => case 5 6 7 8
//	  cout<<e<<" 5678"<<"\n";
	newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e1 );
	newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e2 );
	newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e3 );
	if ( edge_face3(e) == 2 ) {
	  newgrid->putLoc2Glob1234( n4, n7, n5, n8  , e4 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e5 );
	} else {
	  newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e4 );
	  newgrid->putLoc2Glob1234( n3, n5, n8, n7  , e5 );
	}
	break;
      case 6:
	//      default: errorFP("WavesGridFEHier::refine3D",
	//		       "case 5 did not find second node");
      default: cout<<"WavesGridFEHier::refine3D case 5 did not find second node\n";
      }
      break;
//----------------------------------------------------------------
    case 6: // find the edge which not is going to be refined
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e );
      newgrid->putParent( e2 = incCurrElmNo(), e );
      newgrid->putParent( e3 = incCurrElmNo(), e );
      newgrid->putParent( e4 = incCurrElmNo(), e );
      newgrid->putParent( e5 = incCurrElmNo(), e );
      newgrid->putParent( e6 = incCurrElmNo(), e ); 
      if ( n10 )
	newgrid->putParent( e7 = incCurrElmNo(), e );
      for( node2=1; node2<6; node2++ ) 
	if( ! new_nodes( e, node2 ) ) 
	  break;
      switch ( node2+1 ) {
      case 2: // not node 6 => case 5 7 8 9 10
//	  cout<<e<<" 578910"<<"\n";
	newgrid->putLoc2Glob1234( n4, n5, n10, n9  , e1 );
	if( edge_face2(e) == 4 ) {
	  newgrid->putLoc2Glob1234( n2, n3, n5, n9  , e2 );
	  newgrid->putLoc2Glob1234( n3, n5, n9, n10 , e3 );
        } else {
	  newgrid->putLoc2Glob1234( n2, n3, n5, n10 , e2 );
	  newgrid->putLoc2Glob1234( n10, n2, n9, n5 , e3 );
	}
	switch ( edge_face3(e)+1 ) {
	case 3:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e4 );
	  newgrid->putLoc2Glob1234( n4, n5, n7, n10 , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n7, n8  , e6 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e7 );
	  break;
	case 4:
	  newgrid->putLoc2Glob1234( n8, n3, n5, n7  , e4 );
	  newgrid->putLoc2Glob1234( n3, n8, n5, n10 , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e6 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e7 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e4 );
	  newgrid->putLoc2Glob1234( n1, n10, n5, n7 , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e6 );
	  newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e7 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face3");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face3\n";
	}
	break;
      case 3: // not node 7 => case 5 6 8 9 10
//	  cout<<e<<" 568910"<<"\n";
	newgrid->putLoc2Glob1234( n5, n4, n10, n8  , e1 );
	if ( edge_face3(e) == 3 ) {
	  newgrid->putLoc2Glob1234( n3, n1, n5, n8  , e2 );
	  newgrid->putLoc2Glob1234( n5, n3, n8, n10 , e3 );
	} else {
	  newgrid->putLoc2Glob1234( n3, n1, n5, n10 , e2 );
	  newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e3 );
	}
	switch ( edge_face2(e)+1 ) {
	case 2:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e4 );
	  newgrid->putLoc2Glob1234( n4, n5, n6, n9  , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n6, n10 , e6 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e7 );
	  break;
	case 5:
	  newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e4 );
	  newgrid->putLoc2Glob1234( n9, n3, n5, n10 , e5 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e6 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e7 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e4 );
	  newgrid->putLoc2Glob1234( n10, n2, n5, n6 , e5 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e6 );
	  newgrid->putLoc2Glob1234( n2, n10, n5, n9 , e7 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face2");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face2\n";
	}
	break;
      case 4: // not node 8 => case 5 6 7 9 10
//	  cout<<e<<" 567910"<<"\n";
	newgrid->putLoc2Glob1234( n5, n3, n7, n10  , e1 );
	if ( edge_face3(e) == 2 ) {
	  newgrid->putLoc2Glob1234( n1, n4, n5, n7  , e2 );
	  newgrid->putLoc2Glob1234( n7, n4, n5, n10 , e3 );
	} else {
	  newgrid->putLoc2Glob1234( n1, n4, n5, n10  , e2 );
	  newgrid->putLoc2Glob1234( n1, n10, n5, n7  , e3 );
	}
	switch ( edge_face2(e)+1 ) {
	case 2:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e4 );
	  newgrid->putLoc2Glob1234( n6, n4, n5, n9  , e5 );
	  newgrid->putLoc2Glob1234( n4, n6, n5, n10 , e6 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e7 );
	  break;
	case 5:
	  newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e4 );
	  newgrid->putLoc2Glob1234( n9, n3, n5, n10 , e5 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e6 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e7 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e4 );
	  newgrid->putLoc2Glob1234( n10, n2, n5, n6 , e5 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e6 );
	  newgrid->putLoc2Glob1234( n2, n10, n5, n9 , e7 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face2");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face2\n";
	}
	break;
      case 5:   // not node 9 => case 5 6 7 8 10
//	  cout<<e<<" 567810"<<"\n";
	newgrid->putLoc2Glob1234( n5, n3, n10, n6  , e1 );
	if ( edge_face2(e) == 1 ) {
	  newgrid->putLoc2Glob1234( n4, n2, n5, n6  , e2 );
	  newgrid->putLoc2Glob1234( n5, n4, n6, n10 , e3 );
	} else {
	  newgrid->putLoc2Glob1234( n4, n2, n5, n10  , e2 );
	  newgrid->putLoc2Glob1234( n10, n2, n5, n6  , e3 );
	}
	switch ( edge_face3(e)+1 ) {
	case 3:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e4 );
	  newgrid->putLoc2Glob1234( n7, n4, n5, n10 , e5 );
	  newgrid->putLoc2Glob1234( n4, n7, n5, n8  , e6 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e7 );
	  break;
	case 4:
	  newgrid->putLoc2Glob1234( n3, n5, n8, n7  , e4 );
	  newgrid->putLoc2Glob1234( n5, n3, n8, n10 , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e6 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e7 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e4 );
	  newgrid->putLoc2Glob1234( n1, n10, n5, n7 , e5 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e6 );
	  newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e7 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face3");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face3\n";
	}
	break;
      case 6: // not node 10 => case 5 6 7 8 9
//	  cout<<e<<" 56789"<<"\n";
	newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e1 );
	if ( edge_face2(e) == 1 ) {
	  newgrid->putLoc2Glob1234( n3, n4, n5, n6  , e2 );
	  newgrid->putLoc2Glob1234( n6, n4, n5, n9  , e3 );
	} else {
	  newgrid->putLoc2Glob1234( n3, n9, n5, n6  , e2 );
	  newgrid->putLoc2Glob1234( n3, n4, n5, n9  , e3 );
	}
	newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e4 );
	if ( edge_face3(e) == 2 ) {
	  newgrid->putLoc2Glob1234( n4, n7, n5, n8  , e5 );
	  newgrid->putLoc2Glob1234( n4, n3, n5, n7  , e6 );
	} else {
	  newgrid->putLoc2Glob1234( n4, n3, n5, n8  , e5 );
	  newgrid->putLoc2Glob1234( n8, n3, n5, n7  , e6 );
	}
	break;
	//      default: errorFP("WavesGridFEHier::refine3D",
	//		       "case 6 did not find second node");
      default: cout<<"WavesGridFEHier::refine3D  case 6 did not find second node\n";
      }
      break;
//----------------------------------------------------------------
    case 7 : // refine the element into 8 elements
//----------------------------------------------------------------
      newgrid->putParent( e1 = incCurrElmNo(), e );
      newgrid->putParent( e2 = incCurrElmNo(), e );
      newgrid->putParent( e3 = incCurrElmNo(), e );
      newgrid->putParent( e4 = incCurrElmNo(), e );
      newgrid->putParent( e5 = incCurrElmNo(), e );
      newgrid->putParent( e6 = incCurrElmNo(), e );
      newgrid->putParent( e7 = incCurrElmNo(), e );
      newgrid->putParent( e8 = incCurrElmNo(), e );
      
      if ( getRegularUniformRef() ) {
	// regular method cut of four tetrahedra, choice of diagonal
	// only used for uniform refinement with method 6
	newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e1 );
	newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e2 );
	newgrid->putLoc2Glob1234( n7, n6, n3, n10 , e3 );
	newgrid->putLoc2Glob1234( n8, n9, n10, n4 , e4 );
	switch ( getShortestDiagonal(e,newgrid->coord) ) {
	case 1: //	case five2ten:
	  newgrid->putLoc2Glob1234( n6, n7, n5, n10 , e5 );
	  newgrid->putLoc2Glob1234( n5, n8, n10, n7 , e6 );
	  newgrid->putLoc2Glob1234( n8, n9, n5, n10 , e7 );
	  newgrid->putLoc2Glob1234( n5, n6, n10, n9 , e8 );
//	    cout<<e<<" 5-10"<<"\n";
	  break;
	case 2: //	case six2eight:
	  newgrid->putLoc2Glob1234( n9, n8, n6, n5  , e5 );
	  newgrid->putLoc2Glob1234( n9, n8, n10, n6 , e6 );
	  newgrid->putLoc2Glob1234( n7, n6, n10, n8 , e7 );
	  newgrid->putLoc2Glob1234( n7, n6, n8, n5  , e8 );
//	    cout<<e<<" 6-8"<<"\n";
	  break;
	case 3: //	case seven2nine:
	  newgrid->putLoc2Glob1234( n7, n6, n9, n5  , e5 );
	  newgrid->putLoc2Glob1234( n7, n6, n10, n9 , e6 );
	  newgrid->putLoc2Glob1234( n9, n8, n10, n7 , e7 );
	  newgrid->putLoc2Glob1234( n9, n8, n7, n5  , e8 );
//	    cout<<e<<" 7-9"<<"\n";
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","wrong diagonal");
	default: cout<<"WavesGridFEHier:: refine3D wrong diagonal\n";
	}
	//       	if(e==1466){
	//	  printf(" shortestdiag=%i\n",getShortestDiagonal(e,newgrid->coord));
	//	  printf(" n5=%i,n6=%i,n7=%i,n8=%i,n9=%i,n10=%i\n",n5,n6,n7,n8,n9,n10);
	//	  printf(" e5=%i,e6=%i,e7=%i,e8=%i\n",e5,e6,e7,e8);
	//	}
      } else { // bisection method
	switch ( edge_face2(e)+1 ) {
	case 2:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e1 );
	  newgrid->putLoc2Glob1234( n4, n5, n6, n9  , e2 );
	  newgrid->putLoc2Glob1234( n5, n4, n6, n10 , e3 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e4 );
	  break;
	case 5:
	  newgrid->putLoc2Glob1234( n5, n3, n9, n6  , e1 );
	  newgrid->putLoc2Glob1234( n3, n5, n9, n10 , e2 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e3 );
	  newgrid->putLoc2Glob1234( n5, n2, n6, n9  , e4 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n3, n5, n6, n10 , e1 );
	  newgrid->putLoc2Glob1234( n10, n2, n5, n6 , e2 );
	  newgrid->putLoc2Glob1234( n4, n5, n10, n9 , e3 );
	  newgrid->putLoc2Glob1234( n2, n10, n5, n9 , e4 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face2");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face2\n";
	}
	switch ( edge_face3(e)+1 ) {
	case 3:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e5 );
	  newgrid->putLoc2Glob1234( n4, n5, n7, n10 , e6 );
	  newgrid->putLoc2Glob1234( n5, n4, n7, n8  , e7 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e8 );
	  break;
	case 4:
	  newgrid->putLoc2Glob1234( n3, n5, n8, n7  , e5 );
	  newgrid->putLoc2Glob1234( n5, n3, n8, n10 , e6 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e7 );
	  newgrid->putLoc2Glob1234( n1, n5, n7, n8  , e8 );
	  break;
	case 6:
	  newgrid->putLoc2Glob1234( n5, n3, n7, n10 , e5 );
	  newgrid->putLoc2Glob1234( n1, n10, n5, n7 , e6 );
	  newgrid->putLoc2Glob1234( n5, n4, n10, n8 , e7 );
	  newgrid->putLoc2Glob1234( n10, n1, n5, n8 , e8 );
	  break;
	  //	default: errorFP("WavesGridFEHier:: refine3D","illegal no in edge_face3");
	default: cout<<"WavesGridFEHier:: refine3D illegal no in edge_face2\n";
	}
      }
      break; 
//------------------------------------------------------
    default : 
//------------------------------------------------------
      //      fatalerrorFP("WavesGridFEHier:: refine3D",
      //		   "\nMust be 1-7 new elements, now there are %d.",refcase);
      cout<<"WavesGridFEHier:: refine3D \nMust be 1-7 new elements, now there?\n";
    } //end of switch
  }

  //  pos_child( nel + 1 ) = getCurrElmNo() + 1; 
  pos_child( nel ) = getCurrElmNo() + 1; 

// copy the subdomain types for the new grid from the parent elements:
// subdomain_type is newsizeed i GridFE.redim if onesbd false
  if( ! onemat ){
    //    for(e=1;e<=new_nel;e++)
    //      newgrid->subdomain_type(e) = subdomain_type(newgrid->getParent_eff(e));
    cout<<"not onematerial\n"<<flush;
    for(e=0;e<new_nel;e++)
      newgrid->setMaterialType(e,getMaterialType(newgrid->getParent_eff(e)));
  }
  else
    cout<<"onematerial\n"<<flush;

  real x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4;
  if ( T10 ) {
    for( e=0; e<new_nel; e++ ) {
// compute positions and boundary indicators for the new midpoint nodes
      if( getNoChildren( newgrid->getParent_eff( e ) ) > 1 ){
	// the element has at least one new node
	n1 = newgrid->loc2glob( e, 0 );
	n2 = newgrid->loc2glob( e, 1 );
	n3 = newgrid->loc2glob( e, 2 );
	n4 = newgrid->loc2glob( e, 3 );
	n5 = newgrid->loc2glob( e, 4 );
	n6 = newgrid->loc2glob( e, 5 );
	n7 = newgrid->loc2glob( e, 6 );
	n8 = newgrid->loc2glob( e, 7 );
	n9 = newgrid->loc2glob( e, 8 );
	n10= newgrid->loc2glob( e, 9 );
	
	x1 = newgrid->coord(n1,0); 
	y1 = newgrid->coord(n1,1);
	z1 = newgrid->coord(n1,2);
	x2 = newgrid->coord(n2,0);
	y2 = newgrid->coord(n2,1);
	z2 = newgrid->coord(n2,2);
	x3 = newgrid->coord(n3,0);
	y3 = newgrid->coord(n3,1);
	z3 = newgrid->coord(n3,2);
	x4 = newgrid->coord(n4,0);
	y4 = newgrid->coord(n4,1); 
	z4 = newgrid->coord(n4,2);

	if( ! newgrid->coord(n5,0) ) { // to avoid computing things twice
	  newgrid->coord(n5,0) = (x1 + x2)/2;
	  newgrid->coord(n5,1) = (y1 + y2)/2;
	  newgrid->coord(n5,2) = (z1 + z2)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n1,b) == '1' && nind(n2,b) == '1' ) // ends set..
	      nind(n5,b) = '1'; // then set midpoint
	      */
	}
	if( ! newgrid->coord(n6,0) ) { 
	  newgrid->coord(n6,0) = (x2 + x3)/2;
	  newgrid->coord(n6,1) = (y2 + y3)/2;
	  newgrid->coord(n6,2) = (z2 + z3)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n2,b) == '1' && nind(n3,b) == '1' )
	      nind(n6,b) = '1';
	      */
	}
	if( ! newgrid->coord(n7,0) ) {
	  newgrid->coord(n7,0) = (x3 + x1)/2;
	  newgrid->coord(n7,1) = (y3 + y1)/2;
	  newgrid->coord(n7,2) = (z3 + z1)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n3,b) == '1' && nind(n1,b) == '1' )
	      nind(n7,b) = '1';
	      */
	}
	if( ! newgrid->coord(n8,0) ) {
	  newgrid->coord(n8,0) = (x4 + x1)/2;
	  newgrid->coord(n8,1) = (y4 + y1)/2;
	  newgrid->coord(n8,2) = (z4 + z1)/2;
	  /*      for( b=1; b<=nbind; b++ )
	    if( nind(n4,b) == '1' && nind(n1,b) == '1' )
	      nind(n8,b) = '1';
	      */
	}
	if( ! newgrid->coord(n9,0) ) {
	  newgrid->coord(n9,0) = (x4 + x2)/2;
	  newgrid->coord(n9,1) = (y4 + y2)/2;
	  newgrid->coord(n9,2) = (z4 + z2)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n4,b) == '1' && nind(n2,b) == '1' )
	      nind(n9,b) = '1';
	      */
	}
	if( ! newgrid->coord(n10,0) ) {
	  newgrid->coord(n10,0) = (x4 + x3)/2;
	  newgrid->coord(n10,1) = (y4 + y3)/2;
	  newgrid->coord(n10,2) = (z4 + z3)/2;
	  /*	  for( b=1; b<=nbind; b++ )
	    if( nind(n4,b) == '1' && nind(n3,b) == '1' )
	      nind(n10,b) = '1';
	      */
	}
      }
    }
  }

  removeNewNodeInfo(); // remove the data showing the marked edges
  //  edge_face2.newsize(0);
  //  edge_face3.newsize(0);
  newgrid->setParentInfo( true );
  setChildrenInfo( true );

//  cout<<"end1 refine3D\n"<<flush;
//#ifdef ARRAY_RANGECHECK
  if( ! newgrid->checkElementOrientation() ){
    cout<<"WavesGridFEHier:: refine3D\n";
    exit(-1);
  }
//  cout<<"end2 refine3D\n"<<flush;
  if( !checkChildrenInfo() ) 
    cout<<"children info failed."<<"\n";
//  cout<<"end3 refine3D\n"<<flush;
  if( !newgrid->checkParentInfo() ) 
    cout<<"parent info failed."<<"\n";
//#endif
//  cout<<"end4 refine3D\n"<<flush;
  //  for(e=0;e<newgrid->nel;e++)
  //    cout<<"parent "<<e<<" "<<newgrid->parent(e)<<"\n"<<flush;

  //  newgrid->fixCurvesAndSurfaces();
  //  newgrid->projectNodesOnCurvesAndSurfaces();

  //  newgrid->fixCurvesAndSurfacesParameters ();
  //  newgrid->projectNodesWithParameters ();

  //  newgrid->printBoundary("movedrefsrf.mtv",1);
  //  newgrid->printBoundary("movedrefcrv.mtv",2);

  return true; // everything went fine
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeInteriorIndicators()
//-----------------------------------------------------------------------------
{
  if( ! isNeighborsComputed() ) // Is the element-element info computed?
    initNeighbor(); 

  //  VecSimple(int) inner( nno );
  //  inner.fill( 1 ); // set boundary nodes to zero in the inner vector
  MV_Vector<int> inner( nno,1 );
  //  inner.fill( 1 ); // set boundary nodes to zero in the inner vector

  int i,b,e;
  switch ( nsd ) {
  case 1:
    for( e=0; e<nel; e++ ) {
      if( getElmNeighbor(e,0)==-1 ) inner(loc2glob(e,0)) = 0;
      if( getElmNeighbor(e,1)==-1 ) inner(loc2glob(e,1)) = 0;
    }
    break;
  case 2:
    for( e=0; e<nel; e++ ) {
      if( getElmNeighbor(e,0)==-1 )
	inner(loc2glob(e,1)) = inner(loc2glob(e,0)) = 0;
      if( getElmNeighbor(e,1)==-1 )
	inner(loc2glob(e,2)) = inner(loc2glob(e,1)) = 0;
      if( getElmNeighbor(e,2)==-1 )
	inner(loc2glob(e,2)) = inner(loc2glob(e,0)) = 0;
    }
    //    if( "ElmT6n2D" == elm_tp )
    if( maxnne == 6 )
      for( e=0; e<nel; e++ ) {
	if( getElmNeighbor(e,0)==-1 ) inner(loc2glob(e,3)) = 0;
	if( getElmNeighbor(e,1)==-1 ) inner(loc2glob(e,4)) = 0;
	if( getElmNeighbor(e,2)==-1 ) inner(loc2glob(e,5)) = 0;
      }
    break;
  case 3:
    for( e=0; e<nel; e++ ) {
      if( getElmNeighbor(e,0)==-1 )
	inner(loc2glob(e,3)) = inner(loc2glob(e,1)) = inner(loc2glob(e,0)) = 0;
      if( getElmNeighbor(e,1)==-1 )
	inner(loc2glob(e,3)) = inner(loc2glob(e,2)) = inner(loc2glob(e,1)) = 0;
      if( getElmNeighbor(e,2)==-1 )
	inner(loc2glob(e,3)) = inner(loc2glob(e,2)) = inner(loc2glob(e,0)) = 0;
      if( getElmNeighbor(e,3)==-1 )
	inner(loc2glob(e,2)) = inner(loc2glob(e,1)) = inner(loc2glob(e,0)) = 0;
    }
    //    if( "ElmT10n3D" == elm_tp ) 
    if( maxnne == 10 ) 
      for( e=0; e<nel; e++ ) {
	if( getElmNeighbor(e,0)==-1 )
	  inner(loc2glob(e,8)) = inner(loc2glob(e,7)) = inner(loc2glob(e,4))=0;
	if( getElmNeighbor(e,1)==-1 )
	  inner(loc2glob(e,9)) = inner(loc2glob(e,8)) = inner(loc2glob(e,5))=0;
	if( getElmNeighbor(e,2)==-1 )
	  inner(loc2glob(e,9)) = inner(loc2glob(e,7)) = inner(loc2glob(e,6))=0;
	if( getElmNeighbor(e,3)==-1 )
	  inner(loc2glob(e,6)) = inner(loc2glob(e,5)) = inner(loc2glob(e,4))=0;
      }
    break;
  default: 
    //    errorFP("WavesGridFEHier:: removeInteriorIndicators","nsd=%d",nsd);
    cout<<"WavesGridFEHier:: removeInteriorIndicators nsd>3 \n";
  }
  /*
  MatSimple(char)& indicators = boundaryData().indicatorAccess();
// check that size of indicators really is (nno,nbind) ?

  for ( i=1; i <= nno; i++ )
    if( inner(i) )
      for ( b = 1; b <= nbind; b++)
	indicators(i,b) = '0';
	*/
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: isElmNeighborConsistent () const
//-----------------------------------------------------------------------------
// Check if neighbor element numbers are consistent, i.e., a neighbor to
// an element must have the element as its neighbor. return true if consistent 
{
// ksa: could extend this routine to check that the nodes situated on a 
// face/edge also exists in the neighbor element. Not yet done.
  if( ! isNeighborsComputed() ){ // is the neighbor info created?
    //    errorFP("WavesGridFEHier:: isElmNeighborConsistent",
    //	      "Neighbor info is not computed");
    cout<<"WavesGridFEHier:: isElmNeighborConsistent Neighbor info is not computed\n";
    return false;
  }
  int e,n,elm;
  int nfaces = nsd+1; // number of sides 2 in 1D, 3 in 2D and 4 in 3D
  for( e=0; e<nel; e++ )
    for( n=0; n<nfaces; n++ ){
      //      if((elm = elm_neigh( e, n ))) // elm is zero if boundary i.e. no neighbor
      elm=getElmNeighbor( e, n );
      if(elm >= 0) // elm is zero if boundary i.e. no neighbor
	if( ! isElmNeighbor( elm, e ) ) {
	  //	  errorFP("WavesGridFEHier:: isElmNeighborConsistent",
	  //		  "\nElement %d and face/edge/node %d, neighbor is %d\n"
	  //		  "parents are %d and %d."
	  //		  ,e,n,elm,getParent(e),getParent(elm));
	  cout<<"WavesGridFEHier:: isElmNeighborConsistent\n";
	  return false; 
	}
    }
  return true;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier::initNeighbor (bool alwaysSortEdgesIn2D)
//-----------------------------------------------------------------------------
// Initiate the elm_neigh structure so that it contains the neighbors elements.
{
  elm_neigh.newsize( nel, nsd+1 ); // dimensioning of neighbor data structure    
  //  elm_neigh.fill(0); // 0 means boundary edge
  elm_neigh=-1;

  switch ( nsd ) {
  case 1: initNeighbor1D(); break;
  case 2: initNeighbor2D( alwaysSortEdgesIn2D ); break;
  case 3: initNeighbor3D(); break;
    //  default: errorFP("WavesGridFEHier::initNeighbor","nsd = %d not implemented", nsd);
  default: cout<<"WavesGridFEHier::initNeighbornsd >3 not implemented\n";
  }
  setNeighborsComputed( true );

  //#ifdef ARRAY_RANGECHECK
  if( ! isElmNeighborConsistent() )
    //    fatalerrorFP("WavesGridFEHier:: initNeighbor", "Not consistent.");
    cout<<"WavesGridFEHier:: initNeighbor  Not consistent.\n";
  if( !checkElementOrientation() )
    //    fatalerrorFP("WavesGridFEHier:: initNeighbor3D", "Wrong orientation.");
    cout<<"WavesGridFEHier:: initNeighbor3D Wrong orientation\n";
  //#endif
}

//-----------------------------------------------------------------------------
void WavesGridFEHier::initNeighbor1D ()
//-----------------------------------------------------------------------------
// Initiate the elm_neigh structure so that it contains the neighbors elements.
{
  int e,k,ne;

// Find the neighbor element which have one node in common
// with a given side of element e. If no element is found this
// is a boundary edge.


// should check a standard case 

// elm_neigh  1 2 3 ...  nel-1 nel
// --------------------------------   
//            0 1 2 ... nel-2 nel-1
//            2 3 4 ...  nel    0 

// Very efficient way to compute element-element info.
// The n2e vector is initially zero. In the loop over
// the elements each edge is accessed at most twice.
// The first time, the number of the neighbor element
// is put in n2e, the second time the elm_neigh is set
// correctly. The edges which only are accessed once are 
// left as zeros in elm_neigh as they should.

  //  VecSimple(int) n2e(nno);
  //  n2e.fill( 0 );
  MV_Vector<int> n2e(nno,-1);

    // element-element neighbor info.:

  //  int secondnode = ("ElmB2n1D" == elm_tp) ? 2:3;
  int secondnode = (maxnne==2) ? 1:2;

  for( e=0; e<nel; e++ ){
    ne = n2e( k = loc2glob( e, 0 ) );
    if( ne>=0 ) {
      putElmNeighbor( e, 1, ne );
      putElmNeighbor( ne, 0, e );
    }
    else
      n2e( k ) = e;
    ne = n2e( k = loc2glob( e, secondnode ) );
    if( ne>=0 ) {
      putElmNeighbor( e, 0, ne );
      putElmNeighbor( ne, 1, e );
    }
    else
      n2e( k ) = e;
  }
}


//-----------------------------------------------------------------------------
void WavesGridFEHier::initNeighbor2D ( bool alwaysSortEdgesIn2D )
//-----------------------------------------------------------------------------
// Initiate the elm_neigh structure so that it contains the neighbors elements.
{
  //  CPUclock cpu;
  //  cpu.initTime();

  int n,e,j,k,n1,n2,ne;

// Reorder the element-node info so that the first edge in the
// element is longest, (the edge between node 1 and 2.
// warning: the ordering may not be correct after that
// the nodes are moved, e.g., using the move or smoothing functions
// Normally if no nodes has been moved in a refined grid it is not
// necessary to sort the order, they are already correct. Thereof
// the parameter alwaysSortEdgesIn2D

  if( isThisCoarsestGrid () || alwaysSortEdgesIn2D ) { 
    //    if ( "ElmT3n2D" == elm_tp ) {
    if ( maxnne==3 ) {
      for ( e=0; e<nel; e++ ) {
	n = getLongestEdgeNo( e );
	if( n == 0 ) // if true then correct order
	  continue;
	if( n == 1 ) 
	  putLoc2Glob123( loc2glob(e,1), loc2glob(e,2), loc2glob(e,0), e);
	else     
	  putLoc2Glob123( loc2glob(e,2), loc2glob(e,0), loc2glob(e,1), e);
      }
    } else { // Quadratic elements
      for ( e=0; e<nel; e++ ){
	n = getLongestEdgeNo( e );
	if( n == 0 ) // if true then correct order
	  continue;
	if( n == 1 ) {
	  putLoc2Glob123( loc2glob(e,1), loc2glob(e,2), loc2glob(e,0), e);
	  putLoc2Glob456( loc2glob(e,4), loc2glob(e,5), loc2glob(e,3), e);
	} else {
	  putLoc2Glob123( loc2glob(e,2), loc2glob(e,0), loc2glob(e,1), e);
	  putLoc2Glob456( loc2glob(e,5), loc2glob(e,3), loc2glob(e,4), e);
	}
      }
    }
  }
//cout<<"\ninit  after edge ordering, cpu-time = " << cpu.getInterval()<<"\n";

  //  if( "ElmT3n2D" == elm_tp ) {
  if( maxnne==3 ) {

// Find the neighbor elements which have two nodes in common with a given 
// side of element e. If no element is found this is a boundary edge.

    // first compute node-to-element information
    // Key observation: An element can be added to a given node at most once
    // hence no sort is needed for the node-to-element (n2e) case.
    //    VecSimple(int) countn2e( nno );
    //    VecSimple(int) n2ejcol(nel*3); // each element is added to exactly 3 nodes
    //    VecSimple(int) n2eirow(nno+1);
    //    countn2e.fill(0);

    MV_Vector<int> countn2e( nno,0 );
    MV_Vector<int> n2ejcol(nel*3); // each element is added to exactly 3 nodes
    MV_Vector<int> n2eirow(nno+1);

    for ( e=0; e<nel; e++ ) { // count the number of elements a node is in
      countn2e( loc2glob(e,0) )++;
      countn2e( loc2glob(e,1) )++;
      countn2e( loc2glob(e,2) )++;
    }
// set up the n2eirow which tells us where next node begins
    /*    int prev = n2eirow(1) = 1;
    for (n=1; n<=nno; n++ )  prev = n2eirow( n+1 ) = countn2e( n ) + prev;
    countn2e.fill(0);
    for ( e=1; e<=nel; e++ ){
      ne = loc2glob(e,1);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
      ne = loc2glob(e,2);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
      ne = loc2glob(e,3);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
    }
    */
    int prev = n2eirow(0) = 0;
    for (n=0; n<nno; n++ )  prev = n2eirow( n+1 ) = countn2e( n ) + prev;
    countn2e=0;
    for ( e=0; e<nel; e++ ){
      ne = loc2glob(e,0);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
      ne = loc2glob(e,1);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
      ne = loc2glob(e,2);
      n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
    }
//cout << "\ninit after computing n2e, cpu-time = " << cpu.getInterval()<<"\n";

    // compute element-element neighbor info.:

// it is faster to go backwards due to the fact that the elements in n2e
// are ordered starting with lowest element number and increases.
// Therefore the neighbor element will be found earlier.
    /*
    for( e=nel; e>=1; e-- ) 
      for( j=1; j<=3; j++ ){
	if( getElmNeighbor(e,j)==-1 ) { // if not yet set
	  n1 = loc2glob( e, j ); // first node
	  n2 = loc2glob( e, j==3 ? 1 : j+1 );
	  int kstop = n2eirow(n1+1);
	  for( k=n2eirow(n1); k<kstop; k++ ) { // search for neighbor element
	    ne = n2ejcol(k);
	    if ( ne != e ) { // find second common node
	      if ( loc2glob (ne, 1) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 1,  e );
		k = kstop; // to break the k-loop
	      }
	      else if ( loc2glob (ne, 2) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 2,  e );
		k = kstop; // to break the k-loop
	      }
	      else if ( loc2glob (ne, 3) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 3,  e );
		k = kstop; // to break the k-loop
	      }
	    }
	  }
	}
      }
    */
    for( e=nel-1; e>=0; e-- ) 
      for( j=0; j<3; j++ ){
	if( getElmNeighbor(e,j)==-1 ) { // if not yet set
	  n1 = loc2glob( e, j ); // first node
	  n2 = loc2glob( e, j==2 ? 0 : j+1 );
	  int kstop = n2eirow(n1+1);
	  for( k=n2eirow(n1); k<kstop; k++ ) { // search for neighbor element
	    ne = n2ejcol(k);
	    if ( ne != e ) { // find second common node
	      if ( loc2glob (ne, 0) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 0,  e );
		k = kstop; // to break the k-loop
	      }
	      else if ( loc2glob (ne, 1) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 1,  e );
		k = kstop; // to break the k-loop
	      }
	      else if ( loc2glob (ne, 2) == n2 ) {
		putElmNeighbor(  e, j, ne );
		putElmNeighbor( ne, 2,  e );
		k = kstop; // to break the k-loop
	      }
	    }
	  }
	}
      }

//  cout << "\ninit finished e2e, cpu-time = " << cpu.getInterval()<<"\n";
  } else {  // if( "ElmT6n2D" == elm_tp )

// Very efficient way to compute element-element info for quadratic
// elements, based on that the nodes on mid of edges is only shared
// by two elements.
// The n2e vector is initially zero. In the loop over the elements each 
// edge is accessed at most twice. The first time, the number of the 
// neighbor element is put in n2e, the second time the elm_neigh is set
// correctly. The edges accessed only once are left as zeros in elm_neigh 
// as they should, (they are boundary edges).

    MV_Vector<int> n2e(nno,-1);
    for( e=0; e<nel; e++ )
      for( j=3; j<6; j++ ){
	ne = n2e( k = loc2glob( e, j ));
	if( ne>=0 ) {
	  putElmNeighbor( e, j-3, ne );
	  if( loc2glob( ne, 3 ) == k )       putElmNeighbor( ne, 0, e );
	  else if( loc2glob( ne, 4 ) == k )  putElmNeighbor( ne, 1, e );
	  else if( loc2glob( ne, 5 ) == k )  putElmNeighbor( ne, 2, e );
	}
	else
	  n2e( k ) = e;
      }
    /*
    VecSimple(int) n2e(nno);
    n2e.fill( 0 );
    for( e=1; e<=nel; e++ )
      for( j=4; j<=6; j++ )
	if( (ne = n2e( k = loc2glob( e, j ) )) ) {
	  putElmNeighbor( e, j-3, ne );
	  if( loc2glob( ne, 4 ) == k )       putElmNeighbor( ne, 1, e );
	  else if( loc2glob( ne, 5 ) == k )  putElmNeighbor( ne, 2, e );
	  else if( loc2glob( ne, 6 ) == k )  putElmNeighbor( ne, 3, e );
	}
	else
	  n2e( k ) = e;
*/
  }
}


//-----------------------------------------------------------------------------
int WavesGridFEHier:: firstOf256 ( MV_Vector<int>& order )
//-----------------------------------------------------------------------------
{ 
// order contains order of longest edges where edge 1 is longest
// find which of edges 2, 5 and 6 are longest
  int ed = order(1);
  if( ed==1 || ed==4 || ed==5 ) return ed;
  ed = order(2);
  if( ed==1 || ed==4 || ed==5 ) return ed;
  ed = order(3);
  if( ed==1 || ed==4 || ed==5 ) return ed;
  //  errorFP("WavesGridFEHier:: firstOf256","should have returned");
  //  order.print(s_o);
  cout<<"WavesGridFEHier:: firstOf256 should have returned\n";
  for(int i=0;i<6;i++)
    cout<<" "<<order(i);
  cout<<"\n"<<flush;
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: firstOf346( MV_Vector<int>& order )
//-----------------------------------------------------------------------------
{
// order contains order of longest edges where edge 1 is longest
// find which of edges 3, 4 and 6 are longest
  int ed = order(1);
  if( ed==2 || ed==3 || ed==5 ) return ed;
  ed = order(2);
  if( ed==2 || ed==3 || ed==5 ) return ed;
  ed = order(3);
  if( ed==2 || ed==3 || ed==5 ) return ed;
  //  errorFP("WavesGridFEHier:: firstOf346","should have returned");
  //  order.print(s_o);
  cout<<"WavesGridFEHier:: firstOf346 should have returned\n";
  for(int i=0;i<6;i++)
    cout<<" "<<order(i);
  cout<<"\n"<<flush;
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: firstOf256( MV_Vector<int>& order,
			     const int ed2, const int ed5, const int ed6)
//-----------------------------------------------------------------------------
{
  int ed = order(1);
  if( ed==ed2 ) return 1;
  if( ed==ed5 ) return 4;
  if( ed==ed6 ) return 5;
  ed = order(2);
  if( ed==ed2 ) return 1;
  if( ed==ed5 ) return 4;
  if( ed==ed6 ) return 5;
  ed = order(3);
  if( ed==ed2 ) return 1;
  if( ed==ed5 ) return 4;
  if( ed==ed6 ) return 5;

  //  errorFP("WavesGridFEHier:: firstOf256","should have returned");
  //  order.print(s_o);
  cout<<"WavesGridFEHier:: firstOf256 should have returned\n";
  for(int i=0;i<6;i++)
    cout<<" "<<order(i);
  cout<<"\n"<<flush;
  return 0;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: firstOf346( MV_Vector<int>& order,
			     const int ed3, const int ed4, const int ed6)
//-----------------------------------------------------------------------------
{
  int ed = order(1);
  if( ed==ed3 ) return 2;
  if( ed==ed4 ) return 3;
  if( ed==ed6 ) return 5;
  ed = order(2);
  if( ed==ed3 ) return 2;
  if( ed==ed4 ) return 3;
  if( ed==ed6 ) return 5;
  ed = order(3);
  if( ed==ed3 ) return 2;
  if( ed==ed4 ) return 3;
  if( ed==ed6 ) return 5;
  //  errorFP("WavesGridFEHier:: firstOf346","should have returned");
  //  order.print(s_o);
  cout<<"WavesGridFEHier:: firstOf346 should have returned\n";
  for(int i=0;i<6;i++)
    cout<<" "<<order(i);
  cout<<"\n"<<flush;
  return 0;
}

//------------------------------------------------------------------------
void WavesGridFEHier:: initNeighbor3D ()
//------------------------------------------------------------------------
// Initiate the elm_neigh structure so that it contains the numbers of
// the neighbors elements on the other side of the element faces.
{
// Reorder the element-node info so that the first edge is longest,
// warning: this could be changed by the move function

  //  CPUclock cpu;
  //  cpu.initTime();

  int n,e,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,ne;

  //  bool T10 = ( "ElmT10n3D" == elm_tp ) ? true : false;
  bool T10 = ( maxnne==10 ) ? true : false;

  edge_face2.newsize(nel);
  edge_face3.newsize(nel);
  lengthsscratch.newsize(6); // global arrays used in getEdgeLengthOrder
  orderscratch.newsize(6);

  if(maxnne<=10)
  for ( e=0; e<nel; e++ ) {
    getEdgeLengthOrder( orderscratch, e );
    n = orderscratch(0);
    if( n != 0) { // if n is 1 then correct order else reorder
      n1 = loc2glob( e, 0 );  n2 = loc2glob( e, 1 );
      n3 = loc2glob( e, 2 );  n4 = loc2glob( e, 3 );
      if ( T10 ) {
	n5 = loc2glob( e, 4 );	n6 = loc2glob( e, 5 );
	n7 = loc2glob( e, 6 );	n8 = loc2glob( e, 7 );
	n9 = loc2glob( e, 8 );	n10= loc2glob( e, 9 );
      }
    }
    switch ( n ) { // two alternatives, which to choose?
    case 0:
      edge_face2(e) = firstOf256(orderscratch);
      edge_face3(e) = firstOf346(orderscratch);
      break;
    case 1:
      putLoc2Glob1234 ( n2, n3, n1, n4, e );   // now prev edge 2 is first
      if(T10)
	putLoc2Glob5678910 (n6, n7, n5, n9, n10, n8, e );
      edge_face2(e) = firstOf256(orderscratch,2,5,3);
      edge_face3(e) = firstOf346(orderscratch,0,4,3);
      break;
    case 2: 
      putLoc2Glob1234 ( n3, n1, n2, n4, e );	
      if(T10)
	putLoc2Glob5678910 ( n7, n5, n6, n10, n8, n9, e );
      edge_face2(e) = firstOf256(orderscratch,0,3,4);
      edge_face3(e) = firstOf346(orderscratch,1,5,4);
      break;
    case 3: 
      putLoc2Glob1234 ( n1, n4, n2, n3, e );	
      if(T10)
	putLoc2Glob5678910 ( n8, n9, n5, n7, n10, n6, e );
      edge_face2(e) = firstOf256(orderscratch,4,5,1);
      edge_face3(e) = firstOf346(orderscratch,0,2,1);
      break; 
    case 4: 
      putLoc2Glob1234 ( n2, n4, n3, n1, e );	
      if(T10)
	putLoc2Glob5678910 ( n9, n10, n6, n5, n8, n7, e );
      edge_face2(e) = firstOf256(orderscratch,5,3,2);
      edge_face3(e) = firstOf346(orderscratch,1,0,2);
      break; 
    case 5: 
      putLoc2Glob1234 ( n3, n4, n1, n2, e );	
      if(T10)
	putLoc2Glob5678910 ( n10, n8, n7, n6, n9, n5, e );
      edge_face2(e) = firstOf256(orderscratch,3,4,0);
      edge_face3(e) = firstOf346(orderscratch,2,1,0);
      break;
      //    default: errorFP("WavesGridFEHier::initNeighbor3D","illegal value %d",n);
    default: cout<<"WavesGridFEHier::initNeighbor3D illegal value\n";
    }
  }
//  cout << "\ninit0   cpu-time = " << cpu.getInterval() <<"\n";

// Given face of element e, find the neighbor elements which have 
// three nodes in common with the face. If no such element is found this
// is a boundary face.

  //  additional_info.initNeighbor(*this,true,false,false);
  //  WavesNeighborFE& neighbor = additional_info.neighbor;
  if( neighbor.nodeSize() == 0 )
    neighbor.init(*this,true,false,false);

//  cout << "\ninit1   cpu-time = " << cpu.getInterval() <<"\n";

  //  if( "ElmT4n3D" == elm_tp ){
  if( true || maxnne==4 || maxnne==10){

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3
    bool safe = true;
    int ne1,ne2;
    //    for( e=nel; e>=1; e-- ){ // more efficient with decreasing elements
    for( e=nel-1; e>=0; e-- ){ // more efficient with decreasing elements
      n1 = loc2glob( e, 0 );
      n2 = loc2glob( e, 1 );
      n3 = loc2glob( e, 2 );
      n4 = loc2glob( e, 3 );
      if( getElmNeighbor(e,0)==-1 )
	if( getElmNeighbor(e,3)==-1 ){
	  getElmsAtFaces(n1,n2,n4,n3,e,ne1,ne2);// case face 1 and 4 not done
	  if( ne1>=0 ){
	    n = getFaceNo(n1,n2,n4,ne1,safe);
	    putElmNeighbor(  e, 0, ne1 );
	    putElmNeighbor( ne1, n,  e );
	  }
	  if( ne2>=0 ){
	    n = getFaceNo(n1,n2,n3,ne2,safe);
	    putElmNeighbor(  e, 3, ne2 );
	    putElmNeighbor( ne2, n,  e );
	  }
	}
	else {
	  ne = getElmAtFace(n1,n2,n4,e); // case face 1
	  if( ne>=0 ){
	    n = getFaceNo(n1,n2,n4,ne,safe);
	    putElmNeighbor(  e, 0, ne );
	    putElmNeighbor( ne, n,  e );
	  }
	}
      else if( getElmNeighbor(e,3)==-1 ) {
	ne = getElmAtFace(n1,n2,n3,e);    // case face 4
	if( ne>=0 ){
	  n = getFaceNo(n1,n2,n3,ne,safe);
	  putElmNeighbor(  e, 3, ne );
	  putElmNeighbor( ne, n,  e );
	}
      }

      if( getElmNeighbor(e,1)==-1 )
	if( getElmNeighbor(e,2)==-1 ){
	  getElmsAtFaces(n3,n4,n2,n1,e,ne1,ne2);
	  if( ne1>=0 ){
	    n = getFaceNo(n3,n4,n2,ne1,safe);
	    putElmNeighbor(  e, 1, ne1 );
	    putElmNeighbor( ne1, n,  e );
	  }
	  if( ne2>=0 ){
	    n = getFaceNo(n3,n4,n1,ne2,safe);
	    putElmNeighbor(  e, 2, ne2 );
	    putElmNeighbor( ne2, n,  e );
	  }
	}
	else {
	  ne = getElmAtFace(n2,n3,n4,e);
	  if( ne>=0 ){
	    n = getFaceNo(n2,n3,n4,ne,safe);
	    putElmNeighbor(  e, 1, ne );
	    putElmNeighbor( ne, n,  e );
	  }
	}
      else if( getElmNeighbor(e,2)==-1 ){
	ne = getElmAtFace(n1,n3,n4,e);
	if( ne>=0 ){
	  n = getFaceNo(n1,n3,n4,ne,safe);
	  putElmNeighbor(  e, 2, ne );
	  putElmNeighbor( ne, n,  e );
	}
      }
    }
  }
  /*
//  cout << "\ninit2   cpu-time = " << cpu.getInterval() <<"\n";

  //  if( "ElmT10n3D" == elm_tp ){
  if( maxnne==10 ){

// face1 has nodes 1 2 4 and midnodes 5 8 9
// face2 has nodes 2 3 4 and midnodes 6 9 10
// face3 has nodes 1 3 4 and midnodes 7 8 10
// face4 has nodes 1 2 3 and midnodes 5 6 7

    for( e=nel-1; e>=0; e-- ){ 
      n1 = loc2glob( e, 0 );
      n2 = loc2glob( e, 1 );
      n3 = loc2glob( e, 2 );
      n4 = loc2glob( e, 3 );
      bool safe=false;
      if( getElmNeighbor(e,0)==-1 ) {     // case face 1 not yet done
	ne = getElmAtFace(loc2glob(e,4),loc2glob(e,7),e); // find  element
	if( ne>=0 ){                      // ne is 0 if it is a boundary face
	  n = getFaceNo(n1,n2,n4,ne,safe);  // find the face number
	  putElmNeighbor(  e, 0, ne ); // put found elm in neighbor array
	  putElmNeighbor( ne, n,  e );
	}
      }
      if( getElmNeighbor(e,1)==-1 ) {     // case face 2 not yet done
	ne = getElmAtFace(loc2glob(e,5),loc2glob(e,8),e); 
	if( ne>=0 ){
	  n = getFaceNo(n2,n3,n4,ne,safe); 
	  putElmNeighbor(  e, 1, ne );
	  putElmNeighbor( ne, n,  e );
	}
      }
      if( getElmNeighbor(e,2)==-1 ) {     // case face 3 not yet done
	ne = getElmAtFace(loc2glob(e,6),loc2glob(e,7),e);
	if( ne>=0 ){
	  n = getFaceNo(n1,n3,n4,ne,safe); 
	  putElmNeighbor(  e, 2, ne );
	  putElmNeighbor( ne, n,  e );
	}
      }
      if( getElmNeighbor(e,3)==-1 ) {     // case face 4 not yet done
	ne = getElmAtFace(loc2glob(e,4),loc2glob(e,5),e);   
	if( ne>=0 ){
	  n = getFaceNo(n1,n2,n3,ne,safe); 
	  putElmNeighbor(  e, 3, ne );
	  putElmNeighbor( ne, n,  e );
	}
      }
    }
  }
  */
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: isElmNeighbor( const int e, const int elnei ) const
//-----------------------------------------------------------------------------
// Does element e have element elnei as a neighbor?
{
  if( !isNeighborsComputed() ) { // Is the element-element info computed?
    //    warningFP("WavesGridFEHier::isElmNeighbor","Element neighbors not computed");
    cout<<"WavesGridFEHier::isElmNeighbor  Element neighbors not computed\n";
    return false; 
  }
  if( e==-1 ) {
    //    warningFP("WavesGridFEHier::isElmNeighbor","Element number is 0");
    cout<<"WavesGridFEHier::isElmNeighbor Element number is -1\n";
    return false; 
  }
  int nof = nsd+1; 
  for( int i=0; i<nof; i++ )
    if( getElmNeighbor(e,i) == elnei ) return true;
  return false; // no hit was found
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initNumbering() 
//-----------------------------------------------------------------------------
{
  if( !isNeighborsComputed() ) // Is the element-element info computed?
    initNeighbor();  // if not compute it
  //  if( "ElmT4n3D" == elm_tp )
  //    new_nodes.newsize( nel, 6 );
  //  else if( "ElmT3n2D" == elm_tp)
  //    new_nodes.newsize( nel, 3 );
  //  else if( "ElmT6n2D" == elm_tp )
  //    new_nodes.newsize( nel, 6 );
  //  else if( "ElmT10n3D" == elm_tp )
  //    new_nodes.newsize( nel, 31 ); // 6 old 12 on edges 12 on faces 1 internal
  //  else if( "ElmB2n1D" == elm_tp || "ElmB3n1D" == elm_tp)
  //    new_nodes.newsize( nel, 1 ); 
  //  else 
  //    fatalerrorFP("WavesGridFEHier:: initNumbering"," %s not implemented", 
  //		 elm_tp.chars());
  if( nsd==3&&maxnne==4 )
    new_nodes.newsize( nel, 6 );
  else if( nsd==2&&maxnne==3 )
    new_nodes.newsize( nel, 3 );
  else if( nsd==2&&maxnne==6 )
    new_nodes.newsize( nel, 6 );
  else if( nsd==3&&maxnne==10 )
    new_nodes.newsize( nel, 31 ); // 6 old 12 on edges 12 on faces 1 internal
  else if( nsd==1&&maxnne==2 || nsd==1&&maxnne==3)
    new_nodes.newsize( nel, 1 ); 
  else 
    //    fatalerrorFP("WavesGridFEHier:: initNumbering"," %s not implemented", 
    //		 elm_tp.chars());
    cout<<"WavesGridFEHier:: initNumbering element not implemented\n";
  new_nodes=0;
  setNewNodesSet( true );
  //  setCurrNodeNo(nno); 
  setCurrNodeNo(nno-1); 
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: markElementsBool( const MV_Vector<bool>& ele )
//-----------------------------------------------------------------------------
// ele is true for the elements to be refined
{
  initNumbering();
  int e;
  if ( nsd==1 ) { //in 1D node numbering is simple, all new nodes inside elm
    int refm = getRefinementMethod();
    for( e=0; e<nel; e++ ) if( ele(e) ) new_nodes(e,0) = refm;
    return true;
  }

  stackcounter=0;
  if( ele.size() != nel )
    //    errorFP("WavesGridFEHier:: markElementsBool",
    //	    "wrong size of vector (%d)",ele.size());
    cout<<"WavesGridFEHier:: markElementsBool wrong size of vector\n";
  else // in 2D and 3D use recursive method
    for( e=0; e<nel; e++ ) if( ele( e ) ) markAllEdgesRec( e );

  elmstack.newsize(0);
  ednostack.newsize(0);

  extraNumberingForT10 (); // add nodes on faces

  return true;
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: markAllElements()
//-----------------------------------------------------------------------------
// ele is true for the elements to be refined
{
  cout<<"markAllElements begin\n"<<flush;
  initNumbering();

  if ( nsd==1 ) {
    //    new_nodes.fill( getRefinementMethod() );
    new_nodes=getRefinementMethod();
    return true;
  }

  cout<<"markAllElements after initnumbering\n"<<flush;
  int e;
  stackcounter=0;
  if( nsd==2 && getRefinementMethod() == 3 )// all edges refined
    for( e=0; e<nel; e++ ){
      if( ! new_nodes(e,0) ) numberingNewNodes( e, 0 );
      if( ! new_nodes(e,1) ) numberingNewNodes( e, 1 );
      if( ! new_nodes(e,2) ) numberingNewNodes( e, 2 );
    }
  else if( nsd==3 && getRefinementMethod() == 6 ) { // all edges refined
    //  cout<<"markAllElements 2\n"<<flush;
    for( e=0; e<nel; e++ ){
      //      cout<<e<<" element markAllElements\n"<<flush;
      if( ! new_nodes(e,0) ) { numberingNewNodes( e, 0 ); stackcounter=0; }
      if( ! new_nodes(e,1) ) { numberingNewNodes( e, 1 ); stackcounter=0; }
      if( ! new_nodes(e,2) ) { numberingNewNodes( e, 2 ); stackcounter=0; }
      if( ! new_nodes(e,3) ) { numberingNewNodes( e, 3 ); stackcounter=0; }
      if( ! new_nodes(e,4) ) { numberingNewNodes( e, 4 ); stackcounter=0; }
      if( ! new_nodes(e,5) ) { numberingNewNodes( e, 5 ); stackcounter=0; }
    }
  }
  else
    for( e=0; e<nel; e++ )
      markAllEdgesRec( e );

  elmstack.newsize(0);
  ednostack.newsize(0);

  extraNumberingForT10();

  return true;
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: markElements( const MV_Vector<int>& ele )
//-----------------------------------------------------------------------------
// If size of ele is nel, ele is nonzero for the elements to be refined
// otherwise ele contains numbers of the elements to be refined.
{
  initNumbering();
  int e;
  if ( nsd==1 ) { 
// in 1D new_nodes contain the number of new elements caused by
// refinement of the element
    for( e=0; e<nel; e++ )
      new_nodes(e,0) = ele(e);
    return true;
  }

  stackcounter=0;
  if( ele.size() != nel )   
    //    errorFP("WavesGridFEHier:: markElements","wrong size of vector");
    cout<<"WavesGridFEHier:: markElements  wrong size of vector\n";
  else              // ele is nonzero for the elements to be marked 
    for( e=0; e<nel; e++ )
      markAllEdgesRec( e, ele( e ) );

  elmstack.newsize(0);
  ednostack.newsize(0);

  extraNumberingForT10();

  return true;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: extraNumberingForT10()
//-----------------------------------------------------------------------------
{
  if(maxnne != 10) 
    return; // if not ElmT10n3D

  //  errorFP("WavesGridFEHier:: extraNumberingForT10","not implemented yet");
  cout<<"WavesGridFEHier:: extraNumberingForT10  not implemented yet\n";
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: markEdgeRec( int e, int edge_no )
//-----------------------------------------------------------------------------
// numbering of node at edge with edge_no in element e, note the recursion
{ 
  if( ! new_nodes( e, edge_no ) ) { // if not marked before
    int enei;
    if ( nsd==2 ) {
      numberingNewNodes( e, edge_no );
      enei = getElmNeighbor( e, edge_no );
      if( enei>=0 ) // check if neighbor
	markEdgeRec( enei, 0 ); // is longest edge marked recursive call  
    }
    else { // nsd == 3 
//  find the tetrahedra which share the same edge
      int start = stackcounter;
      numberingNewNodes( e, edge_no ); // stackcounter is increased
      for( int i=start; i<stackcounter; i++ ) {
	enei = elmstack(i); // elm number of element neighbor to edge
	switch ( ednostack(i) ) {
	case 0: // already numbered
	  break;
	case 1: 
	case 4: 
	  //	  printf(" enei=%d, edge_face2=%d\n",enei,edge_face2(enei));
	  markEdgeRec( enei, 0 );                // number longest in element,
	  markEdgeRec( enei, edge_face2(enei) ); // and longest of 2, 5 and 6
	  break;
	case 2:
	case 3:
	  markEdgeRec( enei, 0 );                // number longest in element,
	  markEdgeRec( enei, edge_face3(enei) ); // and longest of 3, 4 and 6
	  break;
	case 5:
	  markEdgeRec( enei, 0 );                // number longest in element,
	  markEdgeRec( enei, edge_face2(enei) ); // and longest of 2, 5 and 6,
	  markEdgeRec( enei, edge_face3(enei) ); // and longest of 3, 4 and 6
	  break;
	default:
	  //	  errorFP("WavesGridFEHier:: markEdgeRec",
	  //		  "illegal edge number %d in element %d",
	  //		  ednostack(i),elmstack(i));
	  cout<<"WavesGridFEHier:: markEdgeRec illegal edge number\n";
	} 
      }
      stackcounter = start; // set stackcounter to the previous value
    }
  }
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: markAllEdgesRec( int e )
//-----------------------------------------------------------------------------
// recursive numbering of edges of element e
{
  if ( nsd==2 ) {
    switch ( getRefinementMethod() ){
    case 1 : // mark only the longest edge of the element
      markEdgeRec( e, 0 );
      break;
    case 2 : // mark the two longest edges of the element
      markEdgeRec( e, 0 );
      markEdgeRec( e, getSecondLongestEdgeNo1( e ) );
      break;
    case 3 : // mark recursively all three edges of the element
      markEdgeRec( e, 0 );
      markEdgeRec( e, 1 );
      markEdgeRec( e, 2 );
      break;
    default:
      //      errorFP("WavesGridFEHier:: markAllEdgesRec","Refinement method %d"
      //	      "not implemented, must be 1-3.", getRefinementMethod() );
      cout<<"WavesGridFEHier:: markAllEdgesRec Refinement method %d"<<
	      "not implemented, must be 1-3.\n";
    }
  } else { // nsd==3
    int i;
    int no_mark_edges = getRefinementMethod();
    switch ( no_mark_edges ) {
    case 1 : // mark recursively the longest edge?
      markEdgeRec( e, 0 );  // is longest edge marked?
      break;
    case 2 : 
    case 3 : 
    case 4 : 
    case 5 :
      getEdgeLengthOrder(orderscratch,e);
      for(i=0;i<no_mark_edges;i++) markEdgeRec( e, orderscratch(i) );
      break;
    case 6 : // mark recursively all six edges of the element
      for ( i=0; i<6; i++ ) markEdgeRec( e, i );
      break;
    default:
      //      errorFP("WavesGridFEHier:: markAllEdgesRec","Refinement method %d"
      //	      "not implemented, must be 1-6.", getRefinementMethod() );
      cout<<"WavesGridFEHier:: markAllEdgesRec Refinement method %d"<<
	      "not implemented, must be 1-6.\n";
    }
  }
}


//-----------------------------------------------------------------------------
void WavesGridFEHier:: markAllEdgesRec( int e, int no_mark_edges )
//-----------------------------------------------------------------------------
// recursive marking of edges of element e
// (not used in 1D)
{
  if ( nsd==2 ) {
    switch ( no_mark_edges ){
    case 0 : break; // mark no edge
    case 1 : // mark only the longest edge of the element
      markEdgeRec( e, 0 );
      break;
    case 2 : // mark the two longest edges of the element
      markEdgeRec( e, 0 );
      markEdgeRec( e, getSecondLongestEdgeNo1( e ) );// assumes 1 is longest
      break;
    case 3 : // mark recursively all three edges of the element
      markEdgeRec( e, 0 );
      markEdgeRec( e, 1 );
      markEdgeRec( e, 2 );
      break;
    default:
      //      errorFP("WavesGridFEHier:: markAllEdgesRec","Mixed refinement method %d"
      //	      "not implemented, must be 0-3.", no_mark_edges);
      cout<<"WavesGridFEHier:: markAllEdgesRec Mixed refinement method %d"<<
	      "not implemented, must be 0-3.\n";
    }
  } else { // nsd==3
    int i;
    switch ( no_mark_edges ) {
    case 0 : break;
    case 1 : // mark recursively the longest edge?
      markEdgeRec( e, 0 );
      break;
    case 2 : 
    case 3 : 
    case 4 : 
    case 5 : 
      getEdgeLengthOrder(orderscratch,e);
      for(i=0;i<no_mark_edges;i++)  markEdgeRec( e, orderscratch(i) );
      break;
    case 6 : // mark recursively all six edges of the element
      for(i=0;i<6;i++)  markEdgeRec( e, i );
      break;
    default:
      //      errorFP("WavesGridFEHier:: markAllEdgesRec","Mixed refinement method %d"
      //	      "not implemented, must be 0-6.", no_mark_edges );
      cout<<"WavesGridFEHier:: markAllEdgesRec Mixed refinement method %d"<<
	      "not implemented, must be 0-6.\n";
    }
  }
}

bool WavesGridFEHier:: lessThan( const int ed1, const int ed2, const int element)
{ // unique ordering of edges, which of edge ed1 and ed2 is longest?
// Routine to order edges that are of equal length
// return true if ed1 is longer than ed2
  int n1,n2,n3,n4;
  getNodesOnEdge ( element, ed1, n1, n2 );
  getNodesOnEdge ( element, ed2, n3, n4 );

  if ( test_equal_with_geometry ) { // test by geometry if equal length
    real x1,x2;
    int d;
    for(d=0;d<nsd;d++) {
      x1 = coord(n1,d) + coord(n2,d);
      x2 = coord(n3,d) + coord(n4,d);
      if( x1 < x2 ) return true;
      if( x1 > x2 ) return false;
    }
    return true; // it seems like the points are identical, should not happen
  }
  else { // test by global node numbers
// first check which edge contains the smallest node value
    //  int min12 = min(n1,n2), min34 = min(n3,n4);
  int min12 = (n1<n2)?n1:n2, min34 = (n3<n4)?n3:n4;
  if ( min12 < min34 ) return true;
  if ( min12 > min34 ) return false;
// next if min12 and min34 equal check the sum
  if ( n1 + n2 < n3 + n4 ) return true;
  }
  return false;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: getEdgeLengthOrder( MV_Vector<int>& order, int e )
//-----------------------------------------------------------------------------
// Get the edge number of the longest edge of element e in position 1
// the second longest in position 2 and so on
{
  int n = loc2glob( e, 0 );
  real x1 = coord(n,0), y1 = coord(n,1), z1 = coord(n,2);
  n = loc2glob( e, 1 );
  real x2 = coord(n,0), y2 = coord(n,1), z2 = coord(n,2);
  n = loc2glob( e, 2 );
  real x3 = coord(n,0), y3 = coord(n,1), z3 = coord(n,2);
  n = loc2glob( e, 3 );
  real x4 = coord(n,0), y4 = coord(n,1), z4 = coord(n,2);

  lengthsscratch(0) = sqr( x1 - x2 ) + sqr( y1 - y2 ) + sqr( z1 - z2 );
  lengthsscratch(1) = sqr( x2 - x3 ) + sqr( y2 - y3 ) + sqr( z2 - z3 );
  lengthsscratch(2) = sqr( x3 - x1 ) + sqr( y3 - y1 ) + sqr( z3 - z1 );
  lengthsscratch(3) = sqr( x1 - x4 ) + sqr( y1 - y4 ) + sqr( z1 - z4 );
  lengthsscratch(4) = sqr( x2 - x4 ) + sqr( y2 - y4 ) + sqr( z2 - z4 );
  lengthsscratch(5) = sqr( x3 - x4 ) + sqr( y3 - y4 ) + sqr( z3 - z4 ); 
//  order.newsize(6); // assume that order and lengthsscratch has been newsizeed
  order(0)=0; order(1)=1; order(2)=2; order(3)=3; order(4)=4; order(5)=5;

// sort the lengths if equal length then "lessThan" is used
  int i,j,idx,idx2;
  real v;
  for (j=1;j<6;j++) {
    v=lengthsscratch(idx=order(i=j));
    while( lengthsscratch(idx2=order(i-1)) < v ||
	   ( lengthsscratch(idx2) == v && lessThan(idx2,idx, e ) ) ) 
      { order(i)=idx2; i--; if(i==0) break; }
    order(i)=idx;
  }
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getLongestEdgeNo( int e )
//-----------------------------------------------------------------------------
// Get the local edge number of the longest edge in element e
{
  if ( nsd == 2 ) {
    int n = loc2glob( e, 0 );
    real x1 = coord(n,0), y1 = coord(n,1);
    n = loc2glob( e, 1 );
    real x2 = coord(n,0), y2 = coord(n,1);
    n = loc2glob( e, 2 );
    real x3 = coord(n,0), y3 = coord(n,1);
    real l1 = sqr( x1 - x2 ) + sqr( y1 - y2 );
    real l2 = sqr( x2 - x3 ) + sqr( y2 - y3 );
    real l3 = sqr( x3 - x1 ) + sqr( y3 - y1 );
    if( l1 >= l2 && l1 >= l3 ) return 0;
    return l2 >= l3 ? 1 : 2;
  }
  else if (nsd == 3 ) {
    getEdgeLengthOrder( orderscratch, e );
    return orderscratch(0);
  }
  //  errorFP("WavesGridFEHier:: getLongestEdgeNo",
  //	  "only implemented in 2D and 3D nsd=%d",nsd);
  cout<<"WavesGridFEHier:: getLongestEdgeNo"<<
	  "only implemented in 2D and 3D nsd=%d\n";
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getSecondLongestEdgeNo1( int e )
//-----------------------------------------------------------------------------
// Get the local edge number of the second longest
// edge in element e, in 2D it is assumed that edge number 1 is longest
{
  if ( nsd==2 ) {
    int n = loc2glob( e, 1 );
    real x2 = coord(n,0), y2 = coord(n,1);
    n = loc2glob( e, 2 );
    real x3 = coord(n,0), y3 = coord(n,1);
    real l2 = sqr( x2 - x3 ) + sqr( y2 - y3 );
    real l3 = sqr( x3 - coord(n,0) ) + sqr( y3 - coord(n,1) );
    return l2 < l3 ? 2 : 1;
  }
  else if ( nsd==3 ) {
    getEdgeLengthOrder( orderscratch, e );
    return orderscratch(1);
  }
  //  errorFP("WavesGridFEHier:: getSecondLongestEdgeNo1",
  //	  "only implemented in 2D and 3D nsd=%d",nsd);
  cout<<"WavesGridFEHier:: getSecondLongestEdgeNo1"<<
	  "only implemented in 2D and 3D nsd=%d\n";
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getShortestEdgeNo1( int e )
//-----------------------------------------------------------------------------
// Get the local edge number of the shortest edge in element e, in 2D it is 
// assumed edge 1 is longest
{
  if ( nsd==2 ) {
    int n = loc2glob( e, 1 );
    real x2 = coord(n,0), y2 = coord(n,1);
    n = loc2glob( e, 2 );
    real x3 = coord(n,0), y3 = coord(n,1);
    real l2 = sqr( x2 - x3 ) + sqr( y2 - y3 );
    real l3 = sqr( x3 - coord(n,0) ) + sqr( y3 - coord(n,1) );
    return l2 < l3 ? 1 : 2;
  }
  else if ( nsd==3 ) {
    getEdgeLengthOrder( orderscratch, e );
    return orderscratch(5);
  }
  //  errorFP("WavesGridFEHier:: getShortestEdgeNo1",
  //	  "only implemented in 2D and 3D nsd=%d",nsd);
  cout<<"WavesGridFEHier:: getShortestEdgeNo1"<<
	  "only implemented in 2D and 3D nsd=%d\n";
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getShortestEdgeNo( int e )
//-----------------------------------------------------------------------------
// Get the local edge number of the shortest edge in element e
{
  if ( nsd==2 ) {
    int n = loc2glob( e, 0 );
    real x1 = coord(n,0), y1 = coord(n,1);
    n = loc2glob( e, 1 );
    real x2 = coord(n,0), y2 = coord(n,1);
    n = loc2glob( e, 2 );
    real x3 = coord(n,0), y3 = coord(n,1);
    real l1 = sqr( x1 - x2 ) + sqr( y1 - y2 );
    real l2 = sqr( x2 - x3 ) + sqr( y2 - y3 );
    real l3 = sqr( x3 - x1 ) + sqr( y3 - y1 );
    if( l1 < l2 && l1 < l3 ) return 0;
    return l2 < l3 ? 1 : 2;
  }
  else if ( nsd==3 ) { 
    getEdgeLengthOrder( orderscratch, e );
    return orderscratch(5);
  }
  //  errorFP("WavesGridFEHier:: getShortestEdgeNo",
  //	  "only implemented in 2D and 3D nsd=%d",nsd);
  cout<<"WavesGridFEHier:: getShortestEdgeNo"<<
	  "only implemented in 2D and 3D nsd=%d\n";
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getShortestDiagonal( int e, MV_ColMat<real>& c )
//-----------------------------------------------------------------------------
// used for finding shortest diagonal in tetrahedron
// find the shortest internal diagonal of 5-10, 6-8, and 7-9
// it is not so important to avoid the arbitrariness of the 
// equal case, since the choice is an internal matter for an element.
// c contains the nodal coorddinates (x,y,z)
{
  int n  = new_nodes( e, 0 ),  m  = new_nodes( e, 5 );
  real l1 = sqr(c(n,0)-c(m,0)) + sqr(c(n,1)-c(m,1)) + sqr(c(n,2)-c(m,2));
  n  = new_nodes( e, 1 );  m  = new_nodes( e, 3 );
  real l2 = sqr(c(n,0)-c(m,0)) + sqr(c(n,1)-c(m,1)) + sqr(c(n,2)-c(m,2));
  n  = new_nodes( e, 2 );  m  = new_nodes( e, 4 );
  real l3 = sqr(c(n,0)-c(m,0)) + sqr(c(n,1)-c(m,1)) + sqr(c(n,2)-c(m,2));

  //  if( l1 < l2 && l1 < l3 ) return 0;
  //  return l2 < l3 ? 1 : 2;
  if( l1 < l2 && l1 < l3 ) return 1;
  return l2 < l3 ? 2 : 3;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: getNodesOnEdge( int e, int edge_no, int& n1, int& n2 )
//-----------------------------------------------------------------------------
// get the global node numbers on local edge "edge_no" of a tetrahedron
// (only for ElmT4n3D) (no ordering of nodes on edge defined)
{
  switch ( edge_no ) {
  case 0: n1 = loc2glob(e,0); n2 = loc2glob(e,1); break;
  case 1: n1 = loc2glob(e,1); n2 = loc2glob(e,2); break;
  case 2: n1 = loc2glob(e,0); n2 = loc2glob(e,2); break;
  case 3: n1 = loc2glob(e,0); n2 = loc2glob(e,3); break;
  case 4: n1 = loc2glob(e,1); n2 = loc2glob(e,3); break;
  case 5: n1 = loc2glob(e,2); n2 = loc2glob(e,3); break;
  default: 
    //    errorFP("WavesGridFEHier:: getNodesOnEdge","illegal edge number %d",edge_no);
    cout<<"WavesGridFEHier:: getNodesOnEdge illegal edge number %d\n";
  }
}

//-----------------------------------------------------------------------------
int WavesGridFEHier::getElmAtFace(const int n1, const int n2, const int n3,
			     const int e)
//-----------------------------------------------------------------------------
{ 
// Given an element number e and the three global node numbers on the corners
// of an element face, compute the element number of the element 
// of the other side of the face.
// It is assumed that node-to-element (n2e) info for the grid is computed
// and that the elements are sorted in increasing order
// (which is the case for the present implementation of WavesNeighborFE).
// 0 is returned if no other element is found, i.e., is the face
// is a boundary face.

  //  WavesNeighborFE& n2e = additional_info.neighbor;
  WavesNeighborFE& n2e = (WavesNeighborFE&) neighbor;
  int n1stop = n2e.nodeIrow(n1+1);
  int n2stop = n2e.nodeIrow(n2+1);
  int n3stop = n2e.nodeIrow(n3+1);
  int i = n2e.nodeIrow(n1);
  int j = n2e.nodeIrow(n2);
  int k = n2e.nodeIrow(n3);
  for( ; i<n1stop; i++ ) {
    int v = n2e.nodeJcol(i);
    if( v == e ) continue; // already know this common element
    while( j < n2stop && v > n2e.nodeJcol(j) ) j++; 
    if( j >= n2stop ) return -1; // outside range of elms to node 2
    if( v == n2e.nodeJcol(j) ) { // two elements equal
      while( k < n3stop && v > n2e.nodeJcol(k) ) k++;
      if( k >= n3stop ) return -1; // outside range of elms to node 3
      if( v == n2e.nodeJcol(k) ) return v;// if the last elm also equal, bingo!
    }
  }
  return -1;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getElmAtFace(const int n1, const int n2, const int e)
//-----------------------------------------------------------------------------
{ 
// Version for ElmT10n3D element
// Same as getElmAtFace(const int n1, const int n2, const int n3, const int e),
// but now it is sufficient with two midpoint nodes to identify face.

  //  WavesNeighborFE& n2e = additional_info.neighbor;
  WavesNeighborFE& n2e = (WavesNeighborFE&) neighbor;
  int n1stop = n2e.nodeIrow(n1+1);
  int n2stop = n2e.nodeIrow(n2+1);
  int i = n2e.nodeIrow(n1);
  int j = n2e.nodeIrow(n2);
  for( ; i<n1stop; i++ ) {
    int v = n2e.nodeJcol(i);
    if( v == e ) continue; // already know this common element
    while( j < n2stop && v > n2e.nodeJcol(j) ) j++; 
    if( j >= n2stop ) return -1; // outside range of elms to node 2
    if( v == n2e.nodeJcol(j) ) return v; // if the last elm also equal, bingo!
  }
  return -1;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier::getElmsAtFaces(const int n1, const int n2, const int n3,
				const int n4, const int e, int& e3, int& e4)
//-----------------------------------------------------------------------------
{ 
// given an element number e and the four global node numbers on the corners
// of the element.
// Compute the element numbers of the elements of the other side of the two 
// faces which share the edge n1-n2.
// It is assumed that node-to-element (n2e) info for the grid is computed
// and that the elements their are sorted.
// 0 is returned if the neighbor element was not found, which means that
// a boundary face is found.

// e3 is the neighbor element to face n1 n2 n3 
// e4 is the neighbor element to face n1 n2 n4

  //  WavesNeighborFE& n2e = additional_info.neighbor;
  WavesNeighborFE& n2e = (WavesNeighborFE&) neighbor;
  int n1stop = n2e.nodeIrow(n1+1);
  int n2stop = n2e.nodeIrow(n2+1);
  int n3stop = n2e.nodeIrow(n3+1);
  int n4stop = n2e.nodeIrow(n4+1);
  int i = n2e.nodeIrow(n1); // start for first element n1 is part of
  int j = n2e.nodeIrow(n2);
  int k = n2e.nodeIrow(n3);
  int l = n2e.nodeIrow(n4);
  e3 = e4 = -1;
  for( ; i<n1stop; i++ ) {
    int v = n2e.nodeJcol(i);
    if( v == e ) continue; // already know this common element
    while( j < n2stop && v > n2e.nodeJcol(j) ) j++; 
    if( j >= n2stop ) return; // no more elms to node 2
    if( v == n2e.nodeJcol(j) ) { // both elements equal
      while( k < n3stop && v > n2e.nodeJcol(k) ) k++;
      if( k < n3stop && v == n2e.nodeJcol(k) ) { e3 = v; if(e4>-1) return;}
      while( l < n4stop && v > n2e.nodeJcol(l) ) l++;
      if( l < n4stop && v == n2e.nodeJcol(l) ) { e4 = v; if(e3>-1) return;}
    }
  }
}

//-----------------------------------------------------------------------------
void WavesGridFEHier::getElmsAtEdge( int e, int edge_no )
//-----------------------------------------------------------------------------
{ 
// efficient computation of all the elements which share the edge
// which has local edge number edge_no in element e.
// The local edge no of the elements are also computed.
// The result is stored last in two stacks: elmstack and ednostack,
// where stackcounter gives the position of the first position in the
// stack not used

  int n1,n2;
  getNodesOnEdge(e,edge_no,n1,n2);
  //  cout<<"getElmsAtEdge after getNodes "<<e<<" "<<edge_no<<"\n"<<flush;
  
  //  WavesNeighborFE& n2e = additional_info.neighbor;// assumed n2e data computed
  WavesNeighborFE& n2e = (WavesNeighborFE&) neighbor;
  //  cout<<"getElmsAtEdge neigh "<<e<<" "<<edge_no<<"\n"<<flush;
  int n1start = n2e.nodeIrow(n1);
  int n1stop = n2e.nodeIrow(n1+1);
  int n2start = n2e.nodeIrow(n2);
  int n2stop = n2e.nodeIrow(n2+1);
  //  int maxsize = max(n1stop-n1start,n2stop-n2start);
  int maxsize = (n1stop-n1start>n2stop-n2start)?(n1stop-n1start):
    (n2stop-n2start);
  //  cout<<"getElmsAtEdge max "<<e<<" "<<edge_no<<"\n"<<flush;
  int i;
  if ( elmstack.size() < stackcounter + maxsize ){ // dynamic size of stacks
// need to increase the size of the stacks
    //  cout<<"getElmsAtEdge nsize "<<maxsize<<" "<<elmstack.size()<<"\n"<<flush;
  
    int nsize = 2*elmstack.size() + maxsize; // new safe size for the stacks
// need scratch vector to move the data from the stack while the size of
// the stack is increased
    //    VecSimple(int) scratch(stackcounter-1); 
    MV_Vector<int> scratch(stackcounter); 
    for ( i=0; i<stackcounter; i++ ) scratch(i) = elmstack(i);
    elmstack.newsize(nsize);
    for ( i=0; i<stackcounter; i++ ) elmstack(i) = scratch(i);
    for ( i=0; i<stackcounter; i++ ) scratch(i) = ednostack(i);
    ednostack.newsize(nsize);
    for ( i=0; i<stackcounter; i++ ) ednostack(i) = scratch(i);
    //  cout<<"getElmsAtEdge nsize end"<<e<<" "<<edge_no<<"\n"<<flush;
  }
  int j=n2start;
  for( i=n1start; i<n1stop; i++ ) {
    int v = n2e.nodeJcol(i);
    while( j < n2stop && v > n2e.nodeJcol(j) ) j++;
    if(j >= n2stop) return;
    if(v == n2e.nodeJcol(j)){
      ednostack(stackcounter) = getLocEdgeNo(n1,n2,v); 
      elmstack(stackcounter++) = v;
    }
  }
  //  cout<<"getElmsAtEdge end "<<e<<" "<<edge_no<<"\n"<<flush;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getLocEdgeNo( int n1, int n2, int e )
//-----------------------------------------------------------------------------
// find the local edge number of element e which has nodes n1 and n2 
// used for ElmT4n3D
{ 
  /*
#ifdef ARRAY_RANGECHECK
// this version checks that n1 and n2 are nodes of the element
  int m1 = loc2glob(e,1);
  int m2 = loc2glob(e,2);
  if( ( n1 == m1 && n2 == m2 ) || ( n1 == m2 && n2 == m1 ) ) return 1;
  int m3 = loc2glob(e,3);
  if( n1 == m3 ) {
    if( n2 == m2 ) return 2;
    if( n2 == m1 ) return 3;    
  }
  if( n2 == m3 ) {
    if( n1 == m2 ) return 2;
    if( n1 == m1 ) return 3;    
  }
  int m4 = loc2glob(e,4);
  if( n1 == m4 ) {
    if( n2 == m1 ) return 4;
    if( n2 == m2 ) return 5;
    if( n2 == m3 ) return 6;
  }
  if( n2 == m4 ) {
    if( n1 == m1 ) return 4;
    if( n1 == m2 ) return 5;
    if( n1 == m3 ) return 6;
  }
  errorFP("WavesGridFEHier:: getLocEdgeNo",
	  "did not find edge in element %d with nodes %d and %d",e,n1,n2);
  return 0;
#else
*/
  // in this version it is assumed that n1 and n2 are nodes of element e
  // if not so the result will be wrong
  if( n1 == loc2glob(e,0) ) {
    if( n2 == loc2glob(e,1) ) return 0;
    if( n2 == loc2glob(e,2) ) return 2;
    return 3;
  }
  if( n1 == loc2glob(e,1) ) {
    if( n2 == loc2glob(e,0) ) return 0;
    if( n2 == loc2glob(e,2) ) return 1;
    return 4;
  }
  if( n1 == loc2glob(e,2) ) {
    if( n2 == loc2glob(e,0) ) return 2;
    if( n2 == loc2glob(e,1) ) return 1;
    return 5;
  }
  if( n2 == loc2glob(e,0) ) return 3;
  if( n2 == loc2glob(e,1) ) return 4;
  return 5;
  //#endif
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: numberingNewNodes( int e, int edge_no )
//-----------------------------------------------------------------------------
// Go through the list of elements and their 
// marked edges and assign each of them a unique node number.
// A marked edge has a no value before it is assigned a (positive)
// node number.
{
  switch ( maxnne ) { // could use enumtype element
  case 3: { //      if( "ElmT3n2D" == elm_tp)
    int midptno = incCurrNodeNo();// a new node at the midpoint
    new_nodes( e, edge_no ) = midptno; // assign the primary position
    int elnei = getElmNeighbor( e, edge_no );
    if( elnei>=0 ){ // if nonzero means that the edge is internal
      if( getElmNeighbor( elnei, 0 ) == e ) // find edge no for neighbor elm
	new_nodes( elnei, 0 ) = midptno;    // assign the secondary position
      else if( getElmNeighbor( elnei, 1 ) == e ) 
	new_nodes( elnei, 1 ) = midptno;
      else 
	new_nodes( elnei, 2 ) = midptno;
    }
  } break;
  case 4: { //  if( "ElmT4n3D" == elm_tp)
    //    cout<<"numberingNewNodes "<<e<<" "<<edge_no<<"\n"<<flush;
    int midptno = incCurrNodeNo();// a new node at the midpoint
    int i = stackcounter;
// note: stackcounter is a global variable which is increased in getElmsAtEdge
    getElmsAtEdge( e, edge_no ); // find tetrahedra which share the same edge
    //    cout<<"after getElmsnumberingNewNodes \n"<<flush;
    for( ; i<stackcounter; i++){
      //      printf(" %d i getElmsnumberingNewNodes \n",i);
      //printf(" %d elmstack(i) getElmsnumberingNewNodes \n",elmstack(i));
      //printf(" %d ednostack(i) getElmsnumberingNewNodes \n",ednostack(i));
      //printf(" %d midptno getElmsnumberingNewNodes \n",midptno);
      new_nodes(elmstack(i),ednostack(i)) = midptno;
    }
  } break;
  case 6: { //      if( "ElmT6n2D" == elm_tp)
    int no1 = incCurrNodeNo(); // prior to mdptno
    int no2 = incCurrNodeNo(); // after in counterclockwise order
    new_nodes( e, edge_no     ) = no1; 
    new_nodes( e, edge_no + 3 ) = no2; 
    int elnei = getElmNeighbor( e, edge_no );
    if( elnei>=0 ){ // if nonzero means that the edge is internal
      if( getElmNeighbor( elnei, 0 ) == e ) { // find edge no for neighbor elm
	new_nodes( elnei, 0 ) = no2;          //  reverse order on other side
	new_nodes( elnei, 3 ) = no1;
      }
      else if( getElmNeighbor( elnei, 2 ) == e ) {
	new_nodes( elnei, 1 ) = no2;
	new_nodes( elnei, 4 ) = no1;
      }
      else {
	new_nodes( elnei, 2 ) = no2;
	new_nodes( elnei, 5 ) = no1;
      }
    }
  } break;
  case 10: { // if( "ElmT10n3D" == elm_tp)
    //    fatalerrorFP("WavesGridFEHier:: numberingNewNodes",
    //		 "ElmT10n3D not implemented (yet)");
    cout<<"WavesGridFEHier:: numberingNewNodes"<<
		 "ElmT10n3D not implemented (yet)\n";
// find the elements and local edgeno and put on stack(easier in this case)
    int midptno = loc2glob(e,4+edge_no); // numbering of edges and quad nodes
    int ed1 = incCurrNodeNo();
    int ed2 = incCurrNodeNo();
    int i = stackcounter;
// find better way for T10 using midpoint nodes
    getElmsAtEdge( e, edge_no ); // find tetrahedra which share the same edge
    for( ; i<stackcounter; i++) {
      int ne = elmstack(i);
      int nedno = ednostack(i);
      new_nodes(ne,nedno)    = midptno;
      new_nodes(ne,6+nedno)  = ed1; // check that the ordering is correct
      new_nodes(ne,12+nedno) = ed2; // could maybe change later
    }
  } break;
  default:
    //  fatalerrorFP("WavesGridFEHier:: numberingNewNodes",
    //	       "What element? wrong nne=%d",nne);
    cout<<"WavesGridFEHier:: numberingNewNodes"<<
	       "What element? wrong nne=%d\n";
  }
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: getGlobNodesOnFace( int& n1, int& n2, int& n3, int f, int e )
//-----------------------------------------------------------------------------
// Get the global node numbers on face f of element e
{
  switch ( f ) {
  case 0: n1 = loc2glob(e,0); n2 = loc2glob(e,1); n3 = loc2glob(e,3);  break;
  case 1: n1 = loc2glob(e,1); n2 = loc2glob(e,2); n3 = loc2glob(e,3);  break;
  case 2: n1 = loc2glob(e,0); n2 = loc2glob(e,3); n3 = loc2glob(e,2);  break;
  case 3: n1 = loc2glob(e,0); n2 = loc2glob(e,2); n3 = loc2glob(e,1);  break;
  default:
    //    errorFP("WavesGridFEHier:: getNodesOnFace","illegal face number %d",f);
    cout<<"  WavesGridFEHier:: getNodesOnFace illegal face number %d\n";
  }
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getFaceNo( int n1, int n2, int n3, int e, bool safe )
//-----------------------------------------------------------------------------
// Get the face number of element e which consists of the nodes n1, n2 and n3
// if safe is true the function will return 0 if no match is found.
{
  if(safe){
// a version which is safe in the case that n1,n2,n3 are not nodes of element e
    int m1 = loc2glob(e,0), m2 = loc2glob(e,1);
    int m3 = loc2glob(e,2), m4 = loc2glob(e,3);
    int f1=0, f2=0, f3=0, f4=0;

    if( n1 == m1 )      { f1++; f3++; f4++; } // m1 is a part of face 1 3 and 4
    else if( n1 == m2 ) { f1++; f2++; f4++; } // m2 is a part of face 1 2 and 4
    else if( n1 == m3 ) { f2++; f3++; f4++; } // m3 is a part of face 2 3 and 4
    else if( n1 == m4 ) { f1++; f2++; f3++; } // m4 is a part of face 1 2 and 3

    if( n2 == m1 )      { f1++; f3++; f4++; }
    else if( n2 == m2 ) { f1++; f2++; f4++; }
    else if( n2 == m3 ) { f2++; f3++; f4++; }
    else if( n2 == m4 ) { f1++; f2++; f3++; }

    if( n3 == m1 )      { f1++; f3++; f4++; }
    else if( n3 == m2 ) { f1++; f2++; f4++; }
    else if( n3 == m3 ) { f2++; f3++; f4++; }
    else if( n3 == m4 ) { f1++; f2++; f3++; }

    if( f1 == 3 ) return 0; // if 3 hits for a face then correct face is found
    if( f2 == 3 ) return 1;
    if( f3 == 3 ) return 2;
    if( f4 == 3 ) return 3;

    cout<<"WavesGridFEHier:: getFaceNo, Looking for local face for nodes\n";
    cout<<n1<<" "<<n2<<" "<<n3<<" in element "<<e<<"\n";
    cout<<" the element contains nodes "
	<<m1<<" "<<m2<<" "<<m3<<" "<<m4<<"\n"<<flush;

    return -1; // did not find any match
  } else {
// This version is not safe if n1,n2,n3 are not in element e
// but is is faster than the safe version ...
    int n = loc2glob(e,0);
    if( n != n1 && n != n2 && n != n3 ) return 1;
    n = loc2glob(e,1);
    if( n != n1 && n != n2 && n != n3 ) return 2;
    n = loc2glob(e,2);
    if( n != n1 && n != n2 && n != n3 ) return 0;
    return 3;
  }
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getParent( const int e ) const
//-----------------------------------------------------------------------------
// Get the element number in prev_grid of the parent triangle of element e.
// Return 0 if there is no parent grid or if parent info is not computed
// use getParent_eff for efficiency if one knows that 
// isThisCoarsestGrid() || !isParentInfo() is false
{
  if( isThisCoarsestGrid() || !isParentInfo() ) return -1;
  return parent( e );
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getParent( const int rel_level, const int e ) const
//-----------------------------------------------------------------------------
// Get the parent element number of element e in grid rel_level levels up from
// from this grid.
// Return 0 if there is no parent grid or if parent info is not computed
// use getParent_eff for efficiency
{
  int parent_ele = e;
  WavesGridFEHier* g = (WavesGridFEHier*) this;
  for( int i=1; i<= rel_level; i++) {
    if( g->isThisCoarsestGrid() || ! g->isParentInfo() ) return -1;
    else {
      parent_ele = g->getParent_eff( parent_ele );
      g = g->getParentGrid();
    }
  }
  return parent_ele;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getParent_eff( const int rel_level, const int e ) const
//-----------------------------------------------------------------------------
// More efficient version of getParent( int rel_level, int e )
{
  int parent_ele = e;
  WavesGridFEHier* g = (WavesGridFEHier*) this;
  for( int i=1; i<= rel_level; i++){
    parent_ele = g->getParent_eff( parent_ele );
    g = g->getParentGrid();
  }
  return parent_ele;
}

//-----------------------------------------------------------------------------
int WavesGridFEHier:: getChild( const int e, const int child_no) const
//-----------------------------------------------------------------------------
// Get the element number in child grid of the child_no children
// of element e.
// Return 0 if there is no child grid or if children not computed or
// if the element does not have child_no children. 
// use getChild_eff for efficiency
{
  if( isThisFinestGrid() || ! isChildrenInfo() || 
      child_no > pos_child(e+1) - pos_child(e) || child_no < 1)
    return -1;
  return pos_child( e ) + child_no - 1;
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: checkParentInfo() const
//-----------------------------------------------------------------------------
// Check for each element in this grid that its parent 
// has the element as one of the children.
{
  if( isThisCoarsestGrid() ) return true; // no parent, nothing to check
  WavesGridFEHier* pgrid = getParentGrid();
  if( pgrid->isThisFinestGrid() )
    //    fatalerrorFP("WavesGridFEHier::checkParentInfo",
    //		 "Inconsistency, the two grids do not agree.");
    cout<<"WavesGridFEHier::checkParentInfo"<<
		 "Inconsistency, the two grids do not agree.\n";
  for( int e=0; e < nel; e++ )
    {
      int elm_no_par = getParent_eff(e);
      int no_of_ch = pgrid->getNoChildren( elm_no_par );
      bool found_child = false;
      for( int i=0; i<no_of_ch; i++ )
	if( pgrid->getChild_eff( elm_no_par, i ) == e ){ 
	  found_child = true;
	  break;
	}
      if( ! found_child ) 
	return false; // Could not find child
    }
  return true; // everything seems ok
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: checkChildrenInfo() const
//-----------------------------------------------------------------------------
// Check for all element in this grid that their children have the
// correct parent in this grid.
{
  //  cout<<"in checkChildrenInfo\n"<<flush;
  if( isThisFinestGrid() ) return true; //has no children, nothing to check
  WavesGridFEHier* cgrid = getChildGrid();
  if( cgrid->isThisCoarsestGrid() )
    //    fatalerrorFP("WavesGridFEHier::checkChildrenInfo",
    //		 "Inconsistency, the two grids do not agree.");
    cout<<"WavesGridFEHier::checkChildrenInfo"<<
		 "Inconsistency, the two grids do not agree.\n";
  for( int e=0; e<nel; e++ ) {
    //    cout<<e<<" in checkChildrenInfo\n"<<flush;
    int no_of_ch = getNoChildren(e);
    //    cout<<e<<" "<<no_of_ch<<" in checkChildrenInfo\n"<<flush;
    for( int i=0; i<no_of_ch; i++ )
      if( cgrid->getParent_eff( getChild_eff( e, i )) != e ) 
	return false; 
  }
  return true; // everything seems ok
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: checkElementOrientation () const
//-----------------------------------------------------------------------------
// Check for all elements in this grid that they are counterclockwise oriented.
// for ElmT6n2D check also the "quadratic nodes" i.e. the midnodes.
{
  if ( nsd==1 ) {
    //    if( "ElmB2n1D" == elm_tp ){
    if( maxnne==2 ){
      bool test = true;
      for( int e=0; e<nel; e++ )
	if( coord(loc2glob( e, 0 ),0) > coord(loc2glob( e, 1 ),0) ) {
	  //	  warningFP("WavesGridFEHier:: checkElementOrientation",
	  //		    "Wrong orientation in element %d",e);
	  cout<<"WavesGridFEHier:: checkElementOrientation"<<
		    "Wrong orientation in element %d\n";
	  test = false;
	}
      return test;
    }
    //    else if( "ElmB3n1D" == elm_tp ){
    else if( maxnne==3 ){
      bool test = true;
      for( int e=0; e<nel; e++ )
	if( coord(loc2glob( e, 0 ),0) > coord(loc2glob( e, 1 ),0) ||
	    coord(loc2glob( e, 1 ),0) > coord(loc2glob( e, 2 ),0) ) {
	  //	  warningFP("WavesGridFEHier:: checkElementOrientation",
	  //		    "Wrong orientation in element %d",e);
	  cout<<"WavesGridFEHier:: checkElementOrientation"<<
		    "Wrong orientation in element %d\n";
	  test = false;
	}
      return test;
    }
    else
      //      fatalerrorFP("WavesGridFEHier:: checkElementOrientation",
      //		   "Not implemented yet, element type = %s",elm_tp.chars() );
      cout<<"WavesGridFEHier:: checkElementOrientation"<<
		   "Not implemented yet, element type = %s\n";
  }
  else if ( nsd==2 ) {
    //    if( "ElmT3n2D" != elm_tp && "ElmT6n2D" != elm_tp )
    //      fatalerrorFP("WavesGridFEHier:: checkElementOrientation",
    //		   "Not implemented yet, element type = %s",elm_tp.chars() );
    if( maxnne!=3 && maxnne!=6 )
      cout<<"WavesGridFEHier:: checkElementOrientation"<<
		   "Not implemented yet, element type = %s\n";
    int n1,n2,n3;
    bool test = true;
    
    //    int times = ( "ElmT3n2D" == elm_tp ) ? 1 : 2;
    int times = ( maxnne==3 ) ? 1 : 2;
    
    for ( int i=1; i<=times; i++ ) {
      for( int e=0; e<nel; e++ ) {
	if(i==1){
	  n1 = loc2glob( e, 0 ); n2 = loc2glob( e, 1 ); n3 = loc2glob( e, 2 );
	} else {
	  n1 = loc2glob( e, 3 ); n2 = loc2glob( e, 4 ); n3 = loc2glob( e, 5 );
	}
	real x1 = coord(n1,0),  y1 = coord(n1,1);
	real x2 = coord(n2,0),  y2 = coord(n2,1);
	real x3 = coord(n3,0),  y3 = coord(n3,1);
	real det = (x1-x3)*(y2-y3) - (x2-x3)*(y1-y3);
	if( det < 1e-20 ) {
	  //	  warningFP("WavesGridFEHier:: checkElementOrientation",
	  //		    "Wrong orientation in element %d",e);
	  cout<<"WavesGridFEHier:: checkElementOrientation2D "<<
	    "Wrong orientation in element "<<e<<"\n"<<flush;
	  test = false;
	}
      }
    }
    return test;
  }
  else if ( nsd==3 ) {
    //    if( "ElmT4n3D" != elm_tp && "ElmT10n3D" != elm_tp )
    //      fatalerrorFP("WavesGridFEHier:: checkElementOrientation",
    //		   "Not implemented yet, element type = %s",elm_tp.chars() );
    if( maxnne!=4 && maxnne!=10 )
      cout<<"WavesGridFEHier:: checkElementOrientation"<<
	"Not implemented yet, element type = %s\n"<<flush;

    bool test = true;
    FET4n3D T4n3D(this);
    for( int e=0; e<nel; e++ ) {
      T4n3D.refill(e);
      if( T4n3D.det() < 1e-30 ) {
	//	cout<<"WavesGridFEHier:: checkElementOrientationWrong orientation in element %d\n";
	//	errorFP("WavesGridFEHier:: checkElementOrientation",
	//		"\nWrong orientation in element %d"
	//		" n1=%d, n2=%d, n3=%d, n4=%d",e,
	//		T4n3D.n1(),T4n3D.n2(),T4n3D.n3(),T4n3D.n4());
	printf("WavesGridFEHier:: checkElementOrientation\n Wrong orientation in element %d n1=%d, n2=%d, n3=%d, n4=%d\n\n",e,T4n3D.n1(),T4n3D.n2(),T4n3D.n3(),T4n3D.n4());
	printf("x1=%lf,y1=%lf,z1=%lf\n",T4n3D.x1(),T4n3D.y1(),T4n3D.z1());
	printf("x2=%lf,y2=%lf,z2=%lf\n",T4n3D.x2(),T4n3D.y2(),T4n3D.z2());
	printf("x3=%lf,y3=%lf,z3=%lf\n",T4n3D.x3(),T4n3D.y3(),T4n3D.z3());
	printf("x4=%lf,y4=%lf,z4=%lf\n",T4n3D.x4(),T4n3D.y4(),T4n3D.z4());
	printf("determinant=%lf\n\n",T4n3D.det());
	//	putLoc2glob(e,0,T4n3D.n2()); // does _not_ correct here
	//	putLoc2glob(e,1,T4n3D.n1());
	test = false;
	//	exit(1);
      }
    }
    return test;
  }
  return false; // should have returned before
}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridL2Q( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// construct a grid with piecewise quadratic elements from a linear grid
{
  switch (nsd) {
  case 1: return changeGridB2toB3( newgrid );
  case 2: return changeGridT3toT6( newgrid );
  case 3: return changeGridT4toT10( newgrid );
  default:
    //    fatalerrorFP("WavesGridFEHier::changeGridL2Q","nsd=%d not implemented yet",nsd);
    cout<<"WavesGridFEHier::changeGridL2Q"<<"nsd=%d not implemented yet\n";
  }
  return NULL;
}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridQ2L( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// construct a grid with piecewise linear elements from a quadratic grid
{ 
  switch (nsd) {
  case 1: return changeGridB3toB2( newgrid );
  case 2: return changeGridT6toT3( newgrid );
  case 3: return changeGridT10toT4( newgrid );
  default: 
    //    fatalerrorFP("WavesGridFEHier::changeGridQ2L","nsd=%d not implemented yet",nsd);
    cout<<"WavesGridFEHier::changeGridQ2L"<<"nsd=%d not implemented yet\n";
  }
  return NULL;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initL2QorQ2L( WavesGridFEHier* newgrid, const int new_nno )
//-----------------------------------------------------------------------------
// initialize changed grid. Linear to quadratic or quadratic to linear
{
  //  newgrid->setNonUniformMesh();

// copy the WavesGridFEHier specific information
  /*
  if( ! isThisFinestGrid() ){                       // there is finer grid...
    newgrid->setFinestGrid( false );
    newgrid->setChildGrid ( getChildGrid() ); // set pointer to child
    getChildGrid()->setParentGrid( newgrid ); // and back again
// Note that the old grid (this) will be incorrect since the grid it
// points to does not point back.
    if( isChildrenInfo() ){                   // children info exists
      newgrid->setChildrenInfo( true );
      newgrid->pos_child.newsize(nel+1);
      newgrid->pos_child = pos_child;
    }
  }

  if( ! isThisCoarsestGrid() ){ // there is a coarser grid...
    newgrid->setCoarsestGrid( false );
    newgrid->setParentGrid ( getParentGrid() );// set pointer to parent
    getParentGrid()->setChildGrid( newgrid);   // and back again
// Note that the old grid (this) will be incorrect since the grid it
// points to does not point back.
    if( isParentInfo() ){                // parent info exists
      newgrid->setParentInfo( true );
      newgrid->parent.newsize(nel);
      newgrid->parent = parent;
    }
  }
	
  removeHierInfo();  // maybe drastic?, but the info in old grid is now correct
  */
  if ( isNeighborsComputed() ) {
    newgrid->setNeighborsComputed( true );
    newgrid->elm_neigh.newsize(nel,nsd+1);
    newgrid->elm_neigh = elm_neigh;
  }

// will not copy new_nodes, (the matrix used to mark and number new nodes)

  newgrid->setRefinementMethod( getRefinementMethod() );

// copy the name of the boundary indicators
  int i,b;
  /*
  for ( i=1; i<=nbind; i++) 
    newgrid->putBoIndName( getBoIndName(i), i );

// copy the old nodes which will be present in the new grid
// copy boundary indicators    

  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
  if ( nsd==1 )
    for( i=0; i<new_nno; i++ )
      newgrid->coord(i,0) = coord(i,0); 
  else if ( nsd==2 )
    for( i=0; i<new_nno; i++ ){
      newgrid->coord(i,0) = coord(i,0); 
      newgrid->coord(i,1) = coord(i,1); 
    }
  else if ( nsd==3 )
    for( i=0; i<new_nno; i++ ){
      newgrid->coord(i,0) = coord(i,0); 
      newgrid->coord(i,1) = coord(i,1); 
      newgrid->coord(i,2) = coord(i,2); 
    }
  else
    //    fatalerrorFP("WavesGridFEHier::initL2QorQ2L","Wrong nsd=%d",nsd);
    cout<<"WavesGridFEHier::initL2QorQ2L"<<"Wrong nsd=%d\n";
  /*
  if ( nbind )
    for( i=1; i<=new_nno; i++ )
      for( b=1; b<=nbind; b++ )
	if( ind(i,b) == '1' ) 
	  nind(i,b) = '1';
	  */
  if( !onemat ) 
    newgrid->elid=elid;
}


// commented: should be written in WavesGridB.C elements with type
// ELMLINE2 
//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridB2toB3( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a grid with quadratic elements ElmB3n1d
// from linear elements ElmB2n1d 
{
  /*
  if( "ElmB3n1D" == elm_tp ) {
    warningFP("WavesGridFEHier:: changeGridB2toB3", "Already a ElmB3n1D grid");
    return this;
  }

  if( "ElmB2n1D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: changeGridB2toB3",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
  if( !isNeighborsComputed() ) // Is the element-element info computed?
    initNeighbor();

  int b,e;

  int oldrefmethod = getRefinementMethod();
  setRefinementMethod(1); // set a new node on each edge
  markAllElements();
  setRefinementMethod(oldrefmethod);

  int new_nno = nno + nel;  // one new node for each changed element

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  //  newgrid->redim( nsd, new_nno, nel, 3, nbind, 3, onesbd );
  newgrid->redim( nsd, new_nno, nel, 3, nbind, ELMLINE2 );
  //  newgrid->setElmType( "ElmB3n1D" );

  initL2QorQ2L(newgrid,nno);
  /*
  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
// Go through all elements and compute position of the new node
// at the midpoint

  for( e=0; e<nel; e++ ) {
    int n1 =  loc2glob( e, 0 );
    int n2 =  loc2glob( e, 1 );
    int m1 = incCurrNodeNo(); 
    newgrid->putLoc2Glob123 ( n1, m1, n2, e );
    newgrid->coord(m1,0) = (coord(n1,0)+coord(n2,0))/2;
    /*    for( b=1; b<=nbind; b++ )
      if( ind(n1,b) == '1' && ind(n2,b) == '1' ) // if both ends are set...
 	nind(m1,b) = '1'; // then set midpoint 
	*/
  }

  //#ifdef ARRAY_RANGECHECK
  if( !checkElementOrientation() )
    //    errorFP("WavesGridFEHier:: changeGridB2toB3", "Wrong element orientation.");
    cout<<"WavesGridFEHier:: changeGridB2toB3  Wrong element orientation\n";
  //#endif

  return newgrid;
}


//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridB3toB2( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a linear triangle elements grid ElmT3n2d from a grid 
// with quadratic elements ElmT6n2d.
{
  //  if( "ElmB2n1D" == elm_tp ) {
  //    warningFP("WavesGridFEHier:: changeGridT6toT3", "Already a ElmB2n1D grid");
  if( maxnne==2 ) {
    //    warningFP("WavesGridFEHier:: changeGridT6toT3", "Already a ElmB2n1D grid");
    cout<<"WavesGridFEHier:: changeGridT6toT3  Already a ElmB2n1D grid\n";
    return this;
  }

  //  if( "ElmB3n1D" != elm_tp )
  //    fatalerrorFP("WavesGridFEHier:: changeGridB3toB2",
  //		 "Not implemented yet, element type = %s",elm_tp.chars() );
  if( maxnne==3 )
    cout<<"WavesGridFEHier:: changeGridB3toB2"<<
		 "Not implemented yet, element type = %s\n";

  int e;

// find largest node number for the nodes to get size for node vector
  int new_nno = -1;
  for( e=nel-1; e>=0; e-- ) {
    if( loc2glob(e,0) > new_nno )  new_nno = loc2glob(e,0);
    if( loc2glob(e,2) > new_nno )  new_nno = loc2glob(e,2);
  }

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  //  newgrid->redim( nsd, new_nno, nel, 2, nbind, 2, onesbd );
  //  newgrid->setElmType( "ElmB2n1D" );
  newgrid->redim( nsd, new_nno, nel, 2, nbind, ELMLINE1 );
  //  newgrid->setElmType( "ElmB2n1D" );

  initL2QorQ2L(newgrid,new_nno);

// Go through all elements and copy elements

  for( e=0; e<nel; e++ ) 
    newgrid->putLoc2Glob12 ( loc2glob( e, 0 ), loc2glob( e, 2 ), e );

  return newgrid;
}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridT3toT6( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a grid with quadratic elements ElmT6n2d
// from alinear triangle elements ElmT3n2d  grid
{
  /*
  if( "ElmT6n2D" == elm_tp ) {
    warningFP("WavesGridFEHier:: changeGridT3toT6", "Already a ElmT6n2D grid");
    return this;
  }

  if( "ElmT3n2D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: changeGridT3toT6",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
  if( !isNeighborsComputed() )// Is the element-element info computed?
    initNeighbor();

  int oldrefmethod = getRefinementMethod();
  setRefinementMethod(3); // set a new node on each edge
  markAllElements();
  setRefinementMethod(oldrefmethod);

  // the numbering of the new nodes has incremented the current_node_number

  int new_nno = getCurrNodeNo()+1; 

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  // newgrid->redim( nsd, new_nno, nel, 6, nbind, 6, onesbd );
  newgrid->redim( nsd, new_nno, nel, 6, nbind, ELMTRI2 );
  //  newgrid->setElmType( "ElmT6n2D" );

  initL2QorQ2L(newgrid,nno);
  /*
  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
// Go through all elements and compute position of the new nodes
// at the midpoints of edges which have been assigned a node number.

  int b,e,n1,n2,n3,n4,n5,n6;

  for( e=0; e<nel; e++ ) {

    n1 =  loc2glob( e, 0 );  n2 =  loc2glob( e, 1 );  n3 =  loc2glob( e, 2 );
    n4 = new_nodes( e, 0 );  n5 = new_nodes( e, 1 );  n6 = new_nodes( e, 2 );

    newgrid->putLoc2Glob123 ( n1, n2, n3, e );
    newgrid->putLoc2Glob456 ( n4, n5, n6, e );
    
    if( ! newgrid->coord(n4,0) ){ // to avoid computing things twice
      newgrid->coord(n4,0) = (coord(n1,0)+coord(n2,0))/2;
      newgrid->coord(n4,1) = (coord(n1,1)+coord(n2,1))/2;
      /*      for( b=1; b<=nbind; b++ )
	if( ind(n1,b)=='1' && ind(n2,b)=='1' ) //if both ends are set...
	  nind(n4,b) = '1'; // set midpoint
	  */
    }
    if( ! newgrid->coord(n5,0) ){ 
      newgrid->coord(n5,0) = (coord(n2,0)+coord(n3,0))/2;
      newgrid->coord(n5,1) = (coord(n2,1)+coord(n3,1))/2;
      /*      for( b=1; b<=nbind; b++ )
	if( ind(n2,b)=='1' && ind(n3,b)=='1' )
	  nind(n5,b) = '1';
	  */
    }
    if( ! newgrid->coord(n6,0) ){ 
      newgrid->coord(n6,0) = (coord(n3,0)+coord(n1,0))/2;
      newgrid->coord(n6,1) = (coord(n3,1)+coord(n1,1))/2;
      /*      for( b=1; b<=nbind; b++ )
	if( ind(n3,b)=='1' && ind(n1,b)=='1' )
	  nind(n6,b) = '1';
	  */
    }
  }

  //#ifdef ARRAY_RANGECHECK
  if( !checkElementOrientation() )
    //    errorFP("WavesGridFEHier:: changeGridT3toT6", "Wrong orientation.");
    cout<<"WavesGridFEHier:: changeGridT3toT6 Wrong orientation.\n";
  //#endif

  return newgrid;
}


//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridT6toT3( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a linear triangle elements grid ElmT3n2d from a grid 
// with quadratic elements ElmT6n2d.
{
  /*
  if( "ElmT3n2D" == elm_tp ) {
    warningFP("WavesGridFEHier:: changeGridT6toT3", "Already a ElmT3n2D grid");
    return this;
  }

  if( "ElmT6n2D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: changeGridT6toT3",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
// find largest node number for the corner nodes to get size for node vector
  int e;
  int new_nno = -1;
  for( e=nel-1; e>=0; e-- ) {
    if( loc2glob(e,0) > new_nno )  new_nno = loc2glob(e,0);
    if( loc2glob(e,1) > new_nno )  new_nno = loc2glob(e,1);
    if( loc2glob(e,2) > new_nno )  new_nno = loc2glob(e,2);
  }

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  //  newgrid->redim( nsd, new_nno, nel, 3, nbind, 3, onesbd );
  newgrid->redim( nsd, new_nno, nel, 3, nbind, ELMTRI1 );
  //  newgrid->setElmType( "ElmT3n2D" );

  initL2QorQ2L(newgrid,new_nno);

// Go through all elements and copy elements

  for( e=0; e<nel; e++ ) 
    newgrid->putLoc2Glob123 ( loc2glob( e, 0 ),
			      loc2glob( e, 1 ),
			      loc2glob( e, 2 ), e );

  return newgrid;
}

//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridT4toT10( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a grid with quadratic tetrahedral elements ElmT10n3D
// from a grid with linear tetrahedral elements ElmT4n3D.
{
  /*
  if( "ElmT10n3D" == elm_tp ) {
    warningFP("WavesGridFEHier:: changeGridT4toT10","Already a ElmT10n3D grid");
    return this;
  }

  if( "ElmT4n3D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: changeGridT4toT10",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
  if( !isNeighborsComputed() )// Is the element-element info computed?
    initNeighbor();
  //  cout<<"DEBUG afterinit changeGridT4toT10\n"<<flush;

  int oldrefmethod = getRefinementMethod();
  setRefinementMethod(6);
  markAllElements();
  //  cout<<"DEBUG aftermarkAll changeGridT4toT10\n"<<flush;
  setRefinementMethod(oldrefmethod); // reset the old value

  // the numbering of the new nodes has incremented the current_node_number

  int new_nno = getCurrNodeNo()+1; 

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  //  newgrid->redim( nsd, new_nno, nel, 10, nbind, 10, onesbd );
  int new_maxnne = 10;
  newgrid->redim( nsd, new_nno, nel,  new_maxnne, nbind, ELMTET2 );
  //  newgrid->setElmType( "ElmT10n3D" );

  initL2QorQ2L(newgrid,nno);
  //  cout<<"DEBUG after initL2Qor changeGridT4toT10\n"<<flush;
  /*
  MatSimple(char)&  ind = boundaryData().indicatorAccess();
  MatSimple(char)&  nind = newgrid->boundaryData().indicatorAccess();
  */
// Go through all elements and compute position of the new nodes
// at the midpoints of edges which have been assigned a node number.

  int b,i,e,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,mid,m1,m2;

  for( e=0; e<nel; e++ ) {

    //    cout<<e<<" elem DEBUG changeGridT4toT10\n"<<flush;
    //    cout<<"size newnodes"<<new_nodes.size(0)<<" "<<new_nodes.size(1)<<"\n"<<flush;
    n1 =  loc2glob( e, 0 );  n2 =  loc2glob( e, 1 );
    n3 =  loc2glob( e, 2 );  n4 =  loc2glob( e, 3 );
    n5 = new_nodes( e, 0 );  n6 = new_nodes( e, 1 );
    n7 = new_nodes( e, 2 );  n8 = new_nodes( e, 3 );
    n9 = new_nodes( e, 4 );  n10= new_nodes( e, 5 );
    //    cout<<e<<"DEBUG 0 changeGridT4toT10\n"<<flush;
    newgrid->putLoc2Glob1234 ( n1, n2, n3, n4, e );
    newgrid->putLoc2Glob5678910 ( n5, n6, n7, n8, n9, n10, e );
    //    cout<<e<<"DEBUG 1 changeGridT4toT10\n"<<flush;

    for(i=1;i<=6;i++){
      switch (i) {
      case 1: mid=n5;  m1=n1;  m2=n2; break;
      case 2: mid=n6;  m1=n2;  m2=n3; break;
      case 3: mid=n7;  m1=n3;  m2=n1; break;
      case 4: mid=n8;  m1=n1;  m2=n4; break;
      case 5: mid=n9;  m1=n2;  m2=n4; break;
      case 6: mid=n10; m1=n3;  m2=n4; break;
      }
      //    cout<<e<<"DEBUG 2 changeGridT4toT10\n"<<flush;
      if( mid && ! newgrid->coord(mid,0) ){ // to avoid recomputing things
	newgrid->coord(mid,0) = (coord(m1,0)+coord(m2,0))/2;
	newgrid->coord(mid,1) = (coord(m1,1)+coord(m2,1))/2;
	newgrid->coord(mid,2) = (coord(m1,2)+coord(m2,2))/2;
	/*	for( b=1; b<=nbind; b++ )
	  if( ind(m1,b)=='1' && ind(m2,b)=='1' ) // if both ends are set...
	    nind( mid, b )='1'; // then set midpoint
	    */
      }
    }
  }
  //  cout<<"DEBUG after loop changeGridT4toT10\n"<<flush;

  //#ifdef ARRAY_RANGECHECK
  if( !checkElementOrientation() )
    //    errorFP("WavesGridFEHier:: changeGridT4toT10", "Wrong orientation.");
    cout<<"WavesGridFEHier:: changeGridT4toT10"<< "Wrong orientation.\n";
  //#endif

  return newgrid;
}


//-----------------------------------------------------------------------------
WavesGridFEHier* WavesGridFEHier:: changeGridT10toT4( WavesGridFEHier* newgrid )
//-----------------------------------------------------------------------------
// Create a grid with linear tetrahedral grid ElmT4n3D
// from a grid with quadratic elements ElmT10n3D.
{
  /*
  if( "ElmT4n3D" == elm_tp )
    {
      warningFP("WavesGridFEHier:: changeGridT10toT4",
		"Already a ElmT4n3D grid");
      return this;
    }

  if( "ElmT10n3D" != elm_tp )
    fatalerrorFP("WavesGridFEHier:: changeGridT10toT4",
		 "Not implemented yet, element type = %s",elm_tp.chars() );
		 */
// find largest node number for the corner nodes to get size for node vector

  int e;
  int new_nno = -1;
  for( e=nel-1; e>=0; e-- ) {
    if( loc2glob(e,0) > new_nno )  new_nno = loc2glob(e,0);
    if( loc2glob(e,1) > new_nno )  new_nno = loc2glob(e,1);
    if( loc2glob(e,2) > new_nno )  new_nno = loc2glob(e,2);
    if( loc2glob(e,3) > new_nno )  new_nno = loc2glob(e,3);
  }

  //  WavesGridFEHier* newgrid = new WavesGridFEHier;
  //  newgrid->redim( nsd, new_nno, nel, 4, nbind, 4, onesbd );
  newgrid->redim( nsd, new_nno, nel, 4, nbind, ELMTET1 );
  //  newgrid->setElmType( "ElmT4n3D" );

  initL2QorQ2L(newgrid,new_nno);

// Go through all elements and copy elements

  for( e=0; e<nel; e++ )
    newgrid->putLoc2Glob1234 ( loc2glob( e, 0 ),
			       loc2glob( e, 1 ),
			       loc2glob( e, 2 ),
			       loc2glob( e, 3 ), e);

  return newgrid;
}


//-----------------------------------------------------------------------------
real WavesGridFEHier:: getMinEta( MV_Vector<int>& element_list, int qualcrit )
//-----------------------------------------------------------------------------
{
  real mineta = 99999;
  int length = element_list.size();
  if(length==0) 
    printf("WavesGridFEHier:: getMinEta, length is 0");
  int e;
  real eta;
  if ( nsd==2 ) {
    for( e=0; e<length; e++ ) {
      FET3.refill(element_list(e));
      if( FET3.det() <= 0.0 ) 
	eta = -1; // illegal grid, could happen in smoothing 
      else 
	eta = FET3.qualityMeasure(qualcrit); 
      if( eta < mineta ) mineta = eta; 
    } 
  }
  else if( nsd==3 ) { 
    for( e=0; e<length; e++ ) {
      FET4.refill(element_list(e));
      if( FET4.det() <= 0.0 ) 
	eta = -1; // illegal grid, could happen in smoothing 
      else 
	eta = FET4.qualityMeasure(qualcrit); 
      if( eta < mineta ) mineta = eta; 
    } 
  }
  else 
    printf("WavesGridFEHier:: getMinEta , Wrong dimension %d",nsd); 
  return mineta; 
}


//-----------------------------------------------------------------------------
void WavesGridFEHier::makeStatistics ()
//-----------------------------------------------------------------------------
{ // ksa: this one is not finished add more statistics
  printf("Begin makeStatistics\n");
  real maxlength = 0;
  real minlength = 99999;
  real minangle = 99999;
  real maxangle = -99999;
  real maxeta = 0;
  real aveeta=0.0;
  real mineta = 99999;
  real maxlen,minlen,maxang,minang;
  
  int e;
  if ( nsd==2 ) { 
    for ( e=0; e<nel; e++ ) {
      FET3.refill(e); 
      if( FET3.det() <= 0.0 ) { 
	printf("WavesGridFEHier:: makeStatistics Wrong orientation in element %d, n1=%d,n2=%d,n3=%d,",e, FET3.n1(),FET3.n2(),FET3.n3());
      }
      maxang = FET3.maxAngle(); 
      minang = FET3.minAngle(); 
      real eta = FET3.qualityMeasure(1);
      aveeta += eta;
      maxlen = FET3.maxLength(); 
      minlen = FET3.minLength(); 
      if( eta > maxeta ) maxeta = eta; 
      if( eta < mineta ) mineta = eta; 
      if( maxlen > maxlength ) maxlength = maxlen; 
      if( minlen < minlength ) minlength = minlen; 
      if( maxang > maxangle ) maxangle = maxang; 
      if( minang < minangle ) minangle = minang; 
    }
  }
  else if ( nsd==3 ) {
    FET4.attach(this);
    for ( e=0; e<nel; e++ ) {
      //      printf(" ele %i\n",e);
      FET4.refill(e); 
      if( FET4.det() <= 0.0 ) { 
	printf("WavesGridFEHier:: makeStatistics");
	printf("Wrong orientation in element %d, n1=%d,n2=%d,n3=%d,n4=%d", 
	       e,FET4.n1(),FET4.n2(),FET4.n3(),FET4.n4());
      } 
      maxang = FET4.maxSolidAngle(); 
      minang = FET4.minSolidAngle(); 
      real eta = FET4.qualityMeasure(1); 
      aveeta += eta;
      maxlen = FET4.maxLength();
      minlen = FET4.minLength();  
      if( eta > maxeta ) maxeta = eta; 
      if( eta < mineta ) mineta = eta; 
      if( maxlen > maxlength ) maxlength = maxlen; 
      if( minlen < minlength ) minlength = minlen; 
      if( maxang > maxangle ) maxangle = maxang; 
      if( minang < minangle ) minangle = minang; 
    } 
  } 
  aveeta /= real(nel);
  cout<<"\n";
  cout<<"maxlength = "<<maxlength<<"\n";
  cout<<"minlength = "<<minlength<<"\n";
  cout<<"max angle = "<<maxangle<<"\n";
  cout<<"min angle = "<<minangle<<"\n";
  cout<<"max eta (crit 1)= "<<maxeta<<"\n";
  cout<<"ave eta (crit 1)= "<<aveeta<<"\n";
  cout<<"min eta (crit 1)= "<<mineta<<"\n";
  cout<<"\n"<<flush;
}

/*
//-----------------------------------------------------------------------------
int WavesGridFEHier:: findElm( Ptv(real)& x, int level0, int elm0, int level1 )
//-----------------------------------------------------------------------------
// given a point (x) which is in elm0 on level level0, find the element
// which the point is in in level1
{
  if(level1<level0) 
    warningFP("WavesGridFEHier:: findElm",
	      "getParent is more efficient if level1<level0");
  WavesGridFEHier* grid = getCoarsestGrid();
  int i;
  for(i=1;i<level0;i++) grid = grid->getChildGrid();
  int element = elm0;
  if( ! grid->isInsideElm( x, element ) )
    fatalerrorFP("WavesGridFEHier:: findElm","The point must be in element!");
  for( i = level0; i<level1; i++) {
    if( grid->isChildrenInfo() ){
      int firstchildele = grid->getFirstChild(element);
      int lastchildele = grid->getFirstChild(element+1);
      grid = grid->getChildGrid();
// note that last element can also be tested
      for( element=firstchildele; element<lastchildele; element++ )
	if( grid->isInsideElm( x, element ) )
	  break;
      if( element == lastchildele ) 
	fatalerrorFP("WavesGridFEHier:: findElm","failed to find child!");
    }
    else
      fatalerrorFP("WavesGridFEHier:: findElm", 
		   "No children info exists, level=%d",i);
  }
  return element;
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: findElm( Ptv(real)& x, int level, int& element)
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = getCoarsestGrid();
  if( ! grid->simpleSearch(x,element) )
    return false;
  for( int level_counter = 1; level_counter<level; level_counter++) {
    if( grid->isChildrenInfo() ){
      int firstchildele = grid->getFirstChild(element);
      int lastchildele = grid->getFirstChild(element+1) - 1;
      grid =  grid->getChildGrid();
// note that last element is not tested
      for( element=firstchildele; element<lastchildele; element++ )
	if( grid->isInsideElm( x, element ) ) 
	  break;
    }
    else
      fatalerrorFP("WavesGridFEHier:: findElm", "No children info exists");
  }
  return true;
}


//-----------------------------------------------------------------------------
bool WavesGridFEHier:: simpleSearch( Ptv(real)& x, int& element)
//-----------------------------------------------------------------------------
{
  for( element=1; element<=nel; element++ ) 
    if( isInsideElm( x, element ) ) return true;
  return false;
}

//-----------------------------------------------------------------------------
bool WavesGridFEHier:: isInsideElm( Ptv(real)& gpt, int element )
//-----------------------------------------------------------------------------
{
  if ( nsd==1 ) {
    real x = gpt(1);
    int n1,n2;
    real x1,x2;
    real df; // safety factor to test if inside element
  
    if( isBoundingBoxesComputed() )
      { 
        box_counter++;
        if( x < bounding_box(element,1) || x > bounding_box(element,2) )
	  return false;
        df = 1e-10*exp(getGridLevelNo()*0.742);
      }
    else
      {
        box_counter++;
        df = 2e-10*exp(getGridLevelNo()*0.742);
      
        n1 = loc2glob( element, 1 );
        n2 = loc2glob( element, getNoNodesInElm(element) ); 
        // works for both linear and quadratic
      
        x1 = getCoor(n1,1);
        x2 = getCoor(n2,1);
      
        real xmin,xmax;
        if ( x1<x2 ){ xmin=x1; xmax=x2; }
        else   {      xmin=x2; xmax=x1; }
        real xdiff = (xmax-xmin)*df;
        if( x < xmin - xdiff || x > xmax + xdiff ) 
	  return false;
        df *= 0.5; //using a smaller value of safety factor for the next test
      }
    in_elm_counter++;
  }
  else if ( nsd==2 ) {
    real x = gpt(1), y = gpt(2);
    int n;
    real x1,x2,x3,y1,y2,y3;
    real df; // safety factor to test if inside element
  
    if( isBoundingBoxesComputed() )
      { 
        box_counter++;
        if( x < bounding_box(element,1) || x > bounding_box(element,2) ||
	    y < bounding_box(element,3) || y > bounding_box(element,4) )
	  return false;
	// prepare for element test
        df = 1e-10*exp(getGridLevelNo()*0.742);
        n = loc2glob( element, 1 );
        x1 = getCoor(n,1); y1 = getCoor(n,2);
        n = loc2glob( element, 2 );
        x2 = getCoor(n,1); y2 = getCoor(n,2);
        n = loc2glob( element, 3 );
        x3 = getCoor(n,1); y3 = getCoor(n,2);
      }
    else
      {
        box_counter++;
        df = 2e-10*exp(getGridLevelNo()*0.742);
        
        n = loc2glob( element, 1 );
        x1 = getCoor(n,1); y1 = getCoor(n,2);
        n = loc2glob( element, 2 );
        x2 = getCoor(n,1); y2 = getCoor(n,2);
        n = loc2glob( element, 3 );
        x3 = getCoor(n,1); y3 = getCoor(n,2);
        
        real xmin,xmax;
        if ( x1<x2 ){ xmin=x1; xmax=x2; }
        else   {      xmin=x2; xmax=x1; }
        if (x3<xmin) xmin=x3;
        else if (x3>xmax) xmax=x3;
        real xdiff = (xmax-xmin)*df;
        if( x < xmin - xdiff || x > xmax + xdiff ) 
  	  return false;
        
        real ymin,ymax;
        if ( y1<y2 ){ ymin=y1; ymax=y2; }
        else   {      ymin=y2; ymax=y1; }
        if (y3<ymin) ymin=y3;
        else if (y3>ymax) ymax=y3;
        real ydiff = (ymax-ymin)*df;
        if( y < ymin - ydiff || y > ymax + ydiff ) 
  	  return false;
        df *= 0.5; //using a smaller value of safety factor for the next test
    }
  
    in_elm_counter++;
  
    real x13 = x1 - x3,   y31 = y3 - y1,   x32 = x3 - x2,   y23 = y2 - y3;
    real det = x13 * y23 - x32 * y31;
    
    real meps = -df*det;
    real x03 = x - x3, y03 = y - y3;
    real detphi1 = y23 * x03 + x32 * y03; 
    if( detphi1 < meps )  
      return false;
    real detphi2 = y31 * x03 + x13 * y03;
    if( detphi2 < meps )  
      return false;
    real detphi3 = det - detphi1 - detphi2;
    if( detphi3 < meps )  
      return false;
  }
  else if ( nsd==3 ) {
    real x = gpt(1), y = gpt(2), z = gpt(3);
    real df; // safety factor to test if inside element
    
    if( isBoundingBoxesComputed() )
      {
        box_counter++;
        if( x < bounding_box(element,1) || x > bounding_box(element,2) ||
  	    y < bounding_box(element,3) || y > bounding_box(element,4) ||
  	    z < bounding_box(element,5) || z > bounding_box(element,6) )
  	return false;
        df = 1e-10*exp(getGridLevelNo()*0.742);
      }
    else
      {
        box_counter++;
	int n1,n2,n3,n4;
	real x1,x2,x3,x4, pt,min,max,diff;
        df = 2e-10*exp(getGridLevelNo()*0.742);
        n1 = loc2glob( element, 1 );
        n2 = loc2glob( element, 2 );
        n3 = loc2glob( element, 3 );
        n4 = loc2glob( element, 4 );
  
        for(int d=1;d<=nsd;d++){
  	  x1 = getCoor(n1,d);
  	  x2 = getCoor(n2,d); 
  	  x3 = getCoor(n3,d);
  	  x4 = getCoor(n4,d);
  	  if ( x1<x2 ){ min=x1; max=x2; }
  	  else   {      min=x2; max=x1; }
  	  if (x3<min) min=x3;
  	  else if (x3>max) max=x3;
  	  if (x4<min) min=x4;
  	  else if (x4>max) max=x4;
  	  diff = (max-min)*df;
  	  pt = gpt(d);
  	  if( pt < min - diff || pt > max + diff ) 
	    return false;
        }
        df*=0.5; // use a stricter rule for inclusion
      } 
    in_elm_counter++;

    FET4.refill(element);
    real u,r,s,t;
    FET4.barycentric(x,y,z,u,r,s,t);
    real mdf = -df;
    if( u < mdf || r < mdf || s < mdf || t < mdf )
      return false; // the point is outside the tetrahedron
  }
// if not returned yet, the point (x) is inside the element
  return true;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initSearchStatAllGrids ()
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = getFinestGrid();
  while( ! grid->isThisCoarsestGrid() ){
    grid->initSearchStat();
    grid = grid->getParentGrid();
  }
  grid->initSearchStat();
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: printSearchStatAllGrids (Os os)
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = getFinestGrid();
  int no_of_levels = grid->getGridLevelNo();
  int total_box=0, total_elm_tests=0;
  while( no_of_levels-- ){
    int boxtests = grid->getNoBoxTests ();
    int eletests = grid->getNoElmTests ();
    total_box += boxtests;
    total_elm_tests += eletests;
    os <<" level = "<< no_of_levels+1
       <<" box tests = "<< boxtests 
       <<" ele tests = "<< eletests <<"\n";
    grid = grid->getParentGrid();
  }
  os <<" total on all levels: "
     <<" total box tests = "<< total_box
     <<" total ele tests = "<< total_elm_tests <<"\n";
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeBoundingBoxes()
//-----------------------------------------------------------------------------
{
  bounding_box.newsize(0,0);
  setBoundingBoxesComputed(false);
}


//-----------------------------------------------------------------------------
void WavesGridFEHier:: computeBoundingBoxes()
//-----------------------------------------------------------------------------
{
  bounding_box.newsize(nel,2*nsd);
  setBoundingBoxesComputed(true);
// taking a safety factor twice the one used in the inside test.
  const real df = 2*1e-10*exp(getGridLevelNo()*0.742);
  if ( nsd==1 ) {
    real xmin,xmax;
    
    for( int e=1; e<=nel; e++ ) {
      int n1 = loc2glob( e, 1 );
      int n2 = loc2glob( e, getNoNodesInElm(e) );// works for ElmB2n and ElmB3n
      
      real x1 = getCoor(n1,1);
      real x2 = getCoor(n2,1);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      
      real xdiff = (xmax-xmin)*df;
      
      bounding_box(e,1)=xmin-xdiff;
      bounding_box(e,2)=xmax+xdiff;
    }
  }
  else if ( nsd==2 ) {
    
    real xmin,xmax,ymin,ymax;
    
    for( int e=1; e<=nel; e++ ) {
      int n = loc2glob( e, 1 );
      real x1 = getCoor(n,1),  y1 = getCoor(n,2);
      n = loc2glob( e, 2 );
      real x2 = getCoor(n,1),  y2 = getCoor(n,2);
      n = loc2glob( e, 3 );
      real x3 = getCoor(n,1),  y3 = getCoor(n,2);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      if ( x3<xmin )         xmin=x3;
      else if  ( x3>xmax )   xmax=x3;
      
      if ( y1<y2 ){ ymin=y1; ymax=y2; }
      else   {      ymin=y2; ymax=y1; }
      if ( y3<ymin )         ymin=y3;
      else if  ( y3>ymax )   ymax=y3;
      
      real xdiff = (xmax-xmin)*df;
      real ydiff = (ymax-ymin)*df;
      
      bounding_box(e,1)=xmin-xdiff;
      bounding_box(e,2)=xmax+xdiff;
      bounding_box(e,3)=ymin-ydiff;
      bounding_box(e,4)=ymax+ydiff;
    }
  }
  else {  // if(nsd==3) {
    
    real xmin,xmax,ymin,ymax,zmin,zmax;
    
    for( int e=1; e<=nel; e++ ) {
      
      int n = loc2glob( e, 1 );
      real x1 = coor(n,1), y1 = coord(n,2), z1 = coord(n,3);
      n = loc2glob( e, 2 );
      real x2 = coord(n,1), y2 = coord(n,2), z2 = coord(n,3);
      n = loc2glob( e, 3 );
      real x3 = coord(n,1), y3 = coord(n,2), z3 = coord(n,3);
      n = loc2glob( e, 4 );
      real x4 = coord(n,1), y4 = coord(n,2), z4 = coord(n,3);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      if ( x3<xmin )         xmin=x3;
      else if  ( x3>xmax )   xmax=x3;
      if ( x4<xmin )         xmin=x4;
      else if  ( x4>xmax )   xmax=x4;
      
      if ( y1<y2 ){ ymin=y1; ymax=y2; }
      else   {      ymin=y2; ymax=y1; }
      if ( y3<ymin )         ymin=y3;
      else if  ( y3>ymax )   ymax=y3;
      if ( y4<ymin )         ymin=y4;
      else if  ( y4>ymax )   ymax=y4;
      
      if ( z1<z2 ){ zmin=z1; zmax=z2; }
      else   {      zmin=z2; zmax=z1; }
      if ( z3<zmin )         zmin=z3;
      else if  ( z3>zmax )   zmax=z3;
      if ( z4<zmin )         zmin=z4;
      else if  ( z4>zmax )   zmax=z4;
      
      real xdiff = (xmax-xmin)*df;
      real ydiff = (ymax-ymin)*df;
      real zdiff = (zmax-zmin)*df;
      
      bounding_box(e,1)=xmin-xdiff;
      bounding_box(e,2)=xmax+xdiff;
      bounding_box(e,3)=ymin-ydiff;
      bounding_box(e,4)=ymax+ydiff;
      bounding_box(e,5)=zmin-zdiff;
      bounding_box(e,6)=zmax+zdiff;
    }
  }
}
*/

//-----------------------------------------------------------------------------
// int WavesGridFEHier:: findElm( Ptv(real)& x, int level0, int elm0, int level1 )
int WavesGridFEHier:: findElm( MV_Vector<real>& x, int level0, int elm0, int level1 )
//-----------------------------------------------------------------------------
// given a point (x) which is in elm0 on level level0, find the element
// which the point is in in level1
{
  if(level1<level0) 
    //    warningFP("WavesGridFEHier:: findElm",
    //	      "getParent is more efficient if level1<level0");
    printf("WavesGridFEHier:: findElm getParent is more efficient if level1<level0\n");
  WavesGridFEHier* grid = (WavesGridFEHier*) getCoarsestGrid();
  int i;
  for(i=1;i<level0;i++) grid = grid->getChildGrid();
  int element = elm0;
  if( ! grid->isInsideElm( x, element ) )
    //    fatalerrorFP("WavesGridFEHier:: findElm","The point must be in element!");
    printf("WavesGridFEHier:: findElm , The point must be in element!\n");
  for( i = level0; i<level1; i++) {
    if( grid->isChildrenInfo() ){
      int firstchildele = grid->getFirstChild(element);
      int lastchildele = grid->getFirstChild(element+1);
      grid = grid->getChildGrid();
// note that last element can also be tested
      for( element=firstchildele; element<lastchildele; element++ )
	if( grid->isInsideElm( x, element ) )
	  break;
      if( element == lastchildele ) 
	//	fatalerrorFP("WavesGridFEHier:: findElm","failed to find child!");
	printf("WavesGridFEHier:: findElm, failed to find child!\n");
    }
    else
      //      fatalerrorFP("WavesGridFEHier:: findElm", 
      //		   "No children info exists, level=%d",i);
      printf("WavesGridFEHier:: findElm , No children info exists, level=%d\n",i);
  }
  return element;
}

//-----------------------------------------------------------------------------
//bool WavesGridFEHier:: findElm( Ptv(real)& x, int level, int& element)
bool WavesGridFEHier:: findElm( MV_Vector<real>& x, int level, int& element)
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = (WavesGridFEHier*) getCoarsestGrid();
  if( ! grid->simpleSearch(x,element) )
    return false;
  for( int level_counter = 1; level_counter<level; level_counter++) {
    if( grid->isChildrenInfo() ){
      int firstchildele = grid->getFirstChild(element);
      int lastchildele = grid->getFirstChild(element+1) - 1;
      grid =  grid->getChildGrid();
// note that last element is not tested
      for( element=firstchildele; element<lastchildele; element++ )
	if( grid->isInsideElm( x, element ) ) 
	  break;
    }
    else
      //      fatalerrorFP("WavesGridFEHier:: findElm", "No children info exists");
      printf("WavesGridFEHier:: findElm, No children info exists\n");
  }
  return true;
}


//-----------------------------------------------------------------------------
//bool WavesGridFEHier:: simpleSearch( Ptv(real)& x, int& element)
bool WavesGridFEHier:: simpleSearch( MV_Vector<real>& x, int& element)
//-----------------------------------------------------------------------------
{
  for( element=0; element<nel; element++ ) 
    if( isInsideElm( x, element ) ) return true;
  return false;
}

//-----------------------------------------------------------------------------
//bool WavesGridFEHier:: isInsideElm( Ptv(real)& gpt, int element )
bool WavesGridFEHier:: isInsideElm(  MV_Vector<real>& gpt, int element )
//-----------------------------------------------------------------------------
{
  if ( nsd==1 ) {
    //    real x = gpt(1);
    real x = gpt(0);
    int n1,n2;
    real x1,x2;
    real df; // safety factor to test if inside element
  
    if( isBoundingBoxesComputed() )
      { 
        box_counter++;
        if( x < bounding_box(element,0) || x > bounding_box(element,1) )
	  return false;
        df = 1e-10*exp(getGridLevelNo()*0.742);
      }
    else
      {
        box_counter++;
        df = 2e-10*exp(getGridLevelNo()*0.742);
      
	//        n1 = loc2glob( element, 1 );
        n1 = loc2glob( element, 0 );
	//        n2 = loc2glob( element, getNoNodesInElm(element) ); 
        n2 = loc2glob( element, getMaxNoNodesInElm()-1 ); 
        // works for both linear and quadratic
      
	//        x1 = getCoor(n1,1);
	//        x2 = getCoor(n2,1);
        x1 = getCoor(n1,0);
        x2 = getCoor(n2,0);
      
        real xmin,xmax;
        if ( x1<x2 ){ xmin=x1; xmax=x2; }
        else   {      xmin=x2; xmax=x1; }
        real xdiff = (xmax-xmin)*df;
        if( x < xmin - xdiff || x > xmax + xdiff ) 
	  return false;
        df *= 0.5; //using a smaller value of safety factor for the next test
      }
    in_elm_counter++;
  }
  else if ( nsd==2 ) {
    //    real x = gpt(1), y = gpt(2);
    real x = gpt(0), y = gpt(1);
    int n;
    real x1,x2,x3,y1,y2,y3;
    real df; // safety factor to test if inside element
  
    if( isBoundingBoxesComputed() )
      { 
        box_counter++;
	//        if( x < bounding_box(element,1) || x > bounding_box(element,2) ||
	//	    y < bounding_box(element,3) || y > bounding_box(element,4) )
        if( x < bounding_box(element,0) || x > bounding_box(element,1) ||
	    y < bounding_box(element,2) || y > bounding_box(element,3) )
	  return false;
	// prepare for element test
        df = 1e-10*exp(getGridLevelNo()*0.742);
	//        n = loc2glob( element, 1 );
//        x1 = getCoor(n,1); y1 = getCoor(n,2);
//        n = loc2glob( element, 2 );
//        x2 = getCoor(n,1); y2 = getCoor(n,2);
//        n = loc2glob( element, 3 );
//        x3 = getCoor(n,1); y3 = getCoor(n,2);

        n = loc2glob( element, 0 );
        x1 = getCoor(n,0); y1 = getCoor(n,1);
        n = loc2glob( element, 1 );
        x2 = getCoor(n,0); y2 = getCoor(n,1);
        n = loc2glob( element, 2 );
        x3 = getCoor(n,0); y3 = getCoor(n,1);
      }
    else
      {
        box_counter++;
        df = 2e-10*exp(getGridLevelNo()*0.742);
        
	//        n = loc2glob( element, 1 );
	//        x1 = getCoor(n,1); y1 = getCoor(n,2);
	//        n = loc2glob( element, 2 );
	//        x2 = getCoor(n,1); y2 = getCoor(n,2);
	//        n = loc2glob( element, 3 );
	//        x3 = getCoor(n,1); y3 = getCoor(n,2);

        n = loc2glob( element, 0 );
        x1 = getCoor(n,0); y1 = getCoor(n,1);
        n = loc2glob( element, 1 );
        x2 = getCoor(n,0); y2 = getCoor(n,1);
        n = loc2glob( element, 2 );
        x3 = getCoor(n,0); y3 = getCoor(n,1);
        
        real xmin,xmax;
        if ( x1<x2 ){ xmin=x1; xmax=x2; }
        else   {      xmin=x2; xmax=x1; }
        if (x3<xmin) xmin=x3;
        else if (x3>xmax) xmax=x3;
        real xdiff = (xmax-xmin)*df;
        if( x < xmin - xdiff || x > xmax + xdiff ) 
  	  return false;
        
        real ymin,ymax;
        if ( y1<y2 ){ ymin=y1; ymax=y2; }
        else   {      ymin=y2; ymax=y1; }
        if (y3<ymin) ymin=y3;
        else if (y3>ymax) ymax=y3;
        real ydiff = (ymax-ymin)*df;
        if( y < ymin - ydiff || y > ymax + ydiff ) 
  	  return false;
        df *= 0.5; //using a smaller value of safety factor for the next test
    }
  
    in_elm_counter++;
  
    real x13 = x1 - x3,   y31 = y3 - y1,   x32 = x3 - x2,   y23 = y2 - y3;
    real det = x13 * y23 - x32 * y31;
    
    real meps = -df*det;
    real x03 = x - x3, y03 = y - y3;
    real detphi1 = y23 * x03 + x32 * y03; 
    if( detphi1 < meps )  
      return false;
    real detphi2 = y31 * x03 + x13 * y03;
    if( detphi2 < meps )  
      return false;
    real detphi3 = det - detphi1 - detphi2;
    if( detphi3 < meps )  
      return false;
  }
  else if ( nsd==3 ) {
    //    real x = gpt(1), y = gpt(2), z = gpt(3);
    real x = gpt(0), y = gpt(1), z = gpt(2);
    real df; // safety factor to test if inside element
    
    if( isBoundingBoxesComputed() )
      {
        box_counter++;
	//        if( x < bounding_box(element,1) || x > bounding_box(element,2) ||
	//  	    y < bounding_box(element,3) || y > bounding_box(element,4) ||
	//  	    z < bounding_box(element,5) || z > bounding_box(element,6) )
        if( x < bounding_box(element,0) || x > bounding_box(element,1) ||
  	    y < bounding_box(element,2) || y > bounding_box(element,3) ||
  	    z < bounding_box(element,4) || z > bounding_box(element,5) )
  	return false;
        df = 1e-10*exp(getGridLevelNo()*0.742);
      }
    else
      {
        box_counter++;
	int n1,n2,n3,n4;
	real x1,x2,x3,x4, pt,min,max,diff;
        df = 2e-10*exp(getGridLevelNo()*0.742);
	//        n1 = loc2glob( element, 1 );
	//        n2 = loc2glob( element, 2 );
	//        n3 = loc2glob( element, 3 );
	//        n4 = loc2glob( element, 4 );

        n1 = loc2glob( element, 0 );
        n2 = loc2glob( element, 1 );
        n3 = loc2glob( element, 2 );
        n4 = loc2glob( element, 3 );
  
	//        for(int d=1;d<=nsd;d++){
        for(int d=0;d<nsd;d++){
  	  x1 = getCoor(n1,d);
  	  x2 = getCoor(n2,d); 
  	  x3 = getCoor(n3,d);
  	  x4 = getCoor(n4,d);
  	  if ( x1<x2 ){ min=x1; max=x2; }
  	  else   {      min=x2; max=x1; }
  	  if (x3<min) min=x3;
  	  else if (x3>max) max=x3;
  	  if (x4<min) min=x4;
  	  else if (x4>max) max=x4;
  	  diff = (max-min)*df;
  	  pt = gpt(d);
  	  if( pt < min - diff || pt > max + diff ) 
	    return false;
        }
        df*=0.5; // use a stricter rule for inclusion
      } 
    in_elm_counter++;

    FET4.refill(element);
    real u,r,s,t;
    FET4.barycentric(x,y,z,u,r,s,t);
    real mdf = -df;
    if( u < mdf || r < mdf || s < mdf || t < mdf )
      return false; // the point is outside the tetrahedron
  }
// if not returned yet, the point (x) is inside the element
  return true;
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: initSearchStatAllGrids ()
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = (WavesGridFEHier*) getFinestGrid();
  while( ! grid->isThisCoarsestGrid() ){
    grid->initSearchStat();
    grid = grid->getParentGrid();
  }
  grid->initSearchStat();
}

//-----------------------------------------------------------------------------
//void WavesGridFEHier:: printSearchStatAllGrids (Os os)
void WavesGridFEHier:: printSearchStatAllGrids ()
//-----------------------------------------------------------------------------
{
  WavesGridFEHier* grid = (WavesGridFEHier*) getFinestGrid();
  int no_of_levels = grid->getGridLevelNo();
  int total_box=0, total_elm_tests=0;
  while( no_of_levels-- ){
    int boxtests = grid->getNoBoxTests ();
    int eletests = grid->getNoElmTests ();
    total_box += boxtests;
    total_elm_tests += eletests;
    printf(" level = %d  box tests = %d ele tests = %d\n",no_of_levels+1,boxtests,eletests);
    grid = grid->getParentGrid();
  }
  printf(" total on all levels: total box tests = %d total ele tests = %d\n",
	 total_box,total_elm_tests);
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: removeBoundingBoxes()
//-----------------------------------------------------------------------------
{
  bounding_box.newsize(0,0);
  setBoundingBoxesComputed(false);
}


//-----------------------------------------------------------------------------
void WavesGridFEHier:: computeBoundingBoxes()
//-----------------------------------------------------------------------------
{
  bounding_box.newsize(nel,2*nsd);
  setBoundingBoxesComputed(true);
// taking a safety factor twice the one used in the inside test.
  const real df = 2*1e-10*exp(getGridLevelNo()*0.742);
  if ( nsd==1 ) {
    real xmin,xmax;
    
    //    for( int e=1; e<=nel; e++ ) {
    for( int e=0; e<nel; e++ ) {
      //      int n1 = loc2glob( e, 1 );
      //      int n2 = loc2glob( e, getNoNodesInElm(e) );// works for ElmB2n and ElmB3n
      int n1 = loc2glob( e, 0 );
      int n2 = loc2glob( e, getMaxNoNodesInElm()-1 );// works for ElmB2n and ElmB3n
      
      //      real x1 = getCoor(n1,1);
      //      real x2 = getCoor(n2,1);
      real x1 = getCoor(n1,0);
      real x2 = getCoor(n2,0);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      
      real xdiff = (xmax-xmin)*df;
      
      //      bounding_box(e,1)=xmin-xdiff;
      //      bounding_box(e,2)=xmax+xdiff;
      bounding_box(e,0)=xmin-xdiff;
      bounding_box(e,1)=xmax+xdiff;
    }
  }
  else if ( nsd==2 ) {
    
    real xmin,xmax,ymin,ymax;
    
    //    for( int e=1; e<=nel; e++ ) {
    //      int n = loc2glob( e, 1 );
    //      real x1 = getCoor(n,1),  y1 = getCoor(n,2);
    //      n = loc2glob( e, 2 );
    //      real x2 = getCoor(n,1),  y2 = getCoor(n,2);
    //      n = loc2glob( e, 3 );
    //      real x3 = getCoor(n,1),  y3 = getCoor(n,2);

    for( int e=0; e<nel; e++ ) {
      int n = loc2glob( e, 0 );
      real x1 = getCoor(n,0),  y1 = getCoor(n,1);
      n = loc2glob( e, 1 );
      real x2 = getCoor(n,0),  y2 = getCoor(n,1);
      n = loc2glob( e, 2 );
      real x3 = getCoor(n,0),  y3 = getCoor(n,1);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      if ( x3<xmin )         xmin=x3;
      else if  ( x3>xmax )   xmax=x3;
      
      if ( y1<y2 ){ ymin=y1; ymax=y2; }
      else   {      ymin=y2; ymax=y1; }
      if ( y3<ymin )         ymin=y3;
      else if  ( y3>ymax )   ymax=y3;
      
      real xdiff = (xmax-xmin)*df;
      real ydiff = (ymax-ymin)*df;
      
      //      bounding_box(e,1)=xmin-xdiff;
      //      bounding_box(e,2)=xmax+xdiff;
      //      bounding_box(e,3)=ymin-ydiff;
      //      bounding_box(e,4)=ymax+ydiff;

      bounding_box(e,0)=xmin-xdiff;
      bounding_box(e,1)=xmax+xdiff;
      bounding_box(e,2)=ymin-ydiff;
      bounding_box(e,3)=ymax+ydiff;
    }
  }
  else {  // if(nsd==3) {
    
    real xmin,xmax,ymin,ymax,zmin,zmax;
    
    //    for( int e=1; e<=nel; e++ ) {
    //      
    //      int n = loc2glob( e, 1 );
    //      real x1 = coor(n,1), y1 = coord(n,2), z1 = coord(n,3);
    //      n = loc2glob( e, 2 );
    //      real x2 = coord(n,1), y2 = coord(n,2), z2 = coord(n,3);
    //      n = loc2glob( e, 3 );
    //      real x3 = coord(n,1), y3 = coord(n,2), z3 = coord(n,3);
    //      n = loc2glob( e, 4 );
    //      real x4 = coord(n,1), y4 = coord(n,2), z4 = coord(n,3);

    for( int e=0; e<nel; e++ ) {
      
      int n = loc2glob( e, 0 );
      real x1 = coord(n,0), y1 = coord(n,1), z1 = coord(n,2);
      n = loc2glob( e, 1 );
      real x2 = coord(n,0), y2 = coord(n,1), z2 = coord(n,2);
      n = loc2glob( e, 2 );
      real x3 = coord(n,0), y3 = coord(n,1), z3 = coord(n,2);
      n = loc2glob( e, 3 );
      real x4 = coord(n,0), y4 = coord(n,1), z4 = coord(n,2);
      
      if ( x1<x2 ){ xmin=x1; xmax=x2; }
      else   {      xmin=x2; xmax=x1; }
      if ( x3<xmin )         xmin=x3;
      else if  ( x3>xmax )   xmax=x3;
      if ( x4<xmin )         xmin=x4;
      else if  ( x4>xmax )   xmax=x4;
      
      if ( y1<y2 ){ ymin=y1; ymax=y2; }
      else   {      ymin=y2; ymax=y1; }
      if ( y3<ymin )         ymin=y3;
      else if  ( y3>ymax )   ymax=y3;
      if ( y4<ymin )         ymin=y4;
      else if  ( y4>ymax )   ymax=y4;
      
      if ( z1<z2 ){ zmin=z1; zmax=z2; }
      else   {      zmin=z2; zmax=z1; }
      if ( z3<zmin )         zmin=z3;
      else if  ( z3>zmax )   zmax=z3;
      if ( z4<zmin )         zmin=z4;
      else if  ( z4>zmax )   zmax=z4;
      
      real xdiff = (xmax-xmin)*df;
      real ydiff = (ymax-ymin)*df;
      real zdiff = (zmax-zmin)*df;
            
      //      bounding_box(e,1)=xmin-xdiff;
      //      bounding_box(e,2)=xmax+xdiff;
      //      bounding_box(e,3)=ymin-ydiff;
      //      bounding_box(e,4)=ymax+ydiff;
      //      bounding_box(e,5)=zmin-zdiff;
      //      bounding_box(e,6)=zmax+zdiff;

      bounding_box(e,0)=xmin-xdiff;
      bounding_box(e,1)=xmax+xdiff;
      bounding_box(e,2)=ymin-ydiff;
      bounding_box(e,3)=ymax+ydiff;
      bounding_box(e,4)=zmin-zdiff;
      bounding_box(e,5)=zmax+zdiff;
    }
  }
}

//-----------------------------------------------------------------------------
void WavesGridFEHier:: presentGridHier ()
//-----------------------------------------------------------------------------
{
  int total_number_of_grids = getFinestGridLevelNo();
  cout << "\n";
  cout << " Hierarchical finite element mesh (WavesGridFEHier): "<< "\n";
  cout << "\n";
  cout << " This is mesh number " << getGridLevelNo() <<" out of "
     << total_number_of_grids << " meshes." << "\n";
  cout << "\n";
  
  WavesGridFEHier* grid = (WavesGridFEHier*) getCoarsestGrid();
  
  for( int i=0; i<total_number_of_grids; i++ ){
    cout << " Level " << i << " is of type ";
    if(nsd==3){
      if(grid->getElementType(0)==ELMTET1)
	cout<<"Linear Tetra";
      else if(grid->getElementType(0)==ELMTET2)
	cout<<"Quadratic Tetra";
      /*
      else if(grid->getElementType(0)==ELMTET3)
	cout<<"Cubic Tetra";
      else if(grid->getElementType(0)==ELMTET4)
	cout<<"fourth order Tetra";
      else if(grid->getElementType(0)==ELMTET5)
	cout<<"fifth order Tetra";
      else if(grid->getElementType(0)==ELMTET6)
	cout<<"sixth order Tetra"; 
	*/
    }
    else if(nsd==2) {
      if(grid->getElementType(0)==ELMTRI1)
	cout<<"Linear Triangle";
      else if(grid->getElementType(0)==ELMTRI2)
	cout<<"Quadratic Triangle";
    }
    cout << ", has " << grid->getNoNodes() << " nodes and "
	 << grid->getNoElms() << " elements." << "\n";
    grid = grid->getChildGrid();
  }
  cout << "\n"<<flush;  
}


//-----------------------------------------------------------------------------
void WavesGridFEHier:: printGrid (char *file) 
//-----------------------------------------------------------------------------
// Print out all edges in the grid, assumes that the node numbers are
// ordered, works for linear and quadratic triangles in 
{
  FILE *fp;
  fp = fopen(file, "w");

  if(nsd==1){
    int n;
    real minx = getCoor( 0,0 );
    real maxx = getCoor( 0,0 );
    for ( n=0; n<nno; n++ ){
      real p = getCoor( n,0 );
      if(p<minx) minx=p;
      else if(p>maxx) maxx=p;
    }
    real diff = (maxx-minx)/200;
    for ( n=0; n<nno; n++ ){
      real p = getCoor( n,0 );
// make a small vertical line at the node
      fprintf(fp, "%lf 0 0\n %lf %lf 0\n\n", p,p,diff);
    }
    return;
  }
  
  if(nsd==3) fprintf(fp,"$data=curve3d\n");
  
  if ( nsd==2&&maxnne==6 ) {
//  additional_info.initNeighbor( *this );
    //    additional_info.initNeighbor(*this,true,true,false);
    //    WavesNeighborFE& neighbor = additional_info.neighbor;
    neighbor.init(*this,true,true,false);
    
    real p[2],pp[2];
    int nstart,nstop;
    for ( int n=0; n<nno; n++ ) {
      nstart=neighbor.couplingsIrow(n);
      nstop=neighbor.couplingsIrow(n+1);
      p[0] = getCoor( n,0 );
      p[1] = getCoor( n,1 );
      for(int k=nstart;k<nstop; k++){
	int nn = neighbor.couplingsJcol(k);
	if ( nn >= n ) break;
	pp[0]=getCoor(nn,0);
	pp[1]=getCoor(nn,1);
	fprintf(fp, "%lf %lf 0\n %lf %lf 0\n\n", p[0],p[1],pp[0],pp[1]);
      }
    }
  }
  if(nsd==3){
    int d;
    real p1[3],p2[3],p3[3],p4[3];
    for( int e=0; e<nel; e++ ){
      for(d=0;d<nsd;d++){
	p1[d] = getCoor( loc2glob( e, 0 ),d );
	p2[d] = getCoor( loc2glob( e, 1 ),d );
	p3[d] = getCoor( loc2glob( e, 2 ),d );
	p4[d] = getCoor( loc2glob( e, 3 ),d );
      }
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p2[0],p2[1],p2[2],p3[0],p3[1],p3[2]);
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p3[0],p3[1],p3[2],p1[0],p1[1],p1[2]);
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p1[0],p1[1],p1[2],p4[0],p4[1],p4[2]);
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p2[0],p2[1],p2[2],p4[0],p4[1],p4[2]);
      fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p3[0],p3[1],p3[2],p4[0],p4[1],p4[2]);
    }
  }
  fprintf(fp,"$end\n");
  fclose(fp);
}

void WavesGridFEHier::print_gid_mesh_FEM( WavesGridFEHier& gg,char *file) 
  //-----------------------------------------------------------------------------
{
  FILE *fp;
  int i, j;
  int el;
  int d;
  double tmax;
  double dd;
  fp = fopen(file, "r");
  printf("Writing data to: %s\n", file);

  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();
  int nsd = gg.getNoSpaceDim();
  
  
  int nnode = nno;
  int elnum = nel;
  fp = fopen(file, "w");
 
  

  if(nsd==3){
 
 
  fprintf(fp,"MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
  fprintf(fp,"Coordinates\n");

   for(i = 0; i < nnode; i++)
      fprintf(fp, "%i %f %f %f\n ", i + 1,
	     gg.getCoor(i,0),gg.getCoor(i,1),gg.getCoor(i,2));
  }
  else if(nsd==2)
    {


  fprintf(fp,"MESH   dimension  %i  ElemType Triangle  Nnode  3\n", nsd);
  fprintf(fp,"Coordinates\n");


    for(i = 0; i < nnode; i++)
      fprintf(fp, "%i %f %f %f\n ", i + 1,
	      gg.getCoor(i,0),gg.getCoor(i,1),0.0);
  }
 else if(nsd==1)
    {

  fprintf(fp,"MESH   dimension  %i  ElemType Line  Nnode  2\n", nsd);
  fprintf(fp,"Coordinates\n");


    for(i = 0; i < nnode; i++)
      fprintf(fp, "%i %f %f  %f\n ", i + 1,
	      gg.getCoor(i,0),0.0,0.0);
  }
    
    
 fprintf(fp,"end coordinates\n");
 fprintf(fp,"\n");
 fprintf(fp,"Elements\n");

  for(el = 0; el < elnum; el++) {
    switch ( gg.getElementType(el) ) {
      case ELMTET1:
	fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, 
	       gg.loc2glob(el,1)+1,
	       gg.loc2glob(el,2)+1, 
	       gg.loc2glob(el,3)+1,
	       gg.loc2glob(el,0)+1,
               gg.getMaterialType(el)+1);
       break;
     case ELMPYR1:
       fprintf(fp, "%i %i   %i %i %i %i %i\n", el + 1,
	       gg.loc2glob(el,1) + 1,
	       gg.loc2glob(el,2) + 1, 
	       gg.loc2glob(el,3) + 1,
	       gg.loc2glob(el,4) + 1,
	       gg.loc2glob(el,0) + 1,
	       gg.getMaterialType(el)+1);
       break;
     case ELMTRI1:  
       fprintf(fp, "%i %i  %i %i %i\n", el + 1, 
	       gg.loc2glob(el,0) + 1,
	       gg.loc2glob(el,1) + 1, 
	       gg.loc2glob(el,2) + 1,
	       gg.getMaterialType(el)+1);
       break;
     case ELMQUAD1: 
       fprintf(fp, "%i %i   %i %i %i %i\n", el + 1,
	       gg.loc2glob(el,0) + 1,
	       gg.loc2glob(el,1) + 1,
	       gg.loc2glob(el,2) + 1, 
	       gg.loc2glob(el,3) + 1,
	       gg.getMaterialType(el)+1);
       break;
    case ELMLINE1:  
       fprintf(fp, "%i %i   %i %i\n", el + 1, 
	       gg.loc2glob(el,0) + 1,
	       gg.loc2glob(el,1) + 1, 
	       gg.getMaterialType(el)+1);
       break;
     case NO_ELEMENT:
       printf("error print wrong type of element\n");
    }	
  } 
 fprintf(fp,"end elements\n");
  fclose(fp);
}



//-----------------------------------------------------------------------------
void WavesGridFEHier:: printGrid( char *file, const int sbdtype)
//-----------------------------------------------------------------------------
// Print out all edges in the grid, assumes that the node numbers are
// ordered, works for linear and quadratic triangles in 
{
  if( ! (nsd==2||nsd==3) ){
    printf("WavesGridFEHier:: printGrid not implemented for nsd=%d",nsd);
    return;
  }
  FILE *fp;
  fp = fopen(file,"w");

  if(nsd==3) fprintf(fp,"$data=curve3d\n");

  int d;
  real p1[3],p2[3],p3[3],p4[3];
  for( int e=0; e<nel; e++ ){
    if( getMaterialType (e) == sbdtype ) {
      if(nsd==3){
	for(d=0;d<nsd;d++){
	  p1[d] = getCoor( loc2glob( e, 0 ),d );
	  p2[d] = getCoor( loc2glob( e, 1 ),d );
	  p3[d] = getCoor( loc2glob( e, 2 ),d );
	  p4[d] = getCoor( loc2glob( e, 3 ),d );
	}
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p2[0],p2[1],p2[2],p3[0],p3[1],p3[2]);
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p3[0],p3[1],p3[2],p1[0],p1[1],p1[2]);
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p1[0],p1[1],p1[2],p4[0],p4[1],p4[2]);
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p2[0],p2[1],p2[2],p4[0],p4[1],p4[2]);
	fprintf(fp, "%lf %lf %lf\n%lf %lf %lf\n\n", p3[0],p3[1],p3[2],p4[0],p4[1],p4[2]);
      }
      else if(nsd==2){
	for(d=0;d<nsd;d++){
	  p1[d] = getCoor( loc2glob( e, 0 ),d );
	  p2[d] = getCoor( loc2glob( e, 1 ),d );
	  p3[d] = getCoor( loc2glob( e, 2 ),d );
	}
	fprintf(fp, "%lf %lf 0\n%lf %lf 0\n\n", p1[0],p1[1],p2[0],p2[1]);
	fprintf(fp, "%lf %lf 0\n%lf %lf 0\n\n", p2[0],p2[1],p3[0],p3[1]);
	fprintf(fp, "%lf %lf 0\n%lf %lf 0\n\n", p3[0],p3[1],p1[0],p1[1]);
      }
    }
  }
  fprintf(fp,"$end\n");
  fclose(fp);
}



//-----------------------------------------------------------------------------
void WavesGridFEHier:: printBoundary(char* file) 
//-----------------------------------------------------------------------------
// Print out the edges in the grid, which lies on the outer boundary
// or between different subdomaintypes
{
  printf("Enter printboundary, write to file %s\n",file);

  if( ! (nsd==2||nsd==3) ){
    printf("WavesGridFEHier:: printBoundary,  not implemented for nsd=%d",nsd);
    return;
  }

  FILE *r;
  r = fopen(file, "w");

  printf("Writing mtv data to: %s\n", file);


  if( !isNeighborsComputed() )// Is the element-element info computed?
    initNeighbor();
  //  printf("After initNeighbor\n");
  //  cout << "\n"<<flush;  
  if(nsd==2) {
    for(int e=0;e<nel;e++){
      for( int i1=0; i1<3; i1++ ){
	int sbde = getMaterialType (e);
	int nei = getElmNeighbor(e,i1);
	if( nei == -1 || getMaterialType (nei) != sbde ) {
	  if(nei>=0) 
	    //	    fprintf(r,"%%lc=%d\n",max(sbde,getMaterialType (nei)));
	    fprintf(r,"%%lc=%d\n",sbde>getMaterialType (nei)?sbde:getMaterialType (nei));
	  else 
	    fprintf(r,"%%lc=%d\n",sbde);
	  int i2 = (i1==2) ? 0 : i1+1;
	  int n1 = loc2glob( e, i1 );
	  int n2 = loc2glob( e, i2 );
	  fprintf(r,"%12.5e %12.5e 0\n",getCoor(n1,0),getCoor(n1,1));
	  fprintf(r,"%12.5e %12.5e 0\n\n",getCoor(n2,0),getCoor(n2,1));
	}
      }
    } 
  }
  
  if(nsd==3) {
    //    fprintf(r,"$data=curve3d\n");
    fprintf(r,"$data=curve3d\n");

    for(int e=0;e<nel;e++){
      for( int i1=0; i1<4; i1++ ){
	int sbde = getMaterialType (e);
	int nei = getElmNeighbor(e,i1);
	if( nei == -1 || getMaterialType (nei) != sbde ) {
	  if(nei>=0) 
	    //	    fprintf(r,"%%lc=%i\n",max(sbde,getMaterialType (nei)));
	    fprintf(r,"%%lc=%i\n",sbde>getMaterialType (nei)?sbde:getMaterialType (nei));
	  else 
	    fprintf(r,"%%lc=%i\n",sbde);
	  int n1,n2,n3;
	  getGlobNodesOnFace(n1, n2, n3, i1, e );
	  fprintf(r,"%lf %lf %lf\n",
		  getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
	  fprintf(r,"%lf %lf %lf\n",
		  getCoor(n2,0),getCoor(n2,1),getCoor(n2,2));
	  fprintf(r,"%lf %lf %lf\n",
		  //	  fprintf(r,"%12.5e %12.5e %12.5e\n",
		  getCoor(n3,0),getCoor(n3,1),getCoor(n3,2));
	  fprintf(r,"%lf %lf %lf\n\n",
		  //	  fprintf(r,"%12.5e %12.5e %12.5e\n\n",
		  getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
	}
      }
    } 

    int i,j;
    //    fprintf(r,"\n%%lw=%i\n\n",5);
    /*
    for(i=0;i<nocurves;i++){
      int nonodes = curves[i].size();
      for(j=1;j<nonodes;j++){
	int n1 = curves[i](j-1);
	int n2 = curves[i](j);
	fprintf(r,"%%lc=%i\n%lf %lf %lf\n",i+1,
		getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
	fprintf(r,"%lf %lf %lf\n\n",
		getCoor(n2,0),getCoor(n2,1),getCoor(n2,2));
      }
    }
    */
    /*
    int nocurves = getSurfaceTriangulation().getNoCurves();
    for(i=0;i<nocurves;i++){
      int noedges = getSurfaceTriangulation().getNoCurveEdges(i);
      //      int noedges = curves_elm_nodes[i].size(0);
      for(j=0;j<noedges;j++){
	//	int n1 = curves_elm_nodes[i](j,1);
	//	int n2 = curves_elm_nodes[i](j,2);
	int n1 = getSurfaceTriangulation().getCurveNode(i,j,0);
	int n2 = getSurfaceTriangulation().getCurveNode(i,j,1);
	fprintf(r,"%%lc=%i\n%lf %lf %lf\n",i+1,
		getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
	fprintf(r,"%lf %lf %lf\n\n",
		getCoor(n2,0),getCoor(n2,1),getCoor(n2,2));
      }
    }
    int nosurfaces = getSurfaceTriangulation().getNoSurfaces();
    for(i=0;i<nosurfaces;i++){
      //      int nooffaces = surfaces_elm_nodes[i].size(0);
      int nooffaces = getSurfaceTriangulation().getNoSurfaceFaces(i);
      for( j=0; j<nooffaces; j++ ){
	fprintf(r,"%%lc=%i\n",i+1);
	//	int n1 = surfaces_elm_nodes[i](j,1);
	//	int n2 = surfaces_elm_nodes[i](j,2);
	//	int n3 = surfaces_elm_nodes[i](j,3);
	int n1 = getSurfaceTriangulation().getSurfaceNode(i,j,0);
	int n2 = getSurfaceTriangulation().getSurfaceNode(i,j,1);
	int n3 = getSurfaceTriangulation().getSurfaceNode(i,j,2);
	  //	  getGlobNodesOnFace(n1, n2, n3, i1, e );
	fprintf(r,"%lf %lf %lf\n",
		getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
	fprintf(r,"%lf %lf %lf\n",
		getCoor(n2,0),getCoor(n2,1),getCoor(n2,2));
	fprintf(r,"%lf %lf %lf\n",
		//	  fprintf(r,"%12.5e %12.5e %12.5e\n",
		getCoor(n3,0),getCoor(n3,1),getCoor(n3,2));
	fprintf(r,"%lf %lf %lf\n\n",
		//	  fprintf(r,"%12.5e %12.5e %12.5e\n\n",
		getCoor(n1,0),getCoor(n1,1),getCoor(n1,2));
      }
    }
    */
  }

  fprintf(r,"$end\n");
  fclose(r);
  return;

  /*  
  int i,e,n;
//  additional_info.initNeighbor(*this);
  additional_info.initNeighbor(*this,true,true,false);
  WavesNeighborFE& neighbor = additional_info.neighbor;
  if( neighbor.couplingsSize () == 0 || neighbor.nodeSize() == 0 )
    neighbor.init(*this,true,true,false);
  SparseDS ds;
  int coupsize = neighbor.couplingsSize ();
  ds.newsize( nno, coupsize );
  int nno1 = nno+1;
  for ( n=1; n<=nno1; n++ )
    ds.irow(n) = neighbor.couplingsIrow(n);
  for ( n=1; n<=coupsize; n++ )
    ds.jcol(n) = neighbor.couplingsJcol(n);
  
  MatSparse(real) plotted(ds);
  
  plotted.fill(1);
  
  os<<"$data=curve3d"<<"\n";
  int n1,n2,n3;
  for( e=1; e<=nel; e++ ){
    int sbde = getMaterialType (e);
    for(i=1;i<=4;i++){
      int nei = getElmNeighbor(e,i);
      if( nei == -1 || getMaterialType (nei) != sbde ) {
//      if( nei != 0 && getMaterialType (nei) != sbde ) {
	getGlobNodesOnFace( n1, n2, n3, i, e );
	if( plotted(n1,n2) ) {
	  if(nei>=0) os<<"%lc="<<max(sbde,getMaterialType (nei))<<"\n";
	  else os<<"%lc="<<sbde<<"\n";
	  os<<getCoor(n1,1)<<" "<<getCoor(n1,2)<<" "<<getCoor(n1,3)<<"\n"
	    <<getCoor(n2,1)<<" "<<getCoor(n2,2)<<" "<<getCoor(n2,3)<<"\n"
	    <<"\n";
	  plotted(n1,n2) = plotted(n2,n1) = 0;
	}
	if( plotted(n2,n3) ) {
	  if(nei>=0) os<<"%lc="<<max(sbde,getMaterialType (nei))<<"\n";
	  else os<<"%lc="<<sbde<<"\n";
	  os<<getCoor(n2,1)<<" "<<getCoor(n2,2)<<" "<<getCoor(n2,3)<<"\n"
	    <<getCoor(n3,1)<<" "<<getCoor(n3,2)<<" "<<getCoor(n3,3)<<"\n"
	    <<"\n";
	  plotted(n2,n3) = plotted(n3,n2) = 0;
	}
	if( plotted(n3,n1) ) {
	  if(nei>=0) os<<"%lc="<<max(sbde,getMaterialType (nei))<<"\n";
	  else os<<"%lc="<<sbde<<"\n";
	  os<<getCoor(n1,1)<<" "<<getCoor(n1,2)<<" "<<getCoor(n1,3)<<"\n"
	    <<getCoor(n3,1)<<" "<<getCoor(n3,2)<<" "<<getCoor(n3,3)<<"\n"
	    <<"\n";
	  plotted(n1,n3) = plotted(n3,n1) = 0;
	}
      }
    }
  }
  os<<"$end"<<"\n";
  */
}
