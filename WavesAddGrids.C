#include <include/wavesAddGrids.h>


bool 
AddGrids:: addGrids( WavesGridB& g1,
		     WavesGridB& g2,
		     WavesGridB& gadd )
{
// Assumption: two grids which are disjoint except for common nodes
// at the boundary
// are used to construct a new grid. The numbering of nodes from grid 1
// are kept. The nodes in grid 2 which are common with grid 1 are given
// the grid 1 nodal numbers. The remaining nodes are given increasing
// node numbers with start where the grid 1 numbering stopped.
// The element numbering of grid 1 will be kept in the added grid,
// and the elements which was in grid 2 gets the following element
// numbers in the same order.

  cout<<"entering addGrids\n"<<flush;

  int nno1 = g1.getNoNodes();
  int nel1 = g1.getNoElms();
  int nno2 = g2.getNoNodes();
  int nel2 = g2.getNoElms();

  int maxnne1 = g1.getMaxNoNodesInElm ();
  int maxnne2 = g2.getMaxNoNodesInElm ();
  //  int nbind1 = g1.getNoBoInds ();
  //  bool onemat1 = g1.oneMaterialType();
  bool onemat1 = false;
  //  int nbind2 = g2.getNoBoInds ();
  //  bool onemat2 = g2.oneMaterialType();
  bool onemat2 = false;

  int nnoadd = nno1; // last current node number in the new grid

  int nsd = g1.getNoSpaceDim();
  int i,n,e,d;
  MV_Vector<int> old2new(nno2);

// identify which nodes in grid 2 which are also nodes in grid 1
// this is an O(N^2) algorithm, might be more efficient
// to change place of grid 1 and grid 2 if grid2 large and grid1 small

  real dist;
  bool exact;
  int commonnodecounter=0; // counter for found common nodes
  MV_Vector<real> globpt(nsd);
  for(n=0;n<nno2; n++){
    for(d=0;d<nsd;d++)
      globpt(d) = g2.getCoor(n,d);
    int node = g1.nearestPoint(globpt,dist,exact);
    if ( exact || dist < tolerance ) { // which tolerance?
      commonnodecounter++;
      old2new(n) = node;
    }
    else
      old2new(n) = nnoadd++;
  }

  cout<<" found "<<commonnodecounter<<" common nodes in the grids\n"<<flush;
  cout.flush();
// compute the values of nne, onemat and nbind for the new grid

  int neladd = nel1+nel2; // new number of elements
  int maxnneadd = maxnne1>maxnne2?maxnne1:maxnne2; // maximum number of elms
  bool hybridGrid = false;
  if ( maxnne1 != maxnne2 ) // not completely safe method
    hybridGrid = true;
  int newnneadd = maxnneadd; // assume homogenous grids
  if(hybridGrid)
    newnneadd = 0;
// will not differ indicators in the two grids from each other
  //  int nbindadd = nbind1>nbind2?nbind1:nbind2; 
  int nbindadd=0;
// only one material type?
  bool onematadd = ( !onemat1 || !onemat2 ) ? false : true;

  gadd.redim(nsd,nnoadd,neladd,maxnneadd,nbindadd);

  if( g1.oneElementTypeInGrid() && 
      g2.oneElementTypeInGrid() &&
      g1.getElementType(0) == g2.getElementType(0))
    gadd.setElementType( g1.getElementType(0) );
  else
    printf("AddGrids:: addGrids, Grids of different type added\n");

// 0) set indicator names

// copy the names from grid 1
  /*
  int b;
  for(b=1;b<=nbind2;b++)
    gadd.putBoIndName(g2.getBoIndName(b),b);
// the indicator names from grid 1 overwrites the names of grid 2
  for(b=1;b<=nbind1;b++)
    gadd.putBoIndName(g1.getBoIndName(b),b);
    */
// 1) copy the first nne1 elements in grid1, 

  if ( ! hybridGrid )
    for(e=0;e<nel1;e++){
      gadd.setMaterialType (e,g1.getMaterialType(e));
      for(i=0;i<maxnne1;i++)
	gadd.putLoc2glob(e,i,g1.loc2glob(e,i));
    }
  else
    for(e=0;e<nel1;e++){
      gadd.setMaterialType (e,g1.getMaterialType(e));
      gadd.setElementType(e,g1.getElementType(e));
      int nonodes = g1.getNoNodesInElm(e);
      //      gadd.setNoNodesInElm(e,nonodes);
      for(i=0;i<nonodes;i++)
	gadd.putLoc2glob(e,i,g1.loc2glob(e,i));
    }


// 2) copy the nno1 nodes from grid1, keep indicators

  for( n=0; n<nno1; n++ ){
    for( d=0; d<nsd; d++ )
      gadd.putCoor(n,d,g1.getCoor(n,d));
    //    for(b=0;b<nbind1; b++)
    //      if(g1.boNode(n,b)){
    //	gadd.setBoInd(n,b);
    //      }
  }

// 3) copy the nno1+1 to neladd nodes from grid2
  for( n=0; n<nno2; n++ ){
    int node = old2new(n); // get the new node number
    for( d=0; d<nsd; d++ )
      gadd.putCoor(node,d,g2.getCoor(n,d));
    //    for(b=0;b<nbind2; b++)
    //      if(g2.boNode(n,b)){
    //	gadd.setBoInd(node,b);
    //      }
  }

// 4) copy all elements from grid2
// observe that no "overlap" between grid 1 and grid 2 is allowed

  if ( ! hybridGrid )
    for(e=0;e<nel2;e++){
      gadd.setMaterialType (nel1+e,g2.getMaterialType(e));
      for(i=0;i<maxnne2;i++)
	gadd.putLoc2glob(nel1+e,i,old2new(g2.loc2glob(e,i)));
    }
  else
    for(e=0;e<nel2;e++){
      gadd.setMaterialType (nel1+e,g2.getMaterialType(e));
      gadd.setElementType(nel1+e,g2.getElementType(e));
      int nonodes = g2.getNoNodesInElm(e);
      //      gadd.setNoNodesInElm(nel1+e,nonodes);
      for(i=0;i<nonodes;i++)
	gadd.putLoc2glob(nel1+e,i,old2new(g2.loc2glob(e,i)));
    }

  //  gadd.setNonUniformMesh(); // addgrid is unstructured.

  cout<<"leaving addGrids\n"<<flush;
  return true;

}

bool 
AddGrids:: extractGrid( WavesGridB& g1,
			WavesGridB& g2,
			MV_Vector<int>& partition)
{
// create a new grid g2 which only contains the elements of g1 
// which are numbered in the partition vector
// the elements in the new grid are numbered after the order of
// the elements come in the partition vector
// the nodes are ordered after ...
  cout<<"entering extractGrid\n"<<flush;

// get properties of the input grid
  int nno1 = g1.getNoNodes();
  int nel1 = g1.getNoElms();
  int maxnne1 = g1.getMaxNoNodesInElm ();
  //  int nbind1 = g1.getNoBoInds ();
  //  bool onemat1 = g1.oneMaterialType();
  bool onemat1 = false;
  int nsd = g1.getNoSpaceDim();
  int i,b,n,e,d;
  MV_Vector<int> old2new(nno1); // node numbering in the old grid
  old2new = -1;

  int nel2 = partition.size(); // number of elements in new grid

  int nodecounter=0;
  for(e=0;e<nel2;e++) {
    int ee = partition(e);
    //    cout<<"e="<<e<<" ee="<<ee<<endl<<flush;
    int nne = g1.getNoNodesInElm(ee);
    for(i=0;i<nne;i++) {
      int node = g1.loc2glob(ee,i);
      if( old2new(node)==-1 ) { // has this node been used before?
	old2new(node) = nodecounter;
        nodecounter++;
      }
    }
  }
    
  cout<<" found "<<nodecounter<<" nodes in the new grid\n"<<flush;

// compute the values of nne, onemat and nbind for new grid

  int nel = nel2; // new number of elements
  int maxnne = maxnne1; // maximum number of elms
  int nno = nodecounter;
  int newnne = maxnne; // assume homogenous grids
  //  int nbind = nbind1;
  int nbind = 0;

  //n  g2.redim(nsd,nno,nel,maxnne,nbind,newnne,onemat1);
  g2.redim(nsd,nno,nel,maxnne,nbind);

// should change so that it also works for heterogenous grids
  if ( g1.oneElementTypeInGrid() )
    g2.setElementType( g1.getElementType(0) );
  else
    printf("extractGrid, should be a homogeneous grid\n");

// 0) set boundary indicator names
  /*
// copy the names from grid 1
  for ( b=1; b<=nbind1; b++ )
    g2.putBoIndName(g1.getBoIndName(b),b);
    */
// 1) copy the elements which are in the partition vector
  
  for ( e=0; e<nel2; e++ ) {
    int ee = partition(e);
    g2.setMaterialType (e,g1.getMaterialType(ee));
    int nne = g1.getNoNodesInElm(ee);
    for ( i=0; i<nne; i++ )
      g2.putLoc2glob(e,i,old2new(g1.loc2glob(ee,i)));
  }

// 2) copy the nodes (and indicators) which are also part of the new grid

  for ( n=0; n<nno1; n++ ) {
    int node = old2new(n);
    if ( node != -1 ) {
      for ( d=0; d<nsd; d++ )
	g2.putCoor(node,d,g1.getCoor(n,d));
      //      for ( b=0; b<nbind1; b++ )
      //	if ( g1.boNode(n,b) )
      //	  g2.setBoInd(node,b);
    }
  }

  //  g2.setNonUniformMesh(); // theextracted grid is unstructured.

  cout<<"leaving extractGrid\n"<<flush;
  return true;

}
/*
//----------------------------------
void 
AddGrids:: symmetrizeAndRotateGrid( WavesGridB& grid, 
				    WavesGridB& outgrid,
				    int rotations)
//----------------------------------
{
// Assumptions: A 2D or 3D grid is in a sector with center
// at origin. One part of the boundary of the grid 
// has its nodes at the y=0 plane.
// Also at the boundary at the sector bust be straight
// The angle size of the sector is M_PI/rotations.
// If this assumption fulfilled then the outgrid will be constructed
// by first mirroring the initial grid in the y=0 line/plane. The resulting
// symmetric grid is rotated and added to itself. The outgrid
// is of the form of an annulus or torus or cylinder shape, depending
// on the shape of the initial grid.

  DBP("enter symmetrizeAndRotateGrid");
  Handle(WavesGridB) mirrorgrid,rotgrid;
  mirrorgrid.rebind( new WavesGridB() );
  const int nsd = grid.getNoSpaceDim();
  Ptv(real) plane(nsd+1);
  plane = 0.0; // initilizing all components
  plane(2) = 1.0;
//in 3D  plane(1) = 0.0;  plane(2) = 1.0;  plane(3) = 0.0; plane(4) = 0.0;
//in 2D  plane(1) = 0.0;  plane(2) = 1.0;  plane(3) = 0.0;
// mirror the grid on the line x*0 + y*1 + z*0 + 0 = 0, i.e. y=0
  createMirrorGrid( grid, mirrorgrid(), plane);
  real rotangle = 2*M_PI/rotations; // angle for rotation of grid
  rotgrid.rebind( new WavesGridB( mirrorgrid() ) ); // copy the grid

  for(int i=1;i<rotations-1;i++) {
    rotateGrid( rotgrid(), rotangle );
    Handle(WavesGridB) addedgrids;
    addedgrids.rebind(new WavesGridB);
    addGrids( mirrorgrid(), rotgrid(), addedgrids());
    mirrorgrid.rebind(addedgrids());
  }
// the last rotation and addition is done so that
// the result will end up in the outgrid
  rotateGrid( rotgrid(), rotangle );
  addGrids( mirrorgrid(), rotgrid(), outgrid );

  DBP("leave symmetrizeAndRotateGrid");
}

//----------------------------------
void 
AddGrids:: rotateGrid( WavesGridB& grid, 
		       real alpha ) 
//----------------------------------
{
// rotate the grid around the z-axis (in 3D) and origin (in 2D) alpha radians,
// Even if the grid is three dimensional will not the z-coordinate
// of the nodes be changed
  int nno = grid.getNoNodes();
  int n;
// rotate the nodes of the grid with alpha and radians and then translate 
  real cosa = cos(alpha);
  real sina = sin(alpha);
  for ( n=1; n<=nno; n++ ) {
    real x = grid.getCoor(n,1);
    real y = grid.getCoor(n,2);
// rotate alpha radians around the origin
    grid.putCoor(x*cosa - y*sina,n,1);
    grid.putCoor(x*sina + y*cosa,n,2);
// a possible third coordinate is unchanged
  }
}


//----------------------------------
bool 
AddGrids :: createMirrorGrid( WavesGridB& grid, 
			      WavesGridB& newgrid, 
			      Ptv(real)& plane_ )
//----------------------------------
{
  int nsd = grid.getNoSpaceDim();
  if(plane_.size() != nsd+1)
    errorFP("createMirrorGrid",
	    "hyper plane has dimension %d, it should be %d",
	    plane_.size(),nsd+1);

  real tolerance_ = 1e-10;

  // make a scaled copy of the (hyper) plane

  int n,d;
  Ptv(real) plane(plane_.size());

  real lsum = 0;

  for(d=1;d<=nsd;d++)
    lsum += sqr(plane_(d));

  lsum = sqrt(lsum);

  if(lsum < tolerance_)
    errorFP("createMirrorGrid","incorrect hyper plane given");

  for(d=1;d<=nsd+1;d++)
    plane(d) = plane_(d)/lsum;

  // find the number of nodes which lie in the mirror plane

  int noPtsOnPlane = 0;
  int nno = grid.getNoNodes();
  int nel = grid.getNoElms();
  int maxnne = grid.getMaxNoNodesInElm ();
  int nne;
  if (grid.allElmsOfSameType ())
    nne = maxnne;
  else
    nne = 0;
  int bind = grid.getNoBoInds ();
  //  bool onemattype = grid.oneMaterialType();
  bool onemattype = false;
  MV_Vector<int> foundNodeNumbers(nno); // will not use all
  foundNodeNumbers.fill(0);

  for(n=1;n<=nno;n++) {
    lsum = plane(nsd+1);
    for(d=1;d<=nsd;d++) // A*x+B*y+C*z+D = 0 defines a hyperplane
      lsum += plane(d)*grid.getCoor(n,d);
    if(lsum > -tolerance_ && lsum < tolerance_)
      foundNodeNumbers(++noPtsOnPlane) = n;
  }

  cout<<"Found "<<noPtsOnPlane<<" nodes in the mirror plane\n";

//  if(noPtsOnPlane<2)
//    warningFP("createMirrorGrid","too few nodes in the plane?");
  
  int nnonew,nelnew;
  nnonew = noPtsOnPlane + 2*(nno-noPtsOnPlane);
  nelnew = 2*nel;
  newgrid.redim( nsd, nnonew, nelnew, nne, bind, maxnne, onemattype);
  if(grid.allElmsOfSameType())
    newgrid.setElmType(grid.getElmType(1));

  newgrid.setNonUniformMesh();

// make a new numbering of the nodes in the old grid
// but do not change it.

  MV_Vector<int> newNodeNumbering(nno);
  newNodeNumbering.fill(0);

// first number the nodes in the mirror plane
  for(n=1;n<=noPtsOnPlane;n++)
    newNodeNumbering(foundNodeNumbers(n)) = n;

// then continue with the remaining node numbers
  int ncounter = noPtsOnPlane;
  for ( n=1; n<=nno; n++ )
    if( ! newNodeNumbering(n) )
      newNodeNumbering(n) = ++ncounter;
 
  for(d=1;d<=bind;d++)
    newgrid.putBoIndName(grid.getBoIndName(d),d);

  // now start to copy the node coordinates
  for(n=1;n<=nno;n++) {
    int newnode = newNodeNumbering(n);
    if( newnode <= noPtsOnPlane ) { // nodes in the mirror plane

      for(d=1;d<=nsd;d++) 
	newgrid.putCoor(grid.getCoor(n,d),newnode,d);
      // copy the boundary indicators

      for(d=1;d<=bind;d++) 
	if(grid.boNode(n,d))
	  newgrid.setBoInd(newnode,d);
    }
    else { // case of nodes not in the mirror plane
      for(d=1;d<=nsd;d++) 
	newgrid.putCoor(grid.getCoor(n,d),newnode,d);
      int mirrornode = newnode+nno-noPtsOnPlane;
      // find the position of the mirrored node
      real aaa = plane(nsd+1);
      for(d=1;d<=nsd;d++)
	aaa += plane(d)*grid.getCoor(n,d);
      aaa*=2;
      for(d=1;d<=nsd;d++)
	newgrid.putCoor(grid.getCoor(n,d)-plane(d)*aaa,mirrornode,d);
      for(d=1;d<=bind;d++) 
	if(grid.boNode(n,d)){
	  newgrid.setBoInd(newnode,d);
	  newgrid.setBoInd(mirrornode,d);
	}
    }
  }

// next copy the elements

  int e;
  int nnodes,newnode;

  if( ! grid.allElmsOfSameType() )
    errorFP("createMirrorGrid","elementtype not implemented yet");

  if( grid.getElmType(1) == "ElmT3n2D" || 
      grid.getElmType(1) == "ElmB2n1D" ||
      grid.getElmType(1) == "ElmB3n1D" ) {
    for ( e=1;e<=nel;e++ ) {
      if(!onemattype) {
	newgrid.setMaterialType(e,grid.getMaterialType(e));
	newgrid.setMaterialType(e+nel,grid.getMaterialType(e));
      }

      nnodes = grid.getNoNodesInElm(e);
      
      for(n=1;n<=nnodes;n++){
	newnode = newNodeNumbering(grid.loc2glob(e,n));
	newgrid.putLoc2glob(newnode,e,n);
	if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
	newgrid.putLoc2glob(newnode,nel + e,nnodes + 1 - n);
      }
    }
  }
  else if( grid.getElmType(1) == "ElmT6n2D" ) {
    for ( e=1;e<=nel;e++ ) {
     if(!onemattype) {
       newgrid.setMaterialType(e,grid.getMaterialType(e));
       newgrid.setMaterialType(e+nel,grid.getMaterialType(e));
     }
     newnode = newNodeNumbering(grid.loc2glob(e,1));
     newgrid.putLoc2glob(newnode,e,1);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,2);
     newnode = newNodeNumbering(grid.loc2glob(e,2));
     newgrid.putLoc2glob(newnode,e,2);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,1);
     newnode = newNodeNumbering(grid.loc2glob(e,3));
     newgrid.putLoc2glob(newnode,e,3);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,3);
     newnode = newNodeNumbering(grid.loc2glob(e,4));
     newgrid.putLoc2glob(newnode,e,4);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,4);
     newnode = newNodeNumbering(grid.loc2glob(e,5));
     newgrid.putLoc2glob(newnode,e,5);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,6);
     newnode = newNodeNumbering(grid.loc2glob(e,6));
     newgrid.putLoc2glob(newnode,e,6);
     if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
     newgrid.putLoc2glob(newnode,nel + e,5);
    }
  }
  else if( grid.getElmType(1) == "ElmT4n3D" ||
	   grid.getElmType(1) == "ElmT10n3D" ) {
    for ( e=1;e<=nel;e++ ) {
      if(!onemattype) {
	newgrid.setMaterialType(e,grid.getMaterialType(e));
	newgrid.setMaterialType(e+nel,grid.getMaterialType(e));
      }
      newnode = newNodeNumbering(grid.loc2glob(e,1));
      newgrid.putLoc2glob(newnode,e,1);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,2);

      newnode = newNodeNumbering(grid.loc2glob(e,2));
      newgrid.putLoc2glob(newnode,e,2);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,1);
      newnode = newNodeNumbering(grid.loc2glob(e,3));
      newgrid.putLoc2glob(newnode,e,3);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,3);
      newnode = newNodeNumbering(grid.loc2glob(e,4));
      newgrid.putLoc2glob(newnode,e,4);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,4);
    }
    if( grid.getElmType(1) == "ElmT10n3D" ) {
    for ( e=1;e<=nel;e++ ) {
      newnode = newNodeNumbering(grid.loc2glob(e,5));
      newgrid.putLoc2glob(newnode,e,5);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,5);
      newnode = newNodeNumbering(grid.loc2glob(e,6));
      newgrid.putLoc2glob(newnode,e,6);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,7);
      newnode = newNodeNumbering(grid.loc2glob(e,7));
      newgrid.putLoc2glob(newnode,e,7);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,6);
      newnode = newNodeNumbering(grid.loc2glob(e,8));
      newgrid.putLoc2glob(newnode,e,8);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,9);
      newnode = newNodeNumbering(grid.loc2glob(e,9));
      newgrid.putLoc2glob(newnode,e,9);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,8);
      newnode = newNodeNumbering(grid.loc2glob(e,10));
      newgrid.putLoc2glob(newnode,e,10);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,10);
    }
    }
  }
  else if( grid.getElmType(1) == "ElmB4n2D" || 
	   grid.getElmType(1) == "ElmB8n3D") {
    for ( e=1;e<=nel;e++ ) {
      if(!onemattype) {
	newgrid.setMaterialType(e,grid.getMaterialType(e));
	newgrid.setMaterialType(e+nel,grid.getMaterialType(e));
      }
      newnode = newNodeNumbering(grid.loc2glob(e,1));
      newgrid.putLoc2glob(newnode,e,1);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,2);
      newnode = newNodeNumbering(grid.loc2glob(e,2));
      newgrid.putLoc2glob(newnode,e,2);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,1);
      newnode = newNodeNumbering(grid.loc2glob(e,3));
      newgrid.putLoc2glob(newnode,e,3);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,4);
      newnode = newNodeNumbering(grid.loc2glob(e,4));
      newgrid.putLoc2glob(newnode,e,4);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,3);
    }
    if( grid.getElmType(1) == "ElmB8n3D" ) { 
      for ( e=1;e<=nel;e++ ) {
      newnode = newNodeNumbering(grid.loc2glob(e,5));
      newgrid.putLoc2glob(newnode,e,5);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,6);
      newnode = newNodeNumbering(grid.loc2glob(e,6));
      newgrid.putLoc2glob(newnode,e,6);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,5);
      newnode = newNodeNumbering(grid.loc2glob(e,7));
      newgrid.putLoc2glob(newnode,e,7);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,8);
      newnode = newNodeNumbering(grid.loc2glob(e,8));
      newgrid.putLoc2glob(newnode,e,8);
      if( newnode > noPtsOnPlane) newnode += nno-noPtsOnPlane;
      newgrid.putLoc2glob(newnode,nel + e,7);
      }
    }
  }
  else
    errorFP("createMirrorGrid","elementtype not implemented yet");
  return true;
}
*/
