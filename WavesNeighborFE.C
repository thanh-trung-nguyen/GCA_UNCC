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
  WavesNeighborFE - neighbor information in finite element grids
*/
 

#include "include/wavesNeighborFE.h"
#include "include/wavesGridB.h"
//#include <VecSort_int.h>


int sort00( MV_Vector<int>& values, int sort_range ){
  //  cout<<"begin sort sort_range="<<sort_range<<flush;
  int temp;
  int l=sort_range/2+1;
  int ir=sort_range;
 first_jump:
  if(l > 1) {
    l=l-1;
    temp = values[l-1];}
  else{
    temp = values[ir-1];
    values[ir-1]=values[0];
    ir=ir-1;
    if(ir==1){
      values[0]=temp;
      return temp;}}
  int i=l;
  int j=l+l;
 second_jump:
  if(j<=ir){
    if(j<ir)
      if(values[j-1]<values[j]) j++;
    if(temp<values[j-1]){
      values[i-1]=values[j-1];
      i=j;
      j=j+j;}
    else
      j=ir+1;
    goto second_jump;
  }
  values[i-1]=temp;
  goto first_jump;
}
/*
//-----------------------------------------------------------------------------
void mixedheapsort (MV_Vector<int>& a, int n)
//-----------------------------------------------------------------------------
{
  // -----------------------------------------------------
  // Sort the first n items of the Type vector.
  // three different sort methods are used, where each
  // is most efficient in separate intervals. The values 
  // separating the intervals have been optimized for 
  // sorting integers. But the optimal values may depend 
  // on processor type and use of compiler optimization flags.
  // -----------------------------------------------------
  if(a.size()<n)
    //    fatalerrorFP("VecSort(Type)::  mixedheapsort",
    //	   "\nattempt to sort vector of length %d with arg. %d",length,n);
    cout<<"mixedheapsort \nattempt to sort vector of";
  int i,j;
  int v;

  if(n<13){  // ************* straight insertion **************
    for (j=1;j<n;j++) {
      v=a(j); i=j;
      while ( a(i-1) > v ) { a(i)=a(i-1); i--; if(i==1) break; }
      a(i)=v;
    }
    return;//not necessary
  }
  else if(n<64) { // ************* Shell's method *************
    int inc;
    inc=13;
    do{
      inc *= 3;
      inc++;
    } while ( inc <= n );
    do{
      inc/=3;
      for(j=inc;j<n;j++){
	v=a(j); i=j;
	while(a(i-inc)>v){ a(i)=a(i-inc); i-=inc; if(i<=inc) break; }
	a(i)=v;
      }
    } while ( inc > 4 );
    for (j=1;j<n;j++) {
      v=a(j); i=j;
      while ( a(i-1) > v ) { a(i)=a(i-1); i--; if(i==0) break; }
      a(i)=v;
    }
    return;//not necessary
  }
  else {  // ****************** heap sort *****************
    int l,ir;
    l=(n>>1);
    while(l>0){
      v=a(l);
      i=l--;
      j=i<<1;
      while(j<n){
	if(a(j)<a(j+1)) j++;
	if(v<a(j)){
	  a(i)=a(j);
	  i=j; j <<= 1;
	}
	else  break;
      }
      if(j==n && v<a(j)) {
	a(i)=a(j);
	a(j)=v;
      }else
	a(i)=v;
    }
    v=a(n);
    a(n)=a(1);
    ir=n-1;
    while(ir>1){
      i=1;j=2;
      while(j<ir){
	if(a(j)<a(j+1)) j++;
	if(v<a(j)){
	  a(i)=a(j);
	  i=j; j <<= 1;
	}
	else  break;
      }
      if(j==ir && v<a(j)) {
	a(i)=a(j);
	a(j)=v;
      }else
	a(i)=v;
      v=a(ir);
      a(ir--)=a(1);
    }
    a(1)=v;
  }
}
*/
//------------------------------------------------------------------------
WavesNeighborFE:: WavesNeighborFE	()
//------------------------------------------------------------------------
{
  n2e.redim(0,0);
  n2n.redim(0,0);
  e2e.redim(0,0);
  special_e2e = false;
}

//------------------------------------------------------------------------
void WavesNeighborFE::remove	( bool remove_n2e,
			  bool remove_n2n,
			  bool remove_e2e )
//------------------------------------------------------------------------
{
  if ( remove_n2e )	n2e.redim(0,0);
  if ( remove_n2n )	n2n.redim(0,0);
  if ( remove_e2e )	e2e.redim(0,0);
  if ( remove_e2e ) 
    special_e2e = false;
}

//------------------------------------------------------------------------
void WavesNeighborFE::init (WavesGridB& grid )
//------------------------------------------------------------------------
{
  init( grid, true, true, true );
}

//------------------------------------------------------------------------
void WavesNeighborFE:: init (WavesGridB& grid,
			bool init_n2e,
		        bool init_n2n,
		        bool init_e2e,
		        bool special_e2e_)
//------------------------------------------------------------------------
{
  this->special_e2e = special_e2e_;
  if( ! init_e2e ) 
    this->special_e2e = false;

  // compute on the info which is wanted as given by
  // init_n2e, init_n2n, init_e2e

  if( ! ( init_n2e || init_n2n || init_e2e ) )
    return;  // do not init anything

// create adresses to avoid reference through handle

  WavesSparseDS& n2ea = (WavesSparseDS&) n2e;
  WavesSparseDS& n2na = (WavesSparseDS&) n2n;
  WavesSparseDS& e2ea = (WavesSparseDS&) e2e;

  // initialize neighbor information:
  int n,e,nne,ne,i,j,k,count,prev,jcol_pos;
  int maxcountn2e,nestart,nestop;
  int totsize_n2e = 0;
  int totsize_n2n = 0;
  int totsize_e2e = 0;
  int nno = grid.getNoNodes();
  int nel = grid.getNoElms();
  bool allElmsOfSameType = grid.oneElementTypeInGrid();
  int maxnne = grid.getMaxNoNodesInElm();
  //  allElmsOfSameType = true; //!!!!!!!!!!!!!!!
  // node-element neighbor info.:

  // Key observation: An element can be added to a given node at most once
  // hence no sort is needed for the node-to-element (n2e) case.

  MV_Vector<int> countn2e( nno,0 );

  if(allElmsOfSameType) {
    //cout<<"In WavesNeighborFE using allElmsOfSameType \n";
    totsize_n2e = nel*maxnne; // each element is added to exactly nne nodes
    //cout<<" nel "<<nel<<" maxnne= "<<maxnne<<endl;

    for (e=0; e<nel; e++)  // count the number of elements a node is in
      for (j=0; j<maxnne; j++) {
	countn2e( grid.loc2glob(e,j) )++;
        // cout <<grid.loc2glob(e,j)<<endl;
	}
  } else { // grid is not homogeneous
    //cout<<"In WavesNeighborFE using not allElmsOfSameType \n";
    totsize_n2e = 0;
    for (e=0; e<nel; e++) {
      totsize_n2e += nne = grid.getNoNodesInElm ( e );
      for (j=0; j<nne; j++)
	countn2e( grid.loc2glob(e,j) )++;
    }
  }

  maxcountn2e = 0;

  //cout<<"before  n2ea.redim,  line 243 "<<endl;

  n2ea.redim( nno, nel, totsize_n2e );
  prev = n2ea.irow(0) = 0;
  for (n=0; n<nno; ) {
    k = countn2e( n );
    countn2e( n ) = 0; // set the counters to 0
    prev = n2ea.irow(++n) = k + prev;
    if( k > maxcountn2e ) maxcountn2e = k;
  }
  if( allElmsOfSameType ) {
    for ( e=0; e<nel; e++ )
      for (j=0; j<maxnne; j++) {
	ne = grid.loc2glob(e,j);
	n2ea.jcol ( n2ea.irow( ne ) + countn2e( ne )++ ) = e;
      }
  } else { //the number of nodes in element vary
    for ( e=0; e<nel; e++ ){
      nne = grid.getNoNodesInElm ( e );
      for (j=0; j<nne; j++) {
	ne = grid.loc2glob(e,j);
	n2ea.jcol ( n2ea.irow( ne ) + countn2e( ne )++ ) = e;
      }
    }
  }

  // node-node neighbor info.:
  //cout<<"before init_n2n"<<endl;
  if( init_n2n ){
    int n2nsort_size =  1 + ( maxnne - 1 ) * maxcountn2e; /*beilin*/
    MV_Vector<int> n2nsort( 1 + ( maxnne - 1 ) * maxcountn2e );
// compute the maximal size of the n2njcol array
    int size_of_n2njcol = 0;
    if( allElmsOfSameType )
      size_of_n2njcol = nno + ( maxnne - 1 ) * totsize_n2e;
    else {
      for ( j=0; j<totsize_n2e; j++ )
	size_of_n2njcol += grid.getNoNodesInElm( n2ea.jcol(j) );
      size_of_n2njcol += nno - totsize_n2e;
    }
    MV_Vector<int> n2njcol( size_of_n2njcol );
    //cout<<" horoscho1"<<endl;

    n2na.redimIrow(nno);

    if( allElmsOfSameType ) {
      jcol_pos = -1;
      for ( n=0; n<nno; n++ )  {
	nestart = n2ea.irow(n);
	nestop = n2ea.irow(n+1);
	count = 0;
	//	count = -1;
	for ( j=nestart; j<nestop; j++ ) {
	  e = n2ea.jcol(j);
	  //	  for (i=1; i<=maxnne; i++) {
	  for (i=0; i<maxnne; i++) {
	    k = grid.loc2glob(e,i);
	    if(k != n) // add the node itself later
	      //	      n2nsort( ++count ) = k;
	      n2nsort( count++ ) = k;
	  }
	}   
 
	if (count>= n2nsort_size){ /*beyl*/
	  cerr<<"** n2nsort:  count="<<count<<" >= "<<n2nsort_size<<endl;
	}

	//	n2nsort( ++count ) = n;  // the node is added
	n2nsort( count++ ) = n;  // the node is added
	//	n2nsort.mixedheapsort( count );
	sort00( n2nsort, count );

//	n2nsort.quicksort( count );
	//	prev = n2njcol( ++jcol_pos ) = n2nsort( 1 );
	prev = n2njcol( ++jcol_pos ) = n2nsort( 0 );
	n2na.irow( n ) = jcol_pos; //start of row n
	//	for( j=2; j<=count; j++ ) // remove copies
	for( j=1; j<count; j++ ) // remove copies
	  if( n2nsort( j ) != prev )
	    prev = n2njcol( ++jcol_pos ) = n2nsort( j );

//cout<<" horoscho2"<<endl;
      }
//cout<<" horoscho3"<<endl;
    }
    else {
      jcol_pos = -1;
      for ( n=0; n<nno; n++ )  {
	nestart = n2ea.irow(n);
	nestop = n2ea.irow(n+1); 
	count = 0; 
	for ( j=nestart; j<nestop; j++ ) { 
	  e = n2ea.jcol(j); 
	  nne = grid.getNoNodesInElm(e); 
	  for (i=0; i<nne; i++) {
	    k = grid.loc2glob(e,i); 
	    if(k != n) // add the node itself later 
	      n2nsort( count++ ) = k; 
	  } 
	} 
	//cout<<"after loop, line 332"<<endl;

	//assert (count< n2nsort_size);
	n2nsort( count++ ) = n;  // the node is added 
	//assert (count< n2nsort_size);
	sort00( n2nsort, count ); 
	prev = n2njcol( ++jcol_pos ) = n2nsort( 0 );
	n2na.irow( n ) = jcol_pos; //start of row n
	for( j=1; j<count; j++ ) // remove copies
	  if( n2nsort( j ) != prev )
	    prev = n2njcol( ++jcol_pos ) = n2nsort( j );
      }
    }
    n2na.irow( nno ) = jcol_pos + 1;
    if( !init_n2e && !init_e2e )
      n2e.redim(0,0);
    totsize_n2n = jcol_pos+1;
    n2na.fillJcol( n2njcol, totsize_n2n ); // fill WavesSparseDS with correct size
    n2njcol.newsize(0);  // release the temporary array
  }
  else
    n2n.redim(0,0); // do not want any n2n info

    // element-element neighbor info.:

  if( init_e2e ){
    int e2esize = (1>maxnne * (maxcountn2e - 1))?1:maxnne * (maxcountn2e - 1);
    MV_Vector<int> e2esort( e2esize );
  // compute the maximal size of the e2ejcol array
    int size_of_e2ejcol = 0;
    if( allElmsOfSameType ) {
      for (e=0; e<nel; e++)  {
	for (j=0; j<maxnne; j++)
	  size_of_e2ejcol += countn2e( grid.loc2glob(e,j) );
      }
    }
    else {
      for (e=0; e<nel; e++)  {
	nne = grid.getNoNodesInElm(e);
	for (j=0; j<nne; j++)
	  size_of_e2ejcol += countn2e( grid.loc2glob(e,j) );
      }
    }
    size_of_e2ejcol -= totsize_n2e;
    int e2ejcolsize = (1>size_of_e2ejcol)?1:size_of_e2ejcol;
    MV_Vector<int> e2ejcol( e2ejcolsize );
    e2ea.redimIrow( nel );

    jcol_pos = -1;
    if( allElmsOfSameType ) {
      for (e=0; e<nel; e++)  {
	count = -1;
	for (j=0; j<maxnne; j++) {
	  ne = grid.loc2glob(e,j);
	  nestart = n2ea.irow(ne);
	  nestop = n2ea.irow(ne+1);
	  for( i=nestart; i<nestop; i++ ) {
	    k = n2ea.jcol(i);
	    if( k != e )
	      e2esort( ++count ) = k;
	  }
	}
	sort00(e2esort,count+1);
	prev = e2ejcol( ++jcol_pos ) = e2esort(0);// ok even if count == 0
	e2ea.irow(e) = jcol_pos; // start position of element e
	for( j=1; j<=count; j++ ) // remove copies
	  if( e2esort(j) != prev )
	    prev = e2ejcol( ++jcol_pos ) = e2esort(j);
      }
    }
    else {
      for (e=0; e<nel; e++)  {
	nne = grid.getNoNodesInElm(e);
	count = -1;
	for (j=0; j<nne; j++) {
	  ne = grid.loc2glob(e,j);
	  nestart = n2ea.irow(ne);
	  nestop = n2ea.irow(ne+1);
	  for( i=nestart; i<nestop; i++ ) {
	    k = n2ea.jcol(i);
	    if( k != e )
	      e2esort( ++count ) = k;
	  }
	}
	sort00(e2esort,count+1);
	prev = e2ejcol( ++jcol_pos ) = e2esort(0);// ok even if count == 0
	e2ea.irow(e) = jcol_pos; // start position of element e
	for( j=1; j<=count; j++ ) // remove copies
	  if( e2esort(j) != prev )
	    prev = e2ejcol( ++jcol_pos ) = e2esort(j);
      }
    }

    if( !init_n2e ) // remove n2e info if it is not wanted
      n2e.redim(0,0);
    if ( nel == 1 ) // the case of one single element
      jcol_pos = 0;
    e2ea.irow( nel ) = jcol_pos+1;
    totsize_e2e = jcol_pos+1;

    if(special_e2e){ // use only elements with at least "nsd" common nodes      
      int nsd = grid.getNoSpaceDim ();
      if(nsd >= 2){ // for the nsd==1 the e2e is ok
	jcol_pos = 0;
	int nestart,nestop,nnne;
	if( allElmsOfSameType ) {
	  for (e=0; e<nel; e++)  {
	    nestart = e2ea.irow(e);
	    nestop = e2ea.irow(e+1);
	    e2ea.irow(e) = jcol_pos; // reuse the old arrays, order important!
	    for( i=nestart; i<nestop; i++ ) {
	      k = e2ejcol(i); // a neighbor element
	      count = 0; // check number of common nodes between e and k
	      for (j=0; j<maxnne; j++) {
		ne = grid.loc2glob(e,j); // a node in element e
		for (n=0; n<maxnne; n++) 
		  if( ne == grid.loc2glob(k,n) ) // nodes are equal
		    count++;
	      }
	      if( count >= nsd )  // add to e2e structure
		e2ejcol(jcol_pos++) = k;
	    }
	  }
	} else { // heterogenous grid
	  for (e=0; e<nel; e++)  {
	    nne = grid.getNoNodesInElm(e); 
	    nestart = e2ea.irow(e); 
	    nestop = e2ea.irow(e+1); 
	    e2ea.irow(e) = jcol_pos; // reuse the old arrays, order important!
	    for( i=nestart; i<nestop; i++ ) {
	      k = e2ejcol(i); // a neighbor element
	      count = 0; // check number of common nodes between e and k
	      for (j=0; j<nne; j++) {
		ne = grid.loc2glob(e,j); // a node in element e
		nnne = grid.getNoNodesInElm(k); 
		for (n=0; n<nnne; n++) 
		  if( ne == grid.loc2glob(k,n) ) // nodes are equal
		    count++;
	      }
	      if( count >= nsd )  // add to e2e structure
		e2ejcol(jcol_pos++) = k;
	    }
	  }
	}
	e2ea.irow(nel) = jcol_pos; // last position + 1
	totsize_e2e = jcol_pos;
      }
    }
    e2ea.fillJcol(e2ejcol,totsize_e2e);
    e2ejcol.newsize(0);
  }
  else
    e2e.redim(0,0); // in this case e2e info is not needed

  if(n2e.getNoRows() && !n2e.consistent()) {
    printf("n2e not consistent\n");
    exit(1);
  }
  if(n2n.getNoRows() && !n2n.consistent())  {
    printf("n2n not consistent\n");
    exit(1);
  }
  if(e2e.getNoRows() && !e2e.consistent()) { 
    printf("e2e not consistent\n");
    exit(1);
  }
}


//------------------------------------------------------------------------
bool WavesNeighborFE:: ok() const
//------------------------------------------------------------------------
{
  bool b = true;
  /*  HANDLE_OK("Neighbor::ok", n2e, b);
  HANDLE_OK("Neighbor::ok", n2n, b);
  HANDLE_OK("Neighbor::ok", e2e, b);
  if ( !b )
    errorFP("Neighbor::ok","You have probably forgotten to call init.");
    */
  return b;
}


//-----------------------------------------------------------------------
void WavesNeighborFE::print()
//-----------------------------------------------------------------------
{
  /*
  if (!ok()) {
    errorFP("Neighbor::print",
	    "Cannot print not ok object. You have probably not called\n"
	    "Neighbor::init");
    return;
  }
    n2e.printPattern();
  n2n.printPattern();
  e2e.printPattern();
  */


  int nno, nel, n, e, i,l;
  //  VecSimplest(int) v;

  if ( n2e.ok() )
    nno = n2e.getNoRows();
  else if ( n2n.ok() )
    nno = n2n.getNoRows();

  if ( n2e.ok() )
    nel = n2e.getNoColumns();
  else if( e2e.ok() )
    nel = e2e.getNoColumns();

  cout << "Neighbor information of the finite element mesh:\n\n";
  if( n2e.ok() || n2n.ok() )
    cout << "  Number of nodes = " << nno <<"\n";
  if( n2e.ok() || e2e.ok() )
    cout << "  Number of elements = " << nel << "\n\n";

  if( n2e.ok() ){
    cout << "\n  Nodes neighbor information:\n";
    cout << "    The columns contain:\n";
    cout << "      - node number;\n";
    cout << "      - a set of the numbers of the neighbor elements.\n";
    cout << "n2e#\n";

    for (n=0; n<nno; n++)
      {
	cout << "  (" << n <<"):  ";
	for( i=nodeIrow(n); i<nodeIrow(n+1); i++)
	  cout<<" "<<nodeJcol(i)<<flush;
	cout<<endl;
      }
  }

  if( n2n.ok() ){
    cout << "\n  Nodes neighbor information:\n";
    cout << "    The columns contain:\n";
    cout << "      - node number;\n";
    cout << "      - a set of the numbers of the node and its neighbor nodes.\n";
    cout << "n2n#\n"<<flush;

    for (n=0; n<nno; n++)
      {
	cout << "  (" << n <<"):  ";
	//	couplings(n,v);
	//	l = v.size();
	for( i=couplingsIrow(n); i<couplingsIrow(n+1); i++)
	  cout<<" "<<couplingsJcol(i)<<flush;
	cout<<"\n";
      }
  }
  if( e2e.ok() ){
    cout << "\n  Elements neighbor information: \n";
    cout << "    The columns contain:\n";
    cout << "      - element number;\n";
    cout << "      - a set of the numbers of the neighbor elements.\n";
    if( special_e2e )
      cout << "    The neighbor elements contain nsd common nodes.\n";
    cout << "e2e#\n"<<flush;

    for (e=0; e<nel; e++)
      {
	cout << "  (" << e <<"):  ";
	for( i=elementIrow(e); i<elementIrow(e+1); i++)
	  cout<<" "<<elementJcol(i)<<flush;
	cout<<"\n";
      }
  }
  cout << "\n"<<flush;
}

