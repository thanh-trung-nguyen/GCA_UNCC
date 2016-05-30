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
  Find boundary nodes  at the boundary of a finite element grid
*/



#include "include/wavesfindBoundaryNodes.h"
//-----------------------------------------------------------------------------
void getElmsAtFaces(const WavesNeighborFE& n2e,
		    const int n1, const int n2, const int n3,
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
int getElmAtFace(const WavesNeighborFE& n2e,
		 const int n1, const int n2, const int n3,
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
int getFaceNo(const WavesGridB& grid, 
	      const int n1, const int n2, const int n3, const int e, bool safe)
//-----------------------------------------------------------------------------
// Get the face number of element e which consists of the nodes n1, n2 and n3
// if safe is true the function will return 0 if no match is found.
{
  if(safe){
// a version which is safe in the case that n1,n2,n3 are not nodes of element e
    int m1 = grid.loc2glob(e,0), m2 = grid.loc2glob(e,1);
    int m3 = grid.loc2glob(e,2), m4 = grid.loc2glob(e,3);
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

    return -1; // did not find any match
  } else {
// This version is not safe if n1,n2,n3 are not in element e
// but is is faster than the safe version ...
    int n = grid.loc2glob(e,0);
    if( n != n1 && n != n2 && n != n3 ) return 1;
    n = grid.loc2glob(e,1);
    if( n != n1 && n != n2 && n != n3 ) return 2;
    n = grid.loc2glob(e,2);
    if( n != n1 && n != n2 && n != n3 ) return 0;
    return 3;
  }
}


void findBoundaryNodes(WavesGridB& grid, MV_Vector<int>& boundaryMarkers,
		       MV_ColMat<int>& elmNeighbor)
{

// Initiate the elm_neigh structure so that it contains the numbers of
// the neighbors elements on the other side of the element faces.

  int n,e,n1,n2,n3,n4,ne;
  int nsd = grid.getNoSpaceDim();

  if(nsd==2){  // For a triangulation of linear triangle element
  // find the nodes which lies on the outer boundary
  // The boundary nodes are marked with 1 the interior has the value 0

  int n,e,ne,j,n1,n2,k;
  int nno = grid.getNoNodes();
  int nel = grid.getNoElms();
  MV_Vector<int> countn2e( nno,0 );
  MV_Vector<int> n2ejcol(nel*3); // each element is added to exactly 3 nodes
  MV_Vector<int> n2eirow(nno+1);

  for ( e=0; e<nel; e++ ) { // count the number of elements a node is in
    countn2e( grid.loc2glob(e,0) )++;
    countn2e( grid.loc2glob(e,1) )++;
    countn2e( grid.loc2glob(e,2) )++;
  }
 
  int prev = n2eirow(0) = 0;
  for (n=0; n<nno; n++ )  
    prev = n2eirow( n+1 ) = countn2e( n ) + prev;
  countn2e=0;
  for ( e=0; e<nel; e++ ){
    ne = grid.loc2glob(e,0);
    n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
    ne = grid.loc2glob(e,1);
    n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
    ne = grid.loc2glob(e,2);
    n2ejcol ( n2eirow( ne ) + countn2e( ne )++ ) = e;
  }

// it is faster to go backwards due to the fact that the elements in n2e
// are ordered starting with lowest element number and increases.
// Therefore the neighbor element will be found earlier.

  //  MV_ColMat<int> elmNeighbor(nel,3);
  elmNeighbor.newsize(nel,3);
  elmNeighbor = -1;

  for( e=nel-1; e>=0; e-- ) 
    for( j=0; j<3; j++ ){
      if( elmNeighbor(e,j)==-1 ) { // if not yet set
	n1 = grid.loc2glob( e, j ); // first node
	n2 = grid.loc2glob( e, j==2 ? 0 : j+1 );
	int kstop = n2eirow(n1+1);
	for( k=n2eirow(n1); k<kstop; k++ ) { // search for neighbor element
	  ne = n2ejcol(k);
	  if ( ne != e ) { // find second common node
	    if ( grid.loc2glob (ne, 0) == n2 ) {
	      elmNeighbor(e,j) = ne;
	      elmNeighbor(ne,0) = e;
	      k = kstop; // to break the k-loop
	    }
	    else if ( grid.loc2glob (ne, 1) == n2 ) {
	      elmNeighbor(e,j) = ne;
	      elmNeighbor(ne,1) = e;
	      k = kstop; // to break the k-loop
	    }
	    else if ( grid.loc2glob (ne, 2) == n2 ) {
	      elmNeighbor(e,j) = ne;
	      elmNeighbor(ne,2) = e;
	      k = kstop; // to break the k-loop
	    }
	  }
	}
      }
    }

  boundaryMarkers.newsize(nno);
  boundaryMarkers = 0;

  for ( e=0; e<nel; e++ ) {
    if(elmNeighbor(e,0)==-1){
      boundaryMarkers(grid.loc2glob(e,0))=1;
      boundaryMarkers(grid.loc2glob(e,1))=1;
    }
    if(elmNeighbor(e,1)==-1){
      boundaryMarkers(grid.loc2glob(e,1))=1;
      boundaryMarkers(grid.loc2glob(e,2))=1;
    }
    if(elmNeighbor(e,2)==-1){
      boundaryMarkers(grid.loc2glob(e,2))=1;
      boundaryMarkers(grid.loc2glob(e,0))=1;
    }
  }

  } else if(nsd==3) {
    int n,e,n1,n2,n3,n4,ne;
// Given face of element e, find the neighbor elements which have 
// three nodes in common with the face. If no such element is found this
// is a boundary face.

  //  additional_info.initNeighbor(*this,true,false,false);
  WavesNeighborFE& neighbor = grid.getNeighbor();
  if( neighbor.nodeSize() == 0 )
    neighbor.init(grid,true,false,false);

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3
  int nel = grid.getNoElms();
  int nno = grid.getNoNodes();
  
  //  MV_ColMat<int> elmNeighbor(nel,4);
  elmNeighbor.newsize(nel,4);
  elmNeighbor = -1;

  bool safe = true;
  int ne1,ne2;
  for( e=nel-1; e>=0; e-- ){ // more efficient with decreasing elements
    n1 = grid.loc2glob( e, 0 );
    n2 = grid.loc2glob( e, 1 );
    n3 = grid.loc2glob( e, 2 );
    n4 = grid.loc2glob( e, 3 );
    if( elmNeighbor(e,0)==-1 )
      if( elmNeighbor(e,3)==-1 ){// case face 1 and 4 not done
	getElmsAtFaces(neighbor,n1,n2,n4,n3,e,ne1,ne2);
	if( ne1>=0 ){
	  n = getFaceNo(grid,n1,n2,n4,ne1,safe);
	  elmNeighbor(e,0)=ne1;
	  elmNeighbor(ne1,n)=e;
	}
	if( ne2>=0 ){
	  n = getFaceNo(grid,n1,n2,n3,ne2,safe);
	  elmNeighbor(e,3)=ne2;
	  elmNeighbor(ne2,n)=e;
	}
      }
      else {
	ne = getElmAtFace(neighbor,n1,n2,n4,e); // case face 1
	if( ne>=0 ){
	  n = getFaceNo(grid,n1,n2,n4,ne,safe);
	  elmNeighbor(e,0)=ne;
	  elmNeighbor(ne,n)=e;
	}
      }
    else if( elmNeighbor(e,3)==-1 ) {
      ne = getElmAtFace(neighbor,n1,n2,n3,e);    // case face 4
      if( ne>=0 ){
	n = getFaceNo(grid,n1,n2,n3,ne,safe);
	elmNeighbor(e,3)=ne;
	elmNeighbor(ne,n)=e;
      }
    }

    if( elmNeighbor(e,1)==-1 )
      if( elmNeighbor(e,2)==-1 ){
	getElmsAtFaces(neighbor,n3,n4,n2,n1,e,ne1,ne2);
	if( ne1>=0 ){
	  n = getFaceNo(grid,n3,n4,n2,ne1,safe);
	  elmNeighbor(e,1)=ne1;
	  elmNeighbor(ne1,n)=e;
	}
	if( ne2>=0 ){
	  n = getFaceNo(grid,n3,n4,n1,ne2,safe);
	  elmNeighbor(e,2)=ne2;
	  elmNeighbor(ne2,n)=e;
	}
      }
      else {
	ne = getElmAtFace(neighbor,n2,n3,n4,e);
	if( ne>=0 ){
	  n = getFaceNo(grid,n2,n3,n4,ne,safe);
	  elmNeighbor(e,1)=ne;
	  elmNeighbor(ne,n)=e;
	}
      }
    else if( elmNeighbor(e,2)==-1 ){
      ne = getElmAtFace(neighbor,n1,n3,n4,e);
      if( ne>=0 ){
	n = getFaceNo(grid,n1,n3,n4,ne,safe);
	elmNeighbor(e,2)=ne;
	elmNeighbor(ne,n)=e;
      }
    }
  }
  neighbor.remove();
  boundaryMarkers.newsize(nno);
  boundaryMarkers = 0;

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3

  for ( e=0; e<nel; e++ ) {
    if(elmNeighbor(e,0)==-1){
      boundaryMarkers(grid.loc2glob(e,0))=1;
      boundaryMarkers(grid.loc2glob(e,1))=1;
      boundaryMarkers(grid.loc2glob(e,3))=1;
    }
    if(elmNeighbor(e,1)==-1){
      boundaryMarkers(grid.loc2glob(e,1))=1;
      boundaryMarkers(grid.loc2glob(e,2))=1;
      boundaryMarkers(grid.loc2glob(e,3))=1;
    }
    if(elmNeighbor(e,2)==-1){
      boundaryMarkers(grid.loc2glob(e,0))=1;
      boundaryMarkers(grid.loc2glob(e,2))=1;
      boundaryMarkers(grid.loc2glob(e,3))=1;
    }
    if(elmNeighbor(e,3)==-1){
      boundaryMarkers(grid.loc2glob(e,0))=1;
      boundaryMarkers(grid.loc2glob(e,1))=1;
      boundaryMarkers(grid.loc2glob(e,2))=1;
    }
  }
  }
}
