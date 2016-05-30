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
WavesNeighborFE - neighbor information in finite element grids (GridFE)
 
*/


#ifndef WavesNeighborFE__H
#define WavesNeighborFE__H

#include "wavesSparseDS.h"

class WavesGridB;

//-----------------------------------------------------------------------

class WavesNeighborFE
{

private:

  WavesSparseDS n2e;
  WavesSparseDS n2n;
  WavesSparseDS e2e;

  bool special_e2e;  // parameter for side element-element info

public:

  WavesNeighborFE();
 ~WavesNeighborFE(){}

  void init (WavesGridB& grid);

  void init (WavesGridB& grid, 
	     bool init_n2e, // if true save node-to-element info
	     bool init_n2n, // if true save node-to-node info
	     bool init_e2e, // if true save element-to-element info
	     bool side_e2e = false);//if true only elms with common side

  void remove ( bool remove_n2e = true,
		bool remove_n2n = true, 
		bool remove_e2e = true );

  bool sideElements() { return special_e2e;}

  bool ok () const;

  int nodeSize      () const { return n2e.getNoNonzeroes(); }
  int couplingsSize () const { return n2n.getNoNonzeroes(); }
  int elementSize   () const { return e2e.getNoNonzeroes(); }

  int nodeIrow      (int n) const { return n2e.irow(n); }
  int couplingsIrow (int n) const { return n2n.irow(n); }
  int elementIrow   (int e) const { return e2e.irow(e); }

  int nodeJcol      (int n) const { return n2e.jcol(n); }
  int couplingsJcol (int n) const { return n2n.jcol(n); }
  int elementJcol   (int e) const { return e2e.jcol(e); }

  void print();
};

/*

NAME:  WavesNeighborFE - neighbor information in finite element grids (GridFE)

KEYWORDS:

  finite element grid, neighbor information, neighbor elements


DESCRIPTION:

  The class represents the neighbor information in a finite element grid.
  Three types of neighbor information are computed: (1) the elements that 
  are neighbors of a node, (2) the nodes that are neighbors to a node
  including the node itself, and (3) the elements that are neighbors of 
  an element excluding the element itself.
  See descriptions of the corresponding member functions for a definition 
  of the three types.
  See "GridFEInfo" for the possibility to use the computed neighbor
  information for a grid without reinitialization.


CONSTRUCTORS AND INITIALIZATION:

  Only an empty constructor is offered.
  To initialize the class, one must call the function "init".
  This function requires a "GridFE" object, and
  generates the neighbor information. In one version of "init" it
  is possible to specify which of the neighbor items to be computed.


MEMBER FUNCTIONS:

   Internally the neighbor information is stored in SparseDS structures.
   The most efficient way to access the neighbor information is to use the
   "---Irow" and "---Jcol" functions. 

  "init" - initialize the neighbor information. The overloaded version 
           is used to compute only parts of the information.
	   The fourth argument, which default is false, refers to the 
	   element-to-element computation; if "side_e2e = false" then 
	   only the elements which have at least "nsd" common nodes with 
	   the element are added to the list, where "nsd" is the space 
	   dimension of the grid.

  "remove" - remove all or parts of the neighbor information.

  "node" - given a node number, the function returns a list of all elements
           which the node is part of. 

  "couplings" - same function as "nodes".

  "element" - given an element number, the function returns a list of
              all elements that are neighbors of this element.
	      In the standard version two elements are neighbors if they 
	      have at least one node in common.
              Observe that this definition of neighbors is different from
              other commonly applied definitions, for example, that two
              elements are neighbors if they share a common side.
	      
  "nodes" - given a node number, the function returns a list of all the
            nodes that are coupled to the given node. This includes 
	    constraints if the grid has constraints ("Grid2FE").

  "sideElements" - is true if the element-to-element info only contains
            elements with "nsd", the number of space dimension,
	    common nodes, see "init".

  "nodeSize", "couplingsSize", "nodeIrow"  - total size of node-element, 
               node-node and element-element, respectively

  "nodeIrow", "couplingsIrow", "elementIrow" - give where the row starts, 
                          the pointer refers to the "---Jcol".

  "nodeJcol", "couplingsJcol", "elementJcol" - accessing the items, as given
                       by "---Irow".

  "print" - prints the neighbor information.


EXAMPLES:

  |     GridFE grid;
  |     // generate the grid...
  |     WavesNeighborFE neighbor;
  |     neighbor.init(grid);
  |     neighbor.print(s_o);
  | 
  |     // element-element neighbor info.
  |     int e=1;
  |     VecSimplest(int) neighbor_el;
  |     neighbor.element(e,neighbor_el);
  |     neighbor_el.print(s_o);
  | 
  |     // node-element neighbor info,
  |     int n=1;
  |     neighbor.node(n,neighbor_el);
  |     neighbor_el.print(s_o);
  | 
  |     // Example of access using the Irow- and Jcol-structure:
  |
  |     int ele = 1;
  |     int elestart = neighbor.elementIrow(ele);
  |     int elestop = neighbor.elementIrow(ele+1);
  |     for(n=elestart; n<elestop; n++) {
  |	  int element_neighbor = neighbor.elementJcol(n);
  |       cout << "neighbor element " << element_neighbor << endl;
  |     }
  |


*/
#endif
