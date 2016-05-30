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
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESGRIDB_H
#define __WAVESGRIDB_H

#include "wavesMVvtp.h"
#include "wavesMVvind.h"
#include "wavesMVmtp.h"
#include <stdio.h>
#include "wavesNeighborFE.h"
#include "wavesElementType.h" // Definition of enum variable for element type
#ifdef CC_42
typedef int bool;
#endif
#define true  1
#define false 0
typedef double real;

class WavesGridB
{
	public:

		WavesGridB();
		~WavesGridB();
		int getNoSpaceDim() const
		{
			return nsd;
		}
		int getNoElms() const
		{
			return nel;
		}
		int getNoNodes() const
		{
			return nno;
		}
		int getMaxNoNodesInElm() const
		{
			return maxnne;
		}
		int getMaterialType(const int e) const
		{
			return elid(e);
		}
		int getNoInd() const
		{
			return nbind;
		}
		int getInd(const int n, const int i) const
		{
			return bind(n, i);
		}
		void indOn(const int n, const int i)
		{
			bind(n, i) = 1;
		}
		void indOff(const int n, const int i)
		{
			bind(n, i) = 0;
		}
		void setMaterialType(const int e, const int mat)
		{
			elid(e) = mat;
		}
		real getCoor(const int n, const int d) const
		{
			return coord(n, d);
		}
		void putCoor(const int n, const int d, const real c)
		{
			coord(n, d) = c;
		}
		int loc2glob(const int e, const int i) const
		{
			return nodpek(e, i);
		}
		void putLoc2glob(const int e, const int i, const int nodno)
		{
			nodpek(e, i) = nodno;
		}

		int Coord2Node(double x, double y, double z = 0.0) const;
		bool oneElementTypeInGrid() const
		{
			return (elm_type_for_grid == NO_ELEMENT) ? false : true;
		}
		ElementType getElementType(int e) const
		{
			return (elm_type_for_grid == NO_ELEMENT) ? element_type(e) : elm_type_for_grid;
		}

		int getNoNodesInElm(const int e) const
		{
			return (elm_type_for_grid == NO_ELEMENT) ? getNoNodesForElmType(element_type(e)) : maxnne;
		}
		bool redim(const int nsd, const int nno, const int nel, const int maxnne, const int nbind, const ElementType e_type = NO_ELEMENT);
		bool ok() const;
		void operator =(const WavesGridB& grid);
		void print(const char *file) const;
		void print_inp_amira(const char *file) const;
		void print_vel(const char *file, MV_Vector<int>& velocity) const;
		void scan(const char *file);
		void scan_2Dgid(const char *file, int elnum1);

		WavesNeighborFE& getNeighbor()
		{
			return neighbor;
		}

		void getElmFromNodes(const int n1, const int n2, int& e) const;

		int getNoNodesForElmType(const ElementType etype) const;
		void scan_avs_vel(const char *file, MV_Vector<double>& velocity);
		void scan_avs_amira(const char *file);
		void scan_gid(const char *file, int elnum1);
		void scan_gid_difmat(const char *file);
		void scan_gid(const char *file, bool type_of_mat = 0);

	protected:
		void scan_avs(const char *file);

		bool checkIfOnlyOneElmTypeInGrid();
		int nsd;
		int nno;
		int nel;
		int maxnne;
		int nbind;
		int current_node;

		bool onemat;
		ElementType elm_type_for_grid;

		MV_ColMat<real> coord;
		MV_ColMat<int> nodpek;
		MV_ColMat<int> bind;
		MV_Vector<int> elid;
		MV_Vector<ElementType> element_type;
		WavesNeighborFE neighbor;
};
/*>WavesGridB:*/

/*Class:WavesGridB

 NAME:  WavesGridB 

 SYNTAX:     @WavesGridB


 KEYWORDS:


 DESCRIPTION:


 CONSTRUCTORS AND INITIALIZATION:

 MEMBER FUNCTIONS:


 EXAMPLES:

 SEE ALSO: 

 DEVELOPED BY:   


 AUTHOR:	         
 
 Klas Samuelsson

 End: 
 */

#endif

/* LOG HISTORY of this file:

 * $Log: WavesGridB.h,v $
 * Revision 
 * Version 
 *
 */

