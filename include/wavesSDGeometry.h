//KRISTER

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

#ifndef __WAVESSDGEOMETRY_H
#define __WAVESSDGEOMETRY_H

#include <iostream>
#include "wavesMVvtp.h"
#include "wavesMVvind.h"
#include "wavesMVmtp.h"
#include "wavesSDDefs.h"

/**
 WavesSDGeometry represents a structured, equidistant discretization.

 */

class WavesSDGeometry
{
	public:
		/// Default constructor produces "empty" grid
		WavesSDGeometry();
		/// Produce a 1D grid
		WavesSDGeometry(const int &n_i);
		/// Produce a 2D grid
		WavesSDGeometry(const int &n_i, const int &n_j);
		/// Produce a 3D grid
		WavesSDGeometry(const int &n_i, const int &n_j, const int &n_k);
		/// Produce either a 2d or 3d grid with specific step sizes
		WavesSDGeometry(const int &n_i, const int &n_j, const int &n_k, const int &nsd, const real &dx_, const real dy_, const real &dz_);
		/** Produce either a 2d or 3d grid with specific step sizes
		 * and located with node 0 in specified coordinates */
		WavesSDGeometry(const int &n_i_, const int &n_j_, const int &n_k_, const int &nsd_, const real &dx_, const real dy_, const real &dz_, const real &xmin_, const real ymin_, const real &zmin_);
		/// Copy constructor
		WavesSDGeometry(const WavesSDGeometry &sdg);
		/// Destructor. (Note: Not virtual.)
		~WavesSDGeometry();

		/// Assignment.
		WavesSDGeometry &operator=(const WavesSDGeometry &sdg);

		/// To initialize a grid.
		bool initialize(const int &n_i, const int &n_j, const int &n_k, const int &nsd, const real &dx_, const real dy_, const real &dz_, const real &xmin_, const real ymin_, const real &zmin_);

		bool redim();

		void setMaterialType(const int e, const int mattype) const
		{
			Material_Type[e] = mattype;
		}

		/// Returns true if size, sign of delta and dimensions are consistent.
		bool status() const;

		/// To set specific deltas
		bool setDelta(const real &dx_, const real &dy_ = 1., const real &dz_ = 1.)
		{
			dx = dx_;
			dy = dy_;
			dz = dz_;
			return status();
		}

		/// To set specific coordinates for node 0
		bool setCoordMin(const real &xmin_, const real &ymin_ = 0., const real &zmin_ = 0.)
		{
			xmin = xmin_;
			ymin = ymin_;
			zmin = zmin_;
			return status();
		}

		/// 
		real getDx() const
		{
			return dx;
		}
		///
		real getDy() const
		{
			return dy;
		}
		///
		real getDz() const
		{
			return dz;
		}

		///
		real getXmin() const
		{
			return xmin;
		}
		///
		real getYmin() const
		{
			return ymin;
		}
		///
		real getZmin() const
		{
			return zmin;
		}

		///
		int getNoSpaceDim() const
		{
			return nsd;
		}
		///
		int getN_i() const
		{
			return n_i;
		}
		///
		int getN_j() const
		{
			return n_j;
		}
		///
		int getN_k() const
		{
			return n_k;
		}

		/// Returns x coordinate given i index
		real getX(const int &i) const
		{
			return xmin + i * dx;
		}
		/// Returns y coordinate given j index
		real getY(const int &j) const
		{
			return ymin + j * dy;
		}
		/// Returns z coordinate given k index
		real getZ(const int &k) const
		{
			return zmin + k * dz;
		}

		/// Returns x coordinate given a node
		real node2x(const int &n) const
		{
			return getX(node2i(n));
		}
		/// Returns y coordinate given a node
		real node2y(const int &n) const
		{
			return getY(node2j(n));
		}
		/// Returns z coordinate given a node
		real node2z(const int &n) const
		{
			return getZ(node2k(n));
		}

		/// gives "global" (natural) node number for index i, j. (2d grid.)
		inline int node(const int &i, const int &j = 0) const
		{
			return i + n_i * j;
		}

		/// gives "global" (natural) node number for index i, j, k.
		inline int node(const int &i, const int &j, const int &k) const
		{
			return i + n_i * (j + n_j * k);
		}

		/// gives i for node number with index i, j, k.
		inline int node2i(const int &n) const
		{
			return n % n_i;
		}

		/// gives j for node number with index i, j, k.
		inline int node2j(const int &n) const
		{
			return (nsd > 1) ? (n % (n_i * n_j) / n_i) : 0;
		}

		/// gives k for node number with index i, j, k.
		inline int node2k(const int &n) const
		{
			return (nsd > 2) ? n / (n_i * n_j) : 0;
		}

		int coord2node(const real &x, const real &y, const real &z = 0) const;
		void coord2node3Dz(const double z, MV_Vector<int>& XYnodes); // get the pointer to the array of nodes for each value of z. Thanh added
		void presentWavesSDGeometry() const;

		/**@name Methods with interface conformant with unstructured grids
		 * These methods have been developed in order to be conformant with the
		 * interface for unstructured grids used in Kraftwerk. This makes it possible to use
		 * large structured grids with small memory demands.
		 */
		//@{
		/// Get total no of elements
		int getNoElms() const;
		/// Get total no of nodes
		int getNoNodes() const;
		/// Get no nodes in one element
		int getMaxNoNodesInElm() const;
		/// Get material type for element e
		int getMaterialType(const int e) const;
		/// Get coordinate for node n, coordinate d 
		double getCoor(const int n, const int d) const;
		/// Get global node no from local 
		int loc2glob(const int e, const int i) const;
		//@}

	private:
		/// number of nodes in i, j, and k direction
		int n_i, n_j, n_k;
		/// delta in i, j, and k direction
		real dx, dy, dz;
		/// number of space dimensions
		int nsd;
		/// Location for node 0
		real xmin, ymin, zmin;
		int* Material_Type;

};

#endif
