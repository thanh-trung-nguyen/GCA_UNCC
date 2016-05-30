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

#include "include/wavesSDGeometry.h"

#include <iostream>
#include <math.h>

using namespace std;

//
//--  Constructors etc
//

bool WavesSDGeometry::initialize(const int &n_i_, const int &n_j_, const int &n_k_, const int &nsd_, const real &dx_, const real dy_, const real &dz_, const real &xmin_, const real ymin_, const real &zmin_)
{
	dx = dx_;
	dy = dy_;
	dz = dz_;

	n_i = n_i_;
	n_j = n_j_;
	n_k = n_k_;
	nsd = nsd_;

	xmin = xmin_;
	ymin = ymin_;
	zmin = zmin_;

	return status();
}

WavesSDGeometry::WavesSDGeometry()
{
	initialize(0, 0, 0, 0, 1., 1., 1., 0., 0., 0.);
}

WavesSDGeometry::WavesSDGeometry(const int &n_i_)
{
	// A 1d grid is inited
	initialize(n_i_, 0, 0, 1, 1., 1., 1., 0., 0., 0.);
}

WavesSDGeometry::WavesSDGeometry(const int &n_i_, const int &n_j_)
{
	// A 2d grid is inited
	initialize(n_i_, n_j_, 0, 2, 1., 1., 1., 0., 0., 0.);
}

WavesSDGeometry::WavesSDGeometry(const int &n_i_, const int &n_j_, const int &n_k_)
{
	// A 3d grid is inited
	initialize(n_i_, n_j_, n_k_, 3, 1., 1., 1., 0., 0., 0.);
}

WavesSDGeometry::WavesSDGeometry(const int &n_i_, const int &n_j_, const int &n_k_, const int &nsd_, const real &dx_, const real dy_, const real &dz_)
{
	// All arguments inited
	initialize(n_i_, n_j_, n_k_, nsd_, dx_, dy_, dz_, 0., 0., 0.);
}

WavesSDGeometry::WavesSDGeometry(const int &n_i_, const int &n_j_, const int &n_k_, const int &nsd_, const real &dx_, const real dy_, const real &dz_, const real &xmin_, const real ymin_, const real &zmin_)
{
	// Everything is inited
	initialize(n_i_, n_j_, n_k_, nsd_, dx_, dy_, dz_, xmin_, ymin_, zmin_);
}

WavesSDGeometry::WavesSDGeometry(const WavesSDGeometry &sdg)
{
	// Everything is copied
	initialize(sdg.n_i, sdg.n_j, sdg.n_k, sdg.nsd, sdg.dx, sdg.dy, sdg.dz, sdg.xmin, sdg.ymin, sdg.zmin);
}

WavesSDGeometry::~WavesSDGeometry()
{
}

WavesSDGeometry &WavesSDGeometry::operator=(const WavesSDGeometry &sdg)
{
	if (this == &sdg)
		return *this;
	initialize(sdg.n_i, sdg.n_j, sdg.n_k, sdg.nsd, sdg.dx, sdg.dy, sdg.dz, sdg.xmin, sdg.ymin, sdg.zmin);

	return *this;
}



void WavesSDGeometry::presentWavesSDGeometry() const
{
	cout << "----- Present WavesSDGeometry " << endl;
	cout << "Status is : " << (status() ? 'T' : 'F') << endl;

	cout << "Dimensions(" << nsd << "): " << n_i << " x " << n_j << " x " << n_k << endl;
	cout << "Mesh size: " << dx << ", " << dy << ", " << dz << endl;

	cout << "Min coords: " << xmin << ", " << ymin << ", " << zmin << endl;

	cout << "--" << endl;
}



bool WavesSDGeometry::status() const
{
	bool s = 1;
	switch (nsd)
	{
		case 3:
			s = s && dz > 0. && n_k > 0;
			break;
		case 2:
			s = s && dy > 0. && n_j > 0;
			break;
		case 1:
			s = s && dx > 0. && n_i > 0;
			return s;
		default:
			return 0;
	}

	return 0;
}

// Methods with interface conformant with unstructured grids =================================

int WavesSDGeometry::getNoNodes() const //Get total no of nodes
{
	switch (nsd)
	{
		case 1:
			return n_i;
		case 2:
			return n_i * n_j;
		case 3:
			return n_i * n_j * n_k;
		default:
			return 0;
	}
}

int WavesSDGeometry::getNoElms() const //Get total no of elements
{
	switch (nsd)
	{
		case 1:
			return n_i - 1;
		case 2:
			return (n_i - 1) * (n_j - 1);
		case 3:
			return (n_i - 1) * (n_j - 1) * (n_k - 1);
		default:
			return 0;
	}
}

int WavesSDGeometry::getMaxNoNodesInElm() const //Get no nodes in one element
{
	switch (nsd)
	{
		case 1:
			return 2;
		case 2:
			return 4;
		case 3:
			return 8;
		default:
			return 0;
	}
}

bool WavesSDGeometry::redim()
{
	Material_Type = new int[getNoElms()];
	return true;
}

int WavesSDGeometry::getMaterialType(const int e) const //Get material type for element e
{
	return Material_Type[e];
}

double WavesSDGeometry::getCoor(const int n, const int d) const //Get coordinate for node n, coordinate d 
{
	switch (d)
	{
		case 0:
			return getX(node2i(n));
		case 1:
			return getY(node2j(n));
		case 2:
			return getZ(node2k(n));
		default:
			return 0;
	}
}

int WavesSDGeometry::loc2glob(const int e, const int i) const //Get global node no from local
{
	// 1) Obtain i,j,k for element e. Similar to node2i etc, but use n_i-1 instead of n_i.
	int e_i, e_j, e_k;
	e_i = e % (n_i - 1);
	e_j = e % ((n_i - 1) * (n_j - 1)) / (n_i - 1);
	e_k = e / ((n_i - 1) * (n_j - 1));

	// 2) Compute node number. base is node x_min, y_min, z_min
	int base = e_i + n_i * (e_j + n_j * e_k);
	switch (i)
	{
		case 0:
			return base;
		case 1:
			return base + node(1, 0, 0);
		case 2:
			return base + node(0, 1, 0);
		case 3:
			return base + node(1, 1, 0);
		case 4:
			return base + node(0, 0, 1);
		case 5:
			return base + node(1, 0, 1);
		case 6:
			return base + node(0, 1, 1);
		case 7:
			return base + node(1, 1, 1);
		default:
			return 0;
	}

}

// gives "global" (natural) node number for coord x, y, z.
// Answer is negative unless the given coords are within 
// (dx,dy,dz)*to from a node.  tol = 0.01.
// Answer -1 means failure in x, -2 in y and -3 in z.
// Answer -4 means other failure.

int WavesSDGeometry::coord2node(const real &x, const real &y, const real &z) const
{
	int i, j, k;
	real ir, jr, kr;

	real tol = 0.01;

	if (nsd != 2 && nsd != 3)
		return -4;

	ir = (x - xmin) / dx;
	if (fabs(ir - rint(ir)) < tol)
		i = rint(ir);
	else
		return -1; // failure in x

	jr = (y - ymin) / dy;
	if (fabs(jr - rint(jr)) < tol)
		j = rint(jr);
	else
		return -2; // failure in y

	if (nsd == 2)
		return node(i, j);

	kr = (z - zmin) / dz;
	if (fabs(kr - rint(kr)) < tol)
		k = rint(kr);
	else
		return -3; // failure in z

	return node(i, j, k);
}


// get a pointer to an array of node numbers for a given value of z in 3D. Added by Thanh.
void WavesSDGeometry::coord2node3Dz(const double z, MV_Vector<int>& XYnodes)
{
  int  i, j, k, Nxy, idx, idxz;
        real kr;
        real tol = 0.01;

	kr = (z - zmin) / dz; 	
	k = rint(kr); // index in the z-direction

	Nxy = n_i*n_j;
	idxz= n_i*n_j*k; 
	for (j=0; j < n_j; j++)
	  { for (i =0; i<n_i; i++)
	      {
		idx = i+j*n_i;
		XYnodes(idx) = idx + idxz;
	      }
	  }

}

