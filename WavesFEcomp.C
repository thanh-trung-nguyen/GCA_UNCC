//SINTEF

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

#include "include/wavesFEcomp.h"

using namespace std;


FET3n2D::FET3n2D()

{
	grid_ = NULL;
}


void FET3n2D::attach(const WavesGridB* grid)

{
	this->grid_ = (WavesGridB*) grid;
}


FET3n2D::FET3n2D(const WavesGridB* grid)

{
	attach(grid);
}


FET3n2D::~FET3n2D()

{
}


bool FET3n2D::checkGrid() const

// check that the grid only has ElmT3n2D elements
// (could also be a possibility that allElmsOfSameType says
//  the grid is heterogeneous but it still only has ElmT3n2D elements)
{
	if (!ok() || !(grid().getNoSpaceDim() == 2))
		return false;
	return true;
}


bool FET3n2D::checkGrid(const int e) const

// check that element e is of type ElmT3n2D 
{
	return true;
}


bool FET3n2D::ok() const

{
	bool good;
	if (grid_ == NULL)
		good = false;
	else
		good = true;

	if (!good)
		cout << "errorFP FET3n2D::ok, grid is not attached\n";

	return good;
}


void FET3n2D::refill(const int e)

{
	elm_no_ = e;
	WavesGridB& g = grid();

	n1_ = g.loc2glob(e, 0);
	n2_ = g.loc2glob(e, 1);
	n3_ = g.loc2glob(e, 2);

	x1_ = g.getCoor(n1_, 0);
	y1_ = g.getCoor(n1_, 1);
	x2_ = g.getCoor(n2_, 0);
	y2_ = g.getCoor(n2_, 1);
	x3_ = g.getCoor(n3_, 0);
	y3_ = g.getCoor(n3_, 1);

	g2x = y3_ - y1_;
	g2y = x1_ - x3_;
	g3x = y1_ - y2_;
	g3y = x2_ - x1_;
	det_ = g3y * g2x - g3x * g2y;

	d = 1 / det_;

	g2x *= d;
	g2y *= d;
	g3x *= d;
	g3y *= d;
	g1x = -g2x - g3x;
	g1y = -g2y - g3y;
}


real FET3n2D::angle(const int i) const

{
	switch (i)
	{
		case 1:
			return angle1();
		case 2:
			return angle2();
		case 3:
			return angle3();
		default:
			cout << "errorFP FET3n2D:: angle, wrong angle number\n";
			break;
	}
	return -99999.0;
}


real FET3n2D::edgeLength(const int i) const

{
	switch (i)
	{
		case 1:
			return l1();
		case 2:
			return l2();
		case 3:
			return l3();
		default:
			cout << "errorFP  FET3n2D:: edgeLength wrong edge number\n";
			break;
	}
	return -99999.0;
}


real FET3n2D::qualityMeasure(const int i) const

{
	switch (i)
	{
		case 1:
			return etahH();
		case 2:
			return etaincirc();
		case 3:
			return etaMaxAngle();
		case 4:
			return etapltmg();
		default:
			cout << "errorFP FET3n2D:: qualityMeasure\n";
			break;
	}
	return -9999.0;
}


void FET3n2D::print() const

{
	cout << "FET3n2D:: print (Os os) not finished\n";
}

bool FET3n2D::areQuadNodesAtMidPoints() const

// test if nodes 4, 5 and 6 are at midpoint of edges 1, 2 and 3, respectively.
{
// assume that the element is of type ElmT6n2D
	WavesGridB& g = grid();
	real eps2 = 1e-20 * det(); //  square of allowed distance, size?
	int n = g.loc2glob(elm_no_, 4);
	if (sqr(g.getCoor(n, 1) - midpt1x()) + sqr(g.getCoor(n, 2) - midpt1y()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 5);
	if (sqr(g.getCoor(n, 1) - midpt2x()) + sqr(g.getCoor(n, 2) - midpt2y()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 6);
	if (sqr(g.getCoor(n, 1) - midpt3x()) + sqr(g.getCoor(n, 2) - midpt3y()) > eps2)
		return false;
	return true;
}


// ****************************************************************************
// start of tetrahedral case
// ****************************************************************************



FET4n3D::FET4n3D()

{
	this->grid_ = NULL;
}


void FET4n3D::attach(const WavesGridB* grid)

{
	this->grid_ = (WavesGridB*) grid;
}


FET4n3D::FET4n3D(const WavesGridB* grid)

{
	attach(grid);
}


FET4n3D::~FET4n3D()

{
}


bool FET4n3D::checkGrid() const

// check that the grid only has ElmT4n3D elements
{
	if (!ok() || !(grid().getNoSpaceDim() == 3))
		return false;
	return true;
}


bool FET4n3D::checkGrid(const int /* e*/) const

// check that element e is of type ElmT4n3D elements
{
	if (!ok() || !(grid().getNoSpaceDim() == 3))
		return false;
	return true;
}


bool FET4n3D::ok() const

{
	bool good;
	if (grid_ == NULL)
		good = false;
	else
		good = true;

	if (!good)
		cout << "errorFP FET4n3D::ok grid is not attached\n";

	return good;
}


void FET4n3D::refill(const int e)

{

	elm_no_ = e;
	WavesGridB& g = grid();
	n1_ = g.loc2glob(e, 0);
	n2_ = g.loc2glob(e, 1);
	n3_ = g.loc2glob(e, 2);
	n4_ = g.loc2glob(e, 3);

	x1_ = g.getCoor(n1_, 0);
	y1_ = g.getCoor(n1_, 1);
	z1_ = g.getCoor(n1_, 2);
	x2_ = g.getCoor(n2_, 0);
	y2_ = g.getCoor(n2_, 1);
	z2_ = g.getCoor(n2_, 2);
	x3_ = g.getCoor(n3_, 0);
	y3_ = g.getCoor(n3_, 1);
	z3_ = g.getCoor(n3_, 2);
	x4_ = g.getCoor(n4_, 0);
	y4_ = g.getCoor(n4_, 1);
	z4_ = g.getCoor(n4_, 2);

	j11 = x2_ - x1_;
	j12 = y2_ - y1_;
	j13 = z2_ - z1_;
	j21 = x3_ - x1_;
	j22 = y3_ - y1_;
	j23 = z3_ - z1_;
	j31 = x4_ - x1_;
	j32 = y4_ - y1_;
	j33 = z4_ - z1_;

	g2x = j22 * j33 - j23 * j32;
	g3x = j13 * j32 - j12 * j33;
	g4x = j12 * j23 - j13 * j22;
	g2y = j23 * j31 - j21 * j33;
	g3y = j11 * j33 - j13 * j31;
	g4y = j13 * j21 - j11 * j23;
	g2z = j21 * j32 - j22 * j31;
	g3z = j12 * j31 - j11 * j32;
	g4z = j11 * j22 - j12 * j21;

	det_ = j11 * g2x + j12 * g2y + j13 * g2z;

	d = 1.0 / det_;

	g2x *= d;
	g3x *= d;
	g4x *= d;
	g2y *= d;
	g3y *= d;
	g4y *= d;
	g2z *= d;
	g3z *= d;
	g4z *= d;

	g1x = -g2x - g3x - g4x;
	g1y = -g2y - g3y - g4y;
	g1z = -g2z - g3z - g4z;
}


real FET4n3D::solidAngle(const int i) const

{
	switch (i)
	{
		case 1:
			return solidAngle1();
		case 2:
			return solidAngle2();
		case 3:
			return solidAngle3();
		case 4:
			return solidAngle4();
		default:
			cout << "errorFP FET4n3D:: solidAngle wrong angle number\n";
			break;
	}
	return -9999.0;
}


real FET4n3D::dihedralAngle(const int i) const

{
	switch (i)
	{
		case 1:
			return dihedralAngle1();
		case 2:
			return dihedralAngle2();
		case 3:
			return dihedralAngle3();
		case 4:
			return dihedralAngle4();
		case 5:
			return dihedralAngle5();
		case 6:
			return dihedralAngle6();
		default:
			cout << "errorFP FET4n3D:: sdihedralAngle wrong angle number\n";
			break;
	}
	return -9999.0;
}


real FET4n3D::maxSolidAngle() const

{
	real s1 = 1 / sqrt(a11()), s2 = 1 / sqrt(a22()), s3 = 1 / sqrt(a33()), s4 = 1 / sqrt(a44());
	real d1 = acos(a34() * s3 * s4);
	real d2 = acos(a14() * s1 * s4);
	real d3 = acos(a24() * s2 * s4);
	real d4 = acos(a23() * s2 * s3);
	real d5 = acos(a13() * s1 * s3);
	real d6 = acos(a12() * s1 * s2);
	return 2 * PI - min(min(d1 + d3 + d4, d1 + d2 + d5), min(d2 + d3 + d6, d4 + d5 + d6));
}


real FET4n3D::minSolidAngle() const

{
	real s1 = 1 / sqrt(a11()), s2 = 1 / sqrt(a22()), s3 = 1 / sqrt(a33()), s4 = 1 / sqrt(a44());
	real d1 = acos(a34() * s3 * s4);
	real d2 = acos(a14() * s1 * s4);
	real d3 = acos(a24() * s2 * s4);
	real d4 = acos(a23() * s2 * s3);
	real d5 = acos(a13() * s1 * s3);
	real d6 = acos(a12() * s1 * s2);
	return 2 * PI - max(max(d1 + d3 + d4, d1 + d2 + d5), max(d2 + d3 + d6, d4 + d5 + d6));
}


void FET4n3D::minMaxSolidAngle(real& minsa, real& maxsa) const

{
	real s1 = 1 / sqrt(a11()), s2 = 1 / sqrt(a22()), s3 = 1 / sqrt(a33()), s4 = 1 / sqrt(a44());
	real d1 = acos(a34() * s3 * s4);
	real d2 = acos(a14() * s1 * s4);
	real d3 = acos(a24() * s2 * s4);
	real d4 = acos(a23() * s2 * s3);
	real d5 = acos(a13() * s1 * s3);
	real d6 = acos(a12() * s1 * s2);
	real twopi = 2 * PI;
	s1 = twopi - (d1 + d3 + d4);
	s2 = twopi - (d1 + d2 + d5);
	s3 = twopi - (d2 + d3 + d6);
	s4 = twopi - (d4 + d5 + d6);
	if (s1 < s2)
	{
		minsa = s1;
		maxsa = s2;
	}
	else
	{
		minsa = s2;
		maxsa = s1;
	}
	if (s3 < minsa)
		minsa = s3;
	else if (s3 > maxsa)
		maxsa = s3;
	if (s4 < minsa)
		minsa = s4;
	else if (s4 > maxsa)
		maxsa = s4;
}


void FET4n3D::inCenter(real& cx, real& cy, real& cz) const

// center of the inscribed circle
{
	real s1 = sqrt(a11()), s2 = sqrt(a22()), s3 = sqrt(a33()), s4 = sqrt(a44());
	real Oinv = 1. / (s1 + s2 + s3 + s4);
// weighting with the face areas opposite the nodes
	cx = (x1_ * s1 + x2_ * s2 + x3_ * s3 + x4_ * s4) * Oinv;
	cy = (y1_ * s1 + y2_ * s2 + y3_ * s3 + y4_ * s4) * Oinv;
	cz = (z1_ * s1 + z2_ * s2 + z3_ * s3 + z4_ * s4) * Oinv;
}


real FET4n3D::circumRadius() const

// radius of circumscribed circle
{
//  two alternatives:
	real a = l1sq(), b = l3sq(), c = l4sq();
	return 0.5 * sqrt(sqr(a * g2x + b * g3x + c * g4x) + sqr(a * g2y + b * g3y + c * g4y) + sqr(a * g2z + b * g3z + c * g4z));
}


real FET4n3D::edgeLength(const int i) const

{
	switch (i)
	{
		case 1:
			return l1();
		case 2:
			return l2();
		case 3:
			return l3();
		case 4:
			return l4();
		case 5:
			return l5();
		case 6:
			return l6();
		default:
			cout << "errorFP FET4n3D:: edgeLength wrong edge number\n";
			break;
	}
	return 0.0;
}

real FET4n3D::qualityMeasure(const int i) const

{
	switch (i)
	{
		case 1:
			return eta();
		case 2:
			return etaincirc();
		case 3:
			return etaMinSolidAngle();
		case 4:
			return etaMaxSolidAngle();
		case 5:
			return etahH();
		default:
			cout << "errorFP FET4n3D:: qualityMeasure wrong argument\n";
			break;
	}
	return -9999.0;
}

real FET4n3D::eta() const

{ // from geompack
	return 7.559526299369240 * pow(det_, 0.6666666666666666) / (l1sq() + l2sq() + l3sq() + l4sq() + l5sq() + l6sq());
}

real FET4n3D::etaincirc() const

{
	return 3 * inRadius() / circumRadius();
}

real FET4n3D::etaMinSolidAngle() const

{
	return 1.813941816806565 * minSolidAngle();
}

real FET4n3D::etaMaxSolidAngle() const

{
	return (2 * PI - maxSolidAngle()) * 1.744622290711000e-01;
}

real FET4n3D::etahH() const

{ // normalized by sqrt(3./2.) 
	return 1.224744871391589 * h() / H();
}


void FET4n3D::print() const

{
	cout << "FET4n3D:: print (Os os) not finished\n";
}

bool FET4n3D::areQuadNodesAtMidPoints() const

// test if the nodes 5-10 are at midpoints of edges 1-6,respectively
{
// assume that the element is of type ElmT10n3D
	WavesGridB& g = grid();
	real eps2 = 1e-20 * l1sq(); // square of allowed distance, size?
	int n = g.loc2glob(elm_no_, 5);
	if (sqr(g.getCoor(n, 1) - midpt1x()) + sqr(g.getCoor(n, 2) - midpt1y()) + sqr(g.getCoor(n, 3) - midpt1z()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 6);
	if (sqr(g.getCoor(n, 1) - midpt2x()) + sqr(g.getCoor(n, 2) - midpt2y()) + sqr(g.getCoor(n, 3) - midpt2z()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 7);
	if (sqr(g.getCoor(n, 1) - midpt3x()) + sqr(g.getCoor(n, 2) - midpt3y()) + sqr(g.getCoor(n, 3) - midpt3z()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 8);
	if (sqr(g.getCoor(n, 1) - midpt4x()) + sqr(g.getCoor(n, 2) - midpt4y()) + sqr(g.getCoor(n, 3) - midpt4z()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 9);
	if (sqr(g.getCoor(n, 1) - midpt5x()) + sqr(g.getCoor(n, 2) - midpt5y()) + sqr(g.getCoor(n, 3) - midpt5z()) > eps2)
		return false;
	n = g.loc2glob(elm_no_, 10);
	if (sqr(g.getCoor(n, 1) - midpt6x()) + sqr(g.getCoor(n, 2) - midpt6y()) + sqr(g.getCoor(n, 3) - midpt6z()) > eps2)
		return false;
	return true;
}

/* LOG HISTORY of this file:

 * $Log: FEcomp.C,v $
 *
 */
