//SINITEF

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

#ifndef __WAVESFECOMP_H
#define __WAVESFECOMP_H

#include "wavesGridB.h"
#include <iostream>
// some constants
#define PI 3.141592653589793238462643383279502884197169399375105820974944592308
#define SQUAREROOTTWO 1.4142135623730950488016887242096980785696718753769480732
#define ONETHIRD 0.333333333333333333333333333333333333333333333333333333333333
#define ONESIXTH 0.166666666666666666666666666666666666666666666666666666666666

inline float max(float a, float b)
{
	return (a > b) ? a : b;
}
inline float min(float a, float b)
{
	return (a < b) ? a : b;
}
inline real sqr(float a)
{
	return a * a;
}

/*<FET3n2D:*/
class FET3n2D
{
	private:
		int elm_no_;    // element number
		WavesGridB* grid_; // finite element
		real d;               // help variable,  d = 1.0/det
		real x1_, y1_; // coordinates of nodes of triangle in counter clockwise order
		real x2_, y2_;
		real x3_, y3_;
		real det_;                        // determinant = 2*area
		real g1x, g1y, g2x, g2y, g3x, g3y;   // gradients of shape functions
		int n1_, n2_, n3_;                  // node numbers
	public:
		FET3n2D();
		FET3n2D(const WavesGridB* grid);
		~FET3n2D();
		void attach(const WavesGridB* grid);
		void refill(const int elm_no); // fill for new element         
		bool checkGrid() const;      // check that all elements in grid are ElmT3n2D
		bool checkGrid(const int e) const; // check that element e is of type ElmT3n2D
		WavesGridB& grid() const
		{
			return *grid_;
		}
		int getElmNo() const
		{
			return elm_no_;
		}

		real x1() const
		{
			return x1_;
		} // get the coordinates of nodes
		real y1() const
		{
			return y1_;
		}
		real x2() const
		{
			return x2_;
		}
		real y2() const
		{
			return y2_;
		}
		real x3() const
		{
			return x3_;
		}
		real y3() const
		{
			return y3_;
		}
		real det() const
		{
			return det_;
		} // the determinant
		real dN1x() const
		{
			return g1x;
		} // dN1x is the x-component of the gradient of
		real dN1y() const
		{
			return g1y;
		} // the shape function which is 1 at node 1
		real dN2x() const
		{
			return g2x;
		} // and 0 at nodes 2 and 3.
		real dN2y() const
		{
			return g2y;
		}
		real dN3x() const
		{
			return g3x;
		}
		real dN3y() const
		{
			return g3y;
		}
		int n1() const
		{
			return n1_;
		} // global node numbers
		int n2() const
		{
			return n2_;
		}
		int n3() const
		{
			return n3_;
		}

		real area() const
		{
			return 0.5 * det_;
		}
// angles are measured in radians, if edge 1 (between node 1 and 2)
// is the longest edge then angle3 is the largest angle
		real maxAngle() const
		{
			return max(max(angle1(), angle2()), angle3());
		}
		real minAngle() const
		{
			return min(min(angle1(), angle2()), angle3());
		}
		real angle(const int i) const;
		real angle1() const
		{
			return acos(-a23() / sqrt(a22() * a33()));
		}
		real angle2() const
		{
			return acos(-a13() / sqrt(a11() * a33()));
		}
		real angle3() const
		{
			return acos(-a12() / sqrt(a11() * a22()));
		}
		real inRadius() const
		{
			return det_ / (l1() + l2() + l3());
		}
		real circumRadius() const
		{
			return 0.5 * d * sqrt(l1sq() * l2sq() * l3sq());
		}
		void inCenter(real& cx, real& cy) const
		{
// weighting with the lengths of the edges opposite to the nodes
			real s1 = l1(), s2 = l2(), s3 = l3();
			real Oinv = 1. / (s1 + s2 + s3);
			cx = (x1_ * s2 + x2_ * s3 + x3_ * s1) * Oinv;
			cy = (y1_ * s2 + y2_ * s3 + y3_ * s1) * Oinv;
		}
		void circumCenter(real& cx, real& cy) const
		{
			real a = a23() * sqr(det_);
			cx = 0.5 * (x2_ + x3_ - a * g1x);
			cy = 0.5 * (y2_ + y3_ - a * g1y);
		}

		void centroid(real& cx, real& cy) const
		{
			cx = centroidx();
			cy = centroidy();
		}
		real centroidx() const
		{
			return (x1_ + x2_ + x3_) * ONETHIRD;
		}
		real centroidy() const
		{
			return (y1_ + y2_ + y3_) * ONETHIRD;
		}
		real midpt1x() const
		{
			return (x1_ + x2_) * 0.5;
		}
		real midpt1y() const
		{
			return (y1_ + y2_) * 0.5;
		}
		real midpt2x() const
		{
			return (x2_ + x3_) * 0.5;
		}
		real midpt2y() const
		{
			return (y2_ + y3_) * 0.5;
		}
		real midpt3x() const
		{
			return (x3_ + x1_) * 0.5;
		}
		real midpt3y() const
		{
			return (y3_ + y1_) * 0.5;
		}

		void normal1(real& nx, real& ny) const
		{
			real a = -1.0 / sqrt(a33());
			nx = a * g3x;
			ny = a * g3y;
		}
		void normal2(real& nx, real& ny) const
		{
			real a = -1.0 / sqrt(a11());
			nx = a * g1x;
			ny = a * g1y;
		}
		void normal3(real& nx, real& ny) const
		{
			real a = -1.0 / sqrt(a22());
			nx = a * g2x;
			ny = a * g2y;
		}
		real h() const
		{
			return det_ / maxLength();
		}  //use det()/l1() if edge 1 is longest
		real maxLength() const
		{
			return sqrt(max(l1sq(), max(l2sq(), l3sq())));
		}
		real minLength() const
		{
			return sqrt(min(l1sq(), min(l2sq(), l3sq())));
		}
		real edgeLength(const int i) const;
		real l1() const
		{
			return sqrt(l1sq());
		}

		real l2() const
		{
			return sqrt(l2sq());
		}
		real l3() const
		{
			return sqrt(l3sq());
		}
		real l1sq() const
		{
			return sqr(x1_ - x2_) + sqr(y1_ - y2_);
		}
		real l2sq() const
		{
			return sqr(x2_ - x3_) + sqr(y2_ - y3_);
		}
		real l3sq() const
		{
			return sqr(x3_ - x1_) + sqr(y3_ - y1_);
		}
		real circumference() const
		{
			return l1() + l2() + l3();
		}
		real perimeter() const
		{
			return 0.5 * circumference();
		}
		real xmax() const
		{
			return max(x1_, max(x2_, x3_));
		}
		real xmin() const
		{
			return min(x1_, min(x2_, x3_));
		}
		real ymax() const
		{
			return max(y1_, max(y2_, y3_));
		}
		real ymin() const
		{
			return min(y1_, min(y2_, y3_));
		}
// some common expressions used for computing elemental stiffnesses
		real a11() const
		{
			return sqr(g1x) + sqr(g1y);
		}
		real a12() const
		{
			return g1x * g2x + g1y * g2y;
		}
		real a21() const
		{
			return a12();
		}
		real a13() const
		{
			return g1x * g3x + g1y * g3y;
		}
		real a31() const
		{
			return a13();
		}
		real a22() const
		{
			return sqr(g2x) + sqr(g2y);
		}
		real a23() const
		{
			return g2x * g3x + g2y * g3y;
		}
		real a32() const
		{
			return a23();
		}
		real a33() const
		{
			return sqr(g3x) + sqr(g3y);
		}
// ==================================================
// For 2D assembling of the matrice (grad phi_i, phi_j)
// for forward seismic equation 
//=====================================================
		real glk_a11() const
		{
			return (g1x + g1y);
		}
		real glk_a22() const
		{
			return (g2x + g2y);
		}
		real glk_a33() const
		{
			return (g3x + g3y);
		}
// ==================================================
// ==================================================
// For 2D assembling of the matrice (grad phi_i^2, phi_j)
// for forward seismic equation 
//=====================================================
		real grad_x_a11() const
		{
			return (g1x);
		}
		real grad_x_a22() const
		{
			return (g2x);
		}
		real grad_x_a33() const
		{
			return (g3x);
		}
		real grad_y_a11() const
		{
			return (g1y);
		}
		real grad_y_a22() const
		{
			return (g2y);
		}
		real grad_y_a33() const
		{
			return (g3y);
		}
		//==================================================
		// For 2D assembling of divergence matrix
		// common expressions to compute (div u, div v), where u = (u1,u2)
		//==================================================
		//  for first deformation component u1
		real div_a11() const
		{
			return g1x * g1x;
		}
		real div_a12() const
		{
			return g1x * g2x;
		}
		real div_a13() const
		{
			return g1x * g3x;
		}
		real div_a21() const
		{
			return g2x * g1x;
		}
		real div_a22() const
		{
			return g2x * g2x;
		}
		real div_a23() const
		{
			return g2x * g3x;
		}
		real div_a31() const
		{
			return g3x * g1x;
		}
		real div_a32() const
		{
			return g3x * g2x;
		}
		real div_a33() const
		{
			return g3x * g3x;
		}
		//=====================================================
		// for second deformation component u2
		real div_b11() const
		{
			return g1y * g1x;
		}
		real div_b12() const
		{
			return g1y * g2x;
		}
		real div_b13() const
		{
			return g1y * g3x;
		}
		real div_b21() const
		{
			return g2y * g1x;
		}
		real div_b22() const
		{
			return g2y * g2x;
		}
		real div_b23() const
		{
			return g2y * g3x;
		}
		real div_b31() const
		{
			return g3y * g1x;
		}
		real div_b32() const
		{
			return g3y * g2x;
		}
		real div_b33() const
		{
			return g3y * g3x;
		}
		//=====================================================
		//  for second equation , first def.component u1
		real div_c11() const
		{
			return g1x * g1y;
		}
		real div_c12() const
		{
			return g1x * g2y;
		}
		real div_c13() const
		{
			return g1x * g3y;
		}
		real div_c21() const
		{
			return g2x * g1y;
		}
		real div_c22() const
		{
			return g2x * g2y;
		}
		real div_c23() const
		{
			return g2x * g3y;
		}
		real div_c31() const
		{
			return g3x * g1y;
		}
		real div_c32() const
		{
			return g3x * g2y;
		}
		real div_c33() const
		{
			return g3x * g3y;
		}
		//=====================================================
		// for second equation, second def.component u2
		real div_d11() const
		{
			return g1y * g1y;
		}
		real div_d12() const
		{
			return g1y * g2y;
		}
		real div_d13() const
		{
			return g1y * g3y;
		}
		real div_d21() const
		{
			return g2y * g1y;
		}
		real div_d22() const
		{
			return g2y * g2y;
		}
		real div_d23() const
		{
			return g2y * g3y;
		}
		real div_d31() const
		{
			return g3y * g1y;
		}
		real div_d32() const
		{
			return g3y * g2y;
		}
		real div_d33() const
		{
			return g3y * g3y;
		}
		//=====================================================================
		// For 2D assembling of matrix (phi_i, rot phi_j) = R2D_1 + R2D_2,
		// where R2D_1 = (phi_i, - d_phi_j/dy) and R2D_2 = (phi_i, d_phi_j/dx)
		// with piecewise-linear test functions, located
		// at the vertices of triangle. Then in local assembling matrix will appear
		// only diagonal terms :
		//  R2D_1 = (- d_phi_i/dy),   R2D_2 = (d_phi_j/dx)
		//=====================================================================
		//  for first electric field component E1
		real R2D_a11() const
		{
			return g1y;
		}
		real R2D_a22() const
		{
			return g2y;
		}
		real R2D_a33() const
		{
			return g3y;
		}
		//=====================================================
		// for second electric field component E2
		real R2D_b11() const
		{
			return -g1x;
		}
		real R2D_b22() const
		{
			return -g2x;
		}
		real R2D_b33() const
		{
			return -g3x;
		}
		//=====================================================
		//=====================================================

		real barycentric1(const real x, const real y) const
		{
			return 1 + g1x * (x - x1_) + g1y * (y - y1_);
		}
		real barycentric2(const real x, const real y) const
		{
			return g2x * (x - x1_) + g2y * (y - y1_);
		}
		real barycentric3(const real x, const real y) const
		{
			return g3x * (x - x1_) + g3y * (y - y1_);
		}
		void barycentric(real x, real y, real& phi1, real& phi2, real& phi3) const
		{
			x -= x1_;
			y -= y1_;
			phi2 = g2x * x + g2y * y;
			phi3 = g3x * x + g3y * y;
			phi1 = 1.0 - phi2 - phi3;
		}
		bool inSide(const real x, const real y) const;     // not yet implemented
		bool outSide(const real x, const real y) const;    // not yet implemented
		bool onTriangle(const real x, const real y) const; // not yet implemented
// some quality measures for triangles normalized so that 1 is best and 0 worst
		real qualityMeasure(const int i) const;
		real etahH() const
		{        // qualityMeasure(1), h()/H() normalized by 2/sqrt(3)
			return 1.154700538379251773 * det_ / max(max(l1sq(), l2sq()), l3sq());
		}
		real etaincirc() const
		{   // qualityMeasure(2), inRadius/circumRadius normalized
			real s1 = l1(), s2 = l2(), s3 = l3();
			return 4 * sqr(det_) / ((s1 + s2 + s3) * s1 * s2 * s3);
		}
		real etaMaxAngle() const
		{  // qualityMeasure(3), normalized by 3/(2*PI)
			return (PI - maxAngle()) * 0.4774648292756859185;
		}
		real etapltmg() const
		{     // qualityMeasure(4)
			  // used in pltmg version 7.0, page 29, normalized by 2*sqrt(3)
			return 3.464101615137754422 * det_ / (l1sq() + l2sq() + l3sq());
		}
		//  void testComputations ();
		bool ok() const; // check that a grid is attached
		void print() const;

		void evalGradient(const MV_Vector<real>& u, real gradu[])
		{
			real u1 = u(n1_), u2 = u(n2_), u3 = u(n3_);
			gradu[0] = u1 * dN1x() + u2 * dN2x() + u3 * dN3x();
			gradu[1] = u1 * dN1y() + u2 * dN2y() + u3 * dN3y();
		}
		real evalFunction(const MV_Vector<real>& u, const real loc[])
		{
			return u(n1_) * loc[0] + u(n2_) * loc[1] + u(n3_) * loc[2];
		}

		real getPntx(const real loc[])
		{
			return x1_ * loc[0] + x2_ * loc[1] + x3_ * loc[2];
		}
		real getPnty(const real loc[])
		{
			return y1_ * loc[0] + y2_ * loc[1] + y3_ * loc[2];
		}
		void getPnt(const real loc[], real xy[])
		{
			xy[0] = x1_ * loc[0] + x2_ * loc[1] + x3_ * loc[2];
			xy[1] = y1_ * loc[0] + y2_ * loc[1] + y3_ * loc[2];
		}

// additional functions for ElmT6n2D elements with the nodes 4, 5 and 6 on
// the midpoints of the edges 1, 2, and 3, respectively.
		bool areQuadNodesAtMidPoints() const; // make refill first
// computation of basis functions and gradients given a local coordinate
// use barycentric to find local coordinates

// WARNING, only correct if nodes at midpoints of edges! 
		void quadShape(const real r, const real s, const real t, real& psi1, real& psi2, real& psi3, real& psi4, real& psi5, real& psi6) const
		{
			psi1 = r * (2.0 * r - 1.0);
			psi2 = s * (2.0 * s - 1.0);
			psi3 = t * (2.0 * t - 1.0);
			psi4 = 4.0 * r * s;
			psi5 = 4.0 * s * t;
			psi6 = 4.0 * r * t;
		}
// gradient of the first three shape functions
		void dNq123(const real r, const real s, const real t, real& psi1x, real& psi1y, real& psi2x, real& psi2y, real& psi3x, real& psi3y) const
		{
			real a = 4 * r - 1;
			psi1x = a * g1x;
			psi1y = a * g1y;
			a = 4 * s - 1;
			psi2x = a * g2x;
			psi2y = a * g2y;
			a = 4 * t - 1;
			psi3x = a * g3x;
			psi3y = a * g3y;
		}
// gradient of the shape functions 4, 5 and 6
		void dNq456(real r, real s, real t, real& psi4x, real& psi4y, real& psi5x, real& psi5y, real& psi6x, real& psi6y) const
		{
			r *= 4;
			s *= 4;
			t *= 4;
			psi4x = r * g2x + s * g1x;
			psi4y = r * g2y + s * g1y;
			psi5x = s * g3x + t * g2x;
			psi5y = s * g3y + t * g2y;
			psi6x = t * g1x + r * g3x;
			psi6y = t * g1y + r * g3y;
		}

		// default/dummy functions:

};
/*>FET3n2D:*/

/*Class:FET3n2D

 NAME:  FET3n2D - finite element for ElmT3n2D for efficient programming

 SYNTAX:     @FET3n2D


 KEYWORDS:

 finite elements, triangle, geometry


 DESCRIPTION:

 The purpose of this class is to collect all the computations involving
 the geometry of a triangle in one place. A major consideration is
 efficiency. This motivates the use of inlining where possible
 and maybe the lack of some natural versions.
 Applications of the class include programming of 
 finite element assembly, mesh generation, interpolation and more.
 After call of "refill" it is possible to use the other functions.
 The class is mainly meant for elements of type "ElmT3n2D", but
 if there are nodes at midpoints of edges it may also be used 
 for applicable functions for "ElmT6n2D" computations.

 CONSTRUCTORS AND INITIALIZATION:

 There is one empty constructor and one where the grid is attached.
 "checkGrid" is used to examine if all elements of the attached grid are of
 correct type.

 MEMBER FUNCTIONS:

 Convention for numbering of nodes, edges and angles in triangle. 
 Notation:  n = nodes, e = edges, a = angles.
 Note the counter clockwise order.

 e2
 n3__________n2
 |a3    a2/   
 |      /
 e3|    / e1
 |a1/
 n1|/

 grid: get the grid from which the element is read.
 refill: load in triangle from grid and compute gradients.
 getElmNo: get the current element number.
 x1,y1, x2,y2, x3,y3: the coordinates of the triangle nodes, they must
 must be in counter clockwise direction
 det: the jacobi determinant (twice the area of triangle)
 dN1x,dN1y, dN2x,dN2y, dN3x,dN3y: gradients of the linear shape functions 
 associated with nodes 1, 2 and 3
 n1,n2,n3: the grid node numbers of the vertices.
 checkGrid: checks that 1) a grid is attached, 
 2) the grid only contains "ElmT3n2D" elements.
 area: area of the triangle.
 maxAngle,minAngle,angle,angle1,angle2,angle3: measured in radians. 
 Note that if edge 1 is the longest edge then angle3 is the largest angle.
 See the ascii figure above for the numbering convention.
 inRadius,inCenter: the radius and center of the inscribed circle.
 circumRadius,circumCenter: the radius and center of the circumscribed circle
 centroid,centroidx,centroidy: centroid of the triangle (the gravity point).
 midpt1x,midpt1y, midpt2x,midpt2y, midpt3x,midpt3y: midpoints of the edges
 1, 2 and 3, see ascii figure above.
 normal: the outer normals of length 1 of edges.
 h: small length scale of triangle defined as h = 2*area/diameter,
 where diameter is the longest edge of the triangle. A geometrical 
 interpretation is that "h" isthe height of the normal to the longest edge.
 maxLength,minLength,edgeLength,l1,l2,l3,l1sq,l2sq,l3sq: compute the 
 lengths of edges. The sq-variants returns the square of the lengths to 
 avoid unnecessary computations of square roots.
 circumference, perimeter: self explaining.
 a11,a12,a13,a21,a22,a23,a31,a32,a33: products of gradients often occurring in
 finite element assembly, also useful for computation of other quantities.
 aij = (dNi,dNj), where (,) is the discrete scalar product.
 localCoordinates: given a point, compute its local coordinates.
 barycentric1,barycentric2,barycentric3: compute the barycentric coordinates
 for a given point, this is the same as the value of the basis functions
 at the point.
 
 inSide: is the point inside the triangle? - not yet implemented
 outSide: is the point outside the triangle? - not yet implemented
 onTriangle: is the point on the boundary? - not yet implemented

 qualityMeasure:  some quality measures for triangles normalized 
 so that the output is 1.0 for an equilateral triangle and tends towards
 0.0 for a triangle of degenerate shape. 
 etahH: ratio between small length scale and diameter.
 etaincirc: ratio between radii of incircle and circumcircle.
 etaMaxAngle: based on the maximum angle.
 etapltmg: quality measure used in pltmg.

 testComputations: make a test that the implementation is correct.

 ok: is a grid attached.
 print: not yet implemented.

 -- functions for "ElmT6n2D" elements 
 Warning, "quadShape","dNq123" and "dNq456" are only correct 
 if nodes at midpoint of edges.

 areQuadNodesAtMidPoints: check that there are nodes in the "ElmT6n2D"
 element at midpoint of edges.
 quadShape: compute value of basis functions at a given barycentric position
 dNq123, dNq456: values of gradients at a given barycentric position

 DEVELOPED BY:   

 SINTEF Applied Mathematics, Oslo, Norway, and
 University of Oslo, Dept. of Mathematics, Norway

 AUTHOR:	        
 
 Klas Samuelsson,  Dept. of Informatics,  University of Oslo

 End:
 */

/*<FET4n3D:*/
class FET4n3D
{
	private:
		int elm_no_;    // element number
		WavesGridB* grid_; // the finite element grid
		// positions of nodes 4 is on the left side of the triangle defined
		// by the nodes 1,2, and 3 in counter clockwise direction
		real x1_, y1_, z1_, x2_, y2_, z2_, x3_, y3_, z3_, x4_, y4_, z4_;
		real det_;                         // determinant = 6*volume
		// gradients of shape functions
		real g1x, g1y, g1z, g2x, g2y, g2z, g3x, g3y, g3z, g4x, g4y, g4z;
		int n1_, n2_, n3_, n4_;                     // node numbers

		real d; // help variables, d = 1.0/det_
		real j11, j12, j13, j21, j22, j23, j31, j32, j33;
	public:
		FET4n3D();
		FET4n3D(const WavesGridB* grid);
		~FET4n3D();
		void attach(const WavesGridB* grid); // attach a grid
		void refill(const int elm_no);  // load with element elm_no
		bool checkGrid() const;      // check that all elements in grid are ElmT4n3D
		bool checkGrid(const int e) const; // check that element e is ElmT4n3D
		WavesGridB& grid() const
		{
			return *grid_;
		}
		int getElmNo() const
		{
			return elm_no_;
		} // get element number

		real x1() const
		{
			return x1_;
		} // coordinates of shape function
		real y1() const
		{
			return y1_;
		}
		real z1() const
		{
			return z1_;
		}
		real x2() const
		{
			return x2_;
		}
		real y2() const
		{
			return y2_;
		}
		real z2() const
		{
			return z2_;
		}
		real x3() const
		{
			return x3_;
		}
		real y3() const
		{
			return y3_;
		}
		real z3() const
		{
			return z3_;
		}
		real x4() const
		{
			return x4_;
		}
		real y4() const
		{
			return y4_;
		}
		real z4() const
		{
			return z4_;
		}
		real det() const
		{
			return det_;
		} // jacobian determinant of mapping
		real dN1x() const
		{
			return g1x;
		} // gradients of shape functions
		real dN1y() const
		{
			return g1y;
		}
		real dN1z() const
		{
			return g1z;
		}
		real dN2x() const
		{
			return g2x;
		}
		real dN2y() const
		{
			return g2y;
		}
		real dN2z() const
		{
			return g2z;
		}
		real dN3x() const
		{
			return g3x;
		}
		real dN3y() const
		{
			return g3y;
		}
		real dN3z() const
		{
			return g3z;
		}
		real dN4x() const
		{
			return g4x;
		}
		real dN4y() const
		{
			return g4y;
		}
		real dN4z() const
		{
			return g4z;
		}
		int n1() const
		{
			return n1_;
		}  // get node numbers
		int n2() const
		{
			return n2_;
		}
		int n3() const
		{
			return n3_;
		}
		int n4() const
		{
			return n4_;
		}

		real volume() const
		{
			return det_ * ONESIXTH;
		} // volume of tetrahedron
		real area(const int i) const; // area of face i
		real area1() const
		{
			return det_ * 0.5 * sqrt(a33());
		} // are of face 1
		real area2() const
		{
			return det_ * 0.5 * sqrt(a11());
		}
		real area3() const
		{
			return det_ * 0.5 * sqrt(a22());
		}
		real area4() const
		{
			return det_ * 0.5 * sqrt(a44());
		}
		real maxArea() const
		{
			return det_ * 0.5 * sqrt(max(max(a33(), a11()), max(a22(), a44())));
		}

// computations of angles of tetrahedron
		real maxSolidAngle() const;
		real minSolidAngle() const;
		void minMaxSolidAngle(real& minsa, real& maxsa) const;
		real solidAngle(const int i) const;
		real solidAngle1() const
		{
			return dihedralAngle1() + dihedralAngle3() + dihedralAngle4() - PI;
		}
		real solidAngle2() const
		{
			return dihedralAngle1() + dihedralAngle2() + dihedralAngle5() - PI;
		}
		real solidAngle3() const
		{
			return dihedralAngle2() + dihedralAngle3() + dihedralAngle6() - PI;
		}
		real solidAngle4() const
		{
			return dihedralAngle4() + dihedralAngle5() + dihedralAngle6() - PI;
		}
		real dihedralAngle(const int i) const; // dihedral angle at edge i
		real dihedralAngle1() const
		{
			return PI - acos(a34() / sqrt(a33() * a44()));
		}
		real dihedralAngle2() const
		{
			return PI - acos(a14() / sqrt(a11() * a44()));
		}
		real dihedralAngle3() const
		{
			return PI - acos(a24() / sqrt(a22() * a44()));
		}
		real dihedralAngle4() const
		{
			return PI - acos(a23() / sqrt(a22() * a33()));
		}
		real dihedralAngle5() const
		{
			return PI - acos(a13() / sqrt(a11() * a33()));
		}
		real dihedralAngle6() const
		{
			return PI - acos(a12() / sqrt(a11() * a22()));
		}
		real inRadius() const
		{
			return 1.0 / (sqrt(a11()) + sqrt(a22()) + sqrt(a33()) + sqrt(a44()));
		}
		void inCenter(real& cx, real& cy, real& cz) const;
		real circumRadius() const;
		void circumCenter(real& cx, real& cy, real& cz) const;
		void centroid(real& cx, real& cy, real& cz) const
		{
			cx = centroidx();
			cy = centroidy();
			cz = centroidz();
		}
		void centroid(real c[]) const
		{
			c[0] = centroidx();
			c[1] = centroidy();
			c[2] = centroidz();
		}
		real centroidx() const
		{
			return (x1_ + x2_ + x3_ + x4_) * 0.25;
		}
		real centroidy() const
		{
			return (y1_ + y2_ + y3_ + y4_) * 0.25;
		}
		real centroidz() const
		{
			return (z1_ + z2_ + z3_ + z4_) * 0.25;
		}
		real midptface1x() const
		{
			return (x1_ + x2_ + x4_) * ONETHIRD;
		}
		real midptface1y() const
		{
			return (y1_ + y2_ + y4_) * ONETHIRD;
		}
		real midptface1z() const
		{
			return (z1_ + z2_ + z4_) * ONETHIRD;
		}
		real midptface2x() const
		{
			return (x2_ + x3_ + x4_) * ONETHIRD;
		}
		real midptface2y() const
		{
			return (y2_ + y3_ + y4_) * ONETHIRD;
		}
		real midptface2z() const
		{
			return (z2_ + z3_ + z4_) * ONETHIRD;
		}
		real midptface3x() const
		{
			return (x1_ + x3_ + x4_) * ONETHIRD;
		}
		real midptface3y() const
		{
			return (y1_ + y3_ + y4_) * ONETHIRD;
		}
		real midptface3z() const
		{
			return (z1_ + z3_ + z4_) * ONETHIRD;
		}
		real midptface4x() const
		{
			return (x1_ + x2_ + x3_) * ONETHIRD;
		}
		real midptface4y() const
		{
			return (y1_ + y2_ + y3_) * ONETHIRD;
		}
		real midptface4z() const
		{
			return (z1_ + z2_ + z3_) * ONETHIRD;
		}
		void midptface1(real mpt[])
		{
			mpt[0] = midptface1x();
			mpt[1] = midptface1y();
			mpt[2] = midptface1z();
		}
		void midptface2(real mpt[])
		{
			mpt[0] = midptface2x();
			mpt[1] = midptface2y();
			mpt[2] = midptface2z();
		}
		void midptface3(real mpt[])
		{
			mpt[0] = midptface3x();
			mpt[1] = midptface3y();
			mpt[2] = midptface3z();
		}
		void midptface4(real mpt[])
		{
			mpt[0] = midptface4x();
			mpt[1] = midptface4y();
			mpt[2] = midptface4z();
		}
		real midpt1x() const
		{
			return (x1_ + x2_) * 0.5;
		}
		real midpt1y() const
		{
			return (y1_ + y2_) * 0.5;
		}
		real midpt1z() const
		{
			return (z1_ + z2_) * 0.5;
		}
		real midpt2x() const
		{
			return (x2_ + x3_) * 0.5;
		}
		real midpt2y() const
		{
			return (y2_ + y3_) * 0.5;
		}
		real midpt2z() const
		{
			return (z2_ + z3_) * 0.5;
		}
		real midpt3x() const
		{
			return (x3_ + x1_) * 0.5;
		}
		real midpt3y() const
		{
			return (y3_ + y1_) * 0.5;
		}
		real midpt3z() const
		{
			return (z3_ + z1_) * 0.5;
		}
		real midpt4x() const
		{
			return (x1_ + x4_) * 0.5;
		}
		real midpt4y() const
		{
			return (y1_ + y4_) * 0.5;
		}
		real midpt4z() const
		{
			return (z1_ + z4_) * 0.5;
		}
		real midpt5x() const
		{
			return (x2_ + x4_) * 0.5;
		}
		real midpt5y() const
		{
			return (y2_ + y4_) * 0.5;
		}
		real midpt5z() const
		{
			return (z2_ + z4_) * 0.5;
		}
		real midpt6x() const
		{
			return (x3_ + x4_) * 0.5;
		}
		real midpt6y() const
		{
			return (y3_ + y4_) * 0.5;
		}
		real midpt6z() const
		{
			return (z3_ + z4_) * 0.5;
		}
// outer normals of length 1
		void normal1(real& nx, real& ny, real& nz) const
		{
			real a = -1 / sqrt(a33());
			nx = a * g3x;
			ny = a * g3y;
			nz = a * g3z;
		}
		void normal2(real& nx, real& ny, real& nz) const
		{
			real a = -1 / sqrt(a11());
			nx = a * g1x;
			ny = a * g1y;
			nz = a * g1z;
		}
		void normal3(real& nx, real& ny, real& nz) const
		{
			real a = -1 / sqrt(a22());
			nx = a * g2x;
			ny = a * g2y;
			nz = a * g2z;
		}
		void normal4(real& nx, real& ny, real& nz) const
		{
			real a = -1 / sqrt(a44());
			nx = a * g4x;
			ny = a * g4y;
			nz = a * g4z;
		}

		void normal1(real n[]) const
		{
			real a = -1 / sqrt(a33());
			n[0] = a * g3x;
			n[1] = a * g3y;
			n[2] = a * g3z;
		}
		void normal2(real n[]) const
		{
			real a = -1 / sqrt(a11());
			n[0] = a * g1x;
			n[1] = a * g1y;
			n[2] = a * g1z;
		}
		void normal3(real n[]) const
		{
			real a = -1 / sqrt(a22());
			n[0] = a * g2x;
			n[1] = a * g2y;
			n[2] = a * g2z;
		}
		void normal4(real n[]) const
		{
			real a = -1 / sqrt(a44());
			n[0] = a * g4x;
			n[1] = a * g4y;
			n[2] = a * g4z;
		}

		real h() const
		{
			return 1.0 / sqrt(max(max(a33(), a11()), max(a22(), a44())));
		}
		real H() const
		{
			return maxLength();
		}
		real maxLength() const
		{
			return sqrt(max(max(max(l1sq(), l2sq()), max(l3sq(), l4sq())), max(l5sq(), l6sq())));
		}
		real minLength() const
		{
			return sqrt(min(min(min(l1sq(), l2sq()), min(l3sq(), l4sq())), min(l5sq(), l6sq())));
		}
		real edgeLength(const int i) const;
		real l1() const
		{
			return sqrt(l1sq());
		}
		real l2() const
		{
			return sqrt(l2sq());
		}
		real l3() const
		{
			return sqrt(l3sq());
		}
		real l4() const
		{
			return sqrt(l4sq());
		}
		real l5() const
		{
			return sqrt(l5sq());
		}
		real l6() const
		{
			return sqrt(l6sq());
		}
		real l1sq() const
		{
			return sqr(j11) + sqr(j12) + sqr(j13);
		}
		real l2sq() const
		{
			return sqr(x2_ - x3_) + sqr(y2_ - y3_) + sqr(z2_ - z3_);
		}
		real l3sq() const
		{
			return sqr(j21) + sqr(j22) + sqr(j23);
		}
		real l4sq() const
		{
			return sqr(j31) + sqr(j32) + sqr(j33);
		}
		real l5sq() const
		{
			return sqr(x2_ - x4_) + sqr(y2_ - y4_) + sqr(z2_ - z4_);
		}
		real l6sq() const
		{
			return sqr(x3_ - x4_) + sqr(y3_ - y4_) + sqr(z3_ - z4_);
		}
		real xmax() const
		{
			return max(max(x1_, x2_), max(x3_, x4_));
		}
		real xmin() const
		{
			return min(min(x1_, x2_), min(x3_, x4_));
		}
		real ymax() const
		{
			return max(max(y1_, y2_), max(y3_, y4_));
		}
		real ymin() const
		{
			return min(min(y1_, y2_), min(y3_, y4_));
		}
		real zmax() const
		{
			return max(max(z1_, z2_), max(z3_, z4_));
		}
		real zmin() const
		{
			return min(min(z1_, z2_), min(z3_, z4_));
		}

		real hface1() const
		{
			return det_ * sqrt(a33() / max(l1sq(), max(l4sq(), l5sq())));
		}
		real hface2() const
		{
			return det_ * sqrt(a11() / max(l2sq(), max(l5sq(), l6sq())));
		}
		real hface3() const
		{
			return det_ * sqrt(a22() / max(l3sq(), max(l4sq(), l6sq())));
		}
		real hface4() const
		{
			return det_ * sqrt(a44() / max(l1sq(), max(l2sq(), l3sq())));
		}
// some common expressions used for computing elemental stiffness
		real a11() const
		{
			return sqr(g1x) + sqr(g1y) + sqr(g1z);
		}
		real a12() const
		{
			return g1x * g2x + g1y * g2y + g1z * g2z;
		}
		real a13() const
		{
			return g1x * g3x + g1y * g3y + g1z * g3z;
		}
		real a21() const
		{
			return a12();
		}
		real a22() const
		{
			return sqr(g2x) + sqr(g2y) + sqr(g2z);
		}
		real a23() const
		{
			return g2x * g3x + g2y * g3y + g2z * g3z;
		}
		real a24() const
		{
			return g2x * g4x + g2y * g4y + g2z * g4z;
		}
		real a31() const
		{
			return a13();
		}
		real a32() const
		{
			return a23();
		}
		real a33() const
		{
			return sqr(g3x) + sqr(g3y) + sqr(g3z);
		}
		real a34() const
		{
			return g3x * g4x + g3y * g4y + g3z * g4z;
		}
		real a14() const
		{
			return g1x * g4x + g1y * g4y + g1z * g4z;
		}
		real a41() const
		{
			return a14();
		}
		real a42() const
		{
			return a24();
		}
		real a43() const
		{
			return a34();
		}
		real a44() const
		{
			return sqr(g4x) + sqr(g4y) + sqr(g4z);
		}

// some common expressions used for computing elemental stiffness
// for Maxwell's equation: case when we write  
// rot rot E = grad div E  - div grad E
// then in assembling matrices remains only following terms
//*************************************************************** 

// for first equation in system, for E1:

		real rot1_E11() const
		{
			return sqr(g1y) + sqr(g1z);
		}
		real rot1_E12() const
		{
			return g1y * g2y + g1z * g2z;
		}
		real rot1_E13() const
		{
			return g1y * g3y + g1z * g3z;
		}
		real rot1_E21() const
		{
			return a12();
		}
		real rot1_E22() const
		{
			return sqr(g2y) + sqr(g2z);
		}
		real rot1_E23() const
		{
			return g2y * g3y + g2z * g3z;
		}
		real rot1_E24() const
		{
			return g2y * g4y + g2z * g4z;
		}
		real rot1_E31() const
		{
			return a13();
		}
		real rot1_E32() const
		{
			return a23();
		}
		real rot1_E33() const
		{
			return sqr(g3y) + sqr(g3z);
		}
		real rot1_E34() const
		{
			return g3y * g4y + g3z * g4z;
		}
		real rot1_E14() const
		{
			return g1y * g4y + g1z * g4z;
		}
		real rot1_E41() const
		{
			return a14();
		}
		real rot1_E42() const
		{
			return a24();
		}
		real rot1_E43() const
		{
			return a34();
		}
		real rot1_E44() const
		{
			return sqr(g4y) + sqr(g4z);
		}

		// for second equation in system, for E2

		real rot2_E11() const
		{
			return sqr(g1x) + sqr(g1z);
		}
		real rot2_E12() const
		{
			return g1x * g2x + g1z * g2z;
		}
		real rot2_E13() const
		{
			return g1x * g3x + g1z * g3z;
		}
		real rot2_E21() const
		{
			return a12();
		}
		real rot2_E22() const
		{
			return sqr(g2x) + sqr(g2z);
		}
		real rot2_E23() const
		{
			return g2x * g3x + g2z * g3z;
		}
		real rot2_E24() const
		{
			return g2x * g4x + g2z * g4z;
		}
		real rot2_E31() const
		{
			return a13();
		}
		real rot2_E32() const
		{
			return a23();
		}
		real rot2_E33() const
		{
			return sqr(g3x) + sqr(g3z);
		}
		real rot2_E34() const
		{
			return g3x * g4x + g3z * g4z;
		}
		real rot2_E14() const
		{
			return g1x * g4x + g1z * g4z;
		}
		real rot2_E41() const
		{
			return a14();
		}
		real rot2_E42() const
		{
			return a24();
		}
		real rot2_E43() const
		{
			return a34();
		}
		real rot2_E44() const
		{
			return sqr(g4x) + sqr(g4z);
		}

		// for third equation, component E3:
		real rot3_E11() const
		{
			return sqr(g1x) + sqr(g1y);
		}
		real rot3_E12() const
		{
			return g1x * g2x + g1y * g2y;
		}
		real rot3_E13() const
		{
			return g1x * g3x + g1y * g3y;
		}
		real rot3_E21() const
		{
			return a12();
		}
		real rot3_E22() const
		{
			return sqr(g2x) + sqr(g2y);
		}
		real rot3_E23() const
		{
			return g2x * g3x + g2y * g3y;
		}
		real rot3_E24() const
		{
			return g2x * g4x + g2y * g4y;
		}
		real rot3_E31() const
		{
			return a13();
		}
		real rot3_E32() const
		{
			return a23();
		}
		real rot3_E33() const
		{
			return sqr(g3x) + sqr(g3y);
		}
		real rot3_E34() const
		{
			return g3x * g4x + g3y * g4y;
		}
		real rot3_E14() const
		{
			return g1x * g4x + g1y * g4y;
		}
		real rot3_E41() const
		{
			return a14();
		}
		real rot3_E42() const
		{
			return a24();
		}
		real rot3_E43() const
		{
			return a34();
		}
		real rot3_E44() const
		{
			return sqr(g4x) + sqr(g4y);
		}

		//=======================================================
		// For 3D assembling of divergence matrix
		// common expressions to compute (div u, div v),
		//  where u = (u1, u2, u3)
		//========================================================
		//  for first deformation component u1 , first equation 
		//________________________________________________________

		real div_a11() const
		{
			return g1x * g1x;
		}
		real div_a12() const
		{
			return g1x * g2x;
		}
		real div_a13() const
		{
			return g1x * g3x;
		}
		real div_a14() const
		{
			return g1x * g4x;
		}

		real div_a21() const
		{
			return g2x * g1x;
		}
		real div_a22() const
		{
			return g2x * g2x;
		}
		real div_a23() const
		{
			return g2x * g3x;
		}
		real div_a24() const
		{
			return g2x * g4x;
		}

		real div_a31() const
		{
			return g3x * g1x;
		}
		real div_a32() const
		{
			return g3x * g2x;
		}
		real div_a33() const
		{
			return g3x * g3x;
		}
		real div_a34() const
		{
			return g3x * g4x;
		}

		real div_a41() const
		{
			return g4x * g1x;
		}
		real div_a42() const
		{
			return g4x * g2x;
		}
		real div_a43() const
		{
			return g4x * g3x;
		}
		real div_a44() const
		{
			return g4x * g4x;
		}

		//===========================================================
		// for second deformation component u2, first equation
		real div_b11() const
		{
			return g1y * g1x;
		}
		real div_b12() const
		{
			return g1y * g2x;
		}
		real div_b13() const
		{
			return g1y * g3x;
		}
		real div_b14() const
		{
			return g1y * g4x;
		}

		real div_b21() const
		{
			return g2y * g1x;
		}
		real div_b22() const
		{
			return g2y * g2x;
		}
		real div_b23() const
		{
			return g2y * g3x;
		}
		real div_b24() const
		{
			return g2y * g4x;
		}

		real div_b31() const
		{
			return g3y * g1x;
		}
		real div_b32() const
		{
			return g3y * g2x;
		}
		real div_b33() const
		{
			return g3y * g3x;
		}
		real div_b34() const
		{
			return g3y * g4x;
		}

		real div_b41() const
		{
			return g4y * g1x;
		}
		real div_b42() const
		{
			return g4y * g2x;
		}
		real div_b43() const
		{
			return g4y * g3x;
		}
		real div_b44() const
		{
			return g4y * g4x;
		}

		//==========================================================
		// for third deformation component u3, first equation
		//==========================================================
		real div_c11() const
		{
			return g1z * g1x;
		}
		real div_c12() const
		{
			return g1z * g2x;
		}
		real div_c13() const
		{
			return g1z * g3x;
		}
		real div_c14() const
		{
			return g1z * g4x;
		}

		real div_c21() const
		{
			return g2z * g1x;
		}
		real div_c22() const
		{
			return g2z * g2x;
		}
		real div_c23() const
		{
			return g2z * g3x;
		}
		real div_c24() const
		{
			return g2z * g4x;
		}

		real div_c31() const
		{
			return g3z * g1x;
		}
		real div_c32() const
		{
			return g3z * g2x;
		}
		real div_c33() const
		{
			return g3z * g3x;
		}
		real div_c34() const
		{
			return g3z * g4x;
		}

		real div_c41() const
		{
			return g4z * g1x;
		}
		real div_c42() const
		{
			return g4z * g2x;
		}
		real div_c43() const
		{
			return g4z * g3x;
		}
		real div_c44() const
		{
			return g4z * g4x;
		}

		//===========================================================
		//  for second equation , first def.component u1
		//___________________________________________________________ 

		real div_d11() const
		{
			return g1x * g1y;
		}
		real div_d12() const
		{
			return g1x * g2y;
		}
		real div_d13() const
		{
			return g1x * g3y;
		}
		real div_d14() const
		{
			return g1x * g4y;
		}

		real div_d21() const
		{
			return g2x * g1y;
		}
		real div_d22() const
		{
			return g2x * g2y;
		}
		real div_d23() const
		{
			return g2x * g3y;
		}
		real div_d24() const
		{
			return g2x * g4y;
		}

		real div_d31() const
		{
			return g3x * g1y;
		}
		real div_d32() const
		{
			return g3x * g2y;
		}
		real div_d33() const
		{
			return g3x * g3y;
		}
		real div_d34() const
		{
			return g3x * g4y;
		}

		real div_d41() const
		{
			return g4x * g1y;
		}
		real div_d42() const
		{
			return g4x * g2y;
		}
		real div_d43() const
		{
			return g4x * g3y;
		}
		real div_d44() const
		{
			return g4x * g4y;
		}

		//============================================================
		// for second equation, second def.component u2
		//____________________________________________________________

		real div_e11() const
		{
			return g1y * g1y;
		}
		real div_e12() const
		{
			return g1y * g2y;
		}
		real div_e13() const
		{
			return g1y * g3y;
		}
		real div_e14() const
		{
			return g1y * g4y;
		}

		real div_e21() const
		{
			return g2y * g1y;
		}
		real div_e22() const
		{
			return g2y * g2y;
		}
		real div_e23() const
		{
			return g2y * g3y;
		}
		real div_e24() const
		{
			return g2y * g4y;
		}

		real div_e31() const
		{
			return g3y * g1y;
		}
		real div_e32() const
		{
			return g3y * g2y;
		}
		real div_e33() const
		{
			return g3y * g3y;
		}
		real div_e34() const
		{
			return g3y * g4y;
		}

		real div_e41() const
		{
			return g4y * g1y;
		}
		real div_e42() const
		{
			return g4y * g2y;
		}
		real div_e43() const
		{
			return g4y * g3y;
		}
		real div_e44() const
		{
			return g4y * g4y;
		}

		//============================================================
		// for second equation, third  def.component u3
		//____________________________________________________________

		real div_f11() const
		{
			return g1z * g1y;
		}
		real div_f12() const
		{
			return g1z * g2y;
		}
		real div_f13() const
		{
			return g1z * g3y;
		}
		real div_f14() const
		{
			return g1z * g4y;
		}

		real div_f21() const
		{
			return g2z * g1y;
		}
		real div_f22() const
		{
			return g2z * g2y;
		}
		real div_f23() const
		{
			return g2z * g3y;
		}
		real div_f24() const
		{
			return g2z * g4y;
		}

		real div_f31() const
		{
			return g3z * g1y;
		}
		real div_f32() const
		{
			return g3z * g2y;
		}
		real div_f33() const
		{
			return g3z * g3y;
		}
		real div_f34() const
		{
			return g3z * g4y;
		}

		real div_f41() const
		{
			return g4z * g1y;
		}
		real div_f42() const
		{
			return g4z * g2y;
		}
		real div_f43() const
		{
			return g4z * g3y;
		}
		real div_f44() const
		{
			return g4z * g4y;
		}

		//============================================================
		// for third  equation, first def.component u1
		//____________________________________________________________

		real div_g11() const
		{
			return g1x * g1z;
		}
		real div_g12() const
		{
			return g1x * g2z;
		}
		real div_g13() const
		{
			return g1x * g3z;
		}
		real div_g14() const
		{
			return g1x * g4z;
		}

		real div_g21() const
		{
			return g2x * g1z;
		}
		real div_g22() const
		{
			return g2x * g2z;
		}
		real div_g23() const
		{
			return g2x * g3z;
		}
		real div_g24() const
		{
			return g2x * g4z;
		}

		real div_g31() const
		{
			return g3x * g1z;
		}
		real div_g32() const
		{
			return g3x * g2z;
		}
		real div_g33() const
		{
			return g3x * g3z;
		}
		real div_g34() const
		{
			return g3x * g4z;
		}

		real div_g41() const
		{
			return g4x * g1z;
		}
		real div_g42() const
		{
			return g4x * g2z;
		}
		real div_g43() const
		{
			return g4x * g3z;
		}
		real div_g44() const
		{
			return g4x * g4z;
		}

		//============================================================
		// for third  equation, second def.component u2
		//____________________________________________________________

		real div_h11() const
		{
			return g1y * g1z;
		}
		real div_h12() const
		{
			return g1y * g2z;
		}
		real div_h13() const
		{
			return g1y * g3z;
		}
		real div_h14() const
		{
			return g1y * g4z;
		}

		real div_h21() const
		{
			return g2y * g1z;
		}
		real div_h22() const
		{
			return g2y * g2z;
		}
		real div_h23() const
		{
			return g2y * g3z;
		}
		real div_h24() const
		{
			return g2y * g4z;
		}

		real div_h31() const
		{
			return g3y * g1z;
		}
		real div_h32() const
		{
			return g3y * g2z;
		}
		real div_h33() const
		{
			return g3y * g3z;
		}
		real div_h34() const
		{
			return g3y * g4z;
		}

		real div_h41() const
		{
			return g4y * g1z;
		}
		real div_h42() const
		{
			return g4y * g2z;
		}
		real div_h43() const
		{
			return g4y * g3z;
		}
		real div_h44() const
		{
			return g4y * g4z;
		}

		//============================================================
		// for third  equation, second def.component u2
		//____________________________________________________________

		real div_i11() const
		{
			return g1z * g1z;
		}
		real div_i12() const
		{
			return g1z * g2z;
		}
		real div_i13() const
		{
			return g1z * g3z;
		}
		real div_i14() const
		{
			return g1z * g4z;
		}

		real div_i21() const
		{
			return g2z * g1z;
		}
		real div_i22() const
		{
			return g2z * g2z;
		}
		real div_i23() const
		{
			return g2z * g3z;
		}
		real div_i24() const
		{
			return g2z * g4z;
		}

		real div_i31() const
		{
			return g3z * g1z;
		}
		real div_i32() const
		{
			return g3z * g2z;
		}
		real div_i33() const
		{
			return g3z * g3z;
		}
		real div_i34() const
		{
			return g3z * g4z;
		}

		real div_i41() const
		{
			return g4z * g1z;
		}
		real div_i42() const
		{
			return g4z * g2z;
		}
		real div_i43() const
		{
			return g4z * g3z;
		}
		real div_i44() const
		{
			return g4z * g4z;
		}

		//=================================================================
		real evalDiv_x(real* u)
		{
			real u1 = u[n1_], u2 = u[n2_], u3 = u[n3_], u4 = u[n4_];
			return u1 * dN1x() + u2 * dN2x() + u3 * dN3x() + u4 * dN4x();
		}

		real evalDiv_y(real* u)
		{
			real u1 = u[n1_], u2 = u[n2_], u3 = u[n3_], u4 = u[n4_];

			return u1 * dN1y() + u2 * dN2y() + u3 * dN3y() + u4 * dN4y();
		}

		real evalDiv_z(real* u)
		{
			real u1 = u[n1_], u2 = u[n2_], u3 = u[n3_], u4 = u[n4_];
			return u1 * dN1z() + u2 * dN2z() + u3 * dN3z() + u4 * dN4z();
		}

		//==================================================================

		real barycentric1(const real x, const real y, const real z) const
		{
			return g1x * (x - x2_) + g1y * (y - y2_) + g1z * (z - z2_);
		}
		real barycentric2(const real x, const real y, const real z) const
		{
			return g2x * (x - x1_) + g2y * (y - y1_) + g2z * (z - z1_);
		}
		real barycentric3(const real x, const real y, const real z) const
		{
			return g3x * (x - x1_) + g3y * (y - y1_) + g3z * (z - z1_);
		}
		real barycentric4(const real x, const real y, const real z) const
		{
			return g4x * (x - x1_) + g4y * (y - y1_) + g4z * (z - z1_);
		}
		void barycentric(real x, real y, real z, real& phi1, real& phi2, real& phi3, real& phi4) const
		{
			x -= x1_;
			y -= y1_;
			z -= z1_;
			phi2 = g2x * x + g2y * y + g2z * z;
			phi3 = g3x * x + g3y * y + g3z * z;
			phi4 = g4x * x + g4y * y + g4z * z;
			phi1 = 1 - phi2 - phi3 - phi4;
		}
		bool inSide(const real x, const real y, const real z) const; // not yet impl
		bool outSide(const real x, const real y, const real z) const; // not yet impl
		bool onTetra(const real x, const real y, const real z) const; // not yet impl
// quality measures for tetrahedrals normalized so that 1 is best and 0 worst
		real qualityMeasure(const int i) const;
		real eta() const;               // qualityMeasure(1)
		real etaincirc() const;         // qualityMeasure(2)
		real etaMinSolidAngle() const;  // qualityMeasure(3)
		real etaMaxSolidAngle() const;  // qualityMeasure(4)
		real etahH() const;             // qualityMeasure(5)
		bool ok() const;  // ia a grid attached?
		void print() const; // not yet implemented

		void evalGradient(const MV_Vector<real>& u, real gradu[])
		{
			real u1 = u(n1_), u2 = u(n2_), u3 = u(n3_), u4 = u(n4_);
			gradu[0] = u1 * dN1x() + u2 * dN2x() + u3 * dN3x() + u4 * dN4x();
			gradu[1] = u1 * dN1y() + u2 * dN2y() + u3 * dN3y() + u4 * dN4y();
			gradu[2] = u1 * dN1z() + u2 * dN2z() + u3 * dN3z() + u4 * dN4z();
		}
		real evalFunction(const MV_Vector<real>& u, const real loc[])
		{
			return u(n1_) * loc[0] + u(n2_) * loc[1] + u(n3_) * loc[2] + u(n4_) * loc[3];
		}
		void evalGradientLoc(const MV_Vector<real>& u, real gradu[])
		{
			real u1 = u(0), u2 = u(1), u3 = u(2), u4 = u(3);
			gradu[0] = u1 * dN1x() + u2 * dN2x() + u3 * dN3x() + u4 * dN4x();
			gradu[1] = u1 * dN1y() + u2 * dN2y() + u3 * dN3y() + u4 * dN4y();
			gradu[2] = u1 * dN1z() + u2 * dN2z() + u3 * dN3z() + u4 * dN4z();
		}
		real evalFunctionLoc(const MV_Vector<real>& u, const real loc[])
		{
			return u(0) * loc[0] + u(1) * loc[1] + u(2) * loc[2] + u(3) * loc[3];
		}

// additional functions for ElmT10n3D elements with the nodes 4-10 on
// the midpoints of the edges 1-6, respectively.
// WARNING, nodes must be on midpoints of edges!
		bool areQuadNodesAtMidPoints() const; // make refill first
// computation of basis functions and gradients given a local coordinate
// use barycentric to find local coordinates
		void quadShape(real r, real s, real t, real u, real& psi1, real& psi2, real& psi3, real& psi4, real& psi5, real& psi6, real& psi7, real& psi8, real& psi9, real& psi10) const
		{
			psi1 = r * (2 * r - 1);
			psi2 = s * (2 * s - 1);
			psi3 = t * (2 * t - 1);
			psi4 = u * (2 * u - 1);
			r *= 2;
			s *= 2;
			t *= 2;
			u *= 2;
			psi5 = r * s;
			psi6 = s * t;
			psi7 = r * t;
			psi8 = r * u;
			psi9 = s * u;
			psi10 = t * u;
		}
		void dNq1234(const real r, const real s, const real t, const real u, real& psi1x, real& psi1y, real& psi1z, real& psi2x, real& psi2y, real& psi2z, real& psi3x, real& psi3y, real& psi3z, real& psi4x, real& psi4y, real& psi4z) const
		{
			real a = 4.0 * r - 1.0;
			psi1x = a * g1x;
			psi1y = a * g1y;
			psi1z = a * g1z;
			a = 4.0 * s - 1.0;
			psi2x = a * g2x;
			psi2y = a * g2y;
			psi2z = a * g2z;
			a = 4.0 * t - 1.0;
			psi3x = a * g3x;
			psi3y = a * g3y;
			psi3z = a * g3z;
			a = 4.0 * u - 1.0;
			psi4x = a * g4x;
			psi4y = a * g4y;
			psi4z = a * g4z;
		}
		void dNq5678910(real r, real s, real t, real u, real& psi5x, real& psi5y, real& psi5z, real& psi6x, real& psi6y, real& psi6z, real& psi7x, real& psi7y, real& psi7z, real& psi8x, real& psi8y, real& psi8z, real& psi9x, real& psi9y, real& psi9z, real& psi10x,
				real& psi10y, real& psi10z) const
		{
			r *= 4.0;
			s *= 4.0;
			t *= 4.0;
			u *= 4.0;
			psi5x = r * g2x + s * g1x;
			psi5y = r * g2y + s * g1y;
			psi5z = r * g2z + s * g1z;
			psi6x = s * g3x + t * g2x;
			psi6y = s * g3y + t * g2y;
			psi6z = s * g3z + t * g2z;
			psi7x = t * g1x + r * g3x;
			psi7y = t * g1y + r * g3y;
			psi7z = t * g1z + r * g3z;
			psi8x = r * g4x + u * g1x;
			psi8y = r * g4y + u * g1y;
			psi8z = r * g4z + u * g1z;
			psi9x = s * g4x + u * g2x;
			psi9y = s * g4y + u * g2y;
			psi9z = s * g4z + u * g2z;
			psi10x = t * g4x + u * g3x;
			psi10y = t * g4y + u * g3y;
			psi10z = t * g4z + u * g3z;
		}

};
/*>FET4n3D:*/

/*Class:FET4n3D

 NAME:  FET4n3D - efficient finite element programming for ElmT4n3D elements

 SYNTAX:     @FET4n3D


 KEYWORDS:

 finite elements, tetrahedron, geometry


 DESCRIPTION:

 The purpose of this class is to collect all the computations involving
 the geometry of a tetrahedron in one place. A major consideration is
 efficiency. This motivates the use of inlining where possible
 and maybe the lack of some natural versions.
 Applications of the class include programming of 
 finite element assembly, mesh generation, interpolation and more.
 After call of "refill" it is possible to use the other functions.
 The class is mainly meant for elements of type "ElmT4n3D", but
 if there are nodes at midpoints of edges it may also be used for
 applicable functions for "ElmT10n3D" computations.

 CONSTRUCTORS AND INITIALIZATION:
 

 MEMBER FUNCTIONS:

 Convention for numbering of nodes, edges, faces and angles in tetrahedras:

 The 4 nodes are ordered so that node 4 is on the the left side of the 
 triangle formed by the nodes 1, 2, and 3 in counter clockwise direction.
 In other words: if your right hand is put in the triangle formed by
 nodes 1, 2 and 3, then your thumb points to the interior of the
 tetrahedra.  The orientation can be tested by making refill and check 
 that "det" or equivalently "volume" is positive.

 edge 1 is between nodes 1 and 2
 edge 2 is between nodes 2 and 3
 edge 3 is between nodes 3 and 1
 edge 4 is between nodes 1 and 4
 edge 5 is between nodes 2 and 4
 edge 6 is between nodes 3 and 4

 face 1 contains nodes 1, 2 and 4
 face 2 contains nodes 2, 3 and 4
 face 3 contains nodes 1, 3 and 4
 face 4 contains nodes 1, 2 and 3
 
 The solid angles angles are numbered after the nodes, and the dihedral
 angles after the edges. A dihedral angle at an edge (wedge) is defined 
 as the angle between the normals of the two faces defining the wedge.
 The solid angle at a vertex is computed as the sum of the dihedrals
 of the edges creating the vertex minus PI. 

 grid: get the grid from which the element is read.
 refill: load in triangle from grid and initialize gradients.
 getElmNo: get the current element number.
 x1,y1,z1, x2,y2,z2, x3,y3,z3, x4,y4,z4: the coordinates of the tetrahedron,
 see description of orientation above.
 det: the jacobi determinant (six times the volume of the tetrahedron)
 dN1x,dN1y,dN1z, dN2x,dN2y,dN2z, dN3x,dN3y,dN3z, dN4x,dN4y,dN4z: gradients 
 of the linear shape functions associated with nodes 1, 2, 3 and 4
 n1,n2,n3,n4: the grid node numbers of the vertices 
 checkGrid: checks that 1) a grid is attached, 
 2) the grid only contains "ElmT4n3D" elements.
 volume: volume of the tetrahedron
 area,area1,area2,area3,area4,maxArea: area of the faces, see above for 
 numbering convention.
 maxSolidAngle,minSolidAngle,minMaxSolidAngle,solidAngle,
 solidAngle1,solidAngle2,solidAngle3,solidAngle4: measured in radians 
 dihedralAngle:  measured in radians
 inRadius,inCenter: the radius and center of the inscribed circle
 circumRadius,circumCenter: the radius  and center of the circumscribed circle
 centroid,centroidx,centroidy,centroidz: centroid of triangle (gravity point)
 midptface1x,midptface1y,midptface1z, etc.: gravity midpoint of faces 1-4
 midpt1x,midpt1y,midpt1z, etc.: midpoints of the edges 1-6
 normal: the outer normals of length 1 of faces, the face numbering is used.
 h: small length scale of triangle defined as volume=maxArea*h/3,
 a geometrical interpretation is, "h" the height of the tetrahedron with 
 the face with largest area as base.
 maxLength,minLength,edgeLength,l1,l1sq, etc.: compute the lengths of edges,
 the sq-variants returns the square of the lengths to avoid unnecessary
 computation of square roots.
 a11,a12, etc.: scalar products of gradients often occurring in
 finite element assembly, also useful for computation of other quantities.
 aij = (dNi,dNj), where (.,.) is the discrete scalar product.
 localCoordinates: given a point compute its local coordinates.
 barycentric,barycentric1, etc.: compute the barycentric coordinates for a 
 point, this is the same as the value of the basis functions at the point.
 
 inSide: is the point inside the triangle?   // not implemented yet
 outSide: is the point outside the triangle? // not implemented yet
 onTriangle: is the point on the boundary?   // not implemented yet

 qualityMeasure:  some quality measures for tetrahedron,  normalized 
 so that the output is 1.0 for an equilateral tetrahedra and tends towards
 0.0 for a degenerate tetrahedron. 
 1) eta: used in geompack.
 2) etaincirc: ratio between radii of incircle and circumcircle
 3) etaMinSolidAngle: based on the minimum solid angle
 4) etaMaxSolidAngle: based on the maximum solid angle
 5) etahH: ratio between small length scale and longest edge

 testComputations: // test correct implementation of functions

 ok: is a grid attached
 print: // not implemented yet

 -- functions for "ElmT10n3D" elements 
 Warning, "quadShape","dNq1234" and "dNq5678910" are only correct 
 if nodes at midpoint of edges.

 areQuadNodesAtMidPoints: check that there are nodes in the "ElmT10n3D"
 element at midpoint of edges.
 quadShape: compute value of basis functions at a given barycentric position.
 dNq1234, dNq5678910: values of gradients at a given barycentric position.

 DEVELOPED BY:   

 SINTEF Applied Mathematics, Oslo, Norway, and
 University of Oslo, Dept. of Mathematics, Norway

 AUTHOR:	        
 
 Klas Samuelsson, Dept. of Informatics, University of Oslo

 End:
 */

#endif

/* LOG HISTORY of this file:

 * $Log: FEcomp.h,v $
 * Revision 
 * Version 
 *
 */

