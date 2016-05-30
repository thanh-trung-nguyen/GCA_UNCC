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

#include "include/wavesOptional.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include "include/wavesFEcomp.h"
#include <sys/times.h>
#include "include/wavesOutputs.h"

#if PETSC_USE_DEBUG
#define CHKERRA(e) if(e){PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,e,0,0);exit(1);}
#else
#define CHKERRA(e); /* obsolete in petsc 2.1 */
#endif 
#define SETERRA(n,p,s)  SETERRQ(n,s)  /* obsolete in petsc 2.1 */

typedef double real;
typedef MV_ColMat<real> Mat_real;
typedef MV_ColMat<int> Mat_int;

//HRN
#ifdef MatType
#undef MatType
#endif

WavesOptional::WavesOptional(WavesSDGeometry &sdg_, Grid &outer_gg_, Grid &gg_, bool& EXCHANGE_, bool& USE_FEM_, bool &USE_FDM_, double &EVALPT_X_, double &EVALPT_Y_, double &EVALPT_Z_, double &maxtime_) :

		outerWithHole(sdg_), outerBoundary1(sdg_), outerBoundary2(sdg_)
{
	sdg = sdg_;
	outer_gg = outer_gg_;
	gg = gg_;
	cout << gg.getNoNodes() << endl;
	EXCHANGE = EXCHANGE_;
	USE_FEM = USE_FEM_;
	USE_FDM = USE_FDM_;
	EVALPT_X = EVALPT_X_;
	EVALPT_Y = EVALPT_Y_;
	EVALPT_Z = EVALPT_Z_;
	maxtime = maxtime_;
	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::WavesOptional(WavesSDGeometry &sdg_, Grid &outer_gg_, Grid &gg_, bool& EXCHANGE_, bool& USE_FEM_, bool &USE_FDM_, double &maxtime_) :

		outerWithHole(sdg_), outerBoundary1(sdg_), outerBoundary2(sdg_)

{
	sdg = sdg_;
	outer_gg = outer_gg_;
	gg = gg_;
	cout << gg.getNoNodes() << endl;
	EXCHANGE = EXCHANGE_;
	USE_FEM = USE_FEM_;
	USE_FDM = USE_FDM_;
	maxtime = maxtime_;
	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::WavesOptional(Grid &gg_, bool& USE_FEM_, bool& EXCHANGE_, bool &USE_FDM_)
{
	gg = gg_;
	cout << gg.getNoNodes() << endl;

	EXCHANGE = EXCHANGE_;
	USE_FEM = USE_FEM_;
	USE_FDM = USE_FDM_;
	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::WavesOptional(Grid &gg_, bool& USE_FEM_, double &EVALPT_X_, double &EVALPT_Y_, double &EVALPT_Z_)
{
	gg = gg_;
	cout << gg.getNoNodes() << endl;
	USE_FEM = USE_FEM_;
	EVALPT_X = EVALPT_X_;
	EVALPT_Y = EVALPT_Y_;
	EVALPT_Z = EVALPT_Z_;
	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::WavesOptional(Grid &gg_)
{
	gg = gg_;
	nsd = gg.getNoSpaceDim();

	cout << nsd << endl;

	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::WavesOptional()
{

	bvalues_initialised = 0;
	nbn = 0;
}

WavesOptional::~WavesOptional()
{
	if (bvalues_initialised)
	{
		delete[] boundindex;
		delete[] bound_values;
	}

}

void WavesOptional::myOpen(ofstream &ofs_)
{

	/* We compose the name of output file as outOut001.imp,
	 using a static counter i */

	static int i = 1;
	char buf[20];

	/* format %03d stands for 3 digits, padded by zeros */
	sprintf(buf, "outOut%03d.inp", i++);
	cout << "File " << buf << " opened." << endl;
	ofs_.open(buf);
}

//=================================================================
// CopyValues should be inserted !!! vfrom util2.h

void WavesOptional::copyValues(double *x, real *y, const int *IX, const int *IY, const int &n)
{

	for (int i = 0; i < n; i++)
	{
		y[IY[i]] = x[IX[i]];
		cout << " y ( " << IY[i] << " ) = " << y[IY[i]] << ", x ( " << IX[i] << " ) = " << x[IX[i]] << endl;

	}
}

// get length of the side with nodes n1 and n2

double WavesOptional::GetLength2D(int n1, int n2)
{

	double n1x = gg.getCoor(n1, 0);
	double n1y = gg.getCoor(n1, 1);

	double n2x = gg.getCoor(n2, 0);
	double n2y = gg.getCoor(n2, 1);

	double proj_x = (n2y - n1y);
	double proj_y = (n2x - n1x);

	double length = sqrt(proj_x * proj_x + proj_y * proj_y);

	cout << " length " << length << endl;
	return length;

}

// get length of the side with nodes n1 and n2

double WavesOptional::GetArea3D(int n1, int n2, int n3)
{

	double x1 = gg.getCoor(n1, 0);
	double y1 = gg.getCoor(n1, 1);
	double z1 = gg.getCoor(n1, 2);

	double x2 = gg.getCoor(n2, 0);
	double y2 = gg.getCoor(n2, 1);
	double z2 = gg.getCoor(n2, 2);

	double x3 = gg.getCoor(n3, 0);
	double y3 = gg.getCoor(n3, 1);
	double z3 = gg.getCoor(n3, 2);

	// compute 1/2-perimeter p
	double a = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1));
	double b = sqrt((x2 - x3) * (x2 - x3) + (y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3));
	double c = sqrt((x3 - x1) * (x3 - x1) + (y3 - y1) * (y3 - y1) + (z3 - z1) * (z3 - z1));
	double p = (a + b + c) / 2.0;

	// compute area using Herona formula
	double area = sqrt(p * (p - a) * (p - b) * (p - c));

	//  cout<<" area "<<area<<endl;
	return area;

}

// get  normal n=(nx,ny)  for DG(1) method in 2D for nodes (n1,n2)
//=================================================================

void WavesOptional::GetNormal2D(int n1, int n2, MV_Vector<double>& normal)
{

	normal.newsize(2);

	double n1x = gg.getCoor(n1, 0);
	double n1y = gg.getCoor(n1, 1);

	double n2x = gg.getCoor(n2, 0);
	double n2y = gg.getCoor(n2, 1);

	normal(0) = -(n2y - n1y);
	normal(1) = (n2x - n1x);

	double n_normal = sqrt(normal(0) * normal(0) + normal(1) * normal(1));

	normal(0) /= n_normal;
	normal(1) /= n_normal;
}

// normal = AB X AC, where AB and AC are vectors of the common face
//  with nodes A,B,C which represents global nodes numbers in the mesh, 
// normal have coordinates  (nx,ny,nz)
void WavesOptional::GetNormal3D(int A, int B, int C, MV_Vector<double>& normal)
{

	double Ax = gg.getCoor(A, 0);
	double Ay = gg.getCoor(A, 1);
	double Az = gg.getCoor(A, 2);

	double Bx = gg.getCoor(B, 0);
	double By = gg.getCoor(B, 1);
	double Bz = gg.getCoor(B, 2);

	double Cx = gg.getCoor(C, 0);
	double Cy = gg.getCoor(C, 1);
	double Cz = gg.getCoor(C, 2);

	normal(0) = (Ay - By) * (Az - Cz) - (Ay - Cy) * (Az - Bz);
	normal(1) = (Ax - Bx) * (Az - Cz) - (Ax - Cx) * (Az - Bz);
	normal(2) = (Ax - Bx) * (Ay - Cy) - (Ax - Cx) * (Ay - By);

	double n_normal = sqrt(normal(0) * normal(0) + normal(1) * normal(1) + normal(2) * normal(2));

	normal(0) /= n_normal;
	normal(1) /= n_normal;
	normal(2) /= n_normal;

}
//============ normal direction in 2D  ==================================
// for a given global node numbers of the element in mesh -  n1, n2 , n3 - find normal to side n1,n2 and 
// direction to this normal
// Direction is bool variable and have value=false by definition
// Direction normal_minus is direction to neighbour element ne, or K^-
// and  direction normal_plus is direction to element e or K^+. 
// if normal_minus = true, then normal is to K^- or neighb. elm-t
// else if   normal_plus = true, then normal is to e or K^+.
//========================================================================
void WavesOptional::NormalDir2D(int n1, int n2, int n3, MV_Vector<double>& normal, int& normal_minus, int& normal_plus)
{
	normal_minus = 0;
	normal_plus = 0;

	double n2x = gg.getCoor(n2, 0);
	double n2y = gg.getCoor(n2, 1);

	double n3x = gg.getCoor(n3, 0);
	double n3y = gg.getCoor(n3, 1);

	// compute scalar product normal* n2n3
	double Product = normal(0) * (n2x - n3x) + normal(1) * (n2y - n3y);

	if (Product > 0)
		normal_minus = 1;
	else
		normal_plus = 1;
}
//============ normal direction in 3D  ==================================
// for given global node numbers of element n1, n2 , n3 and n4 find normal to face n1,n2, n3 and 
// direction to this normal
// if normal*n1n4 > 0, then normal = n_minus, else normal = n_plus
//========================================================================
void WavesOptional::NormalDir3D(int n1, int n2, int n3, int n4, MV_Vector<double>& normal, int& normal_minus, int& normal_plus)
{
	normal_minus = 0;
	normal_plus = 0;

	double n4x = gg.getCoor(n4, 0);
	double n4y = gg.getCoor(n4, 1);
	double n4z = gg.getCoor(n4, 2);

	double n1x = gg.getCoor(n1, 0);
	double n1y = gg.getCoor(n1, 1);
	double n1z = gg.getCoor(n1, 2);

	// compute scalar product normal* n1n4
	double Product = normal(1) * (n1x - n4x) + normal(2) * (n1y - n4y) + normal(3) * (n1z - n4z);

	if (Product > 0)
		normal_minus = 1;
	else
		normal_plus = 1;

}

//=========================================================================

void WavesOptional::sort(int *ia, const int &n)
{

	// *Simple* sort

	int swapDone;
	int tmp;

	do
	{
		swapDone = 0;
		for (int i = 1; i < n; i++)
			if (ia[i] < ia[i - 1])
			{
				swapDone = 1;
				tmp = ia[i];
				ia[i] = ia[i - 1];
				ia[i - 1] = tmp;
			}
	}
	while (swapDone);
}

//==================================================================
const int DEBUG = 1;

int *WavesOptional::maskaFDM(MV_Vector<int>& bndNodes)
{
	int i, j, k;

	nsd = sdg.getNoSpaceDim();
	n_i = sdg.getN_i();
	n_j = sdg.getN_j();
	n_k = sdg.getN_k();
	double eps = 1e-6;
	int numb_nod;
	double coord_x, coord_y, coord_z, c_x, c_y, c_z;
	int sch = 0;
	int numb_ones = 0;

	for (i = 0; i < bndNodes.size(); i++)
		if (bndNodes(i) == 1)
			numb_ones++;

	cout << "number of the code 1: " << numb_ones << endl;

	numb_nod = sdg.getNoNodes();
	int *maska_;
	maska_ = new int[numb_nod];
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					maska_[i + n_i * (j + n_j * k)] = 3;
				}
			}
		}
	}
	else
	{
		for (j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				maska_[i + n_i * j] = 3;
			}
		}
	}

	cout << " in maskaFDM: all nodes have code 3" << endl;
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					c_x = sdg.getCoor(i + n_i * (j + n_j * k), 0);
					c_y = sdg.getCoor(i + n_i * (j + n_j * k), 1);
					c_z = sdg.getCoor(i + n_i * (j + n_j * k), 2);
					for (int n = 0; n < bndNodes.size(); n++)
					{
						coord_x = gg.getCoor(n, 0);
						coord_y = gg.getCoor(n, 1);
						coord_z = gg.getCoor(n, 2);

						if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps) && (fabs(coord_z - c_z) < eps))
						{
							maska_[i + n_i * (j + n_j * k)] = bndNodes(n);
						}
					}

				}
			}
		}
	}
	else
	{
		for (int j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				c_x = sdg.getCoor(i + n_i * j, 0);
				c_y = sdg.getCoor(i + n_i * j, 1);
				for (int n = 0; n < bndNodes.size(); n++)
				{
					if (bndNodes(n) == 1 || bndNodes(n) == 2)
					{
						coord_x = gg.getCoor(n, 0);
						coord_y = gg.getCoor(n, 1);

						if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps))
						{
							maska_[i + n_i * j] = bndNodes(n);
							// cout<<" maska "<< maska_[i+n_i*j]<<endl;
						}

					}
				}
			}
		}
	}

	cout << "in maskaFDM: are assigned codes 3,2,1 ; starting to assign code 0" << endl;

	//for innermost nodes code 0
	if (nsd == 3)
	{
		for (k = 1; k < n_k; k++)
		{
			for (j = 1; j < n_j; j++)
			{
				for (i = 1; i < n_i; i++)
				{
					if ((maska_[i + n_i * (j + n_j * k)] == 3 && (maska_[i + n_i * ((j - 1) + n_j * k)] == 1 || maska_[(i - 1) + n_i * (j + n_j * k)] == 1 || maska_[i + n_i * (j + n_j * (k - 1))] == 1))
							|| (maska_[i + n_i * (j + n_j * k)] == 3 && (maska_[i + n_i * ((j - 1) + n_j * k)] == 0 || maska_[(i - 1) + n_i * (j + n_j * k)] == 0 || maska_[i + n_i * (j + n_j * (k - 1))] == 0)))
					{
						maska_[i + n_i * (j + n_j * k)] = 0;
					}

				}
			}
		}
	}
	else
	{
		for (int j = 1; j < n_j; j++)
		{
			for (int i = 1; i < n_i; i++)
			{
				if ((maska_[i + n_i * j] == 3 && (maska_[i + n_i * (j - 1)] == 1 || maska_[(i - 1) + n_i * j] == 1 || maska_[(i - 1) + n_i * (j - 1)] == 1))
						|| (maska_[i + n_i * j] == 3 && (maska_[i + n_i * (j - 1)] == 0 || maska_[(i - 1) + n_i * j] == 0 || maska_[(i - 1) + n_i * (j - 1)] == 0)))
				{
					maska_[i + n_i * j] = 0;
				}

			}
		}
	}

	ofstream outp;
	outp.open("maskaFDM.dat");
	if (outp.fail())
	{
		perror("maskaFDM.dat");
		return 0;
	}
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					outp << "maska(" << i << "," << j << "," << k << ")=" << maska_[i + n_i * (j + n_j * k)] << "\n";
				}
			}
		}
	}
	else
	{
		int sch = 0;
		for (j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				sch++;
				outp << sch << "  maska(" << i << "," << j << ")=" << maska_[i + n_i * j] << "\n";
			}
		}
		outp << sch << "\n";
	}
	cerr << "Maska[0]=" << maska_[0] << endl; // debugger confused at this point
	return maska_;
}

// ============== in the case of structured meshes ==========================0

int *WavesOptional::maskaStructFDM(MV_Vector<int>& bndNodes)
{
	int i, j, k;

	nsd = sdg.getNoSpaceDim();
	n_i = sdg.getN_i();
	n_j = sdg.getN_j();
	n_k = sdg.getN_k();
	double eps = 1e-6;
	int numb_nod;
	double coord_x, coord_y, coord_z, c_x, c_y, c_z;
	int sch = 0;
	int numb_ones = 0;

	for (i = 0; i < bndNodes.size(); i++)
		if (bndNodes(i) == 1)
			numb_ones++;

	cout << "number of the code 1: " << numb_ones << endl;

	numb_nod = sdg.getNoNodes();
	int *maska_;
	maska_ = new int[numb_nod];
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					maska_[i + n_i * (j + n_j * k)] = 3;
				}
			}
		}
	}
	else
	{
		for (j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				maska_[i + n_i * j] = 3;
			}
		}
	}

	cout << " in maskaFDM: all nodes have code 3" << endl;
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					c_x = sdg.getCoor(i + n_i * (j + n_j * k), 0);
					c_y = sdg.getCoor(i + n_i * (j + n_j * k), 1);
					c_z = sdg.getCoor(i + n_i * (j + n_j * k), 2);
					for (int n = 0; n < bndNodes.size(); n++)
					{
						coord_x = gg.getCoor(n, 0);
						coord_y = gg.getCoor(n, 1);
						coord_z = gg.getCoor(n, 2);

						if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps) && (fabs(coord_z - c_z) < eps))
						{
							maska_[i + n_i * (j + n_j * k)] = bndNodes(n);
						}
					}

				}
			}
		}
	}
	else
	{
		for (int j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				c_x = sdg.getCoor(i + n_i * j, 0);
				c_y = sdg.getCoor(i + n_i * j, 1);
				for (int n = 0; n < bndNodes.size(); n++)
				{
					coord_x = gg.getCoor(n, 0);
					coord_y = gg.getCoor(n, 1);

					if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps))
					{
						maska_[i + n_i * j] = bndNodes(n);
						// cout<<" maska "<< maska_[i+n_i*j]<<endl;
					}

				}
			}
		}
	}

	cout << "in maskaFDM: are assigned codes 3,2,1 ; starting to assign code 0" << endl;

	ofstream outp;
	outp.open("maskaFDM.dat");
	if (outp.fail())
	{
		perror("maskaFDM.dat");
		return 0;
	}
	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					outp << "maska(" << i << "," << j << "," << k << ")=" << maska_[i + n_i * (j + n_j * k)] << "\n";
				}
			}
		}
	}
	else
	{
		int sch = 0;
		for (j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				sch++;
				outp << sch << "  maska(" << i << "," << j << ")=" << maska_[i + n_i * j] << "\n";
			}
		}
		outp << sch << "\n";
	}
	cerr << "Maska[0]=" << maska_[0] << endl; // debugger confused at this point
	return maska_;

}

//============================================================================
void WavesOptional::getFileName(char *buf)
{
	static int i = 1;
	sprintf(buf, "outInn%03d.inp", ++i);
}

real *WavesOptional::fromFEMtoFDM(Vec& u2)
{
	int i, j, k;

	nsd = sdg.getNoSpaceDim();
	n_i = sdg.getN_i();
	n_j = sdg.getN_j();
	n_k = sdg.getN_k();
	PetscScalar* uu2;
	ierr = VecGetArray(u2,&uu2);
	CHKERRA(ierr);
	double eps = 1e-6;
	int numb_nod;
	double coord_x, coord_y, coord_z, c_x, c_y, c_z;

	numb_nod = sdg.getNoNodes();
	real *solfdm = new real[numb_nod];
	int sch = 0;

	if (nsd == 3)
	{
		for (k = 0; k < n_k; k++)
		{
			for (j = 0; j < n_j; j++)
			{
				for (i = 0; i < n_i; i++)
				{
					c_x = sdg.getCoor(i + n_i * (j + n_j * k), 0);
					c_y = sdg.getCoor(i + n_i * (j + n_j * k), 1);
					c_z = sdg.getCoor(i + n_i * (j + n_j * k), 2);

					for (int n = 0; n < gg.getNoNodes(); n++)
					{
						coord_x = gg.getCoor(n, 0);
						coord_y = gg.getCoor(n, 1);
						coord_z = gg.getCoor(n, 2);

						if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps) && (fabs(coord_z - c_z) < eps))
						{
							solfdm[i + n_i * (j + n_j * k)] = uu2[n];

						}
					}

				}
			}
		}
	}
	else
	{
		for (k = 0; k < n_k; k++)
		{
			for (int j = 0; j < n_j; j++)
			{
				for (int i = 0; i < n_i; i++)
				{
					c_x = sdg.getCoor(i + n_i * (j + n_j * k), 0);
					c_y = sdg.getCoor(i + n_i * (j + n_j * k), 1);
					for (int n = 0; n < gg.getNoNodes(); n++)
					{
						coord_x = gg.getCoor(n, 0);
						coord_y = gg.getCoor(n, 1);

						if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps))
						{
							solfdm[i + n_i * (j + n_j * k)] = uu2[n];

						}
					}
				}
			}
		}
	}

	ierr = VecRestoreArray(u2,&uu2);
	CHKERRA(ierr);
	return solfdm;
}

//===============================================================
//     make array with codes from REFINED fdm MESH to coarse 
//===============================================================
void WavesOptional::extractNodesFDM(WavesSDGeometry& ref, WavesSDGeometry& coarse, MV_Vector<int>& NodesFDM)
{
	int i, n;

	nsd = ref.getNoSpaceDim();

	double eps = 1e-6;
	int numb_nod;
	double coord_x, coord_y, coord_z, c_x, c_y, c_z;

	numb_nod = coarse.getNoNodes();
	NodesFDM.newsize(numb_nod);

	int sch = 0;

	if (nsd == 3)
	{
		for (i = 0; i < ref.getNoNodes(); i++)
		{
			c_x = ref.getCoor(i, 0);
			c_y = ref.getCoor(i, 1);
			c_z = ref.getCoor(i, 2);

			for (n = 0; n < numb_nod; n++)
			{
				coord_x = coarse.getCoor(n, 0);
				coord_y = coarse.getCoor(n, 1);
				coord_z = coarse.getCoor(n, 2);

				if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps) && (fabs(coord_z - c_z) < eps))
				{
					NodesFDM[n] = i;
					// NodesFDM[i+n_i*(j + n_j*k)] = 1.0;
					n = numb_nod; // to break n-loop
				}
			}

		}
	}
	else
	{

		for (i = 0; i < ref.getNoNodes(); i++)
		{
			c_x = ref.getCoor(i, 0);
			c_y = ref.getCoor(i, 1);
			for (n = 0; n < numb_nod; n++)
			{
				coord_x = coarse.getCoor(n, 0);
				coord_y = coarse.getCoor(n, 1);

				if ((fabs(coord_x - c_x) < eps) && (fabs(coord_y - c_y) < eps))
				{
					NodesFDM[n] = i;
					// NodesFDM[i+n_i*(j + n_j*k)] = 1.0;
					n = numb_nod; // to break n-loop
				}
			}
		}
	}

}
//==============================================================
// perform Laplace transform
//===============================================================

void WavesOptional::LaplaceTransform(WavesGridB& gg, Mat_real& array, double omega, double dt, MV_Vector<real>& LaplaceTr)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int k;

	double time, time_, sol = 0.0;

	int nr_timesteps = array.size(1);

	for (n = 0; n < nno; n++)
	{
		for (k = 0; k < nr_timesteps - 1; k++)
		{
			time = k * dt;
			time_ = (k + 1) * dt;

			LaplaceTr(n) += (array(n, k) / omega) * (-exp(-omega * time_) + exp(-omega * time));
		}
	}

}

//===============================================================

void WavesOptional::LaplaceTransform_(WavesGridB& gg, Mat_real& array, double omega, double dt, MV_Vector<real>& LaplaceTr)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int k;

	double time, time_, sol = 0.0;

	int nr_timesteps = array.size(1);

	for (n = 0; n < nno; n++)
	{
		for (k = 0; k < nr_timesteps - 1; k++)
		{
			time = k * dt;
			time_ = (k + 1) * dt;

			LaplaceTr(n) += (array(n, k) * exp(-omega * time) + array(n, k + 1) * exp(-omega * time_)) * 0.5 * dt;
			//cout << "Laplace_tr " << LaplaceTr(n) << "  function " << array(n, k) << endl;
		}
	}

}

//============================= added by Thanh
// each row is the time-domain data at a point. NOTE: data at t =0 is not provided. It must be zero!!!
void WavesOptional::LaplaceTransform(Mat_real time_array, double dt, double omega, MV_Vector<real>& LaplaceTr)
{
	int nr_timesteps = time_array.size(1);
	int nr_spacepoints = time_array.size(0);
	int n, k; 
	double time, value;

	for (n = 0; n < nr_spacepoints; n++)
	  {
	    value = 0;
		for (k = 0; k < nr_timesteps - 1; k++)
		{
			time = (k+1) * dt;
			value  += time_array(n, k) * exp(-omega*time)* dt;
		}
		LaplaceTr(n) = value;		
	}

}

// laplace transform with data: each column is the time-domain data at a point
void WavesOptional::LaplaceTransform_col(Mat_real time_array, double dt, double omega, MV_Vector<real>& LaplaceTr)
{
	int nr_timesteps = time_array.size(0);
	int nr_spacepoints = time_array.size(1);
	int n, k; 
	double time, value;

	for (n = 0; n < nr_spacepoints; n++)
	  {
	    value = 0;
		for (k = 0; k < nr_timesteps - 1; k++)
		{
			time = (k + 1) * dt;
		        value  += time_array(k,n) * exp(-omega*time)* dt;
		}
		LaplaceTr(n) = value;		
	}

}

// laplace transform, input data is a file, e.g. ExactFEM.m, NOTE: the data at t= 0 is zero!!! and not provided in the data file
void WavesOptional::LaplaceTransform(char* datafile, int NoTimeSteps, double dt, double s, MV_Vector<real>& LaplaceTr)
{
	int nr_spacepoints = LaplaceTr.size();
	int n, k; 
	double time, value, dat;

	ifstream inp;
	inp.open(datafile);
	
	LaplaceTr = 0.0; 
	for (k = 0; k < NoTimeSteps - 1; k++)
	{
		time = (k + 1) * dt;
		value = exp(-s*time)* dt;
		for (n = 0; n < nr_spacepoints; n++)
		{	
			inp >> dat;
			LaplaceTr(n) += dat*value;		
		}
	}
	inp.close();
}

// laplace transform, input data is a file, e.g. ExactFEM.m, NOTE: the data at t= 0 is zero!!! and not provided in the data file
void WavesOptional::LaplaceTransform(char* datafile, int NoTimeSteps, double dt, MV_Vector<double> s, Mat_real& LaplaceTr)
{
	int nr_spacepoints = LaplaceTr.size(1);
	int Ns = s.size();
	int n, k, ns; 
	double time, dat;

	ifstream inp;
	inp.open(datafile);

	LaplaceTr = 0.0; 
	for (k = 0; k < NoTimeSteps - 1; k++)
	{
		time = (k + 1) * dt;
		for (n = 0; n < nr_spacepoints; n++)
		{	
			inp >> dat;
			for (ns=0; ns<Ns; ns++)
				LaplaceTr(ns,n) += dat*exp(-s(ns)*time)*dt;		
		}
	}
	inp.close();
}


//==============================================================
//  compute space int. in finite element mesh
//==============================================================

double WavesOptional::Compute_space_int(WavesGridB& gg, MV_Vector<real>& difference, int rank)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	double sol = 0.0;
	double val1 = 0.0;
	double val2 = 0.0;
	double val3 = 0.0;
	double val4 = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();

				val1 = difference(n_1) * difference(n_1);
				val2 = difference(n_2) * difference(n_2);
				val3 = difference(n_3) * difference(n_3);

				sol += ((val1 + val2 + val3) / 3.0) * volume;

				if (n_1 == 1019)
					cout << "  val1 " << val1 << "  val2 " << val2 << "  val3 " << val3 << "  volume " << volume << " val1 + val2 + val3 " << val1 + val2 + val3 << " sol " << sol << endl;
			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				n_4 = F.n4();

				volume *= 0.25;

				val1 = difference(n_1) * difference(n_1);
				val2 = difference(n_2) * difference(n_2);
				val3 = difference(n_3) * difference(n_3);
				val4 = difference(n_4) * difference(n_4);

				sol += (val1 + val2 + val3 + val4) * volume;

			}
		}
		sol = sqrt(sol);

	}

	cout << "sol =" << sol << endl;
	return (sol);

}

//==============================================================
//  compute space int. in finite element mesh
//===============================================================

double WavesOptional::Compute_space_int(MV_Vector<real>& Values_at_Elms)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	double sol = 0.0;
	double val1 = 0.0;
	double val2 = 0.0;
	double val3 = 0.0;
	double val4 = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();

				val1 = Values_at_Elms(el) * Values_at_Elms(el);

				sol += (val1 / 3.0) * volume;

			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());

				volume *= 0.25;

				val1 = Values_at_Elms(el) * Values_at_Elms(el);

				sol += val1 * volume;

			}
		}
	}

	sol = sqrt(sol);
	return (sol);

}


void WavesOptional::getElmsAtFaces(const WavesNeighborFE& n2e, const int n1, const int n2, const int n3, const int n4, const int e, int& e3, int& e4)

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

	int n1stop = n2e.nodeIrow(n1 + 1);
	int n2stop = n2e.nodeIrow(n2 + 1);
	int n3stop = n2e.nodeIrow(n3 + 1);
	int n4stop = n2e.nodeIrow(n4 + 1);
	int i = n2e.nodeIrow(n1); // start for first element n1 is part of
	int j = n2e.nodeIrow(n2);
	int k = n2e.nodeIrow(n3);
	int l = n2e.nodeIrow(n4);
	e3 = e4 = -1;
	for (; i < n1stop; i++)
	{
		int v = n2e.nodeJcol(i);
		if (v == e)
			continue; // already know this common element
		while (j < n2stop && v > n2e.nodeJcol(j))
			j++;
		if (j >= n2stop)
			return; // no more elms to node 2
		if (v == n2e.nodeJcol(j))
		{ // both elements equal
			while (k < n3stop && v > n2e.nodeJcol(k))
				k++;
			if (k < n3stop && v == n2e.nodeJcol(k))
			{
				e3 = v;
				if (e4 > -1)
					return;
			}
			while (l < n4stop && v > n2e.nodeJcol(l))
				l++;
			if (l < n4stop && v == n2e.nodeJcol(l))
			{
				e4 = v;
				if (e3 > -1)
					return;
			}
		}
	}
}

int WavesOptional::getElmAtFace(const WavesNeighborFE& n2e, const int n1, const int n2, const int n3, const int e)

{
// Given an element number e and the three global node numbers on the corners
// of an element face, compute the element number of the element 
// of the other side of the face.
// It is assumed that node-to-element (n2e) info for the grid is computed
// and that the elements are sorted in increasing order
// (which is the case for the present implementation of WavesNeighborFE).
// 0 is returned if no other element is found, i.e., is the face
// is a boundary face.

	int n1stop = n2e.nodeIrow(n1 + 1);
	int n2stop = n2e.nodeIrow(n2 + 1);
	int n3stop = n2e.nodeIrow(n3 + 1);
	int i = n2e.nodeIrow(n1);
	int j = n2e.nodeIrow(n2);
	int k = n2e.nodeIrow(n3);
	for (; i < n1stop; i++)
	{
		int v = n2e.nodeJcol(i);
		if (v == e)
			continue; // already know this common element
		while (j < n2stop && v > n2e.nodeJcol(j))
			j++;
		if (j >= n2stop)
			return -1; // outside range of elms to node 2
		if (v == n2e.nodeJcol(j))
		{ // two elements equal
			while (k < n3stop && v > n2e.nodeJcol(k))
				k++;
			if (k >= n3stop)
				return -1; // outside range of elms to node 3
			if (v == n2e.nodeJcol(k))
				return v; // if the last elm also equal, bingo!
		}
	}
	return -1;
}


int WavesOptional::getFaceNo(const WavesGridB& grid, const int n1, const int n2, const int n3, const int e, bool safe)

// Get the face number of element e which consists of the nodes n1, n2 and n3
// if safe is true the function will return 0 if no match is found.
{
	if (safe)
	{
// a version which is safe in the case that n1,n2,n3 are not nodes of element e
		int m1 = grid.loc2glob(e, 0), m2 = grid.loc2glob(e, 1);
		int m3 = grid.loc2glob(e, 2), m4 = grid.loc2glob(e, 3);
		int f1 = 0, f2 = 0, f3 = 0, f4 = 0;

		if (n1 == m1)
		{
			f1++;
			f3++;
			f4++;
		} // m1 is a part of face 1 3 and 4
		else if (n1 == m2)
		{
			f1++;
			f2++;
			f4++;
		} // m2 is a part of face 1 2 and 4
		else if (n1 == m3)
		{
			f2++;
			f3++;
			f4++;
		} // m3 is a part of face 2 3 and 4
		else if (n1 == m4)
		{
			f1++;
			f2++;
			f3++;
		} // m4 is a part of face 1 2 and 3

		if (n2 == m1)
		{
			f1++;
			f3++;
			f4++;
		}
		else if (n2 == m2)
		{
			f1++;
			f2++;
			f4++;
		}
		else if (n2 == m3)
		{
			f2++;
			f3++;
			f4++;
		}
		else if (n2 == m4)
		{
			f1++;
			f2++;
			f3++;
		}

		if (n3 == m1)
		{
			f1++;
			f3++;
			f4++;
		}
		else if (n3 == m2)
		{
			f1++;
			f2++;
			f4++;
		}
		else if (n3 == m3)
		{
			f2++;
			f3++;
			f4++;
		}
		else if (n3 == m4)
		{
			f1++;
			f2++;
			f3++;
		}

		if (f1 == 3)
			return 0; // if 3 hits for a face then correct face is found
		if (f2 == 3)
			return 1;
		if (f3 == 3)
			return 2;
		if (f4 == 3)
			return 3;

		return -1; // did not find any match
	}
	else
	{
// This version is not safe if n1,n2,n3 are not in element e
// but is is faster than the safe version ...
		int n = grid.loc2glob(e, 0);
		if (n != n1 && n != n2 && n != n3)
			return 1;
		n = grid.loc2glob(e, 1);
		if (n != n1 && n != n2 && n != n3)
			return 2;
		n = grid.loc2glob(e, 2);
		if (n != n1 && n != n2 && n != n3)
			return 0;
		return 3;
	}
}

void WavesOptional::findBoundaryNodes(WavesGridB& grid, MV_Vector<int>& boundaryMarkers)
{

// Initiate the elm_neigh structure so that it contains the numbers of
// the neighbors elements on the other side of the element faces.

	int n, e, n1, n2, n3, n4, ne;
	nsd = grid.getNoSpaceDim();

	if (nsd == 2)
	{ // For a triangulation of linear triangle element
	  // find the nodes which lies on the outer boundary
	  // The boundary nodes are marked with 1 the interior has the value 0

		int n, e, ne, j, n1, n2, k;
		int nno = grid.getNoNodes();
		int nel = grid.getNoElms();
		MV_Vector<int> countn2e(nno, 0);
		MV_Vector<int> n2ejcol(nel * 3); // each element is added to exactly 3 nodes
		MV_Vector<int> n2eirow(nno + 1);

		for (e = 0; e < nel; e++)
		{ // count the number of elements a node is in
			countn2e(grid.loc2glob(e, 0))++;countn2e
			(grid.loc2glob(e, 1))++;countn2e
			(grid.loc2glob(e, 2))++;}

int 		prev = n2eirow(0) = 0;
		for (n = 0; n < nno; n++)
			prev = n2eirow(n + 1) = countn2e(n) + prev;
		countn2e = 0;
		for (e = 0; e < nel; e++)
		{
			ne = grid.loc2glob(e, 0);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 1);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 2);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
		}

// it is faster to go backwards due to the fact that the elements in n2e
// are ordered starting with lowest element number and increases.
// Therefore the neighbor element will be found earlier.

		MV_ColMat<int> elmNeighbor(nel, 3);
		elmNeighbor = -1;

		for (e = nel - 1; e >= 0; e--)
			for (j = 0; j < 3; j++)
			{
				if (elmNeighbor(e, j) == -1)
				{ // if not yet set
					n1 = grid.loc2glob(e, j); // first node
					n2 = grid.loc2glob(e, j == 2 ? 0 : j + 1);
					int kstop = n2eirow(n1 + 1);
					for (k = n2eirow(n1); k < kstop; k++)
					{ // search for neighbor element
						ne = n2ejcol(k);
						if (ne != e)
						{ // find second common node
							if (grid.loc2glob(ne, 0) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 0) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 1) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 1) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 2) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 2) = e;
								k = kstop; // to break the k-loop
							}
						}
					}
				}
			}

		boundaryMarkers.newsize(nno);
		boundaryMarkers = 0;

		for (e = 0; e < nel; e++)
		{
			if (elmNeighbor(e, 0) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
			}
			if (elmNeighbor(e, 1) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
			}
			if (elmNeighbor(e, 2) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
			}
		}
	}
	else if (nsd == 3)
	{
		int n, e, n1, n2, n3, n4, ne;
// Given face of element e, find the neighbor elements which have 
// three nodes in common with the face. If no such element is found this
// is a boundary face.

		//  additional_info.initNeighbor(*this,true,false,false);
		WavesNeighborFE& neighbor = grid.getNeighbor();
		if (neighbor.nodeSize() == 0)
			neighbor.init(grid, true, false, false);

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3
		int nel = grid.getNoElms();
		int nno = grid.getNoNodes();

		MV_ColMat<int> elmNeighbor(nel, 4);
		elmNeighbor = -1;

		bool safe = true;
		int ne1, ne2;
		for (e = nel - 1; e >= 0; e--)
		{ // more efficient with decreasing elements
			n1 = grid.loc2glob(e, 0);
			n2 = grid.loc2glob(e, 1);
			n3 = grid.loc2glob(e, 2);
			n4 = grid.loc2glob(e, 3);
			if (elmNeighbor(e, 0) == -1)
				if (elmNeighbor(e, 3) == -1)
				{ // case face 1 and 4 not done
					getElmsAtFaces(neighbor, n1, n2, n4, n3, e, ne1, ne2);
					if (ne1 >= 0)
					{
						n = getFaceNo(grid, n1, n2, n4, ne1, safe);
						elmNeighbor(e, 0) = ne1;
						elmNeighbor(ne1, n) = e;
					}
					if (ne2 >= 0)
					{
						n = getFaceNo(grid, n1, n2, n3, ne2, safe);
						elmNeighbor(e, 3) = ne2;
						elmNeighbor(ne2, n) = e;
					}
				} // to case face 1 and 4
				else
				{
					ne = getElmAtFace(neighbor, n1, n2, n4, e); // case face 1
					if (ne >= 0)
					{
						n = getFaceNo(grid, n1, n2, n4, ne, safe);
						elmNeighbor(e, 0) = ne;
						elmNeighbor(ne, n) = e;
					}
				}
			else if (elmNeighbor(e, 3) == -1) // to  if( elmNeighbor(e,0)==-1 )
			{
				ne = getElmAtFace(neighbor, n1, n2, n3, e); // case face 4
				if (ne >= 0)
				{
					n = getFaceNo(grid, n1, n2, n3, ne, safe);
					elmNeighbor(e, 3) = ne;
					elmNeighbor(ne, n) = e;
				}
			}

			if (elmNeighbor(e, 1) == -1)
				if (elmNeighbor(e, 2) == -1)
				{
					getElmsAtFaces(neighbor, n3, n4, n2, n1, e, ne1, ne2);
					if (ne1 >= 0)
					{
						n = getFaceNo(grid, n3, n4, n2, ne1, safe);
						elmNeighbor(e, 1) = ne1;
						elmNeighbor(ne1, n) = e;
					}
					if (ne2 >= 0)
					{
						n = getFaceNo(grid, n3, n4, n1, ne2, safe);
						elmNeighbor(e, 2) = ne2;
						elmNeighbor(ne2, n) = e;
					}
				}
				else
				{
					ne = getElmAtFace(neighbor, n2, n3, n4, e);
					if (ne >= 0)
					{
						n = getFaceNo(grid, n2, n3, n4, ne, safe);
						elmNeighbor(e, 1) = ne;
						elmNeighbor(ne, n) = e;
					}
				}
			else if (elmNeighbor(e, 2) == -1)
			{
				ne = getElmAtFace(neighbor, n1, n3, n4, e);
				if (ne >= 0)
				{
					n = getFaceNo(grid, n1, n3, n4, ne, safe);
					elmNeighbor(e, 2) = ne;
					elmNeighbor(ne, n) = e;
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

		for (e = 0; e < nel; e++)
		{
			if (elmNeighbor(e, 0) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
			}
			if (elmNeighbor(e, 1) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
			}
			if (elmNeighbor(e, 2) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
			}
			if (elmNeighbor(e, 3) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
			}
		}
	}
}

//=============== make array with neighbour elements in 3D
//===============================================================
void WavesOptional::findNeighborElms3D(WavesGridB& grid, Mat_int& elmNeighbor)
{

	nsd = grid.getNoSpaceDim();

	int n, e, n1, n2, n3, n4, ne;
// Given face of element e, find the neighbor elements which have 
// three nodes in common with the face. If no such element is found this
// is a boundary face.

	//  additional_info.initNeighbor(*this,true,false,false);
	WavesNeighborFE& neighbor = grid.getNeighbor();
	if (neighbor.nodeSize() == 0)
		neighbor.init(grid, true, false, false);

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3
	int nel = grid.getNoElms();
	int nno = grid.getNoNodes();

	// MV_ColMat<int> elmNeighbor(nel,4);
	elmNeighbor.newsize(nel, 4);
	elmNeighbor = -1;

	bool safe = true;
	int ne1, ne2;
	for (e = nel - 1; e >= 0; e--)
	{ // more efficient with decreasing elements
		n1 = grid.loc2glob(e, 0);
		n2 = grid.loc2glob(e, 1);
		n3 = grid.loc2glob(e, 2);
		n4 = grid.loc2glob(e, 3);
		if (elmNeighbor(e, 0) == -1)
			if (elmNeighbor(e, 3) == -1)
			{ // case face 1 and 4 not done
				getElmsAtFaces(neighbor, n1, n2, n4, n3, e, ne1, ne2);
				if (ne1 >= 0)
				{
					n = getFaceNo(grid, n1, n2, n4, ne1, safe);
					elmNeighbor(e, 0) = ne1;
					elmNeighbor(ne1, n) = e;
				}
				if (ne2 >= 0)
				{
					n = getFaceNo(grid, n1, n2, n3, ne2, safe);
					elmNeighbor(e, 3) = ne2;
					elmNeighbor(ne2, n) = e;
				}
			}
			else
			{
				ne = getElmAtFace(neighbor, n1, n2, n4, e); // case face 1
				if (ne >= 0)
				{
					n = getFaceNo(grid, n1, n2, n4, ne, safe);
					elmNeighbor(e, 0) = ne;
					elmNeighbor(ne, n) = e;
				}
			}
		else if (elmNeighbor(e, 3) == -1)
		{
			ne = getElmAtFace(neighbor, n1, n2, n3, e); // case face 4
			if (ne >= 0)
			{
				n = getFaceNo(grid, n1, n2, n3, ne, safe);
				elmNeighbor(e, 3) = ne;
				elmNeighbor(ne, n) = e;
			}
		}

		if (elmNeighbor(e, 1) == -1)
			if (elmNeighbor(e, 2) == -1)
			{
				getElmsAtFaces(neighbor, n3, n4, n2, n1, e, ne1, ne2);
				if (ne1 >= 0)
				{
					n = getFaceNo(grid, n3, n4, n2, ne1, safe);
					elmNeighbor(e, 1) = ne1;
					elmNeighbor(ne1, n) = e;
				}
				if (ne2 >= 0)
				{
					n = getFaceNo(grid, n3, n4, n1, ne2, safe);
					elmNeighbor(e, 2) = ne2;
					elmNeighbor(ne2, n) = e;
				}
			}
			else
			{
				ne = getElmAtFace(neighbor, n2, n3, n4, e);
				if (ne >= 0)
				{
					n = getFaceNo(grid, n2, n3, n4, ne, safe);
					elmNeighbor(e, 1) = ne;
					elmNeighbor(ne, n) = e;
				}
			}
		else if (elmNeighbor(e, 2) == -1)
		{
			ne = getElmAtFace(neighbor, n1, n3, n4, e);
			if (ne >= 0)
			{
				n = getFaceNo(grid, n1, n3, n4, ne, safe);
				elmNeighbor(e, 2) = ne;
				elmNeighbor(ne, n) = e;
			}
		}
	}

	for (e = 0; e < nel; e++)
		cout << "elmNeighbor(" << e << ") = " << elmNeighbor(e, 0) << "  " << elmNeighbor(e, 1) << "  " << elmNeighbor(e, 2) << "    " << elmNeighbor(e, 3) << endl;

}

// make array El2Nodes with common nodes for neighbours elements if they share common face in 3D
// as also make array El2NotCommonNode for neighbours elements - place not common nodes with neighbour element into array 

void WavesOptional::findElmsCommonNodes3D(WavesGridB& grid, Mat_int& El2Nodes, Mat_int& El2NotCommonNode)
{

	//  additional_info.initNeighbor(*this,true,false,false);
	WavesNeighborFE& neighbor = grid.getNeighbor();
	if (neighbor.nodeSize() == 0)
		neighbor.init(grid, true, false, false);

	int n, e, sch, ne, i, j, n1, n2, n3, n4, nn1, nn2, ne1, ne2, ne3, ne4;
	int nno = grid.getNoNodes();
	int nel = grid.getNoElms();

// Array El2Nodes: in the first column are elements, eatch element repeats 4 times , second column is neighbour to the element from the first column, can be -1 if the are no neigbours, can be maximum 4 neighb. to 1 element in 3D. Third, fourth and fifth columns are global node numbers to the shared face of the tetrahedron in 3D.

	El2Nodes.newsize(nel * 4, 5);
	El2Nodes = -1;
	// El2NotCommonNode array in 3D for neighbours elements and node,
// which not belongs to the  two neighbours elements. We need such array to compute jumps in the solution and direction of the normal to the side. Tested on the conforming grids.

	El2NotCommonNode.newsize(nel * 4, 3);
	El2NotCommonNode = -1;

	for (e = 0; e < nel; e++)
	{

		// count_ne2e(e)++;
		n1 = grid.loc2glob(e, 0); // fix first node
		n2 = grid.loc2glob(e, 1);
		n3 = grid.loc2glob(e, 2);
		n4 = grid.loc2glob(e, 3);

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3

		ne1 = getElmAtFace(neighbor, n1, n2, n4, e);
		ne2 = getElmAtFace(neighbor, n2, n3, n4, e);
		ne3 = getElmAtFace(neighbor, n1, n3, n4, e);
		ne4 = getElmAtFace(neighbor, n1, n2, n3, e);

		El2Nodes(4 * e, 0) = e;
		El2Nodes(4 * e, 1) = ne1;
		El2Nodes(4 * e, 2) = n1;
		El2Nodes(4 * e, 3) = n2;
		El2Nodes(4 * e, 4) = n4;

		El2NotCommonNode(4 * e, 0) = e;
		El2NotCommonNode(4 * e, 1) = ne1;
		El2NotCommonNode(4 * e, 2) = n3;

		El2Nodes(4 * e + 1, 0) = e;
		El2Nodes(4 * e + 1, 1) = ne2;
		El2Nodes(4 * e + 1, 2) = n2;
		El2Nodes(4 * e + 1, 3) = n3;
		El2Nodes(4 * e + 1, 4) = n4;

		El2NotCommonNode(4 * e + 1, 0) = e;
		El2NotCommonNode(4 * e + 1, 1) = ne2;
		El2NotCommonNode(4 * e + 1, 2) = n1;

		El2Nodes(4 * e + 2, 0) = e;
		El2Nodes(4 * e + 2, 1) = ne3;
		El2Nodes(4 * e + 2, 2) = n1;
		El2Nodes(4 * e + 2, 3) = n3;
		El2Nodes(4 * e + 2, 4) = n4;

		El2NotCommonNode(4 * e + 2, 0) = e;
		El2NotCommonNode(4 * e + 2, 1) = ne3;
		El2NotCommonNode(4 * e + 2, 2) = n2;

		El2Nodes(4 * e + 3, 0) = e;
		El2Nodes(4 * e + 3, 1) = ne4;
		El2Nodes(4 * e + 3, 2) = n1;
		El2Nodes(4 * e + 3, 3) = n2;
		El2Nodes(4 * e + 3, 4) = n3;

		El2NotCommonNode(4 * e + 3, 0) = e;
		El2NotCommonNode(4 * e + 3, 1) = ne4;
		El2NotCommonNode(4 * e + 3, 2) = n4;

	} // for e

	cout << " El2Nodes array " << endl;

	for (e = 0; e < 4 * nel; e++)
		cout << El2Nodes(e, 0) << "  " << El2Nodes(e, 1) << "   " << El2Nodes(e, 2) << "     " << El2Nodes(e, 3) << "   " << El2Nodes(e, 4) << " not common node " << El2NotCommonNode(e, 2) << endl;

}

//==================================================================
// ====== make array with neighbor elements in 2D. 
//        Neighbor elements are elements which have common nodes with
//        other element in the mesh.
//==================================================================
void WavesOptional::findNeighborElms(WavesGridB& grid, Mat_int& elmNeighbor)
{

	int n, e, n1, n2, n3, n4, ne;
	nsd = grid.getNoSpaceDim();

	if (nsd == 2)
	{ // For a triangulation of linear triangle element
	  // find the nodes which lies on the outer boundary
	  // The boundary nodes are marked with 1 the interior has the value 0

		int n, e, ne, j, n1, n2, k;
		int nno = grid.getNoNodes();
		int nel = grid.getNoElms();
		MV_Vector<int> countn2e(nno, 0);
		MV_Vector<int> n2ejcol(nel * 3); // each element is added to exactly 3 nodes
		MV_Vector<int> n2eirow(nno + 1);

		elmNeighbor.newsize(nel, 3);
		elmNeighbor = -1;

		for (e = 0; e < nel; e++)
		{ // count the number of elements a node is in
			countn2e(grid.loc2glob(e, 0))++;countn2e
			(grid.loc2glob(e, 1))++;countn2e
			(grid.loc2glob(e, 2))++;}

int 		prev = n2eirow(0) = 0;
		for (n = 0; n < nno; n++)
			prev = n2eirow(n + 1) = countn2e(n) + prev;
		countn2e = 0;
		for (e = 0; e < nel; e++)
		{
			ne = grid.loc2glob(e, 0);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 1);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 2);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
		}

// it is faster to go backwards due to the fact that the elements in n2e
// are ordered starting with lowest element number and increases.
// Therefore the neighbor element will be found earlier.

		for (e = nel - 1; e >= 0; e--)
			for (j = 0; j < 3; j++)
			{
				if (elmNeighbor(e, j) == -1)
				{ // if not yet set
					n1 = grid.loc2glob(e, j); // first node
					n2 = grid.loc2glob(e, j == 2 ? 0 : j + 1);
					int kstop = n2eirow(n1 + 1);
					for (k = n2eirow(n1); k < kstop; k++)
					{ // search for neighbor element
						ne = n2ejcol(k);
						if (ne != e)
						{ // find second common node
							if (grid.loc2glob(ne, 0) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 0) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 1) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 1) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 2) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 2) = e;
								k = kstop; // to break the k-loop
							}
						}
					}
				}
			}
	}
}

// make array El2Nodes with common nodes for neighbours elements if they share common side in 2D
// as also make array El2NotCommonNode for neighbours elements - place not common nodes with neighbour element into array 

void WavesOptional::findElmsCommonNodes(WavesGridB& grid, Mat_int& elmNeighbor, Mat_int& El2Nodes, Mat_int& El2NotCommonNode)
{

	int n, e, ne, i, j, n1, n2, n3, nn1, nn2;
	int nno = grid.getNoNodes();
	int nel = grid.getNoElms();

	MV_Vector<int> count_ne2e(nel);
	count_ne2e = 0;

// Array El2Nodes: in the first column are elements, eatch element repeats 3 times , second column is neighbour to the element in the first column, can be -1 if the are no neigbours, can be maximum 3 neighb. to elemnt in 2D. Third and fourth columns are global node numbers to the shared face ( side of the triangle in 2D).

	El2Nodes.newsize(nel * 3, 4);
	El2Nodes = -1;
	// El2NotCommonNode array in 2D for neighbours elements and node,
// which not belongs to the  two neighbours elements. We need such array to compute jumps in the solution and direction of the normal to the side. Tested on the conforming grids.

	El2NotCommonNode.newsize(nel * 3, 3);
	El2NotCommonNode = -1;

	for (e = 0; e < nel; e++)
	{

		for (j = 0; j < 3; j++) // to eatch element are maximum 3 neighbours elements
		{
			ne = elmNeighbor(e, j); // take neighbour element

			El2Nodes(3 * e + j, 0) = e;
			El2Nodes(3 * e + j, 1) = ne;

			El2NotCommonNode(3 * e + j, 0) = e;
			El2NotCommonNode(3 * e + j, 1) = ne;

			if (ne > -1) //if exists such neighbor element
			{
				count_ne2e(e)++;

				for
(				i=0; i < 3; i++) // for 3 nodes in triangle
				{
					n1 = grid.loc2glob( e,i); // fix first node
					n2 = grid.loc2glob( e, i==2 ? 0 : i+1 );

					if ( ne != e )
					{ // if e and ne are not the same 
						if ( grid.loc2glob (ne, 0) == n1 )
						{
							if ( El2Nodes(3*e+j,2) == -1)
							{
								El2Nodes(3*e+j,2) = n1;} // place them in array
							else
							{
								El2Nodes(3*e+j,3) = n1;}

						}
						else if ( grid.loc2glob (ne, 1) == n1 )
						{
							if ( El2Nodes(3*e+j,2) == -1)
							{
								El2Nodes(3*e+j,2) = n1;} // place them in array
							else
							{
								El2Nodes(3*e+j,3) = n1;}

						}
						else if ( grid.loc2glob (ne, 2) == n1 )
						{
							if ( El2Nodes(3*e+j,2) == -1)
							{
								El2Nodes(3*e+j,2) = n1;} // place them in array
							else
							{
								El2Nodes(3*e+j,3) = n1;}

						}
						else
						{
// this node is not common
							if ( El2NotCommonNode(3*e+j,2) == -1)
							{
								El2NotCommonNode(3*e+j,2) = n1;} // place them in array

						}
					} // if e and ne are not the same
				} // for i
			} //if ne != -1
		} // for j
	} // for e

// numerate also sides of the elements which not have a neighbours
	int l, k, ne1, ne2, ne2_1, ne2_2;

	for (e = 0; e < nel; e++)
	{
		if (count_ne2e(e) == 2)
		{ // if there are 2 neighb. to element e     
			for (i = 0; i < 3; i++)
			{ // for 3 nodes in triangle
				k = (i == 2 ? 0 : i + 1);
				l = (i == 1 ? 0 : k + 1);
				ne1 = elmNeighbor(e, i); // take neighbour element
				ne2 = elmNeighbor(e, k);
				if (ne1 != -1 && ne2 != -1)
				{
					n1 = El2Nodes(3 * e + i, 2); // take numbers of the side
					n2 = El2Nodes(3 * e + i, 3); // for element ne1
					n3 = El2NotCommonNode(3 * e + i, 2);

					ne2_1 = El2Nodes(3 * e + k, 2); // and element ne2
					ne2_2 = El2Nodes(3 * e + k, 3);

					if ((n1 != ne2_1 && n3 != ne2_2) || (n1 != ne2_2 && n3 != ne2_1))
					{
						El2Nodes(3 * e + l, 2) = n1;
						El2Nodes(3 * e + l, 3) = n3;
						El2NotCommonNode(3 * e + l, 2) = n2;
					}
					else if ((n2 != ne2_1 && n3 != ne2_2) || (n2 != ne2_2 && n3 != ne2_1))
					{
						El2Nodes(3 * e + l, 2) = n2;
						El2Nodes(3 * e + l, 3) = n3;
						El2NotCommonNode(3 * e + l, 2) = n1;
					}
				} // for  2 neighbours: ne1 ne2
			} // for 3 nodes in triangle
		} // for 2 neighb. elements in count_ne2e
		else if (count_ne2e(e) == 1)
		{ // if there is 1 neighb. to element e     
			for (i = 0; i < 3; i++)
			{ // for 3 nodes in triangle
				k = (i == 2 ? 0 : i + 1);
				l = (i == 1 ? 0 : k + 1);
				ne1 = elmNeighbor(e, i); // take neighbour element

				if (ne1 != -1)
				{ // if ne1 is neighbor element
					n1 = El2Nodes(3 * e + i, 2); // take numbers of the side
					n2 = El2Nodes(3 * e + i, 3); // for element ne1
					n3 = El2NotCommonNode(3 * e + i, 2);

// then n1n3 & n2n3 will be in array El2Nodes
					El2Nodes(3 * e + l, 2) = n1;
					El2Nodes(3 * e + l, 3) = n3;
					El2NotCommonNode(3 * e + l, 2) = n2;

					El2Nodes(3 * e + k, 2) = n2;
					El2Nodes(3 * e + k, 3) = n3;
					El2NotCommonNode(3 * e + k, 2) = n1;
				}
			} // for 3 nodes in triangle
		} // for 1 neighb. to element e
	} // for e
}

//===========================================================================
// make array NormalArray2D with coordinates of the  normals n^+ =(nx^+,ny^+)
// to elements  K^+. Array NormalArray2D have size (3*nel,2),
// where in the first column are values n_x and in the second n_y.
// Values (n_x,n_y) corresponds to the  element e and it's neighbour ne 
//(see  first and second columns in  array El2Nodes), and computes using
// third and fourth columns of the array El2Nodes, where are common nodes
// of the element e and it's neighbour element.  
//=========================================================================

void WavesOptional::makeNormalArray2D(Mat_real& NormalArray2D)
{

	int n1, n2, n3, i, e, nel;

	nel = gg.getNoElms();

	Mat_int elmNeighbor;
	Mat_int El2Nodes;
	Mat_int El2NotCommonNode;

	MV_Vector<double> Normal_e(2);
	Normal_e = 0.0;

	NormalArray2D.newsize(3 * nel, 2);

	int normal_minus, normal_plus;

	findNeighborElms(gg, elmNeighbor);
	findElmsCommonNodes(gg, elmNeighbor, El2Nodes, El2NotCommonNode);
	for (e = 0; e < nel; e++)
	{
		for (i = 0; i < 3; i++)
		{
			n1 = El2Nodes(3 * e + i, 2);
			n2 = El2Nodes(3 * e + i, 3);
			n3 = El2NotCommonNode(3 * e + i, 2);
			GetNormal2D(n1, n2, Normal_e);
			NormalDir2D(n1, n2, n3, Normal_e, normal_minus, normal_plus);

			if (normal_plus == 1)
			{
				NormalArray2D(3 * e + i, 0) = Normal_e(0);
				NormalArray2D(3 * e + i, 1) = Normal_e(1);
			}
			else if (normal_minus == 1)
			{
				NormalArray2D(3 * e + i, 0) = -Normal_e(0);
				NormalArray2D(3 * e + i, 1) = -Normal_e(1);
			}
		} // for i
	} // for e
}

//===========================================================================
// make array NormalArray3D with coordinates of the  normals n^+ =(nx^+,ny^+,nz^+)
// to elements  K^+. Array NormalArray3D have size (4*nel,3),
// where in the first column are values n_x, in the second n_y and 
// in the third - n_z.
// Values (n_x,n_y,n_z) corresponds to the  element e and it's neighbour ne 
//(see  first and second columns in  array El2Nodes), and computes using
// third and fourth columns of the array El2Nodes, where are common nodes
// of the element e and it's neighbour element.  
//=========================================================================

void WavesOptional::makeNormalArray3D(Mat_real& NormalArray3D)
{
	int n1, n2, n3, n4, i, e, nel;
	nel = gg.getNoElms();

	Mat_int elmNeighbor;
	Mat_int El2Nodes;
	Mat_int El2NotCommonNode;

	MV_Vector<double> Normal_e(3);
	Normal_e = 0.0;
	NormalArray3D.newsize(4 * nel, 3);

	int normal_minus, normal_plus;

	findElmsCommonNodes3D(gg, El2Nodes, El2NotCommonNode);
	for (e = 0; e < nel; e++)
	{
		for (i = 0; i < 4; i++)
		{
			n1 = El2Nodes(4 * e + i, 2);
			n2 = El2Nodes(4 * e + i, 3);
			n3 = El2Nodes(4 * e + i, 4);

			n4 = El2NotCommonNode(4 * e + i, 2);

			GetNormal3D(n1, n2, n3, Normal_e);

			NormalDir3D(n1, n2, n3, n4, Normal_e, normal_minus, normal_plus);

			if (normal_plus == 1)
			{
				NormalArray3D(4 * e + i, 0) = Normal_e(0);
				NormalArray3D(4 * e + i, 1) = Normal_e(1);
				NormalArray3D(4 * e + i, 2) = Normal_e(2);

			}
			else if (normal_minus == 1)
			{
				NormalArray3D(4 * e + i, 0) = -Normal_e(0);
				NormalArray3D(4 * e + i, 1) = -Normal_e(1);
				NormalArray3D(4 * e + i, 2) = -Normal_e(2);

			}
		} // for i
	} // for e
}

//==================================================================

void WavesOptional::findB_Elms_Nodes(WavesGridB& grid, MV_Vector<int>& boundaryElms, MV_Vector<int>& boundaryMarkers)
{

// Initiate the elm_neigh structure so that it contains the numbers of
// the neighbors elements on the other side of the element faces.

	int n, e, n1, n2, n3, n4, ne;
	nsd = grid.getNoSpaceDim();

	if (nsd == 2)
	{ // For a triangulation of linear triangle element
	  // find the nodes which lies on the outer boundary
	  // The boundary nodes are marked with 1 the interior has the value 0

		int n, e, ne, j, n1, n2, k;
		int nno = grid.getNoNodes();
		int nel = grid.getNoElms();
		MV_Vector<int> countn2e(nno, 0);
		MV_Vector<int> n2ejcol(nel * 3); // each element is added to exactly 3 nodes
		MV_Vector<int> n2eirow(nno + 1);

		for (e = 0; e < nel; e++)
		{ // count the number of elements a node is in
			countn2e(grid.loc2glob(e, 0))++;countn2e
			(grid.loc2glob(e, 1))++;countn2e
			(grid.loc2glob(e, 2))++;}

int 		prev = n2eirow(0) = 0;
		for (n = 0; n < nno; n++)
			prev = n2eirow(n + 1) = countn2e(n) + prev;
		countn2e = 0;
		for (e = 0; e < nel; e++)
		{
			ne = grid.loc2glob(e, 0);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 1);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
			ne = grid.loc2glob(e, 2);
			n2ejcol(n2eirow(ne) + countn2e(ne)++) = e;
		}

// it is faster to go backwards due to the fact that the elements in n2e
// are ordered starting with lowest element number and increases.
// Therefore the neighbor element will be found earlier.

		MV_ColMat<int> elmNeighbor(nel, 3);
		elmNeighbor = -1;

		for (e = nel - 1; e >= 0; e--)
			for (j = 0; j < 3; j++)
			{
				if (elmNeighbor(e, j) == -1)
				{ // if not yet set
					n1 = grid.loc2glob(e, j); // first node
					n2 = grid.loc2glob(e, j == 2 ? 0 : j + 1);
					int kstop = n2eirow(n1 + 1);
					for (k = n2eirow(n1); k < kstop; k++)
					{ // search for neighbor element
						ne = n2ejcol(k);
						if (ne != e)
						{ // find second common node
							if (grid.loc2glob(ne, 0) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 0) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 1) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 1) = e;
								k = kstop; // to break the k-loop
							}
							else if (grid.loc2glob(ne, 2) == n2)
							{
								elmNeighbor(e, j) = ne;
								elmNeighbor(ne, 2) = e;
								k = kstop; // to break the k-loop
							}
						}
					}
				}
			}

		boundaryMarkers.newsize(nno);
		boundaryMarkers = 0;

		for (e = 0; e < nel; e++)
		{
			if (elmNeighbor(e, 0) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
			}
			if (elmNeighbor(e, 1) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
			}
			if (elmNeighbor(e, 2) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
			}
		}
	}
	else if (nsd == 3)
	{
		int n, e, n1, n2, n3, n4, ne;
// Given face of element e, find the neighbor elements which have 
// three nodes in common with the face. If no such element is found this
// is a boundary face.

		//  additional_info.initNeighbor(*this,true,false,false);
		WavesNeighborFE& neighbor = grid.getNeighbor();
		if (neighbor.nodeSize() == 0)
			neighbor.init(grid, true, false, false);

// face1 has nodes 1 2 4
// face2 has nodes 2 3 4
// face3 has nodes 1 3 4
// face4 has nodes 1 2 3
		int nel = grid.getNoElms();
		int nno = grid.getNoNodes();

		MV_ColMat<int> elmNeighbor(nel, 4);
		elmNeighbor = -1;

		bool safe = true;
		int ne1, ne2;
		for (e = nel - 1; e >= 0; e--)
		{ // more efficient with decreasing elements
			n1 = grid.loc2glob(e, 0);
			n2 = grid.loc2glob(e, 1);
			n3 = grid.loc2glob(e, 2);
			n4 = grid.loc2glob(e, 3);
			if (elmNeighbor(e, 0) == -1)
				if (elmNeighbor(e, 3) == -1)
				{ // case face 1 and 4 not done
					getElmsAtFaces(neighbor, n1, n2, n4, n3, e, ne1, ne2);
					if (ne1 >= 0)
					{
						n = getFaceNo(grid, n1, n2, n4, ne1, safe);
						elmNeighbor(e, 0) = ne1;
						elmNeighbor(ne1, n) = e;
					}
					if (ne2 >= 0)
					{
						n = getFaceNo(grid, n1, n2, n3, ne2, safe);
						elmNeighbor(e, 3) = ne2;
						elmNeighbor(ne2, n) = e;
					}
				}
				else
				{
					ne = getElmAtFace(neighbor, n1, n2, n4, e); // case face 1
					if (ne >= 0)
					{
						n = getFaceNo(grid, n1, n2, n4, ne, safe);
						elmNeighbor(e, 0) = ne;
						elmNeighbor(ne, n) = e;
					}
				}
			else if (elmNeighbor(e, 3) == -1)
			{
				ne = getElmAtFace(neighbor, n1, n2, n3, e); // case face 4
				if (ne >= 0)
				{
					n = getFaceNo(grid, n1, n2, n3, ne, safe);
					elmNeighbor(e, 3) = ne;
					elmNeighbor(ne, n) = e;
				}
			}

			if (elmNeighbor(e, 1) == -1)
				if (elmNeighbor(e, 2) == -1)
				{
					getElmsAtFaces(neighbor, n3, n4, n2, n1, e, ne1, ne2);
					if (ne1 >= 0)
					{
						n = getFaceNo(grid, n3, n4, n2, ne1, safe);
						elmNeighbor(e, 1) = ne1;
						elmNeighbor(ne1, n) = e;
					}
					if (ne2 >= 0)
					{
						n = getFaceNo(grid, n3, n4, n1, ne2, safe);
						elmNeighbor(e, 2) = ne2;
						elmNeighbor(ne2, n) = e;
					}
				}
				else
				{
					ne = getElmAtFace(neighbor, n2, n3, n4, e);
					if (ne >= 0)
					{
						n = getFaceNo(grid, n2, n3, n4, ne, safe);
						elmNeighbor(e, 1) = ne;
						elmNeighbor(ne, n) = e;
					}
				}
			else if (elmNeighbor(e, 2) == -1)
			{
				ne = getElmAtFace(neighbor, n1, n3, n4, e);
				if (ne >= 0)
				{
					n = getFaceNo(grid, n1, n3, n4, ne, safe);
					elmNeighbor(e, 2) = ne;
					elmNeighbor(ne, n) = e;
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

		for (e = 0; e < nel; e++)
		{
			if (elmNeighbor(e, 0) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
				boundaryElms(e) = 1;

			}
			if (elmNeighbor(e, 1) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
				boundaryElms(e) = 1;
			}
			if (elmNeighbor(e, 2) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryMarkers(grid.loc2glob(e, 3)) = 1;
				boundaryElms(e) = 1;
			}
			if (elmNeighbor(e, 3) == -1)
			{
				boundaryMarkers(grid.loc2glob(e, 0)) = 1;
				boundaryMarkers(grid.loc2glob(e, 1)) = 1;
				boundaryMarkers(grid.loc2glob(e, 2)) = 1;
				boundaryElms(e) = 1;
			}
		}
	}
}

real WavesOptional::sortMakeIndex(MV_Vector<real>& values, MV_Vector<int>& index, int sort_range)
{
	real temp;
	int l = sort_range / 2 + 1;
	int ir = sort_range;
	index.newsize(sort_range);
	int itemp;
	int k;
	for (k = 0; k < sort_range; k++)
		index(k) = k;

	first_jump: if (l > 1)
	{
		l = l - 1;
		temp = values[l - 1];
		itemp = index(l - 1);
	}
	else
	{
		temp = values[ir - 1];
		itemp = index(ir - 1);
		values[ir - 1] = values[0];
		index(ir - 1) = index(0);
		ir = ir - 1;
		if (ir == 1)
		{
			values[0] = temp;
			index(0) = itemp;
			return temp;
		}
	}
	int i = l;
	int j = l + l;
	second_jump: if (j <= ir)
	{
		if (j < ir)
			if (values[j - 1] < values[j])
				j++;
		if (temp < values[j - 1])
		{
			values[i - 1] = values[j - 1];
			index(i - 1) = index(j - 1);
			i = j;
			j = j + j;
		}
		else
			j = ir + 1;
		goto second_jump;
	}
	values[i - 1] = temp;
	index[i - 1] = itemp;
	goto first_jump;

	return 0;
}

//==========================================================================
void WavesOptional::extractNodeNumbers(MV_Vector<int>& exchangeMask, MV_Vector<int>& unodes, int num)
//========================================================================
{
	int n = exchangeMask.size();
	int count = 0;
	int i;
	for (i = 0; i < n; i++)
		if (exchangeMask(i) == num)
			count++;
	unodes.newsize(count);
	count = 0;
	for (i = 0; i < n; i++)
		if (exchangeMask(i) == num)
			unodes(count++) = i;
}
void WavesOptional::computeSortVectorfromGrid(MV_Vector<int>& unodes, MV_Vector<real>& sortkey)
//------
{
	nsd = gg.getNoSpaceDim();
	int nsize = unodes.size();
	sortkey.newsize(nsize);
	int i;
	for (i = 0; i < nsize; i++)
	{
		if (nsd == 2)
		{
			real x = gg.getCoor(unodes(i), 0);
			real y = gg.getCoor(unodes(i), 1);
			sortkey(i) = x * 1e-4 + y;
		}
		else
		{
			real x = gg.getCoor(unodes(i), 0);
			real y = gg.getCoor(unodes(i), 1);
			real z = gg.getCoor(unodes(i), 2);
			sortkey(i) = x * 1e-8 + y * 1e-4 + z;
		}
	}
}

void WavesOptional::makeSortedNodesIndex(MV_Vector<int>& exchangeMask, MV_Vector<int>& unodes, int num)
{
/*	for (int i = 0; i < unodes.size(); i++)
	{
		cout << "i=" << i << " unodes(i)=" << unodes(i) << endl << flush;
	}
*/
	extractNodeNumbers(exchangeMask, unodes, num);
//	cout << "ok1<<" << flush;
	MV_Vector<real> sortkey;
	computeSortVectorfromGrid(unodes, sortkey);
//	cout << "ok2<<" << flush;
	MV_Vector<int> index;
	sortMakeIndex(sortkey, index, sortkey.size());
//	cout << "ok3<<" << flush;
	int nsize = index.size();
	int i;

/*	for (i = 0; i < nsize; i++)
	{
		cout << "i=" << i << " unodes(i)=" << unodes(i) << " index(i)=" << index(i) << endl << flush;
	}
*/
	MV_Vector<int> ucopy(unodes.size());
	for (i = 0; i < nsize; i++)
		ucopy(i) = unodes(i);
	for (i = 0; i < nsize; i++)
		unodes(i) = ucopy(index(i));
}

//============================================================================
void WavesOptional::computeExNodesOuterbord(Grid& grid, MV_Vector<int>& bndNodes, MV_Vector<int>& markNodes, int &outer_bord_code)
{
//============================================================================
	int i, j, k;
	int nno = grid.getNoNodes();
	int nel = grid.getNoElms();

	bndNodes.newsize(nno);
	bndNodes = 0;
	// to fill boundary with code outer_bord_code
	for (i = 0; i < nno; i++)
	{
		if (markNodes(i) == 1)
		{
			markNodes(i) = outer_bord_code;
			bndNodes(i) = outer_bord_code;
		}
		else
		{
			markNodes(i) = 0;
			bndNodes(i) = 0;
		}
	}

}

//=======================================================================
//============================================================================
void WavesOptional::extractBoundNodes(Grid& grid_outer, Grid& grid, MV_Vector<int>& bndNodes)
{
//============================================================================
	int i, j, k;
	int nno = grid.getNoNodes();
	int nel = grid.getNoElms();

	int nno_out = grid_outer.getNoNodes();
	int nel_out = grid_outer.getNoElms();

	double eps = 1e-6;

	MV_Vector<int> outer_markNodes(nno_out);
	outer_markNodes = 1;

	MV_Vector<int> markNodes(nno);
	findBoundaryNodes(grid, markNodes);

	for (i = 0; i < nno; i++)
		cout << "markNodes(" << i << ") = " << markNodes(i) << endl;

	bndNodes.newsize(nno);
	bndNodes = 0;
	int nsd = grid.getNoSpaceDim();

	// to fill boundary with code outer_bord_code
	for (i = 0; i < nno; i++)
	{

		if (nsd == 2)
		{
			double x = grid.getCoor(i, 0);
			double y = grid.getCoor(i, 1);

			for (k = 0; k < nno_out; k++)
			{
				double x_out = grid_outer.getCoor(k, 0);
				double y_out = grid_outer.getCoor(k, 1);

				if (fabs(x - x_out) < eps && fabs(y - y_out) < eps)
				{
					if (markNodes(i) == 1)
					{
						bndNodes(i) = 2;
						break;
					}

					if (markNodes(i) == 0)
					{
						bndNodes(i) = 1;
						break;
					}
				}
			} // loop for grid_outer
		} // for nsd == 2
		else // if nsd == 3
		{
			double x = grid.getCoor(i, 0);
			double y = grid.getCoor(i, 1);
			double z = grid.getCoor(i, 2);

			for (k = 0; k < nno_out; k++)
			{
				double x_out = grid_outer.getCoor(k, 0);
				double y_out = grid_outer.getCoor(k, 1);
				double z_out = grid_outer.getCoor(k, 2);

				if (fabs(x - x_out) < eps && fabs(y - y_out) < eps && fabs(z - z_out) < eps)
				{
					if (markNodes(i) == 1)
					{
						bndNodes(i) = 2;
						break;
					}

					if (markNodes(i) == 0)
					{
						bndNodes(i) = 1;
						break;
					}
				}
			} // for k
		} // else
	} //loop for grid

/* commented by Thanh, Dec 05, 2012.
	ofstream outp;
	int sch = 0;
	outp.open("bndNodes1.dat");
	for (int l = 0; l < nno; l++)
	{
		if (bndNodes(l) == 2)
		{
			sch++;
			outp << "bndNodes(" << l << ")=" << bndNodes(l) << "\n";
		}
	}
	outp << sch << "\n";
*/ 
}

//=========================================================================
// to fill inner boundary with code 1
//==================================================================
void WavesOptional::computeExNodesCommon(Grid& grid, MV_Vector<int>& bndNodes, MV_Vector<int>& markNodes, const int &bord_code1, const int &bord_code2)
{
	int nel = grid.getNoElms();
	int i, j, k;
	// to fill second boundary with code bord_code1
	for (i = 0; i < nel; i++)
	{
		int ne = grid.getNoNodesInElm(i);
		for (j = 0; j < ne; j++)
		{
			int node1 = grid.loc2glob(i, j);
			if (markNodes(node1) == bord_code2)
			{
				for (k = 0; k < ne; k++)
				{
					int node2 = grid.loc2glob(i, k);
					if (markNodes(node2) == 0)
					{
						bndNodes(node2) = bord_code1;
						markNodes(node2) = bord_code1;
					}
				}
			}
		}
	}
}

//===========================================================================  
void WavesOptional::makeExchangeNodes(Grid& grid, const overlap_code &ovc, MV_Vector<int>& bndNodes)
//===========================================================================
{
	int code1, code2, code3, code4, code5;
	int one = 1;
	int nno = grid.getNoNodes();
	MV_Vector<int> markNodes(nno);
	findBoundaryNodes(grid, markNodes);

/* commented by Thanh, Dec 05, 2012.
	int sch = 0;
	ofstream outpfem;
	outpfem.open("bndNodes_fem.dat");
	for (int l = 0; l < nno; l++)
	{
		if (markNodes(l) == 1)
			sch++;

		outpfem << "bndNodes(" << l << ")=" << markNodes(l) << "\n";

	}
*/
	//cout << "number of b.points is " << sch << endl;

	switch (ovc)
	{
		case over5:
		{
			code1 = 22222;
			code2 = 2222;
			code3 = 222;
			code4 = 22;
			code5 = 2;
			computeExNodesOuterbord(grid, bndNodes, markNodes, code1);
			computeExNodesCommon(grid, bndNodes, markNodes, code2, code1);
			computeExNodesCommon(grid, bndNodes, markNodes, code3, code2);
			computeExNodesCommon(grid, bndNodes, markNodes, code4, code3);
			computeExNodesCommon(grid, bndNodes, markNodes, code5, code4);
			computeExNodesCommon(grid, bndNodes, markNodes, one, code5);
			break;
		}
		case over4:
		{
			code2 = 2222;
			code3 = 222;
			code4 = 22;
			code5 = 2;
			computeExNodesOuterbord(grid, bndNodes, markNodes, code2);
			computeExNodesCommon(grid, bndNodes, markNodes, code3, code2);
			computeExNodesCommon(grid, bndNodes, markNodes, code4, code3);
			computeExNodesCommon(grid, bndNodes, markNodes, code5, code4);
			computeExNodesCommon(grid, bndNodes, markNodes, one, code5);
			break;
		}
		case over3:
		{
			code3 = 222;
			code4 = 22;
			code5 = 2;
			computeExNodesOuterbord(grid, bndNodes, markNodes, code3);
			computeExNodesCommon(grid, bndNodes, markNodes, code4, code3);
			computeExNodesCommon(grid, bndNodes, markNodes, code5, code4);
			computeExNodesCommon(grid, bndNodes, markNodes, one, code5);

			//for (int n = 0; n < nno; n++)
			//	cout << "bndNodes=" << bndNodes(n) << endl;
			break;
		}
		case over2:
		{
			code4 = 22;
			code5 = 2;

			computeExNodesOuterbord(grid, bndNodes, markNodes, code4);

			computeExNodesCommon(grid, bndNodes, markNodes, code5, code4);
			computeExNodesCommon(grid, bndNodes, markNodes, one, code5);
			// print to file 
/* commented by Thanh, Dec 05, 2012.
			ofstream outp;
			outp.open("bndNodes1.dat");
			for (int l = 0; l < nno; l++)
			{
				if (bndNodes(l) != 0)
				{
					outp << "bndNodes(" << l << ")=" << bndNodes(l) << "\n";
				}
			}
*/

			break;
		}
		case over1:
		{

			code5 = 2;

			// the nodes which are at the boundary marked with code5=2
			computeExNodesOuterbord(grid, bndNodes, markNodes, code5);

			// in the array bndNodes remains only nodes at the outer boundary
			// of the inner domain

			computeExNodesCommon(grid, bndNodes, markNodes, one, code5);

/* commented by Thanh, Dec 05, 2012.
			// print to file 
			ofstream outp;
			outp.open("bndNodes1.dat");
			for (int l = 0; l < nno; l++)
			{
				outp << "bndNodes(" << l << ")=" << bndNodes(l) << "\n";
			}
*/
			break;
		}
	}

}

int WavesOptional::findInitDisturbPoint(MV_Vector<real>& evalpt)
{
	nsd = evalpt.size();
	int i;
	int nno = gg.getNoNodes();
	// HRN
	int closestnode = 0;

	if (nsd == 2)
	{
		double x0 = evalpt(0);
		double y0 = evalpt(1);
		double dx = gg.getCoor(0, 0) - x0;
		double dy = gg.getCoor(0, 1) - y0;
		double mindist = dx * dx + dy * dy;
		for (i = 0; i < nno; i++)
		{
			dx = gg.getCoor(i, 0) - x0;
			dy = gg.getCoor(i, 1) - y0;
			if (dx * dx + dy * dy < mindist)
			{
				mindist = dx * dx + dy * dy;
				closestnode = i;
			}
		}
	}
	else if (nsd == 3)
	{
		double x0 = evalpt(0);
		double y0 = evalpt(1);
		double z0 = evalpt(2);
		double dx = gg.getCoor(0, 0) - x0;
		double dy = gg.getCoor(0, 1) - y0;
		double dz = gg.getCoor(0, 2) - z0;
		double mindist = dx * dx + dy * dy + dz * dz;
		for (i = 0; i < nno; i++)
		{
			dx = gg.getCoor(i, 0) - x0;
			dy = gg.getCoor(i, 1) - y0;
			dz = gg.getCoor(i, 2) - z0;
			if (dx * dx + dy * dy < mindist)
			{
				mindist = dx * dx + dy * dy + dz * dz;
				closestnode = i;
			}
		}
	}

	// HRN
	return closestnode;
}

void WavesOptional::print_gid_mesh(char *file)

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

	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), gg.getCoor(i, 2));
	}
	else if (nsd == 2)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Triangle  Nnode  3\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");

	for (el = 0; el < elnum; el++)
	{
		switch (gg.getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i   %i %i %i %i %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 4) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el));
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.getMaterialType(el));
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i   %i %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.getMaterialType(el));
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "end elements\n");
	fclose(fp);
}

void WavesOptional::print_gid_mesh_FEM(WavesGridB& gg, char *file)

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

	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), gg.getCoor(i, 2));
	}
	else if (nsd == 2)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Triangle  Nnode  3\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");

	for (el = 0; el < elnum; el++)
	{
		switch (gg.getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i   %i %i %i %i %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 4) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i   %i %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.getMaterialType(el) + 1);
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "end elements\n");
	fclose(fp);
}

void WavesOptional::print_gid_mesh_Amira(WavesGridB& gg, char *file)

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

	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), gg.getCoor(i, 2));
	}
	else if (nsd == 2)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Triangle  Nnode  3\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");

	for (el = 0; el < elnum; el++)
	{
		switch (gg.getElementType(el))
		{

			case ELMTET1:
				fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i   %i %i %i %i %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 4) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i   %i %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.getMaterialType(el) + 1);
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "end elements\n");
	fclose(fp);
}

//==========================================================================

void WavesOptional::print_gid_jump(WavesGridB& gg, char *file, MV_Vector<int>& jump)

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

	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), gg.getCoor(i, 2));
	}
	else if (nsd == 2)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Triangle  Nnode  3\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");

	for (el = 0; el < elnum; el++)
	{
		switch (gg.getElementType(el))
		{

			case ELMTET1:
				fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 0) + 1, jump(el));
				break;
			case ELMPYR1:
				fprintf(fp, "%i %i   %i %i %i %i %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 4) + 1, gg.loc2glob(el, 0) + 1, jump(el));
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, jump(el));
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i   %i %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, jump(el));
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fprintf(fp, "end elements\n");
	fclose(fp);
}

//==========================================================================

void WavesOptional::print_gid_nodes_FEM(WavesGridB& gg, char *file)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;

	fp = fopen(file, "w");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int nsd = gg.getNoSpaceDim();
	int nno_sdg = sdg.getNoNodes();
	int nel_sdg = sdg.getNoElms();

	int nnode = nno;
	int elnum = nel;

	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Tetrahedra  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), gg.getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		fprintf(fp, "MESH   dimension  %i  ElemType Triangles  Nnode  3\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n", i + 1, gg.getCoor(i, 0), gg.getCoor(i, 1), 0.0);
	}

	fclose(fp);
}

void WavesOptional::print_gid_elements_FEM(WavesGridB& gg, char *file)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int nsd = gg.getNoSpaceDim();
	int nno_sdg = sdg.getNoNodes();

	int nnode = nno;
	int elnum = nel;

	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");

	for (el = 0; el < elnum; el++)
	{
		switch (gg.getElementType(el))
		{
			case ELMTET1:
				fprintf(fp, "%i   %i    %i   %i   %i   %i\n", el + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.loc2glob(el, 0) + 1, gg.getMaterialType(el) + 1);
				break;
			case ELMTRI1:
				fprintf(fp, "%i %i  %i %i \n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1);
				break;
			case ELMQUAD1:
				fprintf(fp, "%i %i   %i %i %i %i\n", el + 1, gg.loc2glob(el, 0) + 1, gg.loc2glob(el, 1) + 1, gg.loc2glob(el, 2) + 1, gg.loc2glob(el, 3) + 1, gg.getMaterialType(el));
				break;
			case NO_ELEMENT:
				printf("error print wrong type of element\n");
				break;
		}
	}
	fclose(fp);
}
//===============================================================================

int WavesOptional::print_1Dgid_result(char *file, MV_Vector<double>& results, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{

		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	// here, we print out displacement of the elastic vector in  2D
	// Vector_results is results title
	// 2 - type of analysis (displacement analysis)
	// %i - time step
	// 1 - kind of results - scalar
	// 1 - position of the data - in nodes
	// 0 - no description inside

	fprintf(fp, "DISPLACEMENT   2   %i    1    1    0 \n", k);

	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %f \n", i + 1, results(i));
	}
	fprintf(fp, "\n");

	fclose(fp);

	return 0;

}

int WavesOptional::print_params(char *file, MV_Vector<double>& results, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{
		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	// here, we print out displacement of the elastic vector in  2D
	// Vector_results is results title
	// 2 - type of analysis (displacement analysis)
	// %i - time step
	// 1 - kind of results - scalar
	// 1 - position of the data - in nodes
	// 0 - no description inside

	fprintf(fp, "DISPLACEMENT   2   %i    2   1    0 \n", k);

	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %f %f %f \n", i + 1, results(i), results(i), results(i));
	}
	fprintf(fp, "\n");

	fclose(fp);

	return 0;
}

int WavesOptional::print_params_int(char *file, MV_Vector<int>& results, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{

		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	// here, we print out displacement of the elastic vector in  2D
	// Vector_results is results title
	// 2 - type of analysis (displacement analysis)
	// %i - time step
	// 1 - kind of results - scalar
	// 1 - position of the data - in nodes
	// 0 - no description inside

	fprintf(fp, "DISPLACEMENT   2   %i    1   1    0 \n", k);

	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %i %i %i \n", i + 1, results(i), results(i), results(i));
	}
	fprintf(fp, "\n");

	fclose(fp);

	return 0;
}

//=========================================================================
//===========================================================================

int WavesOptional::print_2Dgid_result(char *file, real* u_1_, real* u_2_, real* u_array, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{
		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

//      fprintf(fp,"DISPLACEMENT   2   %i    2    1    0 \n", k);

		if (k == 1)
			fprintf(fp, "GiD Post Results File 1.0");

		fprintf(fp, "\n");
		fprintf(fp, "Result  \"Displacements\"  \"Wave movement\"  %i    Vector OnNodes   \n", k);
		fprintf(fp, "Values");
		fprintf(fp, "\n");
		for (i = 0; i < nnode; i++)
		{
			fprintf(fp, "%i %f %f %f\n", i + 1, u_1_[i], u_2_[i], u_array[i]);

			cout << " u_1_(" << i << ") = " << u_1_[i] << " u_2_" << u_2_[i] << endl;
		}

		fprintf(fp, "end values");
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 0;

}

//===========================================================================
//===========================================================================

int WavesOptional::print_2Dgid_result(char *file, MV_Vector<real>& u1, MV_Vector<real>& u2, MV_Vector<real>& u3, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	nno = gg.getNoNodes();
	nel = gg.getNoElms();
	nsd = gg.getNoSpaceDim();

	int nnode = nno;
	int elnum = nel;

	// here, we print out displacement of the elastic vector in  2D
	// Vector_results is results title
	// 2 - type of analysis (displacement analysis)
	// %i - time step
	// 2 - kind of results - vector
	// 1 - position of the data - in nodes
	// 0 - no description inside

	if (k == 1)
		fprintf(fp, "GiD Post Results File 1.0");

	fprintf(fp, "\n");
	fprintf(fp, "Result  \"Displacements\"  \"Wave movement\"  %i    Vector OnNodes   \n", k);
	fprintf(fp, "Values");
	fprintf(fp, "\n");
	for (i = 0; i < nnode; i++)
	{
		fprintf(fp, "%i %f %f %f\n", i + 1, u1(i), u2(i), u3(i));

	}

	fprintf(fp, "end values");
	fprintf(fp, "\n");

	fclose(fp);

	return 0;

}

//===========================================================================

int WavesOptional::print_2Dgid_result_common(char *file, real* u_1_, real* u_2_, real* u_array, int k)

{

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{

		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		if (k == 1)
			fprintf(fp, "GiD Post Results File 1.0");

		fprintf(fp, "\n");
		fprintf(fp, "Result  \"Displacements\"  \"Wave movement\"  %i    Vector OnNodes   \n", k);
		fprintf(fp, "Values");
		fprintf(fp, "\n");
		for (i = 0; i < nnode; i++)
		{
			fprintf(fp, "%i %f %f %f\n", i + 1, u_1_[i], u_2_[i], u_array[i]);

			cout << " u_1_(" << i << ") = " << u_1_[i] << " u_2_" << u_2_[i] << endl;
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 0;

}

int WavesOptional::print_3Dgid_result(char *file, real* u_array, real* u_1_, real* u_2_, real* u_3_, int k)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{
		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{
		cout << " Call function print_2Dgid_result" << endl;
	}
	else if (nsd == 3)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f %f\n", i + 1, u_1_[i], u_2_[i], u_3_[i], u_array[i]);

		fprintf(fp, "\n");

	}
	fclose(fp);

	return 0;

}

int WavesOptional::print_3Dgid_result(char *file, MV_Vector<double>& u_array, MV_Vector<double>& u_1_, MV_Vector<double>& u_2_, MV_Vector<double>& u_3_, int k)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	nno = gg.getNoNodes();
	nel = gg.getNoElms();
	nsd = gg.getNoSpaceDim();

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{
		cout << " Call function print_2Dgid_result" << endl;
	}
	else if (nsd == 3)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		cout << "inside print 3D " << endl;

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f %f\n", i + 1, u_1_(i), u_2_(i), u_3_(i), u_array(i));

		fprintf(fp, "\n");

	}
	fclose(fp);

	return 0;

}




int WavesOptional::print_3Dgid_result_common( char *file, 
				  real*  u_array,
				  real* u_1_,
				  real* u_2_,
				  real* u_3_,
				  int k) 
  //-----------------------------------------------------------------------------
{
  FILE *fp;   
  int i, j;
  int el;
  int d;
  double tmax;
  double dd;
  int nno,nel,nsd;  


  fp = fopen(file, "a+");
  
  
  if (!fp) 
    {
      cout<< " cannot open file \""<<file<<"\"\n";
      exit(1);
    }
    

  fseek(fp,1,SEEK_END);



  if (USE_FEM || EXCHANGE)
    {
      nno = gg.getNoNodes();
      nel = gg.getNoElms();
      nsd = gg.getNoSpaceDim();
    }
 
 
  int nnode = nno;
  int elnum = nel;

  if (nsd == 2)
    {cout<<" Call function print_2Dgid_result_common"<<endl;
    }
  else if(nsd==3)
    {

      // here, we print out displacement of the elastic vector in  2D
      // Vector_results is results title
      // 2 - type of analysis (displacement analysis)
      // %i - time step
      // 2 - kind of results - vector
      // 1 - position of the data - in nodes
      // 0 - no description inside
      
     

	if (k==1)
	fprintf(fp,"GiD Post Results File 1.0");
  
	fprintf(fp,"\n");  
        fprintf(fp,"Result  \"Displacements\"  \"Wave movement\"  %i    Vector OnNodes   \n", k);
	fprintf(fp,"Values");
       fprintf(fp,"\n");    
      for(i = 0; i < nnode; i++)
	{
	fprintf(fp, "%i %f %f %f %f\n ", i + 1,
		u_1_[i], u_2_[i], u_3_[i], u_array[i]);
      
	//  cout<<" u_1_("<<i<<") = "<<u_1_[i]<<" u_2_"<<u_2_[i]<<endl;
	}
          
     
        fprintf(fp,"\n");    
	     /*    
      fprintf(fp,"DISPLACEMENT   2   %i    2    1    0 \n", k);
      
      
      for(i = 0; i < nnode; i++)
	fprintf(fp, "%i %f %f %f %f\n ", i + 1,
		u_1_[i], u_2_[i], u_3_[i], u_array[i]);
      
      
      
      fprintf(fp,"\n");          
	     */
    
    
    }
    fclose(fp);
    
  
  return 0;

}

int WavesOptional::print_3Dgid_LENS(char *file, real* u_1_, real* u_2_, real* u_3_, int k)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{
		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{
		cout << " Call function print_2Dgid_result" << endl;
	}
	else if (nsd == 3)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f \n", i + 1, u_1_[i], u_2_[i], u_3_[i]);

		fprintf(fp, "\n");

	}
	fclose(fp);

	return 0;

}
//=====================================================================

int WavesOptional::print_onepoint(char *file, real& u_array, real& u_1_, real& u_2_, int k)

{

	cout << "inside " << endl;

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{

		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}
	else if (USE_FDM)
	{
		nno = sdg.getNoNodes();
		nel = sdg.getNoElms();
		nsd = sdg.getNoSpaceDim();
		cout << " nno = " << nno << " nel = " << nel << "  nsd = " << nsd << endl;
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		if (k == 1)
			fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		fprintf(fp, "%i %f %f %f\n", k, u_1_, u_2_, u_array);

		fprintf(fp, "\n");

	}
	fclose(fp);
	return 0;

}

int WavesOptional::print_3D_onepoint(char *file, real& u_array, real& u_1_, real& u_2_, real& u_3_, int k)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nel, nsd;

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	if (USE_FEM || EXCHANGE)
	{
		nno = gg.getNoNodes();
		nel = gg.getNoElms();
		nsd = gg.getNoSpaceDim();
	}

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{
		cout << " Call function print_2D_onepoint" << endl;
	}
	else if (nsd == 3)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		if (k == 1)
			fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		fprintf(fp, "%i %f %f %f %f\n", k, u_1_, u_2_, u_3_, u_array);

		fprintf(fp, "\n");

	}
	fclose(fp);

	return 0;

}

//====================================================================

int WavesOptional::print_common_2DFEM(char *file, real* u_array, real* u_1_, real* u_2_, int k)

{

	cout << "inside " << endl;

	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nno_sdg, nel, nsd;

	fp = fopen(file, "a+");
	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	nno = gg.getNoNodes();

	nel = gg.getNoElms();
	nsd = gg.getNoSpaceDim();

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (i = 0; i < nnode; i++)
		{
			fprintf(fp, "%i %f %f %f\n", i + 1, u_1_[i], u_2_[i], u_array[i]);

			cout << " u_1_(" << i << ") = " << u_1_[i] << " u_2_" << u_2_[i] << endl;
		}
	}
	fclose(fp);

	return 0;

}

//====================================================================

int WavesOptional::print_common_3DFEM(char *file, real* u_array, real* u_1_, real* u_2_, real* u_3_, int k)

{
	FILE *fp;
	int i, j;
	int el;
	int d;
	double tmax;
	double dd;
	int nno, nno_sdg_, nel, nsd;

	fp = fopen(file, "a+");
	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	nno = gg.getNoNodes();

	nel = gg.getNoElms();
	nsd = gg.getNoSpaceDim();

	int nnode = nno;
	int elnum = nel;

	if (nsd == 2)
	{
		cout << " Call function print_common_2DFEM" << endl;
	}
	else if (nsd == 3)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f %f\n", i + 1, u_1_[i], u_2_[i], u_3_[i], u_array[i]);

	}
	fclose(fp);

	return 0;

}

//====================================================================
void WavesOptional::Write_To_File(char *file_name, MV_Vector<double>& Ar1)
{
	ofstream outp;
	int i, j;

	outp.open(file_name);

	for (i = 0; i < Ar1.size(); i++)
		outp << Ar1(i) << "   ";

	outp.close();

	//print("Solution at one point is written to file: %s\n", file_name);
}

int WavesOptional::writeInp(char *file, Grid *grid, Vec *u, int nsys)
{
	FILE *fp;
	int i, j;
	int el;
	PetscScalar *solArr;
	int NODE_DATA = 1;

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	for (i = 0; i < grid->getNoNodes(); i++)
	{
		fprintf(fp, "%i %f %f %f\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), grid->getCoor(i, 2));
	}
	for (i = 0; i < grid->getNoElms(); i++)
	{
		fprintf(fp, "%i %i tet %i %i %i %i\n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1, grid->loc2glob(i, 3) + 1);
	}

	ierr = VecGetArray(*u, &solArr);

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%f ", solArr[j + i * nsys]);
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < grid->getNoElms(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", grid->getMaterialType(i));
		}
	}
	fclose(fp);

	ierr = VecRestoreArray(*u, &solArr);
	CHKERRQ( ierr);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//========================================================================
// write out mesh in 3D
//=======================================================================
int WavesOptional::writeInpElType(char *file, Grid *grid, Vec *u, int nsys, int type_of_mat, int NODE_DATA)
{
	FILE *fp;
	int i, j;
	int el;
	PetscScalar *solArr;

	int sch_mat = 0;
	int sch_nodes = 0;

	ierr = VecGetArray(*u, &solArr);

	Mat_int Materials(grid->getNoElms(), 5);

	Mat_int Markers(grid->getNoNodes(), 2);
	Markers = 0;

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	for (i = 0; i < grid->getNoElms(); i++)
	{
		if (grid->getMaterialType(i) != type_of_mat)
		{

			Materials(sch_mat, 0) = grid->loc2glob(i, 0);
			Materials(sch_mat, 1) = grid->loc2glob(i, 1);
			Materials(sch_mat, 2) = grid->loc2glob(i, 2);
			Materials(sch_mat, 3) = grid->loc2glob(i, 3);

			Materials(sch_mat, 4) = grid->getMaterialType(i);
			cout << " Materials " << Materials(sch_mat, 4) << "  Materials(sch_mat,0) " << Materials(sch_mat, 0) << endl;
			sch_mat++;

			Markers(grid->loc2glob(i, 0), 0) = 1;
			Markers(grid->loc2glob(i, 1), 0) = 1;
			Markers(grid->loc2glob(i, 2), 0) = 1;
			Markers(grid->loc2glob(i, 3), 0) = 1;

			Markers(grid->loc2glob(i, 0), 1) = grid->getMaterialType(i);
			Markers(grid->loc2glob(i, 1), 1) = grid->getMaterialType(i);
			Markers(grid->loc2glob(i, 2), 1) = grid->getMaterialType(i);
			Markers(grid->loc2glob(i, 3), 1) = grid->getMaterialType(i);

		}
	}

	for (i = 0; i < grid->getNoNodes(); i++)
		if (Markers(i, 0) == 1)
			sch_nodes++;

	Mat_real Coordinates(sch_nodes, 4);
	MV_Vector<int> Loc_Glob(sch_nodes);
	MV_Vector<int> MatType(sch_nodes);

	int n = 0;

	for (i = 0; i < grid->getNoNodes(); i++)
	{

		if (Markers(i, 0) == 1)
		{
			Coordinates(n, 0) = grid->getCoor(i, 0);
			Coordinates(n, 1) = grid->getCoor(i, 1);
			Coordinates(n, 2) = grid->getCoor(i, 2);

			if (NODE_DATA)
				Coordinates(n, 3) = solArr[i];
			else
				MatType(n) = Markers(i, 1);

			Loc_Glob(n) = i;

			n++;
		}

	}

// change global number to local

	for (i = 0; i < sch_mat; i++)
		for (j = 0; j < sch_nodes; j++)
		{
			if (Materials(i, 0) == Loc_Glob(j))
				Materials(i, 0) = j;
			if (Materials(i, 1) == Loc_Glob(j))
				Materials(i, 1) = j;
			if (Materials(i, 2) == Loc_Glob(j))
				Materials(i, 2) = j;
			if (Materials(i, 3) == Loc_Glob(j))
				Materials(i, 3) = j;
		}

//===============================  print out mesh =============

	fprintf(fp, "%i %i %i 0 0\n", sch_nodes, sch_mat, nsys);

	for (i = 0; i < sch_nodes; i++)
	{
		fprintf(fp, "%i %f %f %f\n", i + 1, Coordinates(i, 0), Coordinates(i, 1), Coordinates(i, 2));
	}
	for (i = 0; i < sch_mat; i++)
	{

		fprintf(fp, "%i %i tet %i %i %i %i\n", i + 1, Materials(i, 4), Materials(i, 0) + 1, Materials(i, 1) + 1, Materials(i, 2) + 1, Materials(i, 3) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < sch_nodes; i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%f ", Coordinates(i, 3));
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < sch_nodes; i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", MatType(i));
		}
	}
	fclose(fp);

	ierr = VecRestoreArray(*u, &solArr);
	CHKERRQ( ierr);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//===================================================================0
// write out mesh in 2D
//====================================================================
int WavesOptional::writeInpElType2D(char *file, Grid *grid, Vec *u, int nsys, int type_of_mat, int NODE_DATA)
{
	FILE *fp;
	int i, j;
	int el;
	PetscScalar *solArr;

	int sch_mat = 0;
	int sch_nodes = 0;

	ierr = VecGetArray(*u, &solArr);

	Mat_int Materials(grid->getNoElms(), 4);

	Mat_int Markers(grid->getNoNodes(), 2);
	Markers = 0;

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	for (i = 0; i < grid->getNoElms(); i++)
	{
		if (grid->getMaterialType(i) != type_of_mat)
		{

			Materials(sch_mat, 0) = grid->loc2glob(i, 0);
			Materials(sch_mat, 1) = grid->loc2glob(i, 1);
			Materials(sch_mat, 2) = grid->loc2glob(i, 2);

			Materials(sch_mat, 3) = grid->getMaterialType(i) + 1;
			cout << " Materials " << Materials(sch_mat, 3) << "  Materials(sch_mat,0) " << Materials(sch_mat, 0) << endl;
			sch_mat++;

			Markers(grid->loc2glob(i, 0), 0) = 1;
			Markers(grid->loc2glob(i, 1), 0) = 1;
			Markers(grid->loc2glob(i, 2), 0) = 1;

			Markers(grid->loc2glob(i, 0), 1) = grid->getMaterialType(i);
			Markers(grid->loc2glob(i, 1), 1) = grid->getMaterialType(i);
			Markers(grid->loc2glob(i, 2), 1) = grid->getMaterialType(i);
		}
	}

	for (i = 0; i < grid->getNoNodes(); i++)
		if (Markers(i, 0) == 1)
			sch_nodes++;

// cout<<"after  creating array with materials "<<"sch_nodes "<<sch_nodes<<endl;

	Mat_real Coordinates(sch_nodes, 3);
	MV_Vector<int> Loc_Glob(sch_nodes);
	MV_Vector<int> MatType(sch_nodes);

	int n = 0;

	for (i = 0; i < grid->getNoNodes(); i++)
	{

		if (Markers(i, 0) == 1)
		{
			Coordinates(n, 0) = grid->getCoor(i, 0);
			Coordinates(n, 1) = grid->getCoor(i, 1);

			if (NODE_DATA)
				Coordinates(n, 2) = solArr[i];
			else
				MatType(n) = Markers(i, 1);

			Loc_Glob(n) = i;

			n++;
		}

	}

// change global number to local

	for (i = 0; i < sch_mat; i++)
		for (j = 0; j < sch_nodes; j++)
		{
			if (Materials(i, 0) == Loc_Glob(j))
				Materials(i, 0) = j;
			if (Materials(i, 1) == Loc_Glob(j))
				Materials(i, 1) = j;
			if (Materials(i, 2) == Loc_Glob(j))
				Materials(i, 2) = j;
		}

//===============================  print out mesh =============

	fprintf(fp, "%i %i %i 0 0\n", sch_nodes, sch_mat, nsys);

	for (i = 0; i < sch_nodes; i++)
	{
		fprintf(fp, "%i %f %f 0.0 \n", i + 1, Coordinates(i, 0), Coordinates(i, 1));
	}
	for (i = 0; i < sch_mat; i++)
	{

		fprintf(fp, "%i %i tri %i %i %i \n", i + 1, Materials(i, 3), Materials(i, 0) + 1, Materials(i, 1) + 1, Materials(i, 2) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < sch_nodes; i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%f ", Coordinates(i, 2));
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < sch_nodes; i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", MatType(i));
		}
	}
	fclose(fp);

	ierr = VecRestoreArray(*u, &solArr);
	CHKERRQ( ierr);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//==================================================================

int WavesOptional::writeInp2D(char *file, Grid *grid, Vec *u, int nsys, int p3d)
{
	FILE *fp;
	int i, j;
	int el;
	PetscScalar *solArr;
	int NODE_DATA = 1;
	ierr = VecGetArray(*u, &solArr);

	// p3d ==-1 plot third component 0.0
	// p3d ==i plot as third component field i

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	if (p3d >= 0)
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %15.8e\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), solArr[p3d + i * nsys]);

		}
	}
	else
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %7.5f\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), 0.0 /*solArr[i*nsys]*/);
		}
	}
	for (i = 0; i < grid->getNoElms(); i++)
	{
		ElementType etyp = grid->getElementType(i);
		if (etyp == ELMQUAD1)
			fprintf(fp, "%i %i quad %i %i %i %i\n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1, grid->loc2glob(i, 3) + 1);
		else
			fprintf(fp, "%i %i tri %i %i %i\n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%15.8e ", solArr[i * nsys]);
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < grid->getNoElms(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", grid->getMaterialType(i));
		}
	}
	fclose(fp);

	ierr = VecRestoreArray(*u, &solArr);
	CHKERRQ( ierr);
//	printf("Solution written to file: %s\n", file);

	return 0;
}

//============== another parameters ==================================

//==================================================================

int WavesOptional::writeInp2D(char *file, Grid *grid, MV_Vector<double>& b, int nsys, int p3d)
{
	FILE *fp;
	int i, j;
	int el;
	int NODE_DATA = 1;
	// p3d ==-1 plot third component 0.0
	// p3d ==i plot as third component field i

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	if (p3d >= 0)
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %15.8e\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), b(p3d + i * nsys));

		}
	}
	else
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %7.5f\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), 0.0 /*solArr[i*nsys]*/);
		}
	}
	for (i = 0; i < grid->getNoElms(); i++)
	{
		ElementType etyp = grid->getElementType(i);
		if (etyp == ELMQUAD1)
			fprintf(fp, "%i %i quad %i %i %i %i\n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1, grid->loc2glob(i, 3) + 1);
		else
			fprintf(fp, "%i %i tri %i %i %i\n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%15.8e ", b(i * nsys));
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < grid->getNoElms(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", grid->getMaterialType(i));
		}
	}
	fclose(fp);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//===================== 3D case =======================================
//==================================================================

int WavesOptional::writeInp3D(char *file, Grid *grid, MV_Vector<double>& b, int nsys)
{
	FILE *fp;
	int i, j;
	int el;
	int NODE_DATA = 1;
	// p3d ==-1 plot third component 0.0
	// p3d ==i plot as third component field i

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	for (i = 0; i < grid->getNoNodes(); i++)
	{
		fprintf(fp, "%i %f %f %f \n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), grid->getCoor(i, 2));
	}
	for (i = 0; i < grid->getNoElms(); i++)
	{
		ElementType etyp = grid->getElementType(i);
		if (etyp == ELMTET1)
			fprintf(fp, "%i %i tet %i %i %i %i \n", i + 1, grid->getMaterialType(i), grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 2) + 1, grid->loc2glob(i, 3) + 1);
	}

	if (NODE_DATA)
	{
		fprintf(fp, "%i ", nsys);
		for (i = 0; i < nsys; i++)
			fprintf(fp, "1 ");
		fprintf(fp, "\n");
		for (i = 0; i < nsys; i++)
			fprintf(fp, "u %i, Unit\n", i + 1);
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			for (j = 0; j < nsys; j++)
				fprintf(fp, "%15.8e ", b(i * nsys));
			fprintf(fp, "\n");
		}
	}
	else
	{
		fprintf(fp, "1 1");
		fprintf(fp, "\n");
		fprintf(fp, "u , Unit\n");
		for (i = 0; i < grid->getNoElms(); i++)
		{
			fprintf(fp, "%i ", i + 1);
			fprintf(fp, "%i \n", grid->getMaterialType(i));
		}
	}
	fclose(fp);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//====================================================================

int WavesOptional::writeInpFDM2D(char *file, WavesSDGeometry *grid, double *u, int nsys, int p3d)
{
	FILE *fp;
	int i, j;
	int el;

	// p3d ==-1 plot third component 0.0
	// p3d ==i plot as third component field i

	fp = fopen(file, "w");
	if (!fp)
		SETERRA(1, 1, "Cannot open grid file");

	fprintf(fp, "%i %i %i 0 0\n", grid->getNoNodes(), grid->getNoElms(), nsys);

	if (p3d >= 0)
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %15.8e\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), u[p3d + i * nsys]);

		}
	}
	else
	{
		for (i = 0; i < grid->getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %7.5f\n", i + 1, grid->getCoor(i, 0), grid->getCoor(i, 1), 0.0 /*solArr[i*nsys]*/);
		}
	}

	for (i = 0; i < grid->getNoElms(); i++)
	{
		if (nsd == 3)
			fprintf(fp, "%i %i hex  %i %i %i %i  %i %i %i %i\n", i + 1, 1, grid->loc2glob(i, 4) + 1, grid->loc2glob(i, 5) + 1, grid->loc2glob(i, 7) + 1, grid->loc2glob(i, 6) + 1, grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 3) + 1,
					grid->loc2glob(i, 2) + 1);
		if (nsd == 2)
			fprintf(fp, "%i %i quad %i %i %i %i\n", i + 1, 1, grid->loc2glob(i, 0) + 1, grid->loc2glob(i, 1) + 1, grid->loc2glob(i, 3) + 1, grid->loc2glob(i, 2) + 1);
	}

	fprintf(fp, "1 1");
	fprintf(fp, "\n");
	fprintf(fp, "u , Unit\n");

	for (int n = 0; n < grid->getNoNodes(); n++)
		fprintf(fp, "%i  %15.8e \n", i + 1, u[n]);

	fclose(fp);

//	printf("Solution written to file: %s\n", file);

	return 0;
}

//====================================================================

void WavesOptional::ApplyAvsOut(char* filename, Vec& sol)
{
	int noeq = 1;
	cout << "nsd = " << nsd << endl;

	if (nsd == 2)
	{
		ierr = writeInp2D(filename, &gg, &sol, noeq, 0);
		CHKERRA(ierr);
	}
	else
	{
		ierr = writeInp(filename, &gg, &sol, noeq);
		CHKERRA(ierr);
	}

}

void WavesOptional::getFileNameFDM(char* buf)
{
	// static int i = 100; // for f= e_3rhs0_5_0_3
	//static int i = 10;  // for f = e_3rhs0_5_0_5
	static int i = 1; // for f = e_3rhs0_5_0_5 for inner and f= e_3rhs0_5_0_3 for outer 
	sprintf(buf, "outOut%03d.out", i++);
}

//===================================================================0
void WavesOptional::getFileNameElastic1(char* buf)
{
	// static int i = 100; // for f= e_3rhs0_5_0_3
	//static int i = 10;  // for f = e_3rhs0_5_0_5
	static int i = 1; // for f = e_3rhs0_5_0_5 for inner and f= e_3rhs0_5_0_3 for outer 
	sprintf(buf, "Wave_u1_%03d.inp", i++);
	return;
}

//==============================================================
void WavesOptional::getFileNameElastic2(char* buf)
{
	// static int i = 100; // for f= e_3rhs0_5_0_3
	//static int i = 10;  // for f = e_3rhs0_5_0_5
	static int i = 1; // for f = e_3rhs0_5_0_5 for inner and f= e_3rhs0_5_0_3 for outer 
	sprintf(buf, "Wave_u2_%03d.inp", i++);
}

//=================================================================

void WavesOptional::ApplyAvsOutElastic(int k, Vec& u21, Vec& u22, Vec& ucommon, double* vOuter)
{

	if (USE_FEM || EXCHANGE)
	{
		char tmp1[256];
		char tmp2[256];
		char tmp3[256];
		int noeq = 1;

		getFileNameElastic1(tmp1);
		getFileNameElastic2(tmp2);
		getFileName(tmp3);

		if (nsd == 2)
		{

			ierr = writeInp2D(tmp1, &gg, &u21, noeq, 0);
			CHKERRA(ierr);
			ierr = writeInp2D(tmp2, &gg, &u22, noeq, 0);
			CHKERRA(ierr);
			ierr = writeInp2D(tmp3, &gg, &ucommon, noeq, 0);
			CHKERRA(ierr);
		}
		else
		{
			ierr = writeInp(tmp1, &gg, &u21, noeq);
			CHKERRA(ierr);
		}
	}

	//=======================================

	AVSOutputOp avs2d(&outerWithHole, ofs);

	char tmp[256];
	int noeq = 1;
	ofstream ofs;

	if (USE_FDM || EXCHANGE)
	{
		if (EXCHANGE)
		{
			char tmp_fdm[256];
			int noeq = 1;

			getFileNameFDM(tmp_fdm);

			ierr = writeInpFDM2D(tmp_fdm, &sdg, vOuter, noeq, 0);
		}

	}

}

//==================================================================

void WavesOptional::ApplyAvsOut(int k, Vec& u2)
{

	char tmp[256];
	int noeq = 1;

	if (USE_FEM)
	{
		getFileName(tmp);
		if (nsd == 2)
		{
			ierr = writeInp2D(tmp, &gg, &u2, noeq, 0);
			CHKERRA(ierr);
		}
		else
		{
			ierr = writeInp(tmp, &gg, &u2, noeq);
			CHKERRA(ierr);
		}
	}

	if (EXCHANGE)
	{
		getFileName(tmp);
		if (nsd == 2)
		{
			ierr = writeInp2D(tmp, &gg, &u2, noeq, 0);
			CHKERRA(ierr);
		}
		else
		{
			ierr = writeInp(tmp, &gg, &u2, noeq);
			CHKERRA(ierr);
		}
	}

}

void WavesOptional::ApplyAvsOutElType(int sch, Vec& sol, int EL_TYPE)
{

	char tmp[256];
	int noeq = 1;
	getFileName(tmp);

	// here type of material is 1 - outer mesh 

	if (nsd == 3)
		ierr = writeInpElType(tmp, &gg, &sol, noeq, 1, EL_TYPE);
	else
		ierr = writeInpElType2D(tmp, &gg, &sol, noeq, 1, EL_TYPE);
	CHKERRA(ierr);

}

//====================================================================
// write results file for Gid-postprocessor
// input parametrs: name of the file (example: file.flavia.res)
//                   results, number of iteration k
//  Values of the parameter "results" should be:
//                         1 - if results is scalar ( outputs are
//                             real values at each node of the grid,
//                         format:  i   result(i) )
//                         2  - if results is vector, format:
//                         i   x(i) y(i) z(i), 
//                       where x(i),y(i),z(i) are results at the each
//                       component (x,y,z) at the node i.
//                         3 - if results is matrix , output formats
//                         i  M_xx(i) M_yy(i) M_zz(i) M_xy(i) M_yz(i) M_xz(i),
//                         where i is global node and M_... are values
// of the partial derivatives _xx,_yy, ... at this node. Used for stresses.
//=========================================================================== 
bool WavesOptional::Inp_to_RES_Gid(Vec& u2, char *inp_file, int results, int k)
{
	FILE *fp;

	int i, j;
	int el;

	Mat_real displ;

	PetscScalar *u_array;
	ierr = VecGetArray(u2, &u_array);

	if (results == 1)
		for (i = 0; i < gg.getNoNodes(); i++)
		{
			displ.newsize(gg.getNoNodes(), 1);

			displ(i, 0) = u_array[i];

		}
	else if (results == 2)
	{
		displ.newsize(gg.getNoNodes(), 3);

		displ(i, 0) = u_array[i];
		displ(i, 1) = u_array[i];
		displ(i, 2) = u_array[i];
	}
	else if (results == 3)
	{
		displ.newsize(gg.getNoNodes(), 6);

		displ(i, 0) = u_array[i];
		displ(i, 1) = u_array[i];
		displ(i, 2) = u_array[i];
		displ(i, 3) = u_array[i];
		displ(i, 4) = u_array[i];
		displ(i, 5) = u_array[i];
	}

	fp = fopen(inp_file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << inp_file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	// write 1 line to gid file

	ierr = VecRestoreArray(u2, &u_array);
	CHKERRQ( ierr);

	if (results == 1)
	{
		fprintf(fp, " Scalar_result   1    %i       1   1   0 \n", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f \n", i + 1, displ(i, 0));
		}
	}
	else if (results == 2)
	{
		fprintf(fp, " Vector_result   1    %i       2   1   0 \n", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %f\n", i + 1, displ(i, 0), displ(i, 1), displ(i, 2));
		}
	}
	else if (results == 3)
	{
		fprintf(fp, " Matrix_result   1    %i       3   1   0 \n", k);

		for (i = 0; i < gg.getNoNodes(); i++)
		{
			fprintf(fp, "%i %f %f %f %f %f %f\n", i + 1, displ(i, 0), displ(i, 1), displ(i, 2), displ(i, 3), displ(i, 4), displ(i, 5));
		}
	}

	fprintf(fp, "\n");

	fclose(fp);

//	printf("Solution written to file: %s\n", inp_file);

	return true;

}

//========================================================================

int WavesOptional::writeMTV2D(char *file, Grid *grid, Vec *u, int nsys)
{
	return 0;
}

//=======================================================================
int WavesOptional::writeMTV4D(char *file, Grid *grid, Vec *u, int nsys)
{
	return 0;
}

real WavesOptional::computeMinElmSize(Grid& grid)
{
	int e;
	int nel = grid.getNoElms();
	nsd = grid.getNoSpaceDim();
	if (nsd == 2)
	{
		FET3n2D F(&grid);
		F.refill(0);
		real minh = F.h();
		cout << "minelmsize" << minh << endl;
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			minh = min(minh, F.h());
		}
		return fabs(minh);
	}
	else if (nsd == 3)
	{
		FET4n3D F(&grid);
		F.refill(0);
		real minh = F.h();
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			minh = min(minh, F.h());
		}
		return minh;
	}
	else
		cout << " computeGridSize: dimension " << nsd << " not implemented\n" << flush;
	return -999.0;
}
real WavesOptional::computeMaxElmSize(Grid& grid)
{
	int e;
	int nel = grid.getNoElms();
	nsd = grid.getNoSpaceDim();
	if (nsd == 2)
	{
		FET3n2D F(&grid);
		F.refill(0);
		real maxh = F.h();
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			maxh = max(maxh, F.h());
		}
		return maxh;
	}
	else if (nsd == 3)
	{
		FET4n3D F(&grid);
		F.refill(0);
		real maxh = F.h();
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			maxh = max(maxh, F.h());
		}
		return maxh;
	}
	else
		cout << " computemaxGridSize: dimension " << nsd << " not implemented\n" << flush;
	return -999.0;
}
//==========================================================================  
void WavesOptional::MTVOpen(ofstream &ofs_)
{
	static int i = 1;
	char buf[20];
	sprintf(buf, "MTVOut%3d.mtv", i++);
}

void WavesOptional::writeMatlabDataInPnt(MV_Vector<real>& solFEM, MV_Vector<real>& solFDM, MV_Vector<real>& timedata)
{
	ofstream outp;
	int i;

	outp.open("hyb.m");
	int notimesteps = solFEM.size();
	for (i = 0; i < notimesteps; i++)
	{
		outp << solFEM(i) << "  " << solFDM(i) << "   " << timedata(i) << "\n";
	}
	outp.close();

	//print("Solution at one point is written to file: %s\n", "hyb.m");
}

void WavesOptional::writeSolInMatlabFile(Vec& u1, Vec& u2)
{
	ofstream outp;
	int i;

	outp.open("sol.m");

	int nno = gg.getNoNodes();

	PetscScalar* u2D;
	ierr = VecGetArray(u2,&u2D);
	CHKERRA(ierr);
	PetscScalar* u1D;
	ierr = VecGetArray(u1,&u1D);
	CHKERRA(ierr);

	for (i = 0; i < nno; i++)
	{
		outp << u2D[i] << "  " << u1D[i] << "   " << i << "\n";
	}
	outp.close();

	ierr = VecRestoreArray(u2,&u2D);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1,&u1D);
	CHKERRA(ierr);

	//print("Solution at one point is written to file: %s\n", "sol.m");
}

void WavesOptional::sortFDM(int *ia, const int &n)
{
	int swapDone;
	int tmp;

	do
	{
		swapDone = 0;
		for (int i = 1; i < n; i++)
			if (ia[i] < ia[i - 1])
			{
				swapDone = 1;
				tmp = ia[i];
				ia[i] = ia[i - 1];
				ia[i - 1] = tmp;
			}
	}
	while (swapDone);
}

void WavesOptional::getMTVName(char* buf)
{
	// static int i = 100; // for f= e_3rhs0_5_0_3
	//static int i = 10;  // for f = e_3rhs0_5_0_5
	static int i = 1; // for f = e_3rhs0_5_0_5 for inner and f= e_3rhs0_5_0_3 for outer 
	sprintf(buf, "MTVInn%03d.mtv", i++);
}

double WavesOptional::InitTime()
{

	nsd = sdg.getNoSpaceDim();
	//initialization time 

	if (nsd == 2)
	{
		dt = maxtime / nrSTEPS;
		cout << "dt is " << dt << endl;
	}
	else
	{
		dt = maxtime / nrSTEPS;
		cout << "dt is " << dt << endl;
	}
	return dt;
}

void WavesOptional::InitExchangeFEM(Grid& grid, Grid& grid_outer)
{

	//------ initialization  for exchange information -  
	// Index arrays for exchange
	// exchangeMask is array with codes 0,1,2 in the nodes with global num.
	// code 2 for outer boundary nodes,
	// code 1 for innerst layer of the boundary layer
	// code 0 for inner nodes, which don't belong to boundary and innerst nodes 
	// of the boundary. 

	makeExchangeNodes(grid, over1, exchangeMask);

	WavesOutputs write_mask;

	write_mask.WriteToFile((char *) "BoundNodes.m", exchangeMask);

//	cout << "before  makeSortedNodesIndex(exchangeMask,u1nodes,1)" << endl;

	makeSortedNodesIndex(exchangeMask, u1nodes, 1);

	makeSortedNodesIndex(exchangeMask, u2nodes, 2);

	innerIndexArray1 = new int[u1nodes.size()];
	innerIndexArray2 = new int[u2nodes.size()];
	int k;

//	cout << " FEM:innerIndexArray1 size =" << u1nodes.size() << endl;
	for (k = 0; k < u1nodes.size(); k++)
	{
		innerIndexArray1[k] = u1nodes(k);
//		cout << innerIndexArray1[k] << endl;
	}
//	cout << " FEM: innerIndexArray2 size =" << u2nodes.size() << endl;
	for (k = 0; k < u2nodes.size(); k++)
	{
		innerIndexArray2[k] = u2nodes(k);
//		cout << innerIndexArray2[k] << endl;
	}

}

void WavesOptional::InitExchangeConvergFEM(Grid& grid)
{

	//------ initialization  for exchange information -  
	// Index arrays for exchange
	// exchangeMask is array with codes 0,1,2 in the nodes with global num.
	// code 2 for outer boundary nodes,
	// code 1 for innerst layer of the boundary layer
	// code 0 for inner nodes, which don't belong to boundary and innerst nodes 
	// of the boundary. 

	makeExchangeNodes(grid, over1, exchangeMask);

	cout << "before  makeSortedNodesIndex(exchangeMask,u0nodes,0)" << endl;

	makeSortedNodesIndex(exchangeMask, u0nodes, 0);

	innerIndexArray0 = new int[u0nodes.size()];

	int k;

	cout << " FEM:innerIndexArray0 size =" << u0nodes.size() << endl;
	for (k = 0; k < u0nodes.size(); k++)
	{
		innerIndexArray0[k] = u0nodes(k);
		cout << innerIndexArray0[k] << endl;
	}

}

//=========================================================================

void WavesOptional::InitExchangeFEMfromFile(Grid& grid, Grid& grid_outer)
{

	//------ initialization  for exchange information -  
	// Index arrays for exchange
	// exchangeMask is array with codes 0,1,2 in the nodes with global num.
	// code 2 for outer boundary nodes,
	// code 1 for innerst layer of the boundary layer
	// code 0 for inner nodes, which don't belong to boundary and innerst nodes 
	// of the boundary. 

	exchangeMask.newsize(grid.getNoNodes());

	WavesOutputs read_mask;
	read_mask.ReadSolfromFile((char *) "BoundNodes.m", exchangeMask);
//	cout << "  exchangeMask.size() " << exchangeMask.size() << endl;

//	cout << "before  makeSortedNodesIndex(exchangeMask,u1nodes,1)" << endl;

	makeSortedNodesIndex(exchangeMask, u1nodes, 1);

	makeSortedNodesIndex(exchangeMask, u2nodes, 2);

	innerIndexArray1 = new int[u1nodes.size()];
	innerIndexArray2 = new int[u2nodes.size()];
	int k;

//	cout << " FEM:innerIndexArray1 size =" << u1nodes.size() << endl;
	for (k = 0; k < u1nodes.size(); k++)
	{
		innerIndexArray1[k] = u1nodes(k);
//		cout << innerIndexArray1[k] << endl;
	}
//	cout << " FEM: innerIndexArray2 size =" << u2nodes.size() << endl;
	for (k = 0; k < u2nodes.size(); k++)
	{
		innerIndexArray2[k] = u2nodes(k);
//		cout << innerIndexArray2[k] << endl;
	}

}

//================================
void WavesOptional::initExchangeFEMtoFDM(WavesSDIndexes& innerHole)
{

// to check  convergence of the hybrid  scheme:
	// to do sdindexes for FEM domain 
	// with codes 0,1 corresponding to code_inter
	// and condition9 (inner FEM domain)

	WavesSDMaskIndexes maskFEM(sdg, maska, code_inter);

	int nOfLoops_ = maskFEM.nOfLoopIndex();
	innerHole.copyLI(maskFEM, nOfLoops_);

	innerHole.presentWavesSDIndexes();

}
//================================

void WavesOptional::initExchangeFDM()

{
	int k;
	int curCoor;
	int nrOuterNodes = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	if (EXCHANGE)
	{
		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		WavesSDMaskIndexes maskWithHole(sdg, maska, code12);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

//		cout << "FDM:outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
//		cout << "FDM:outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
//		cout << "FEM: u1nodes.size()" << u1nodes.size() << endl;
//		cout << "FEM: u2nodes.size()" << u2nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());
	} //for EXCHANGE

	if (USE_FEM)
	{
//		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		// for FEM with absorbing b.c.: include boundary
		WavesSDMaskIndexes maskWithHole(sdg, maska, code6);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

		cout << "outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
		cout << "outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
		cout << "u1nodes.size()" << u1nodes.size() << endl;
		cout << "u2nodes.size()" << u2nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());
	} //for USE_FEM

//	cout << "after return outerwithhole" << endl;
}

//======================  exchange for structured FEM mesh ====================

void WavesOptional::initExchangeStructFDM()

{
	int k;
	int curCoor;
	int nrOuterNodes = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();
	MV_Vector<int> out_maska;
	out_maska.newsize(nrOuterNodes);

	if (EXCHANGE)
	{
//		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
//		for (k = 0; k < exchangeMask.size(); k++)
//			cout << exchangeMask(k) << endl;

		//new version for struct. meshes: works faster, then maskaFDM
		maska = maskaStructFDM(exchangeMask);

		for (k = 0; k < sdg.getNoNodes(); k++)
			out_maska(k) = maska[k];

		WavesOutputs write_mask;
		write_mask.WriteToFile((char *) "exchangeMask.m", out_maska);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		WavesSDMaskIndexes maskWithHole(sdg, maska, code12);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

//		cout << " outerWithHole " << endl;
		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

/*		cout << "FDM:outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
		cout << "FDM:outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
		cout << "FEM: u1nodes.size()" << u1nodes.size() << endl;
		cout << "FEM: u2nodes.size()" << u2nodes.size() << endl;
*/
		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());

	} //for EXCHANGE

	if (USE_FEM)
	{
/*		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;
*/
		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		// for FEM with absorbing b.c.: include boundary
		WavesSDMaskIndexes maskWithHole(sdg, maska, code6);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

/*		cout << "outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
		cout << "outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
		cout << "u1nodes.size()" << u1nodes.size() << endl;
		cout << "u2nodes.size()" << u2nodes.size() << endl;
*/
		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());
	} //for USE_FEM

//	cout << "after return outerwithhole" << endl;
}

void WavesOptional::initExchangeMaskFDM()

{
	int k;
	int curCoor;
	int nrOuterNodes = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	MV_Vector<int> out_maska;
	out_maska.newsize(nrOuterNodes);

	maska = new int[nrOuterNodes];

	if (EXCHANGE)
	{

		WavesOutputs read_mask;
		read_mask.ReadSolfromFile((char *) "exchangeMask.m", out_maska);

		for (k = 0; k < sdg.getNoNodes(); k++)
			maska[k] = out_maska(k);

		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		WavesSDMaskIndexes maskWithHole(sdg, maska, code12);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

//		cout << " outerWithHole " << endl;
		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

//		cout << "FDM:outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;//
//		cout << "FDM:outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
//		cout << "FEM: u1nodes.size()" << u1nodes.size() << endl;
//		cout << "FEM: u2nodes.size()" << u2nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());

	} //for EXCHANGE

	if (USE_FEM)
	{
		//cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		//WavesSDMaskIndexes maskWithHole(sdg,maska,code12);

		// for FEM with absorbing b.c.: include boundary
		WavesSDMaskIndexes maskWithHole(sdg, maska, code6);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

		//cout << "outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
		//cout << "outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
		//cout << "u1nodes.size()" << u1nodes.size() << endl;
		//cout << "u2nodes.size()" << u2nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());

	} //for USE_FEM

//	cout << "after return outerwithhole" << endl;
}

//=================================================================0

void WavesOptional::initExchangeFDM(WavesSDIndexes& innerHole)

{
	int k;
	int curCoor;
	int nrOuterNodes = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	if (EXCHANGE)
	{
		//cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

//new version for struct. meshes
		maska = maskaStructFDM(exchangeMask);

		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		WavesSDMaskIndexes maskWithHole(sdg, maska, code12);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

// to check  convergence of the hybrid  scheme:
		// to do sdindexes for FEM domain 
		// with code 0 corresponding to code_inter
		// and condition9 (inner FEM domain)

		WavesSDMaskIndexes maskFEM(sdg, maska, code_inter);

		nOfLoops_ = maskFEM.nOfLoopIndex();
		innerHole.copyLI(maskFEM, nOfLoops_);

		innerHole.presentWavesSDIndexes();

		int size_ar0;

// Integer arrays for exchange
		size_ar0 = innerHole.nOfNodeIndex();

		outerIndexArray0 = new int[size_ar0];

		innerHole.makeIndexArray(outerIndexArray0);
		sortFDM(outerIndexArray0, innerHole.nOfNodeIndex());

//		cout << "FDM:outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
//		cout << "FDM:outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
//		cout << "FDM:innerHole.nOfNodeIndex()" << innerHole.nOfNodeIndex() << endl;
//		cout << "FEM: u1nodes.size()" << u1nodes.size() << endl;
//		cout << "FEM: u2nodes.size()" << u2nodes.size() << endl;
//		cout << "FEM: u0nodes.size()" << u0nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());
		assert( u0nodes.size()==innerHole.nOfNodeIndex());

	}
	if (USE_FEM)
	{
		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain

		WavesSDMaskIndexes maskBoundary1(sdg, maska, codeb1);
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		//WavesSDMaskIndexes maskWithHole(sdg,maska,code12);

		// for FEM with absorbing b.c.: include boundary
		WavesSDMaskIndexes maskWithHole(sdg, maska, code6);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		int nOfLoops = maskBoundary1.nOfLoopIndex();
		outerBoundary1.copyLI(maskBoundary1, nOfLoops);

		nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar1, size_ar2;

// Integer arrays for exchange
		size_ar1 = outerBoundary1.nOfNodeIndex();
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray1 = new int[size_ar1];
		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		outerBoundary1.makeIndexArray(outerIndexArray1);
		sortFDM(outerIndexArray1, outerBoundary1.nOfNodeIndex());

		cout << "outerBoundary1.nOfNodeIndex()" << outerBoundary1.nOfNodeIndex() << endl;
		cout << "outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;
		cout << "u1nodes.size()" << u1nodes.size() << endl;
		cout << "u2nodes.size()" << u2nodes.size() << endl;

		assert( u1nodes.size()==outerBoundary1.nOfNodeIndex());
		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());

	} //for USE_FEM
}

//===================================================================

void WavesOptional::ApplyExchangeFDMtoFEM(Vec& u2, real* uOuter)
{

	PetscScalar* uinner;

	ierr = VecGetArray(u2,&uinner);
	CHKERRA(ierr);

	for (int i = 0; i < outerBoundary2.nOfNodeIndex(); i++)
		cout << "before exchange: uOuter(" << outerIndexArray2[i] << " = " << uOuter[outerIndexArray2[i]] << " , uInner(" << innerIndexArray2[i] << ") = " << uinner[innerIndexArray2[i]] << endl;

	for (int i = 0; i < outerBoundary2.nOfNodeIndex(); i++)
	{
		uinner[innerIndexArray2[i]] = uOuter[outerIndexArray2[i]];
	}

	ierr = VecRestoreArray(u2,&uinner);
	CHKERRA(ierr);

}

void WavesOptional::ApplyExchangeFEMtoFDM(Vec& u2, real* uOuter)
{

	PetscScalar* uinner;

	ierr = VecGetArray(u2,&uinner);
	CHKERRA(ierr);

	for (int i = 0; i < outerBoundary2.nOfNodeIndex(); i++)
	{
		uOuter[outerIndexArray2[i]] = uinner[innerIndexArray2[i]];
	}

	ierr = VecRestoreArray(u2,&uinner);
	CHKERRA(ierr);

}

//=======================================================================
// ======== exchange only in one direction: from FDM to FEM
//========================================================================

void WavesOptional::initExchangeFDMtoFEM()

{
	int k;
	int curCoor;
	int nrOuterNodes = sdg.getNoNodes();
	int nsd = sdg.getNoSpaceDim();

	if (EXCHANGE)
	{
		cout << "before overlapping, size of exchangeMask  = " << exchangeMask.size() << endl;
		for (k = 0; k < exchangeMask.size(); k++)
			cout << exchangeMask(k) << endl;

		maska = maskaFDM(exchangeMask);
		//loopindexes  for outer domain
		WavesSDMaskIndexes maskBoundary2(sdg, maska, codeb2);

		// to make sdindexes for domain with hole 
		// and including boundary, works in 2D
		//  WavesSDMaskIndexes maskWithHole(sdg,maska,code6);

		// to make sdindexes for domain with hole 
		// and without boundary in 2D and 3D

		WavesSDMaskIndexes maskWithHole(sdg, maska, code12);

		int nOfLoops_ = maskWithHole.nOfLoopIndex();
		outerWithHole.copyLI(maskWithHole, nOfLoops_);

		outerWithHole.presentWavesSDIndexes();

		// int nOfLoops = maskBoundary1.nOfLoopIndex();
		// outerBoundary1.copyLI(maskBoundary1,nOfLoops);

		int nOfLoops = maskBoundary2.nOfLoopIndex();
		outerBoundary2.copyLI(maskBoundary2, nOfLoops);

		int size_ar2;

// Integer arrays for exchange
		size_ar2 = outerBoundary2.nOfNodeIndex();

		outerIndexArray2 = new int[size_ar2];

		outerBoundary2.makeIndexArray(outerIndexArray2);
		sortFDM(outerIndexArray2, outerBoundary2.nOfNodeIndex());

		cout << "FDM:outerBoundary2.nOfNodeIndex()" << outerBoundary2.nOfNodeIndex() << endl;

		cout << "FEM: u2nodes.size()" << u2nodes.size() << endl;

		assert( u2nodes.size()==outerBoundary2.nOfNodeIndex());
	} //for EXCHANGE

}
//=======================================================================

void WavesOptional::ApplyExchange(Vec& u2, real* uOuter)
{

//	cout << "exchange if not system" << endl;

	PetscScalar* uinner;

	ierr = VecGetArray(u2,&uinner);
	CHKERRA(ierr);

	for (int i = 0; i < outerBoundary1.nOfNodeIndex(); i++)
	{
		uOuter[outerIndexArray1[i]] = uinner[innerIndexArray1[i]];
	}

	for (int i = 0; i < outerBoundary2.nOfNodeIndex(); i++)
	{
		uinner[innerIndexArray2[i]] = uOuter[outerIndexArray2[i]];
	}

	ierr = VecRestoreArray(u2,&uinner);
	CHKERRA(ierr);

}

void WavesOptional::ApplyExchangeConverg(Vec& u2, real* uOuter, WavesSDIndexes& innerHole)
{

//	cout << "exchange if not system" << endl;

	PetscScalar* uinner;

	ierr = VecGetArray(u2,&uinner);
	CHKERRA(ierr);

	for (int i = 0; i < innerHole.nOfNodeIndex(); i++)
		cout << "code 0 before exchange: uOuter(" << outerIndexArray0[i] << " = " << uOuter[outerIndexArray0[i]] << " , uInner(" << innerIndexArray0[i] << ") = " << uinner[innerIndexArray0[i]] << endl;

	for (int i = 0; i < innerHole.nOfNodeIndex(); i++)
	{
		uOuter[outerIndexArray0[i]] = uinner[innerIndexArray0[i]];
		cout << "code 0 after exchange: uOuter(" << outerIndexArray0[i] << " = " << uOuter[outerIndexArray0[i]] << " , uInner(" << innerIndexArray0[i] << ") = " << uinner[innerIndexArray0[i]] << endl;
	}

	ierr = VecRestoreArray(u2,&uinner);
	CHKERRA(ierr);

}

void WavesOptional::ApplyPlotMTVOut(int k, real* vOuter, Vec& u2)
{
	PlotMTVOutputOp mtv2D(&outerWithHole, ofs);
	PlotMTVOutputOp mtv2DFDM(&sdg, ofs);
	double t = k * dt;
	char tmp[256];
	int noeq = 1;
	cout << "Time is " << t << " nrSteps is " << k << endl;

	if (EXCHANGE)
	{
		MTVOpen(ofs);
		mtv2D.apply(vOuter);
		ofs.close();
	}

	if (USE_FDM)
	{
		MTVOpen(ofs);
		mtv2DFDM.apply(vOuter);
		ofs.close();
	}

	if (EXCHANGE || USE_FEM)
	{
		getMTVName(tmp);
		if (nsd == 2)
		{
			ierr = writeMTV2D(tmp, &gg, &u2, noeq);
			CHKERRA(ierr);
		}
		else
		{
			ierr = writeMTV4D(tmp, &gg, &u2, noeq);
			CHKERRA(ierr);
		}
	}

}

void WavesOptional::PlotMTVFDMGrid()
{
	PlotMTVOutputOp mtv2D(&outerWithHole, ofs);
	MTVOpen(ofs);
	mtv2D.writeMTV2D();
	ofs.close();
}

void WavesOptional::ComputeL2Norm(double* func_to_norm, double& l2norm)
{

	int nel = gg.getNoElms();
	int i, j, node;
	double val1;
	double val_new = 0.0;
	double val_new1 = 0.0;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (i = 0; i < nel; i++)
		{
			if (gg.getElementType(i) == ELMTRI1)
			{
				int ne = gg.getNoNodesInElm(i);
				F.refill(i);
				real volume = F.area();
				val1 = 0.0;

				for (j = 0; j < ne; j++)
				{

					node = gg.loc2glob(i, j);
					double nv1 = func_to_norm[node];
					val1 += nv1 * nv1;

				}
				val_new1 += (val1 / ne) * volume;
			}
		}
		cout << "func*func*S(el) = " << val_new1 << endl;
		l2norm = sqrt(val_new1);
		cout << "L2 norm is  " << l2norm << endl;
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (i = 0; i < nel; i++)
		{
			if (gg.getElementType(i) == ELMTET1)
			{
				int ne = gg.getNoNodesInElm(i);
				F.refill(i);
				real volume = fabs(F.volume());
				val1 = 0.0;

				for (j = 0; j < ne; j++)
				{
					node = gg.loc2glob(i, j);
					double nv1 = func_to_norm[node];
					val1 += nv1 * nv1;
				}

				val_new1 += (val1 / ne) * volume;

			}
		}
		l2norm = sqrt(val_new1);

	}
	cout << "L2 norm is" << l2norm << endl;

	for (int sch = 0; sch < nno; sch++)
	{
		func_to_norm[sch] = func_to_norm[sch] / l2norm;
	}

}

double WavesOptional::L2Norma(Grid& gg, MV_Vector<double>& Gradient)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	double sol = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				double val1 = Gradient(n_1) * Gradient(n_1);
				double val2 = Gradient(n_2) * Gradient(n_2);
				double val3 = Gradient(n_3) * Gradient(n_3);

				cout << " val1 " << val1 << "val2 " << val2 << " val3 " << val3 << endl;

				sol += ((val1 + val2 + val3) / 3.0) * volume;

			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				n_4 = F.n4();

				volume *= 0.25;

				double val1 = Gradient(n_1) * Gradient(n_1);
				double val2 = Gradient(n_2) * Gradient(n_2);
				double val3 = Gradient(n_3) * Gradient(n_3);
				double val4 = Gradient(n_4) * Gradient(n_4);

				sol += (val1 + val2 + val3 + val4) * volume;

			}
		}
	}

	sol = sqrt(sol);
	return (sol);

}

double WavesOptional::Norma(Grid& gg, MV_Vector<double>& Gradient)
{
	int i, n, el;
	int ierr;
	int nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	int n_1, n_2, n_3, n_4;

	double sol = 0.0;

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				double val1 = Gradient(n_1);
				double val2 = Gradient(n_2);
				double val3 = Gradient(n_3);

				cout << " val1 " << val1 << "val2 " << val2 << " val3 " << val3 << endl;

				sol += ((val1 + val2 + val3) / 3.0) * volume;

			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);
		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = fabs(F.volume());
				n_1 = F.n1();
				n_2 = F.n2();
				n_3 = F.n3();
				n_4 = F.n4();

				volume *= 0.25;

				double val1 = Gradient(n_1);
				double val2 = Gradient(n_2);
				double val3 = Gradient(n_3);
				double val4 = Gradient(n_4);

				sol += (val1 + val2 + val3 + val4) * volume;

			}
		}
	}

	return (sol);

}
