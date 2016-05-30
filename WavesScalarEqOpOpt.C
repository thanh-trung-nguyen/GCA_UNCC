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

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "include/wavesFEcomp.h"
#include "include/waveskraftwerk.hh"
#include <sys/times.h>
#include "include/wavesScalarEqOpOpt.h"

#if PETSC_USE_DEBUG
#define CHKERRA(e) if(e){PetscError(__LINE__,__FUNCT__,__FILE__,__SDIR__,e,0,0);exit(1);}
#else
#define CHKERRA(e); /* obsolete in petsc 2.1 */
#endif 
#define SETERRA(n,p,s)  SETERRQ(n,s)  /* obsolete in petsc 2.1 */

typedef double real;
typedef MV_ColMat<real> Mat_real;

WavesScalarEqOpOpt::WavesScalarEqOpOpt(Grid &gg_, Grid &outer_gg_, WavesSDGeometry &sdg_, int &nsd, bool& EXCHANGE_, bool& USE_FEM_, bool& USE_FDM_, bool &USE_RHS_, bool &USE_DIRICHLET_FEM_, bool &USE_DIRICHLET_FDM_, bool &USE_ABSORB_, bool &PRINT_FILES_, int &nrSTEPS_,
		double &maxtime_, double &rhs_, int &type_of_material_, double &velocity_, double &guess_velocity_, bool &doIncludeCorners_) :
		opt(sdg_, outer_gg_, gg_, EXCHANGE_, USE_FEM_, USE_FDM_, maxtime_), outerWithHole(sdg_), outerBoundary1(sdg_), outerBoundary2(sdg_), initialize2DTimeFunction(&sdg_, e_3rhs0_5_0_5), initialize3DTimeFunction(&sdg_, D3rhs0_5_0_5),
				initialize3D_1(&sdg_, D3rhs0_1_0_1), initrhs1(&sdg_, rhs1), initrhs2(&sdg_, rhs2), add_sol_and_function(&sdg_), avs2dFDM(&sdg_, ofs), dirichletBC(&sdg_), leftBoundary(sdg_, SDiLow, doIncludeCorners_), rightBoundary(sdg_, SDiHigh, doIncludeCorners_),
				lowBoundary(sdg_, SDjLow, doIncludeCorners_), topBoundary(sdg_, SDjHigh, doIncludeCorners_), perBoundary(sdg_, SDkLow, doIncludeCorners_), backBoundary(sdg_, SDkHigh, doIncludeCorners_), corners(sdg_, SDCorners, doIncludeCorners_), cornerBC(&corners),
				abs1(&sdg_, SDiLow, dt), abs2(&sdg_, SDjLow, dt), abs3(&sdg_, SDiHigh, dt), abs4(&sdg_, SDjHigh, dt), abs5(&sdg_, SDkLow, dt), abs6(&sdg_, SDkHigh, dt)
{

	sdg = sdg_;
	gg = gg_;

	outer_gg = outer_gg_;
	doIncludeCorners = doIncludeCorners_;
	USE_FEM = USE_FEM_;
	USE_FDM = USE_FDM_;
	USE_RHS = USE_RHS_;

	EXCHANGE = EXCHANGE_;

	USE_DIRICHLET_FEM = USE_DIRICHLET_FEM_;
	USE_DIRICHLET_FDM = USE_DIRICHLET_FDM_;
	USE_ABSORB = USE_ABSORB_;
	PRINT_FILES = PRINT_FILES_;
	nrSTEPS = nrSTEPS_;
	maxtime = maxtime_;
	rhs = rhs_;
	type_of_material = type_of_material_;

	// here velocity = ro in equation !!!

	velocity = velocity_;
	guess_velocity = guess_velocity_;

	bvalues_initialised = 0;

	nbn = 0;

	if (USE_FDM)
		nno = sdg.getNoNodes();

	if (EXCHANGE)
		nno = gg.getNoNodes();

	if (USE_FEM)
		nno = gg.getNoNodes();

}

WavesScalarEqOpOpt::~WavesScalarEqOpOpt()
{
	if (bvalues_initialised)
	{

		delete[] boundindex;
		delete[] bound_values;

	}

	if (FEM_initialised)
	{

//		ierr = MatDestroy(A);
//		CHKERRA(ierr);
//		delete[] boundindex;
//		delete[] bound_values;
//		delete[] bvalues;
	}

	if (EXCHANGE || USE_FDM)
	{
		delete[] u_Outer;
		delete[] v_Outer;

		if (USE_RHS)
		{
			delete[] fn_Outer;

		}

	}

}

double WavesScalarEqOpOpt::InitTime()
{

	nsd = gg.getNoSpaceDim();

	//initialization time 

	if (nsd == 2)
	{
		dt = maxtime / nrSTEPS;
		//cout << "dt is " << dt << endl;
	}
	else
	{
		dt = maxtime / nrSTEPS;
		//cout << "dt is " << dt << endl;

	}
	return dt;
}

const int DEBUG = 1;

int WavesScalarEqOpOpt::dropZeroWavesElements(Mat& Ain, Mat& Aout, double droptol)
{
	int i, j;
	int numbercols;
	int nrows;
	// first count the nonzero elements
	PetscFunctionBegin;
	ierr = MatGetSize(Ain, &nrows, &numbercols);
	CHKERRA(ierr);
	int* number_of_nbs = new int[nrows];
	const int *cols;
	int ncols;
	const PetscScalar *vals;
	//cout << "Size of matrix is " << nrows << " rows and " << numbercols << " columns\n";
	int counter = 0;
	int nocounter = 0;
	for (i = 0; i < nrows; i++)
	{
		ierr = MatGetRow(Ain, i, &ncols, &cols, &vals);
		CHKERRA(ierr);
		int loccounter = 0;
		// ncols - number of nonzero el-s in the row with number i
		// cols - number of the column
		// vals - nonzero values 

		for (j = 0; j < ncols; j++)
		{
			if (fabs(vals[j]) > droptol)
			{
				counter++;
				loccounter++;
			}
			else
			{
				nocounter++;
			}
		}
		number_of_nbs[i] = loccounter;
		ierr = MatRestoreRow(Ain, i, &ncols, &cols, &vals);
		CHKERRA(ierr);
	}
	//cout << "Found " << counter << " nonzero and " << nocounter << " zero elements ( < 1e-12)\n" << flush;

	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nrows, numbercols, 0, number_of_nbs, &Aout);
	CHKERRA(ierr);
	delete[] number_of_nbs;

	Array1dReal resultElMat(50);
	Array1dInt rowIdx(1);
	Array1dInt colIdx(50);

	for (i = 0; i < nrows; i++)
	{
		ierr = MatGetRow(Ain, i, &ncols, &cols, &vals);
		CHKERRA(ierr);
		int loccounter = 0;
		rowIdx(0) = i;
		for (j = 0; j < ncols; j++)
		{
			if (fabs(vals[j]) > droptol)
			{
				colIdx(loccounter) = cols[j];
				resultElMat(loccounter) = vals[j];
				loccounter++;
			}
		}
		ierr = MatRestoreRow(Ain, i, &ncols, &cols, &vals);
		CHKERRA(ierr);

		ierr = MatSetValues(Aout, 1, &rowIdx(0), loccounter, &colIdx(0), &resultElMat(0), ADD_VALUES);
	}
	ierr = MatAssemblyBegin(Aout, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(Aout, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	PetscFunctionReturn(0);

	return 0;
}

void WavesScalarEqOpOpt::getTriangleQuadrature(const int noQpts, MV_ColMat<real>& points, MV_Vector<real>& weights)
{
	// should also check that the size is correct??
	if (noQpts == 3)
	{
		weights = 1.0 / 3.0;

		points(0, 0) = 1.0;
		points(0, 1) = 0.0;
		points(0, 2) = 0.0;

		points(1, 0) = 0.0;
		points(1, 1) = 1.0;
		points(1, 2) = 0.0;

		points(2, 0) = 0.0;
		points(2, 1) = 0.0;
		points(2, 2) = 1.0;

	}
	else if (noQpts == 6)
	{
		weights(0) = 0.223381589678011;
		weights(1) = 0.223381589678011;
		weights(2) = 0.223381589678011;
		weights(3) = 0.109951743655322;
		weights(4) = 0.109951743655322;
		weights(5) = 0.109951743655322;

		points(0, 0) = 0.108103018168070;
		points(0, 1) = 0.445948490915965;
		points(0, 2) = 0.445948490915965;

		points(1, 0) = 0.445948490915965;
		points(1, 1) = 0.445948490915965;
		points(1, 2) = 0.108103018168070;

		points(2, 0) = 0.445948490915965;
		points(2, 1) = 0.108103018168070;
		points(2, 2) = 0.445948490915965;

		points(3, 0) = 0.816847572980459;
		points(3, 1) = 0.091576213509771;
		points(3, 2) = 0.091576213509771;

		points(4, 0) = 0.091576213509771;
		points(4, 1) = 0.091576213509771;
		points(4, 2) = 0.816847572980459;

		points(5, 0) = 0.091576213509771;
		points(5, 1) = 0.816847572980459;
		points(5, 2) = 0.091576213509771;
	}
	else if (noQpts == 12)
	{

		weights(0) = 0.116786275726379;
		weights(1) = 0.116786275726379;
		weights(2) = 0.116786275726379;
		weights(3) = 0.050844906370207;
		weights(4) = 0.050844906370207;
		weights(5) = 0.050844906370207;
		weights(6) = 0.082851075618374;
		weights(7) = 0.082851075618374;
		weights(8) = 0.082851075618374;
		weights(9) = 0.082851075618374;
		weights(10) = 0.082851075618374;
		weights(11) = 0.082851075618374;

		points(0, 0) = 0.501426509658179;
		points(0, 1) = 0.249286745170910;
		points(0, 2) = 0.249286745170910;

		points(1, 0) = 0.249286745170910;
		points(1, 1) = 0.249286745170910;
		points(1, 2) = 0.501426509658179;

		points(2, 0) = 0.249286745170910;
		points(2, 1) = 0.501426509658179;
		points(2, 2) = 0.249286745170910;

		points(3, 0) = 0.873821971016996;
		points(3, 1) = 0.063089014491502;
		points(3, 2) = 0.063089014491502;

		points(4, 0) = 0.063089014491502;
		points(4, 1) = 0.063089014491502;
		points(4, 2) = 0.873821971016996;

		points(5, 0) = 0.063089014491502;
		points(5, 1) = 0.873821971016996;
		points(5, 2) = 0.063089014491502;

		points(6, 0) = 0.053145049844817;
		points(6, 1) = 0.310352451033784;
		points(6, 2) = 0.636502499121399;

		points(7, 0) = 0.310352451033784;
		points(7, 1) = 0.636502499121399;
		points(7, 2) = 0.053145049844817;

		points(8, 0) = 0.636502499121399;
		points(8, 1) = 0.053145049844817;
		points(8, 2) = 0.310352451033784;

		points(9, 0) = 0.636502499121399;
		points(9, 1) = 0.310352451033784;
		points(9, 2) = 0.053145049844817;

		points(10, 0) = 0.053145049844817;
		points(10, 1) = 0.636502499121399;
		points(10, 2) = 0.310352451033784;

		points(11, 0) = 0.310352451033784;
		points(11, 1) = 0.053145049844817;
		points(11, 2) = 0.636502499121399;
	}
}

void WavesScalarEqOpOpt::getTetraQuadrature(const int noQpts, MV_ColMat<real>& points, MV_Vector<real>& weights)

{
	// should also check that the size is correct??
	if (noQpts == 1)
	{
		points = 0.25;
		weights = 1.0 / 6.0;
	}
	else if (noQpts == 4)
	{
		points(0, 0) = 1.0;
		points(0, 1) = 0.0;
		points(0, 2) = 0.0;
		points(0, 3) = 0.0;

		points(1, 0) = 0.0;
		points(1, 1) = 1.0;
		points(1, 2) = 0.0;
		points(1, 3) = 0.0;

		points(2, 0) = 0.0;
		points(2, 1) = 0.0;
		points(2, 2) = 1.0;
		points(2, 3) = 0.0;

		points(3, 0) = 0.0;
		points(3, 1) = 0.0;
		points(3, 2) = 0.0;
		points(3, 3) = 1.0;

		weights = .041666666666666667;
	}
	else if (noQpts == 12)
	{

		points(0, 0) = 8.525816240380812e-02;
		points(0, 1) = 3.517084860955024e-01;
		points(0, 2) = 2.113248654051870e-01;
		points(1, 0) = 3.517084860955024e-01;
		points(1, 1) = 3.517084860955024e-01;
		points(1, 2) = 2.113248654051870e-01;
		points(2, 0) = 3.517084860955024e-01;
		points(2, 1) = 8.525816240380812e-02;
		points(2, 2) = 2.113248654051870e-01;
		points(3, 0) = 6.442273695638098e-01;
		points(3, 1) = 7.222388251550198e-02;
		points(3, 2) = 2.113248654051870e-01;
		points(4, 0) = 7.222388251550198e-02;
		points(4, 1) = 7.222388251550198e-02;
		points(4, 2) = 2.113248654051870e-01;
		points(5, 0) = 7.222388251550198e-02;
		points(5, 1) = 6.442273695638098e-01;
		points(5, 2) = 2.113248654051870e-01;
		points(6, 0) = 2.284485576426188e-02;
		points(6, 1) = 9.424000482046258e-02;
		points(6, 2) = 7.886751345948130e-01;
		points(7, 0) = 9.424000482046258e-02;
		points(7, 1) = 9.424000482046258e-02;
		points(7, 2) = 7.886751345948130e-01;
		points(8, 0) = 9.424000482046258e-02;
		points(8, 1) = 2.284485576426188e-02;
		points(8, 2) = 7.886751345948130e-01;
		points(9, 0) = 1.726202034166492e-01;
		points(9, 1) = 1.935233099426903e-02;
		points(9, 2) = 7.886751345948130e-01;
		points(10, 0) = 1.935233099426903e-02;
		points(10, 1) = 1.935233099426903e-02;
		points(10, 2) = 7.886751345948130e-01;
		points(11, 0) = 1.935233099426903e-02;
		points(11, 1) = 1.726202034166492e-01;
		points(11, 2) = 7.886751345948130e-01;
		int i;
		for (i = 0; i < noQpts; i++)
			points(i, 3) = 1.0 - points(i, 0) - points(i, 1) - points(i, 2);

		weights(0) = 3.473631008974336e-02;
		weights(1) = 3.473631008974336e-02;
		weights(2) = 3.473631008974336e-02;
		weights(3) = 1.709772890426878e-02;
		weights(4) = 1.709772890426878e-02;
		weights(5) = 1.709772890426878e-02;
		weights(6) = 2.493954856591809e-03;
		weights(7) = 2.493954856591809e-03;
		weights(8) = 2.493954856591809e-03;
		weights(9) = 1.227561704951555e-03;
		weights(10) = 1.227561704951555e-03;
		weights(11) = 1.227561704951555e-03;
	}
	else if (noQpts == 36)
	{

		points(0, 0) = 4.449149069543934e-01;
		points(0, 1) = 2.211917138331736e-01;
		points(0, 2) = 1.127016653792585e-01;

		points(1, 0) = 2.211917138331736e-01;
		points(1, 1) = 2.211917138331736e-01;
		points(1, 2) = 1.127016653792585e-01;

		points(2, 0) = 2.211917138331736e-01;
		points(2, 1) = 4.449149069543934e-01;
		points(2, 2) = 1.127016653792585e-01;

		points(3, 0) = 7.753407796383943e-01;
		points(3, 1) = 5.597877749117355e-02;
		points(3, 2) = 1.127016653792585e-01;

		points(4, 0) = 5.597877749117355e-02;
		points(4, 1) = 5.597877749117355e-02;
		points(4, 2) = 1.127016653792585e-01;

		points(5, 0) = 5.597877749117355e-02;
		points(5, 1) = 7.753407796383943e-01;
		points(5, 2) = 1.127016653792585e-01;

		points(6, 0) = 4.715551422064242e-02;
		points(6, 1) = 2.753752129477418e-01;
		points(6, 2) = 1.127016653792585e-01;

		points(7, 0) = 2.753752129477418e-01;
		points(7, 1) = 5.647676074523573e-01;
		points(7, 2) = 1.127016653792585e-01;

		points(8, 0) = 5.647676074523573e-01;
		points(8, 1) = 4.715551422064242e-02;
		points(8, 2) = 1.127016653792585e-01;

		points(9, 0) = 5.647676074523573e-01;
		points(9, 1) = 2.753752129477418e-01;
		points(9, 2) = 1.127016653792585e-01;

		points(10, 0) = 4.715551422064242e-02;
		points(10, 1) = 5.647676074523573e-01;
		points(10, 2) = 1.127016653792585e-01;

		points(11, 0) = 2.753752129477418e-01;
		points(11, 1) = 4.715551422064242e-02;
		points(11, 2) = 1.127016653792585e-01;

		points(12, 0) = 2.507132548290895e-01;
		points(12, 1) = 1.246433725854550e-01;
		points(12, 2) = 5.000000000000000e-01;

		points(13, 0) = 1.246433725854550e-01;
		points(13, 1) = 1.246433725854550e-01;
		points(13, 2) = 5.000000000000000e-01;

		points(14, 0) = 1.246433725854550e-01;
		points(14, 1) = 2.507132548290895e-01;
		points(14, 2) = 5.000000000000000e-01;

		points(15, 0) = 4.369109855084980e-01;
		points(15, 1) = 3.154450724575100e-02;
		points(15, 2) = 5.000000000000000e-01;

		points(16, 0) = 3.154450724575100e-02;
		points(16, 1) = 3.154450724575100e-02;
		points(16, 2) = 5.000000000000000e-01;

		points(17, 0) = 3.154450724575100e-02;
		points(17, 1) = 4.369109855084980e-01;
		points(17, 2) = 5.000000000000000e-01;

		points(18, 0) = 2.657252492240850e-02;
		points(18, 1) = 1.551762255168920e-01;
		points(18, 2) = 5.000000000000000e-01;

		points(19, 0) = 1.551762255168920e-01;
		points(19, 1) = 3.182512495606995e-01;
		points(19, 2) = 5.000000000000000e-01;

		points(20, 0) = 3.182512495606995e-01;
		points(20, 1) = 2.657252492240850e-02;
		points(20, 2) = 5.000000000000000e-01;

		points(21, 0) = 3.182512495606995e-01;
		points(21, 1) = 1.551762255168920e-01;
		points(21, 2) = 5.000000000000000e-01;

		points(22, 0) = 2.657252492240850e-02;
		points(22, 1) = 3.182512495606995e-01;
		points(22, 2) = 5.000000000000000e-01;

		points(23, 0) = 1.551762255168920e-01;
		points(23, 1) = 2.657252492240850e-02;
		points(23, 2) = 5.000000000000000e-01;

		points(24, 0) = 5.651160270378563e-02;
		points(24, 1) = 2.809503133773639e-02;
		points(24, 2) = 8.872983346207415e-01;

		points(25, 0) = 2.809503133773639e-02;
		points(25, 1) = 2.809503133773639e-02;
		points(25, 2) = 8.872983346207415e-01;

		points(26, 0) = 2.809503133773639e-02;
		points(26, 1) = 5.651160270378563e-02;
		points(26, 2) = 8.872983346207415e-01;

		points(27, 0) = 9.848119137860162e-02;
		points(27, 1) = 7.110237000328450e-03;
		points(27, 2) = 8.872983346207415e-01;

		points(28, 0) = 7.110237000328450e-03;
		points(28, 1) = 7.110237000328450e-03;
		points(28, 2) = 8.872983346207415e-01;

		points(29, 0) = 7.110237000328450e-03;
		points(29, 1) = 9.848119137860162e-02;
		points(29, 2) = 8.872983346207415e-01;

		points(30, 0) = 5.989535624174581e-03;
		points(30, 1) = 3.497723808604224e-02;
		points(30, 2) = 8.872983346207415e-01;

		points(31, 0) = 3.497723808604224e-02;
		points(31, 1) = 7.173489166904171e-02;
		points(31, 2) = 8.872983346207415e-01;

		points(32, 0) = 7.173489166904171e-02;
		points(32, 1) = 5.989535624174581e-03;
		points(32, 2) = 8.872983346207415e-01;

		points(33, 0) = 7.173489166904171e-02;
		points(33, 1) = 3.497723808604224e-02;
		points(33, 2) = 8.872983346207415e-01;

		points(34, 0) = 5.989535624174581e-03;
		points(34, 1) = 7.173489166904171e-02;
		points(34, 2) = 8.872983346207415e-01;

		points(35, 0) = 3.497723808604224e-02;
		points(35, 1) = 5.989535624174581e-03;
		points(35, 2) = 8.872983346207415e-01;

		weights(0) = 1.277022783138014e-02;
		weights(1) = 1.277022783138014e-02;
		weights(2) = 1.277022783138014e-02;
		weights(3) = 5.559737515168268e-03;
		weights(4) = 5.559737515168268e-03;
		weights(5) = 5.559737515168268e-03;
		weights(6) = 9.059515813317084e-03;
		weights(7) = 9.059515813317084e-03;
		weights(8) = 9.059515813317084e-03;
		weights(9) = 9.059515813317084e-03;
		weights(10) = 9.059515813317084e-03;
		weights(11) = 9.059515813317084e-03;
		weights(12) = 6.488126429243278e-03;
		weights(13) = 6.488126429243278e-03;
		weights(14) = 6.488126429243278e-03;
		weights(15) = 2.824717020567056e-03;
		weights(16) = 2.824717020567056e-03;
		weights(17) = 2.824717020567056e-03;
		weights(18) = 4.602837534354112e-03;
		weights(19) = 4.602837534354112e-03;
		weights(20) = 4.602837534354112e-03;
		weights(21) = 4.602837534354112e-03;
		weights(22) = 4.602837534354112e-03;
		weights(23) = 4.602837534354112e-03;
		weights(24) = 2.060250271064264e-04;
		weights(25) = 2.060250271064264e-04;
		weights(26) = 2.060250271064264e-04;
		weights(27) = 8.969652596584605e-05;
		weights(28) = 8.969652596584605e-05;
		weights(29) = 8.969652596584605e-05;
		weights(30) = 1.461592553911414e-04;
		weights(31) = 1.461592553911414e-04;
		weights(32) = 1.461592553911414e-04;
		weights(33) = 1.461592553911414e-04;
		weights(34) = 1.461592553911414e-04;
		weights(35) = 1.461592553911414e-04;

		int i;
		for (i = 0; i < noQpts; i++)
			points(i, 3) = 1.0 - points(i, 0) - points(i, 1) - points(i, 2);

	}

	else
		cout << "No quadrature points = " << noQpts << " not available!\n" << flush;
}

//=====================================================================
// assemble rhs like (f,fi_j)
//=====================================================================

void WavesScalarEqOpOpt::evalRHS(Vec& Fn_, Fcn3d rhs_)
{
	int n, e, q;
	nsd = gg.getNoSpaceDim();
	;
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();
	int noWavesBasisFcn = nsd + 1;
	int noQpts = nsd + 1;
	MV_ColMat<real> points(noQpts, nsd + 1);
	MV_Vector<real> weights(noQpts);
	if (nsd == 2)
		getTriangleQuadrature(noQpts, points, weights);
	else
		getTetraQuadrature(noQpts, points, weights);
	PetscScalar zero = 0.0;
	ierr = VecSet(Fn_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( Fn_, &tmpvec);
	CHKERRA(ierr);

	if (nsd == 2)
	{

		FET3n2D F(&gg);
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = F.area();
			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, 0.0);
				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
			}
		}
	}

	else if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = fabs(F.det());

			if (area < 0.0)
			{
				cout << " area is negative !!!" << endl;
				exit(1);
			}

			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2) + F.x4() * points(q, 3);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2) + F.y4() * points(q, 3);
				real z = F.z1() * points(q, 0) + F.z2() * points(q, 1) + F.z3() * points(q, 2) + F.z4() * points(q, 3);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, z);

				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
				tmpvec[F.n4()] += points(q, 3) * areaWeightRhs;
			}
		}
	}

	ierr = VecRestoreArray( Fn_, &tmpvec );
	CHKERRA(ierr);

}

//================================================================0

void WavesScalarEqOpOpt::applyPlaneWaveFEM(MV_Vector<int>& markNodes, double t, double y_fix)
{

	opt.findBoundaryNodes(gg, markNodes);

	int i;
	double eps = 0.5;
	bound_values = new double[nno];

	for (i = 0; i < nno; i++)
	{
		if (markNodes(i) == 1)
		{
			real x = gg.getCoor(i, 0);
			real y = gg.getCoor(i, 1);
			real z = gg.getCoor(i, 2);

			if (fabs(y - y_fix) < eps)
			{
				markNodes(i) = 5;
				bound_values[i] = planeWave(x, y, z, t);

				// code 5 for node, which is at the top boundary
			}
		}

	}
}

//==================== to find nodes at the lens boundary ============

//================================================================0

void WavesScalarEqOpOpt::findNodesLENS(MV_Vector<int>& markNodes)
{

	int i;
	double eps = 0.0000001;
	WavesOutputs out;
	int sch = 0;
	double norma = 0.0;

	for (i = 0; i < nno; i++)
	{

		real x = gg.getCoor(i, 0);
		real y = gg.getCoor(i, 1);
		real z = gg.getCoor(i, 2);

		if ((y <= 10.0 && y >= 9.53939) && (x >= -3.0 && x <= 3.0) && (z >= -3.0 && z <= 3.0) || (x == 0.0 && y == 10.0 && z == 0.0))

		{
			markNodes(i) = 5;
			norma = sqrt(4 * x * x + 4 * y * y + 4 * z * z);
			out.WriteMyPart((char *) "markNodes.dat", i + 1, x, y, z, sch, norma);
			sch++;
			// code 5 for node, which is at the top boundary
		}

	}
}
//==================== to find nodes at the earth boundary - structured nodes 
// ============
//================================================================0

void WavesScalarEqOpOpt::findNodesEarth(MV_Vector<int>& markNodes)
{

	int i;
	double eps = 0.0000001;
	WavesOutputs out;
	int sch = 0;
	double norma = 0.0;

	for (i = 0; i < nno; i++)
	{

		real x = gg.getCoor(i, 0);
		real y = gg.getCoor(i, 1);

		if ((y == 0.62) && (x >= 0.37 && x <= 0.63))

		{
			markNodes(i) = 5;
			norma = sqrt(4 * x * x + 4 * y * y);
			out.WriteMyPart((char *) "markNodes.dat", i + 1, x, y, 0.0, sch, norma);
			sch++;
			// code 5 for node, which is at the top boundary
		}

	}
}

//====== for fdm mesh ============================================

void WavesScalarEqOpOpt::findNodesEarthFDM(MV_Vector<int>& markNodes, int code)
{

	int i;
	double eps = 0.0000001;
	WavesOutputs out;
	int sch = 0;
	double norma = 0.0;

	for (i = 0; i < sdg.getNoNodes(); i++)
	{

		real x = sdg.getCoor(i, 0);
		real y = sdg.getCoor(i, 1);

// next line for tests with GLK method (see ellipse_OPT.dat)
//	 if ( (y ==  0.3 ) && (x >= 0.3 && x <= 0.7)  )
// next if is for convergence tests for wave equation (see testconv2D04.dat)
		if ((y == 0.6) && (x >= 0.2 && x <= 0.8))
		{
			markNodes(i) = 5;
			norma = sqrt(4 * x * x + 4 * y * y);

			if (code == 1)
				out.WriteMyPart((char *) "markNodes.dat", i + 1, x, y, 0.0, sch, norma);
			sch++;
			// code 5 for node, which is at the top boundary
		}

	}

}

//====== for fdm mesh ============================================

void WavesScalarEqOpOpt::findNodesPhotonFDMX(MV_Vector<int>& markNodes, double fix, double fix_min, double fix_max, int code)
{

	int i;
	double eps = 0.0000001;
	WavesOutputs out;
	int sch = 0;
	double norma = 0.0;

	for (i = 0; i < sdg.getNoNodes(); i++)
	{

		real x = sdg.getCoor(i, 0);
		real y = sdg.getCoor(i, 1);

// next line for tests with GLK method (see ellipse_OPT.dat)
//	 if ( (y ==  0.3 ) && (x >= 0.3 && x <= 0.7)  )
// next if is for convergence tests for wave equation (see testconv2D04.dat)
		//  	 if ( (y == -3.0 ) && (x >= -3.0 && x <= 3.0)  )
		if ((y == fix) && (x >= fix_min && x <= fix_max))
		{
			markNodes(i) = 5;
			norma = sqrt(4 * x * x + 4 * y * y);

			if (code == 1)
				out.WriteMyPart((char *) "markNodes.dat", i + 1, x, y, 0.0, sch, norma);
			sch++;
			// code 5 for node, which is at the top boundary
		}

	}

}

void WavesScalarEqOpOpt::findNodesPhotonFDMY(MV_Vector<int>& markNodes, double fix, double fix_min, double fix_max, int code)
{

	int i;
	double eps = 0.0000001;
	WavesOutputs out;
	int sch = 0;
	double norma = 0.0;

	for (i = 0; i < sdg.getNoNodes(); i++)
	{

		real x = sdg.getCoor(i, 0);
		real y = sdg.getCoor(i, 1);

		if ((x == fix) && (y >= fix_min && y <= fix_max))
		{
			markNodes(i) = 5;
			norma = sqrt(4 * x * x + 4 * y * y);

			if (code == 1)
				out.WriteMyPart((char *) "markNodes.dat", i + 1, x, y, 0.0, sch, norma);
			sch++;
			// code 5 for node, which is at the top boundary
		}

	}

}

//======================================================================
void WavesScalarEqOpOpt::evalRHSforDifMat(Vec& Fn_, double t, Fcn3dtime rhs_, Vec& lumpedAmass_, int type_of_material, double velocity)
{
	int n, e, q;
	//  int  ierr;
	nsd = gg.getNoSpaceDim();
	;
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();
	int noWavesBasisFcn = nsd + 1;
	int noQpts = nsd + 1;
	MV_ColMat<real> points(noQpts, nsd + 1);
	MV_Vector<real> weights(noQpts);
	if (nsd == 2)
		getTriangleQuadrature(noQpts, points, weights);
	else
		getTetraQuadrature(noQpts, points, weights);
	PetscScalar zero = 0.0;
	ierr = VecSet(Fn_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( Fn_, &tmpvec);
	CHKERRA(ierr);
	PetscScalar *Amassinv;
	ierr = VecGetArray( lumpedAmass_, &Amassinv);
	CHKERRA(ierr);

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = F.area();
			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, 0.0, t);
				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
			}
		}
	}

	else if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = F.det();

			if (area < 0.0)
			{
				cout << " volume of element " << e << " is negative ! " << area << endl;
				exit(1);
			}

			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2) + F.x4() * points(q, 3);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2) + F.y4() * points(q, 3);
				real z = F.z1() * points(q, 0) + F.z2() * points(q, 1) + F.z3() * points(q, 2) + F.z4() * points(q, 3);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, z, t);
				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
				tmpvec[F.n4()] += points(q, 3) * areaWeightRhs;
			}
		}
	}

	for (n = 0; n < nno; n++)
	{
		// we divide by dt� because in method
		// iterationmatrix  we multiplyed by -dt�

		tmpvec[n] *= -Amassinv[n] / (dt * dt);
	}

	MV_Vector<int> Markers(nno);
	Markers = 0;

	// type of material for element is type of mat in file -1
	// example: in file difmat.inp ellipse have tofmat = 2
	// then gg.getMaterialType(ellipse) = 2 -1 = 1
	// and 1 we should put in pulse2D.dat file

	for (int el = 0; el < nel; el++)
	{
		// for inner small ellipse and for outer ellipse
		if (gg.getMaterialType(el) == type_of_material)
		{
			if (nsd == 2)
			{
				int n1 = gg.loc2glob(el, 0);
				int n2 = gg.loc2glob(el, 1);
				int n3 = gg.loc2glob(el, 2);
				//to avoid repeating multiplications
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
			}
			if (nsd == 3)
			{
				int n1 = gg.loc2glob(el, 0);
				int n2 = gg.loc2glob(el, 1);
				int n3 = gg.loc2glob(el, 2);
				int n4 = gg.loc2glob(el, 3);

				//to avoid repeating multiplications
				Markers(n1) = 1;
				Markers(n2) = 1;
				Markers(n3) = 1;
				Markers(n4) = 1;

			}
		}
	}

	int sch = 0;

	int i;
	for (i = 0; i < nno; i++)
		if (Markers(i) == 1)
		{
			tmpvec[i] = tmpvec[i] / velocity;
			//cout << " tmpvec[i]" << tmpvec[i] << " velocity " << velocity << endl;
			sch++;
		}

	//cout << "chislo markerov" << sch++ << endl;

	ierr = VecRestoreArray( Fn_, &tmpvec );
	CHKERRA(ierr);
	ierr = VecRestoreArray( lumpedAmass_, &Amassinv );
	CHKERRA(ierr);

}

//======================================================================
void WavesScalarEqOpOpt::evalRHSforDifMat(Vec& Fn_, double t, Fcn3dtime rhs_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity)
{
	int n, e, q;
	//  int  ierr;
	nsd = gg.getNoSpaceDim();
	;
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();
	int noWavesBasisFcn = nsd + 1;
	int noQpts = nsd + 1;
	MV_ColMat<real> points(noQpts, nsd + 1);
	MV_Vector<real> weights(noQpts);
	if (nsd == 2)
		getTriangleQuadrature(noQpts, points, weights);
	else
		getTetraQuadrature(noQpts, points, weights);
	PetscScalar zero = 0.0;
	ierr = VecSet(Fn_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( Fn_, &tmpvec);
	CHKERRA(ierr);
	PetscScalar *Amassinv;
	ierr = VecGetArray( lumpedAmass_, &Amassinv);
	CHKERRA(ierr);

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = F.area();
			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, 0.0, t);
				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
			}
		}
	}

	else if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = fabs(F.det());
			for (q = 0; q < noQpts; q++)
			{
				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2) + F.x4() * points(q, 3);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2) + F.y4() * points(q, 3);
				real z = F.z1() * points(q, 0) + F.z2() * points(q, 1) + F.z3() * points(q, 2) + F.z4() * points(q, 3);

				real areaWeightRhs = area * weights(q) * rhs_(x, y, z, t);

				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
				tmpvec[F.n4()] += points(q, 3) * areaWeightRhs;
			}
		}
	}

	for (n = 0; n < nno; n++)
	{
		// we divide by dt� because in method
		// iterationmatrixfordivergence we multiplyed by -dt�

		tmpvec[n] *= -Amassinv[n] / (dt * dt);
	}

	for (int i = 0; i < nno; i++)
		tmpvec[i] = tmpvec[i] / velocity(i);

	ierr = VecRestoreArray( Fn_, &tmpvec );
	CHKERRA(ierr);
	ierr = VecRestoreArray( lumpedAmass_, &Amassinv );
	CHKERRA(ierr);

}

//======================================================================
void WavesScalarEqOpOpt::evalAdjRHSforDifMat(Vec& Fn_, double t, MV_Vector<double>& Dif_sol, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity)
{

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n, e, q;
	real areaWeightRhs;
	nsd = gg.getNoSpaceDim();
	;
	int nel = gg.getNoElms();
	int nno = gg.getNoNodes();
	int noWavesBasisFcn = nsd + 1;
	int noQpts = nsd + 1;
	double Adj_x, Adj_y, Adj_z;
	double eps = 1e-6;

	MV_ColMat<real> points(noQpts, nsd + 1);
	MV_Vector<real> weights(noQpts);
	if (nsd == 2)
		getTriangleQuadrature(noQpts, points, weights);
	else
		getTetraQuadrature(noQpts, points, weights);
	PetscScalar zero = 0.0;
	ierr = VecSet(Fn_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( Fn_, &tmpvec);
	CHKERRA(ierr);
	PetscScalar *Amassinv;
	ierr = VecGetArray( lumpedAmass_, &Amassinv);
	CHKERRA(ierr);

	if (nsd == 2)
	{
		FET3n2D F(&gg);
		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = F.area();
			bool found = false;

			for (q = 0; q < noQpts; q++)
			{

				real x = F.x1() * points(q, 0) + F.x2() * points(q, 1) + F.x3() * points(q, 2);
				real y = F.y1() * points(q, 0) + F.y2() * points(q, 1) + F.y3() * points(q, 2);

				int globnum = gg.loc2glob(e, q);
				Adj_x = gg.getCoor(globnum, 0);
				Adj_y = gg.getCoor(globnum, 1);

				if (fabs(Adj_x - x) < eps && fabs(Adj_y - y) < eps)
				{
					areaWeightRhs = area * weights(q) * Dif_sol(globnum);
					found = true;
					break;
				}assert(found);

				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
			}
		}
	}
	else if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (e = 0; e < nel; e++)
		{
			F.refill(e);
			real area = fabs(F.det());
			bool found = false;

			for (q = 0; q < noQpts; q++)
			{
				int globnum = gg.loc2glob(e, q);

				areaWeightRhs = area * weights(q) * Dif_sol(globnum);

				tmpvec[F.n1()] += points(q, 0) * areaWeightRhs;
				tmpvec[F.n2()] += points(q, 1) * areaWeightRhs;
				tmpvec[F.n3()] += points(q, 2) * areaWeightRhs;
				tmpvec[F.n4()] += points(q, 3) * areaWeightRhs;
			}
		}
	}

	for (n = 0; n < nno; n++)
	{
		// we divide by dt� because in method
		// iterationmatrixfordivergence we multiplyed by -dt�

		tmpvec[n] *= -Amassinv[n] / (dt * dt);
	}

	ierr = VecRestoreArray( Fn_, &tmpvec );
	CHKERRA(ierr);
	ierr = VecRestoreArray( lumpedAmass_, &Amassinv );
	CHKERRA(ierr);

}

//================================================================

// HRN
//#include "Waves/opt/petsc-2.1.0/src/mat/impls/aij/seq/aij.h"
//#include "Wavessrc/mat/impls/aij/seq/aij.h"
int WavesScalarEqOpOpt::MatMult0(Mat A_, Vec xx, Vec yy)
{
	ierr = MatMult(A_, xx, uhelp);
	ierr = VecAYPX(yy, -1.0, uhelp);

	return ierr;
}

void WavesScalarEqOpOpt::copyVal(const real *x, real *y, const int *IX, const int *IY, const int &n)
{

	for (int i = 0; i < n; i++)
		y[IY[i]] = x[IX[i]];

}

//======================================================================
void WavesScalarEqOpOpt::IterationMatrix2D(Mat& A_, Vec& lumpedAmass_, int type_of_material)
{
	//cout << "works iteration matrix for dif mat, line 1357" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_SELF, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_SELF, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 2)
	{
		FET3n2D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();

				resultElMat(0) = volume * F.a11();
				resultElMat(4) = volume * F.a22();
				resultElMat(8) = volume * F.a33();
				resultElMat(3) = resultElMat(1) = volume * F.a12();
				resultElMat(6) = resultElMat(2) = volume * F.a13();
				resultElMat(5) = resultElMat(7) = volume * F.a23();
				//cout << "============ el-t ============" << el << endl;

				//cout << " a_11 = " << F.a11() << " a_12 = " << F.a12() << " a_13 = " << F.a13() << endl;
				//cout << " a_21 = " << F.a21() << " a_22 = " << F.a22() << " a_23 = " << F.a23() << endl;
				//cout << " a_31 = " << F.a31() << " a_32 = " << F.a32() << " a_33 = " << F.a33() << endl;
				//cout << " area for element " << volume << "n1 = " << F.n1() << " n2 = " << F.n2() << " n3= " << F.n3() << endl;

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);

				CHKERRA(ierr);

				//============================================================
				volume *= 0.3333333333333333;

				if (gg.getMaterialType(el) == type_of_material)
				{
					tmpvec[F.n1()] += velocity * volume;
					tmpvec[F.n2()] += velocity * volume;
					tmpvec[F.n3()] += velocity * volume;
				}
				else
				{
					tmpvec[F.n1()] += volume;
					tmpvec[F.n2()] += volume;
					tmpvec[F.n3()] += volume;
				}

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	MV_Vector<int> Markers(nno);
	Markers = 0;

	if (type_of_material > 0)
	{
		for (int el = 0; el < nel; el++)
		{
			// for inner small ellipse and for outer ellipse
			if (gg.getMaterialType(el) == type_of_material)
			{
				if (nsd == 2)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
				}
				if (nsd == 3)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					int n4 = gg.loc2glob(el, 3);

					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
					Markers(n4) = 1;
				}

			}
		}
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}




//==================================================================
//  to compute parameter in the global convergent method in 2D case
//======================================================================
void  WavesScalarEqOpOpt::ParameterMatrix2D( Mat& A_, 
					     Vec& lumpedAmass_,
					     int type_of_material)
{
  //cout<<"works iteration matrix for dif mat, line 1727"<<endl;
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

  //  cout<<"before veccreate"<<endl;

	ierr = VecCreate(PETSC_COMM_SELF, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);
	
  PetscScalar zero=0.0;

  // cout<<"before vecset"<<endl;
  ierr = VecSet(lumpedAmass_,zero); CHKERRA(ierr);
  PetscScalar *tmpvec;
  ierr = VecGetArray(lumpedAmass_, &tmpvec); CHKERRA(ierr);//pointer to lAm

  ///=======================================================================
  
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  //cout<<"after  gridNeighbors.init (gg,false,true,false)"<<endl;
  
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  
  
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &A_);
  CHKERRA(ierr);
 
  // if nodes of con.  does't work, this is more slower variant:
  // ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
  // CHKERRA(ierr);
 
   //cout << "Assembling matrix 01." << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes);  
  Array1dInt rowAndColIdx(nobasenodes); 

  // cout<<"loop"<<endl;
  if(nsd==2)
    {
       FET3n2D F(&gg);
    
 
      for( el = 0; el < nel; el++ )
	{
	  if(gg.getElementType(el)==ELMTRI1)
	    {
	      F.refill(el);
	      real volume = F.area();
	      if (volume < 0.0)
		{
		  cout<<"element "<<el<<"have volume < 0 "<< volume<<endl;
		  exit(1);
		}

	      rowAndColIdx(0) = F.n1();
	      rowAndColIdx(1) = F.n2();
	      rowAndColIdx(2) = F.n3();
	   
          
	     resultElMat(0) = volume*F.a11();
	    resultElMat(4) = volume*F.a22();
	    resultElMat(8) = volume*F.a33();

	    resultElMat(3) = resultElMat(1) = volume*F.a12();
	    resultElMat(6) = resultElMat(2) = volume*F.a13();
	    resultElMat(5) = resultElMat(7) = volume*F.a23();

	   
	
	    /*
	   	    cout<<"============ el-t ============"<<el<<endl;

	    cout<<" a_11 = "<<F.a11()<<" a_12 = "<<F.a12()<<" a_13 = "<<F.a13()<<endl;
	    cout<<" a_21 = "<<F.a21()<<" a_22 = "<<F.a22()<<" a_23 = "<<F.a23()<<endl;
	    cout<<" a_31 = "<<F.a31()<<" a_32 = "<<F.a32()<<" a_33 = "<<F.a33()<<endl;
	    cout<<" area for element "<<volume<<"n1 = "<<F.n1()<<" n2 = "<<F.n2()<<" n3= "<<F.n3()<<endl;

	    */

	      ierr = MatSetValues(A_,nobasenodes,&rowAndColIdx(0),
				  nobasenodes,&rowAndColIdx(0),
				  &resultElMat(0),
				  ADD_VALUES); 

	    

	      CHKERRA(ierr);
	
	//============================================================
	volume *= 0.3333333333333333;


	tmpvec[F.n1()] += volume;
	tmpvec[F.n2()] += volume;
	tmpvec[F.n3()] += volume;
  

	}
      }
    }
 
 
  // cout<<"after loop"<<endl;

  ierr = MatAssemblyBegin(A_,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A_,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  /*

  for( i = 0; i < nno; i++){
    if( tmpvec[i] < +1e-6 ){
      //      cout<<"Warning: too low tmpvec["<<i<<"]=="<<tmpvec[i]<<endl;
    }
    //  here we do (M^L(-1))
    tmpvec[i] = 1.0/tmpvec[i];

  }

  */
  ierr = VecRestoreArray(lumpedAmass_, &tmpvec); CHKERRA(ierr);  

  //  ierr = MatDiagonalScale(A_,lumpedAmass_,PETSC_NULL); CHKERRA(ierr);

  /*
 PetscReal nrm;
	    ierr=MatNorm(A_,NORM_1,&nrm);
	    CHKERRA(ierr);

	    cout<<"Norma of matrice A_ "<<nrm<<endl;

	    exit(1);

  */
}

//========= to compute parameter in globally convergent method in 3D ======
// ===========  assembling of the stiffness matrix ========================
void  WavesScalarEqOpOpt::ParameterMatrix3D(Mat& A_, 
					    Vec& lumpedAmass_,
					    int type_of_material)
{
  int i,n,el;
  int ierr;
  nsd = gg.getNoSpaceDim();
  int nno = gg.getNoNodes();
  int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_SELF, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	// for many processors 
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

  //==========================================================================
  WavesNeighborFE& gridNeighbors = gg.getNeighbor ();
  gridNeighbors.init (gg,false,true,false);
  MV_Vector<int> number_of_nbs(nno);
  int count=0;
  for(i=0;i<nno;i++) // Node connecticivity
    number_of_nbs(count++) = (gridNeighbors.couplingsIrow (i+1) -
			      gridNeighbors.couplingsIrow (i));
  gridNeighbors.remove();
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, &number_of_nbs(0), &A);
  CHKERRA(ierr);

  //cout << "Assembling matrix 01." << endl;
  int nobasenodes = nsd+1;// for triangles and tetras
  Array1dReal resultElMat(nobasenodes*nobasenodes); 
  Array1dInt rowAndColIdx(nobasenodes); 

  //cout<<"loop"<<endl;
  if(nsd==2){
   
   // cout<<"call method  ParameterMatrix2D( Mat& A_,Vec& lumpedAmass_,int type_of_material)"<<endl;

    exit(1);

  } 
  else if(nsd==3){
    FET4n3D F(&gg);
    for( el = 0; el < nel; el++ ){
      if(gg.getElementType(el)==ELMTET1) {
	F.refill(el);
	real volume = fabs(F.volume());

        if (volume < 0.0)
	  {
	    cout<<"negative volume !!! "<<volume<<" in element "<<el<<endl;
            exit(1);

	  }
	rowAndColIdx(0) = F.n1();
	rowAndColIdx(1) = F.n2();
	rowAndColIdx(2) = F.n3();
	rowAndColIdx(3) = F.n4();
	resultElMat(0) = volume*F.a11();
	resultElMat(5) = volume*F.a22();
	resultElMat(10) = volume*F.a33();
	resultElMat(15) = volume*F.a44();
	resultElMat(4) = resultElMat(1) = volume*F.a12();
	resultElMat(8) = resultElMat(2) = volume*F.a13();
	resultElMat(12) = resultElMat(3) = volume*F.a14();
	resultElMat(9) = resultElMat(6) = volume*F.a23();
	resultElMat(13) = resultElMat(7) = volume*F.a24();
	resultElMat(14) = resultElMat(11) = volume*F.a34();
	ierr = MatSetValues(A,nobasenodes,&rowAndColIdx(0),
			    nobasenodes,&rowAndColIdx(0),&resultElMat(0),
			    ADD_VALUES);  
	CHKERRA(ierr);
	//============================================================
	volume *= 0.25;


	tmpvec[F.n1()] += volume;
	tmpvec[F.n2()] += volume;
	tmpvec[F.n3()] += volume;
  	tmpvec[F.n4()] += volume;
      }
    }
  }

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRA(ierr);
  
 ierr = VecRestoreArray(lumpedAmass_, &tmpvec); CHKERRA(ierr);  

}
  
//============= end of assembling in globally convergent method in 3D  ====



//======================================================================
void WavesScalarEqOpOpt::IterationMatrix2DDifMat(Mat& A_,
						 Vec& lumpedAmass_,
						 int type_of_material,
						 MV_Vector<double>& bb)
{
	//cout << "works iteration matrix for dif mat, line 1727" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_SELF, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	// for many processors 
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_SELF, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 2)
	{
		FET3n2D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();

				resultElMat(0) = volume * F.a11();
				resultElMat(4) = volume * F.a22();
				resultElMat(8) = volume * F.a33();
				resultElMat(3) = resultElMat(1) = volume * F.a12();
				resultElMat(6) = resultElMat(2) = volume * F.a13();
				resultElMat(5) = resultElMat(7) = volume * F.a23();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);

				CHKERRA(ierr);

				//============================================================
				volume *= 0.3333333333333333;

				tmpvec[F.n1()] += bb(F.n1()) * volume;
				tmpvec[F.n2()] += bb(F.n2()) * volume;
				tmpvec[F.n3()] += bb(F.n3()) * volume;

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{

		if (tmpvec[i] < +1e-6)
		{
			cout << "Warning: too low tmpvec[" << i << "]==" << tmpvec[i] << endl;
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	MV_Vector<int> Markers(nno);
	Markers = 0;

	if (type_of_material > 0)
	{
		for (int el = 0; el < nel; el++)
		{
			// for inner small ellipse and for outer ellipse
			if (gg.getMaterialType(el) == type_of_material)
			{
				if (nsd == 2)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
				}
				if (nsd == 3)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					int n4 = gg.loc2glob(el, 3);

					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
					Markers(n4) = 1;
				}

			}
		}
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}

//==================================================================
// start 
//======================================================================
void WavesScalarEqOpOpt::IterationMatrix2DGLK(Mat& A_, Vec& C_, Vec& lumpedAmass_, MV_Vector<double>& bb)
{
	//cout << "works iteration matrix for dif mat GLK" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_WORLD, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

//========================================================================
	ierr = VecCreate(PETSC_COMM_WORLD, &C_);
	CHKERRA(ierr);
	ierr = VecSetSizes(C_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(C_, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(C_, zero);
	CHKERRA(ierr);
	PetscScalar *vec_C;
	ierr = VecGetArray(C_, &vec_C);
	CHKERRA(ierr);
	//pointer to C_

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_WORLD, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 2)
	{
		FET3n2D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();

// assembling of stiffness matrix (grad phi_i, grad_phi_j)          

				if (gg.getMaterialType(el) == type_of_material)
				{
					resultElMat(0) = (1.0 / bb(F.n1())) * volume * F.a11();
					resultElMat(4) = (1.0 / bb(F.n2())) * volume * F.a22();
					resultElMat(8) = (1.0 / bb(F.n3())) * volume * F.a33();

					resultElMat(3) = (1.0 / bb(F.n1())) * volume * F.a12();
					resultElMat(6) = (1.0 / bb(F.n1())) * volume * F.a13();
					resultElMat(5) = (1.0 / bb(F.n3())) * volume * F.a23();

					resultElMat(1) = (1.0 / bb(F.n2())) * volume * F.a12();
					resultElMat(2) = (1.0 / bb(F.n3())) * volume * F.a13();
					resultElMat(7) = (1.0 / bb(F.n2())) * volume * F.a23();

				}
				else
				{
					resultElMat(0) = volume * F.a11();
					resultElMat(4) = volume * F.a22();
					resultElMat(8) = volume * F.a33();
					resultElMat(3) = resultElMat(1) = volume * F.a12();
					resultElMat(6) = resultElMat(2) = volume * F.a13();
					resultElMat(5) = resultElMat(7) = volume * F.a23();
				}
				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);

				CHKERRA(ierr);

				volume *= 0.3333333333333333;

				vec_C[F.n1()] += volume * F.glk_a11();
				vec_C[F.n2()] += volume * F.glk_a22();
				vec_C[F.n3()] += volume * F.glk_a33();

				if (gg.getMaterialType(el) == type_of_material)
				{
					tmpvec[F.n1()] += (1.0 / bb(F.n1())) * volume;
					tmpvec[F.n2()] += (1.0 / bb(F.n2())) * volume;
					tmpvec[F.n3()] += (1.0 / bb(F.n3())) * volume;

				}
				else
				{
					tmpvec[F.n1()] += volume;
					tmpvec[F.n2()] += volume;
					tmpvec[F.n3()] += volume;

				}

				tmpmat[F.n1()] += volume;
				tmpmat[F.n2()] += volume;
				tmpmat[F.n3()] += volume;

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = -dt2 / tmpmat[i];
	}

	for (n = 0; n < nno; n++)
	{

		vec_C[n] *= tmpmat[n];

		if (fabs(vec_C[n]) > +1e-20)
			cout << "tmpvec_C[" << n << "]=" << vec_C[n] << " bb " << bb(n) << "  tmpmat " << tmpmat[n] << endl;
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);
	ierr = VecRestoreArray(C_, &vec_C);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, lumpedAmass_, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}
//============================================================
// =======================  end ============

// =============== start new version glk===================0

//======================================================================
void WavesScalarEqOpOpt::IterationMatrix2DGLKnew(Mat& A_, Vec& C_, Vec& lumpedAmass_, MV_Vector<double>& bb)
{
	//cout << "works iteration matrix for dif mat GLK" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_WORLD, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

//========================================================================
	ierr = VecCreate(PETSC_COMM_WORLD, &C_);
	CHKERRA(ierr);
	ierr = VecSetSizes(C_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(C_, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(C_, zero);
	CHKERRA(ierr);
	PetscScalar *vec_C;
	ierr = VecGetArray(C_, &vec_C);
	CHKERRA(ierr);
	//pointer to C_

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_WORLD, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 2)
	{
		FET3n2D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTRI1)
			{
				F.refill(el);
				real volume = F.area();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();

// assembling of stiffness matrix (grad phi_i, grad_phi_j)          

				resultElMat(0) = volume * F.a11();
				resultElMat(4) = volume * F.a22();
				resultElMat(8) = volume * F.a33();
				resultElMat(3) = resultElMat(1) = volume * F.a12();
				resultElMat(6) = resultElMat(2) = volume * F.a13();
				resultElMat(5) = resultElMat(7) = volume * F.a23();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);

				CHKERRA(ierr);

				volume *= 0.3333333333333333;

				vec_C[F.n1()] += volume * F.glk_a11();
				vec_C[F.n2()] += volume * F.glk_a22();
				vec_C[F.n3()] += volume * F.glk_a33();

				tmpvec[F.n1()] += volume;
				tmpvec[F.n2()] += volume;
				tmpvec[F.n3()] += volume;

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	for (n = 0; n < nno; n++)
	{

		vec_C[n] = vec_C[n] * tmpmat[n];

		if (fabs(vec_C[n]) > +1e-20)
			cout << "tmpvec_C[" << n << "]=" << vec_C[n] << " bb " << bb(n) << "  tmpmat " << tmpmat[n] << endl;
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);
	ierr = VecRestoreArray(C_, &vec_C);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, lumpedAmass_, PETSC_NULL);
	CHKERRA(ierr);
	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}

//=========end of new version glk ============================
//======================================================================
void WavesScalarEqOpOpt::IterationMatrix3D(Mat& A_, Vec& lumpedAmass_, int type_of_material)
{
	//cout << "works iteration matrix for dif mat, line 1262" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	ierr = VecCreate(PETSC_COMM_WORLD, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_WORLD, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = F.volume();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();
				rowAndColIdx(3) = F.n4();

				resultElMat(0) = volume * F.a11();
				resultElMat(5) = volume * F.a22();
				resultElMat(10) = volume * F.a33();
				resultElMat(15) = volume * F.a44();
				resultElMat(4) = resultElMat(1) = volume * F.a12();
				resultElMat(8) = resultElMat(2) = volume * F.a13();
				resultElMat(12) = resultElMat(3) = volume * F.a14();
				resultElMat(9) = resultElMat(6) = volume * F.a23();
				resultElMat(13) = resultElMat(7) = volume * F.a24();
				resultElMat(14) = resultElMat(11) = volume * F.a34();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);
				CHKERRA(ierr);

				//============================================================

				volume *= 0.25;
				if (gg.getMaterialType(el) == type_of_material)
				{
					tmpvec[F.n1()] += velocity * volume;
					tmpvec[F.n2()] += velocity * volume;
					tmpvec[F.n3()] += velocity * volume;
					tmpvec[F.n4()] += velocity * volume;
				}
				else
				{
					tmpvec[F.n1()] += volume;
					tmpvec[F.n2()] += volume;
					tmpvec[F.n3()] += volume;
					tmpvec[F.n4()] += volume;
				}

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	MV_Vector<int> Markers(nno);
	Markers = 0;

	if (type_of_material > 0)
	{
		for (int el = 0; el < nel; el++)
		{
			cout << "type of mat for el-t" << el << "  " << gg.getMaterialType(el) << endl;
			// for inner small ellipse and for outer ellipse
			if (gg.getMaterialType(el) == type_of_material)
			{
				if (nsd == 2)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
				}
				if (nsd == 3)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					int n4 = gg.loc2glob(el, 3);

					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
					Markers(n4) = 1;
				}

			}
		}
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}

//======================================================================
void WavesScalarEqOpOpt::IterationMatrix3Dnew(Mat& A_, Vec& lumpedAmass_, int type_of_material)
{
	//cout << "works iteration matrix for dif mat, line 1262" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_WORLD, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_WORLD, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = F.volume();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();
				rowAndColIdx(3) = F.n4();

				resultElMat(0) = volume * F.a11();
				resultElMat(5) = volume * F.a22();
				resultElMat(10) = volume * F.a33();
				resultElMat(15) = volume * F.a44();
				resultElMat(4) = resultElMat(1) = volume * F.a12();
				resultElMat(8) = resultElMat(2) = volume * F.a13();
				resultElMat(12) = resultElMat(3) = volume * F.a14();
				resultElMat(9) = resultElMat(6) = volume * F.a23();
				resultElMat(13) = resultElMat(7) = volume * F.a24();
				resultElMat(14) = resultElMat(11) = volume * F.a34();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);
				CHKERRA(ierr);

				//============================================================

				volume *= 0.25;
				if (gg.getMaterialType(el) == type_of_material)
				{
					tmpvec[F.n1()] += velocity * volume;
					tmpvec[F.n2()] += velocity * volume;
					tmpvec[F.n3()] += velocity * volume;
					tmpvec[F.n4()] += velocity * volume;
				}
				else
				{
					tmpvec[F.n1()] += volume;
					tmpvec[F.n2()] += volume;
					tmpvec[F.n3()] += volume;
					tmpvec[F.n4()] += volume;
				}

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	MV_Vector<int> Markers(nno);
	Markers = 0;

	if (type_of_material > 0)
	{
		for (int el = 0; el < nel; el++)
		{
			//cout << "type of mat for el-t" << el << "  " << gg.getMaterialType(el) << endl;
			// mark all nodes which don't belong to lower cylinder with 1
			// for sphere1_lens2.inp ( lensref2_inc.dat):
			//   if ( gg.getMaterialType(el) == 5 ||  gg.getMaterialType(el) == 2 || gg.getMaterialType(el) == 3)
			// here == 1 - rest of the FEM domain
			//      == 2 - inclusion
			//      == 5 - lower cylinder without inclusion
			//      == 3  - upper cylinder
			// for new_test.inp (lensref3_inc.dat)  we should use 
			//   if ( gg.getMaterialType(el) == 2 ||  gg.getMaterialType(el) == 3 )

			if (gg.getMaterialType(el) == 2 || gg.getMaterialType(el) == 5)
			{
				if (nsd == 2)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
				}
				if (nsd == 3)
				{
					int n1 = gg.loc2glob(el, 0);
					int n2 = gg.loc2glob(el, 1);
					int n3 = gg.loc2glob(el, 2);
					int n4 = gg.loc2glob(el, 3);

					//to avoid repeating multiplications
					Markers(n1) = 1;
					Markers(n2) = 1;
					Markers(n3) = 1;
					Markers(n4) = 1;
				}

			}
		}
	}

	WavesOutputs out;
	out.WriteToFile((char *) "Markers.m", Markers);

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
}
//************************************************************************

//======================================================================
void WavesScalarEqOpOpt::IterationMatrix3DDifMatnew(Mat& A_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity)
{
	//cout << "works iteration matrix for dif mat, line 2039" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();

	ierr = VecCreate(PETSC_COMM_SELF, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;

	ierr = VecCreate(PETSC_COMM_SELF, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	//  ierr = VecSetType(Materials,VECSEQ);    CHKERRA(ierr);
	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	//cout << "after  gridNeighbors.init (gg,false,true,false)" << endl;

	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));

	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);

	// if nodes of con.  does't work, this is more slower variant:
	// ierr = MatCreateSeqAIJ(PETSC_COMM_SELF,nno,nno,0, PETSC_NULL, &A_);
	// CHKERRA(ierr);

	//cout << "Assembling matrix 01." << endl;
	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = F.volume();
				if (volume < 0.0)
				{
					cout << "element " << el << "have volume < 0 " << volume << endl;

				}

				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();
				rowAndColIdx(3) = F.n4();

				resultElMat(0) = volume * F.a11();
				resultElMat(5) = volume * F.a22();
				resultElMat(10) = volume * F.a33();
				resultElMat(15) = volume * F.a44();
				resultElMat(4) = resultElMat(1) = volume * F.a12();
				resultElMat(8) = resultElMat(2) = volume * F.a13();
				resultElMat(12) = resultElMat(3) = volume * F.a14();
				resultElMat(9) = resultElMat(6) = volume * F.a23();
				resultElMat(13) = resultElMat(7) = volume * F.a24();
				resultElMat(14) = resultElMat(11) = volume * F.a34();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);
				CHKERRA(ierr);

				//============================================================

				volume *= 0.25;

				tmpvec[F.n1()] += velocity(F.n1()) * volume;
				tmpvec[F.n2()] += velocity(F.n2()) * volume;
				tmpvec[F.n3()] += velocity(F.n3()) * volume;
				tmpvec[F.n4()] += velocity(F.n4()) * volume;
				if (tmpvec[F.n1()] < 0.0 || tmpvec[F.n2()] < 0.0 || tmpvec[F.n3()] < 0.0 || tmpvec[F.n4()] < 0.0)
					cout << " tmpvec" << tmpvec[F.n1()] << "  " << tmpvec[F.n2()] << "  " << tmpvec[F.n3()] << "  " << tmpvec[F.n4()] << endl;
			}
		}
	}

	//cout << "after loop" << endl;

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}

		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];

	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);
	//cout << " at the end IterationMatrix3DDifMatnew " << endl;
}

//************************************************************************
//======================================================================
void WavesScalarEqOpOpt::IterationMatrix3DDifMat(Mat& A_, Vec& lumpedAmass_, int type_of_material, MV_Vector<double>& velocity)
{
	//cout << "works iteration matrix for dif mat" << endl;
	int i, n, el;
	int ierr;
	nsd = gg.getNoSpaceDim();
	int nno = gg.getNoNodes();
	int nel = gg.getNoElms();
	ierr = VecCreate(PETSC_COMM_WORLD, &lumpedAmass_);
	CHKERRA(ierr);
	ierr = VecSetSizes(lumpedAmass_, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecSetType(lumpedAmass_, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;
	ierr = VecSet(lumpedAmass_, zero);
	CHKERRA(ierr);
	PetscScalar *tmpvec;
	ierr = VecGetArray( lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	//pointer to lAm

	//=======================================================================
	Vec Materials;
	ierr = VecCreate(PETSC_COMM_WORLD, &Materials);
	CHKERRA(ierr);
	ierr = VecSetSizes(Materials, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(Materials, VECMPI);
	CHKERRA(ierr);

	ierr = VecSet(Materials, zero);
	CHKERRA(ierr);
	PetscScalar *tmpmat;
	ierr = VecGetArray(Materials, &tmpmat);
	CHKERRA(ierr);
	//pointer 

	///=======================================================================

	WavesNeighborFE& gridNeighbors = gg.getNeighbor();
	gridNeighbors.init(gg, false, true, false);
	MV_Vector<int> number_of_nbs(nno);
	int count = 0;

	for (i = 0; i < nno; i++) // Node connecticivity
		number_of_nbs(count++) = (gridNeighbors.couplingsIrow(i + 1) - gridNeighbors.couplingsIrow(i));
	gridNeighbors.remove();
	ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, nno, nno, 0, &number_of_nbs(0), &A_);
	CHKERRA(ierr);


	int nobasenodes = nsd + 1; // for triangles and tetras
	Array1dReal resultElMat(nobasenodes * nobasenodes);

	Array1dInt rowAndColIdx(nobasenodes);

	if (nsd == 3)
	{
		FET4n3D F(&gg);

		for (el = 0; el < nel; el++)
		{
			if (gg.getElementType(el) == ELMTET1)
			{
				F.refill(el);
				real volume = F.volume();
				rowAndColIdx(0) = F.n1();
				rowAndColIdx(1) = F.n2();
				rowAndColIdx(2) = F.n3();
				rowAndColIdx(3) = F.n4();

				resultElMat(0) = volume * F.a11();
				resultElMat(5) = volume * F.a22();
				resultElMat(10) = volume * F.a33();
				resultElMat(15) = volume * F.a44();
				resultElMat(4) = resultElMat(1) = volume * F.a12();
				resultElMat(8) = resultElMat(2) = volume * F.a13();
				resultElMat(12) = resultElMat(3) = volume * F.a14();
				resultElMat(9) = resultElMat(6) = volume * F.a23();
				resultElMat(13) = resultElMat(7) = volume * F.a24();
				resultElMat(14) = resultElMat(11) = volume * F.a34();

				ierr = MatSetValues(A_, nobasenodes, &rowAndColIdx(0), nobasenodes, &rowAndColIdx(0), &resultElMat(0), ADD_VALUES);
				CHKERRA(ierr);

				//============================================================

				volume *= 0.25;

				tmpvec[F.n1()] += volume;
				tmpvec[F.n2()] += volume;
				tmpvec[F.n3()] += volume;
				tmpvec[F.n4()] += volume;

			}
		}
	}

	ierr = MatAssemblyBegin(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);
	ierr = MatAssemblyEnd(A_, MAT_FINAL_ASSEMBLY);
	CHKERRA(ierr);

	double dt2 = dt * dt;

	for (i = 0; i < nno; i++)
	{
		if (tmpvec[i] < +1e-6)
		{
		}
		tmpvec[i] = -dt2 / tmpvec[i];

		tmpmat[i] = tmpvec[i];
	}

	for (i = 0; i < nno; i++)
	{

		tmpmat[i] = tmpmat[i] / velocity(i);

		//cout << "tmpmat(" << i << ") = " << tmpmat[i] << " velocity " << velocity(i) << endl;
	}

	ierr = VecRestoreArray(lumpedAmass_, &tmpvec);
	CHKERRA(ierr);
	ierr = VecRestoreArray(Materials, &tmpmat);
	CHKERRA(ierr);

	ierr = MatDiagonalScale(A_, Materials, PETSC_NULL);
	CHKERRA(ierr);

	PetscScalar two = 2.0;
	ierr = MatShift(A_, two);
	CHKERRA(ierr);

}

//==================================================================

void WavesScalarEqOpOpt::ApplyExchangeCommon()
{
	opt.ApplyExchange(u2, u_Outer);

}

void WavesScalarEqOpOpt::ApplyExchange()
{
	opt.ApplyExchangeFDMtoFEM(u2, u_Outer);

}

void WavesScalarEqOpOpt::ApplyExchange1(WavesScalarEqOpOpt& p)
{

	//cout << "exchange from sdindexes" << endl;

	p.opt.ApplyExchange(u2, u_Outer);

}

void WavesScalarEqOpOpt::ApplyExchange_LENS()
{
	opt.ApplyExchange(u2_x, u_x_Outer);
	opt.ApplyExchange(u2_y, u_y_Outer);
	opt.ApplyExchange(u2_z, u_z_Outer);
	opt.ApplyExchange(u2, u_Outer);
}


//======================================================================
//to compute parameter in glob.conv.method; here parameter s is  pseudo frequency
//=================================================================================
void WavesScalarEqOpOpt::Init2DParameter()
{  
  //=================================================
  MV_Vector<int> markNodes(nno);
  opt.findBoundaryNodes(gg,markNodes);
  int sch=0; 
  
  //	ofstream out;
  //	out.open("bndNodes_fem.dat");
  
  for ( int l = 0; l < nno; l++)
    {
      if ( markNodes(l) == 1)
	sch++;
      
      //   out<<"bndNodes("<<l<<")="<<markNodes(l)<<"\n";
//      cout<<"bndNodes("<<l<<")="<<markNodes(l)<<"\n";
    }
  
  //cout<<"program ScalarEqOpOpt.C:l.3317:number of b.points is "<<sch<<endl;
  //	exit(1);
  
  //=====================================================
  
  int noprints = 40;
  
  Mat Ahelp;
  Mat Khelp;



	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &k2);
	CHKERRA(ierr);
	ierr = VecSetSizes(k2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(k2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);

	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

  
   PetscScalar zero = 0.0;
  
   // first component u1  : u1 is solution on k-1 iter.,
   //                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(k2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);


  int nno = gg.getNoNodes();
  int nsd = gg.getNoSpaceDim();
  
  // iteration matrix for system
  tms utimeassemblystart,utimeassemblystop;
  times(&utimeassemblystart);
 // cout<<"before iterationmatrix..."<<endl;
 // cout<<"nsd "<<nsd<<" nno"<<nno<<endl;
 
  if (nsd == 2)
    ParameterMatrix2D(Ahelp, lumpedAmass, type_of_material);
  
  times(&utimeassemblystop); 
  
  //cout<<"##### total time assembly= "
  //    <<(utimeassemblystop.tms_utime-utimeassemblystart.tms_utime)/100.0
  //    <<" #####\n";
  
  
       ierr = dropZeroWavesElements(Ahelp,A,1e-13);
 //      cout<<"after drop zero Ahelp"<<endl;
    
       ierr = MatDestroy(Ahelp); CHKERRA(ierr);
   
	 

   FEM_initialised = 1;
   /*
   if (USE_FEM)
     opt.InitExchangeFEM(gg, outer_gg);
  
   */


}

//======================================================================
//to compute parameter in glob.conv.method in 3D; here parameter s is  pseudo frequency
//=================================================================================
void WavesScalarEqOpOpt::Init3DParameter()
{  
//=================================================
 MV_Vector<int> markNodes(nno);
 opt.findBoundaryNodes(gg,markNodes);
 int sch=0; 
 
 //	ofstream out;
 //	out.open("bndNodes_fem.dat");
 
	for ( int l = 0; l < nno; l++)
	  {
	    if ( markNodes(l) == 1)
	      sch++;

	    //   out<<"bndNodes("<<l<<")="<<markNodes(l)<<"\n";
//	    cout<<"bndNodes("<<l<<")="<<markNodes(l)<<"\n";
	  }

//	cout<<"program ScalarEqOpOpt.C:l.3317:number of b.points is "<<sch<<endl;
	//	exit(1);

		//=====================================================
  
  int noprints = 40;
  
  Mat Ahelp;
  Mat Khelp;


	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &k2);
	CHKERRA(ierr);
	ierr = VecSetSizes(k2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(k2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);

	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

  
  PetscScalar zero = 0.0;
  
   // first component u1  : u1 is solution on k-1 iter.,
   //                       u2 is solution on k iteration.



	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(k2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);
 
  int nno = gg.getNoNodes();
  int nsd = gg.getNoSpaceDim();
  
  // iteration matrix for system
  tms utimeassemblystart,utimeassemblystop;
  times(&utimeassemblystart);
 // cout<<"before iterationmatrix..."<<endl;
 // cout<<"nsd "<<nsd<<" nno"<<nno<<endl;
  if (nsd == 2)
    {
 //     cout<<" call void ScalarEqOpOpt::Init2DParameter() !!!"<<endl;
      exit(1);
    }
  else if (nsd == 3)
    {
      //ParameterMatrix3D(Ahelp, lumpedAmass, type_of_material);
      // when  dropZeroWavesElements(Ahelp,A,1e-13) does not work we choose matrix A
    ParameterMatrix3D(A, lumpedAmass, type_of_material);
  
    }
  times(&utimeassemblystop); 
  
 // cout<<"##### total time assembly= "
//      <<(utimeassemblystop.tms_utime-utimeassemblystart.tms_utime)/100.0
//      <<" #####\n";
  
  /* commented 17.01.2012 - does not work properly MatGetSize   
       ierr = dropZeroElements(Ahelp,A,1e-13);
       cout<<"after drop zero Ahelp"<<endl;
    
  
       ierr = MatDestroy(Ahelp); CHKERRA(ierr);
   
  */	 

   FEM_initialised = 1;
   /*
   if (USE_FEM)
     opt.InitExchangeFEM(gg, outer_gg);
  
   */


}


//======================================================================
// Compute parameter in a globally convergent algorithm in 2D and 3D.
// Description of the algorithm is on the page 183 of our book.
// By default we set value 1 at the  FEM boundary.
//=======================================================================

void WavesScalarEqOpOpt::ParameterFEM(MV_Vector<double>& V_tilde, double& omega,  
				      MV_Vector<double>& param)
{
  
  // double t = k*dt;
  int nno = gg.getNoNodes();
  
  MV_Vector<int> markNodes(nno);
  PetscScalar minusone = -1.0;
  PetscScalar one = 1.0;
  PetscScalar zero = 0.0;
  
  PetscScalar* u1D;
  PetscScalar* u2D;
  PetscScalar* tmpvec;
  ierr =  VecGetArray(u1,&u1D);CHKERRA(ierr);
  
  //here, u1 = V_tilde(i)
  for( int i=0; i < gg.getNoNodes(); i++) 
    {	
      //get function v from v_tilde
      u1D[i] =  V_tilde(i)*omega*omega;
      //get function w from v
      u1D[i] = exp(u1D[i]);
      
    }
  
  ierr =  VecGetArray(lumpedAmass, &tmpvec);CHKERRA(ierr);
  
  for(int i = 0; i <  gg.getNoNodes(); i++){
    if( tmpvec[i] < +1e-6 ){
      //      cout<<"Warning: too low tmpvec["<<i<<"]=="<<tmpvec[i]<<endl;
    }
    //  here we do ((u1D*M^L)^-1))
    tmpvec[i] = 1.0/(tmpvec[i]*u1D[i]);
    
  }
  
  ierr =  VecRestoreArray(u1,&u1D);CHKERRA(ierr);  
  ierr = VecRestoreArray(lumpedAmass, &tmpvec); CHKERRA(ierr);  
	 
  ierr = MatDiagonalScale(A,lumpedAmass,PETSC_NULL); CHKERRA(ierr);


  
  ierr = MatMult(A,u1,u2);         // u2 = A * u1 
  
  ierr =  VecGetArray(u1,&u1D);CHKERRA(ierr);
  ierr =  VecGetArray(u2,&u2D);CHKERRA(ierr);
  
  
  
  for( int i=0; i < gg.getNoNodes(); i++) 
    {   
      
      u1D[i] = -u2D[i]/(omega*omega);  // for visualization of results
      param(i) = -u2D[i]/(omega*omega);
    }
  
  opt.findBoundaryNodes(gg,markNodes); 
  
  
  for( int i=0; i < gg.getNoNodes(); i++) 
    {
      
      
      if (u1D[i] < 0.0 || param(i) < 0.0)
	{
	  u1D[i] = 1.0;  
	  param(i) = 1.0;
	}
      
      if (markNodes(i) ==1)
	{
	  u1D[i] = 1.0;  // for visualization of results
	  param(i) = 1.0;
	}
      
      
      }
 
    
    ierr =  VecRestoreArray(u1,&u1D);CHKERRA(ierr); 
    ierr =  VecRestoreArray(u2,&u2D);CHKERRA(ierr); 
    
}
 
// given the function w, added by Thanh
void WavesScalarEqOpOpt::ParameterFEM2(MV_Vector<double> w, double omega,  
				      MV_Vector<double>& param)
{
  
  int nno = gg.getNoNodes();
  
  MV_Vector<int> markNodes(nno);
  PetscScalar minusone = -1.0;
  PetscScalar one = 1.0;
  PetscScalar zero = 0.0;
  
  PetscScalar* u1D;
  PetscScalar* u2D;
  PetscScalar* tmpvec;
  ierr =  VecGetArray(u1,&u1D);CHKERRA(ierr);
  
  //here, u1 = V_tilde(i)
  for( int i=0; i < gg.getNoNodes(); i++) 
    {	
 	u1D[i] = w(i);      
    }
  
  ierr =  VecGetArray(lumpedAmass, &tmpvec);CHKERRA(ierr);
  
  for(int i = 0; i <  gg.getNoNodes(); i++){
    if( tmpvec[i] < +1e-8 ){
            cout<<"Warning: too low tmpvec["<<i<<"]=="<<tmpvec[i]<<endl;
    }
    //  here we do ((u1D*M^L)^-1))
    tmpvec[i] = 1.0/(tmpvec[i]*u1D[i]);
    if (tmpvec[i] < +1e-10)
	cout << "Too small value " << tmpvec[i] << endl;
  }
  
  ierr =  VecRestoreArray(u1,&u1D);CHKERRA(ierr);  
  ierr = VecRestoreArray(lumpedAmass, &tmpvec); CHKERRA(ierr);  
	 
  ierr = MatDiagonalScale(A,lumpedAmass,PETSC_NULL); CHKERRA(ierr);


  
  ierr = MatMult(A,u1,u2);         // u2 = A * u1 
  
  ierr =  VecGetArray(u1,&u1D);CHKERRA(ierr);
  ierr =  VecGetArray(u2,&u2D);CHKERRA(ierr);
  
  
  
  for( int i=0; i < gg.getNoNodes(); i++) 
    {   
      
      u1D[i] = -u2D[i]/(omega*omega);  // for visualization of results
      param(i) = -u2D[i]/(omega*omega);
    }
  
  opt.findBoundaryNodes(gg,markNodes); 
  
  
  for( int i=0; i < gg.getNoNodes(); i++) 
    {
      
      
      if (u1D[i] < 0.0 || param(i) < 0.0)
	{
	  u1D[i] = 1.0;  
	  param(i) = 1.0;
	}
      
      if (markNodes(i) ==1)
	{
	  u1D[i] = 1.0;  // for visualization of results
	  param(i) = 1.0;
	}
      
      
      }
 
    
    ierr =  VecRestoreArray(u1,&u1D);CHKERRA(ierr); 
    ierr =  VecRestoreArray(u2,&u2D);CHKERRA(ierr); 
    
}

//=====================================================================
void WavesScalarEqOpOpt::Init2DFEM(MV_Vector<double>& b)
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	for (int l = 0; l < nno; l++)
	{
		if (markNodes(l) == 1)
			sch++;
	}

	//cout << "program WavesScalarEqOpOpt.C:l.3317:number of b.points is " << sch << endl;

	//=====================================================

	int noprints = 40;

	Mat Ahelp;

	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_SELF, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);

	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
	//cout << "before iterationmatrix..." << endl;
	//cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 2)
	{
		IterationMatrix2DDifMat(Ahelp, lumpedAmass, type_of_material, b);

	}

	times(&utimeassemblystop);

	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
	//cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;
}

//======================================================================

void WavesScalarEqOpOpt::Init3DFEM()
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	for (int l = 0; l < nno; l++)
	{
		if (markNodes(l) == 1)
			sch++;
	//	cout << "bndNodes(" << l << ")=" << markNodes(l) << "\n";
	}

//	cout << "program WavesScalarEqOpOpt.C:l.1822:number of b.points is " << sch << endl;

	//=====================================================

	int noprints = 40;

	Mat Ahelp;

	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
//	cout << "before iterationmatrix..." << endl;
//	cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 3)
	{
		IterationMatrix3D(Ahelp, lumpedAmass, type_of_material);

	}
	else if (nsd == 2)
	{
		IterationMatrix2D(Ahelp, lumpedAmass, type_of_material);

	}

	times(&utimeassemblystop);

//	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
//	cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;
}

//=======================================================================
//==========  note: PETSC_COMM_WORLD should be replaced by PETSC_COMM_SELF
//=========================================================================
void WavesScalarEqOpOpt::Init3DFEM(MV_Vector<double>& b)
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	for (int l = 0; l < nno; l++)
	{
		if (markNodes(l) == 1)
			sch++;
	}

//	cout << "program WavesScalarEqOpOpt.C:Init3DFEM(MV_Vector<double>&  b)number of b.points is " << sch << endl;

	//=====================================================

	int noprints = 40;

	Mat Ahelp;

	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_SELF, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
//	cout << "before iterationmatrix..." << endl;
//	cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 3)
	{
		IterationMatrix3DDifMatnew(Ahelp, lumpedAmass, type_of_material, b);

	}
	else if (nsd == 2)
	{
//		cout << "call Init2DFEM(MV_Vector<double>&  b) " << endl;
		exit(1);
	}

	times(&utimeassemblystop);

//	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
//	cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;
}

//=====================================================================

void WavesScalarEqOpOpt::Init2DGLK(MV_Vector<double>& b)
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	int noprints = 40;

	Mat Ahelp;
	Mat Chelp;
	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &uhelpC);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelpC, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelpC, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelpC, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
	//cout << "before iterationmatrix..." << endl;
	//cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 2)
	{
		IterationMatrix2DGLK(Ahelp, C, lumpedAmass, b);

	}

	PetscScalar* u_C;
	ierr = VecGetArray(C,&u_C);
	CHKERRA(ierr);

	for (int ii = 0; ii < nno; ii++)
	{
		if (fabs(u_C[ii]) > +1e-20)
			cout << "  u_C " << u_C[ii] << " b " << b(ii) << endl;

	}

	ierr = VecRestoreArray(C,&u_C);
	CHKERRA(ierr);

	times(&utimeassemblystop);

	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);

	//cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;
}
//=====================================================================

void WavesScalarEqOpOpt::InitExchangeCommon()
{

	opt.InitExchangeFEM(gg, outer_gg);

	opt.initExchangeFDM();

	opt.outerWithHole.presentWavesSDIndexes();

}

// init for exchange in the case when we compute in one program different wave equation solvers, to avoid computing
// maskFDM 
void WavesScalarEqOpOpt::InitExchangeCommon(WavesScalarEqOpOpt& p)
{

	opt.InitExchangeFEM(gg, outer_gg);

	int nOfLoops_ = p.opt.outerWithHole.nOfLoopIndex();
	outerWithHole.copyLI(p.opt.outerWithHole, nOfLoops_);

	nOfLoops_ = p.opt.outerBoundary1.nOfLoopIndex();
	outerBoundary1.copyLI(p.opt.outerBoundary1, nOfLoops_);

	nOfLoops_ = p.opt.outerBoundary2.nOfLoopIndex();
	outerBoundary2.copyLI(p.opt.outerBoundary2, nOfLoops_);
}

void WavesScalarEqOpOpt::InitExchangeReadMask()
{

	opt.InitExchangeFEM(gg, outer_gg);
	opt.initExchangeMaskFDM();

//	cout << "present   opt.outerWithHole " << endl;
	opt.outerWithHole.presentWavesSDIndexes();

}

void WavesScalarEqOpOpt::InitExchangeStructCommon()
{

	opt.InitExchangeFEM(gg, outer_gg);

	opt.initExchangeStructFDM();

//	cout << "present   opt.outerWithHole " << endl;
	opt.outerWithHole.presentWavesSDIndexes();

}

//======================================================================

void WavesScalarEqOpOpt::Init3DFEM_LENS()
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	for (int l = 0; l < nno; l++)
	{
		if (markNodes(l) == 1)
			sch++;
	//	cout << "bndNodes(" << l << ")=" << markNodes(l) << "\n";
	}

	//cout << "program WavesScalarEqOpOpt.C:l.1822:number of b.points is " << sch << endl;

	//=====================================================

	int noprints = 40;

	Mat Ahelp;

	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1_x);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_x, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u1_y);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_y, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u1_z);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_z, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u2_x);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_x, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2_y);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_y, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2_z);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_z, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_x, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_y, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_z, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_x, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_y, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_z, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);

	ierr = VecSet(u2_x, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2_y, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2_z, zero);
	CHKERRA(ierr);

	ierr = VecSet(u1_x, zero);
	CHKERRA(ierr);
	ierr = VecSet(u1_y, zero);
	CHKERRA(ierr);
	ierr = VecSet(u1_z, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
	//cout << "before iterationmatrix..." << endl;
	//cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 3)
	{

		IterationMatrix3Dnew(Ahelp, lumpedAmass, type_of_material);

	}

	times(&utimeassemblystop);

	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
	//cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;
}

//***********************************************************************

void WavesScalarEqOpOpt::Init3DFEMDifMat_LENS(MV_Vector<double>& new_velocity)
{

	//=================================================
	MV_Vector<int> markNodes(nno);
	opt.findBoundaryNodes(gg, markNodes);

	int sch = 0;

	//cout << "program WavesScalarEqOpOpt.C:l.1822:number of b.points is " << sch << endl;

	//=====================================================

	int noprints = 40;

	Mat Ahelp;

	// needs two solutions on the field : old and new one

	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1_x);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_x, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u1_y);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_y, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u1_z);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1_z, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u2_x);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_x, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2_y);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_y, PETSC_DECIDE, nno);
	CHKERRA(ierr);
	ierr = VecCreate(PETSC_COMM_WORLD, &u2_z);
	CHKERRA(ierr);
	ierr = VecSetSizes(u2_z, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &uhelp);
	CHKERRA(ierr);
	ierr = VecSetSizes(uhelp, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_x, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_y, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u2_z, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_x, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_y, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(u1_z, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(uhelp, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u1 is solution on k-1 iter.,
	//                       u2 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(uhelp, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);

	ierr = VecSet(u2_x, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2_y, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2_z, zero);
	CHKERRA(ierr);

	ierr = VecSet(u1_x, zero);
	CHKERRA(ierr);
	ierr = VecSet(u1_y, zero);
	CHKERRA(ierr);
	ierr = VecSet(u1_z, zero);
	CHKERRA(ierr);
	ierr = VecSet(Fn, zero);
	CHKERRA(ierr);

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
	//cout << "before iterationmatrix..." << endl;
	//cout << "nsd " << nsd << " nno" << nno << endl;

	if (nsd == 3)
	{
		IterationMatrix3DDifMatnew(Ahelp, lumpedAmass, type_of_material, new_velocity);
	}

	times(&utimeassemblystop);

	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
	//cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	FEM_initialised = 1;

}

//================ exchange in one direction:FDM to FEM ===============

void WavesScalarEqOpOpt::ExFDMtoFEM()
{

	opt.InitExchangeFEM(gg, outer_gg);

	opt.initExchangeFDMtoFEM();

}

//======================================================================

void WavesScalarEqOpOpt::InitFEMDifMat(int type_of_material, MV_Vector<double>& velocity)
{

	int noprints = 40;

	Mat Ahelp;

	int noeq = 1;
	int nno = noeq * gg.getNoNodes();

	// iteration matrix for system
	tms utimeassemblystart, utimeassemblystop;
	times(&utimeassemblystart);
	//cout << "before iterationmatrix..." << endl;
	//cout << "nsd " << nsd << " nno" << nno << endl;

	// here,we get Ahelp - assembled stiffness mass matrix,
	// lumpedAmass = -tau�(M^L)^�1

	IterationMatrix3DDifMat(Ahelp, lumpedAmass, type_of_material, velocity);

	times(&utimeassemblystop);

	cout << "##### total time assembly= " << (utimeassemblystop.tms_utime - utimeassemblystart.tms_utime) / 100.0 << " #####\n";

	ierr = dropZeroWavesElements(Ahelp, A, 1e-13);
	//cout << "after drop zero Ahelp" << endl;

	ierr = MatDestroy(Ahelp);
	CHKERRA(ierr);

	// needs two solutions on the field : old and new one

	Vec uhelp;
	ierr = VecCreate(PETSC_COMM_SELF, &Fn);
	CHKERRA(ierr);
	ierr = VecSetSizes(Fn, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecCreate(PETSC_COMM_WORLD, &u1);
	CHKERRA(ierr);
	ierr = VecSetSizes(u1, PETSC_DECIDE, nno);
	CHKERRA(ierr);

	ierr = VecSetType(u1, VECMPI);
	CHKERRA(ierr);
	ierr = VecSetType(Fn, VECMPI);
	CHKERRA(ierr);

	ierr = VecDuplicate(u1, &u2);
	CHKERRA( ierr);
	ierr = VecDuplicate(u1, &uhelp);
	CHKERRA( ierr);

	PetscScalar zero = 0.0;

	// first component u1  : u11 is solution on k-1 iter.,
	//                       u21 is solution on k iteration.

	ierr = VecSet(u1, zero);
	CHKERRA(ierr);
	ierr = VecSet(u2, zero);
	CHKERRA(ierr);

	FEM_initialised = 1;

}

//=================================================================
void WavesScalarEqOpOpt::InitDirichletFEM()
{

	// create bndNodes  on the boundary 
	MV_Vector<int> markNodes;

	opt.findBoundaryNodes(gg, markNodes);

	MV_Vector<int> bndNodes;

	opt.makeSortedNodesIndex(markNodes, bndNodes, 1);
	int nbn_ = bndNodes.size();

	bvalues_initialised = 1;

	boundindex = new int[nbn_];
	bound_values = new double[nbn_];

//initialize initial values = 0.0 in the all points
	for (int i = 0; i < nbn_; i++)
	{
		boundindex[i] = bndNodes(i);

		bound_values[i] = 0.0;
		//cout << "boundindex " << boundindex[i] << endl;
	}
	//cout << "before applyDirichlet, nbn = " << nbn_ << endl;

	nbn = nbn_;

	MV_Vector<int> boundary_Nodes(gg.getNoNodes());
	boundary_Nodes = 0;

	for (int i = 0; i < nbn_; i++)
		boundary_Nodes(boundindex[i]) = 2;
}

//====================================================================
void WavesScalarEqOpOpt::InitDirichletFEM(Mat_real& BoundaryCoord)
{
	// create bndNodes  on the boundary 
	MV_Vector<int> markNodes(BoundaryCoord.size(0));
	double eps = 1e-6;
	int i, j;
	int count = 0;

	for (i = 0; i < gg.getNoNodes(); i++)
	{
		double x = gg.getCoor(i, 0);
		double y = gg.getCoor(i, 1);
		for (j = 0; j < BoundaryCoord.size(0); j++)
		{
			double x_b = BoundaryCoord(j, 0);
			double y_b = BoundaryCoord(j, 1);

			if (fabs(x - x_b) < eps && fabs(y - y_b) < eps)
			{
				markNodes(count) = i;
				count++;
			}
		}
	}

	nbn = markNodes.size();

	PetscMalloc(nbn*sizeof(int), boundindex);

	PetscMalloc(nbn*sizeof(PetscScalar), bvalues);

	bvalues_initialised = 1;

	//initialize initial values = 0.0 in the all points
	for (int i = 0; i < nbn; i++)
	{
		boundindex[i] = markNodes(i);
		bvalues[i] = 0.0;
	}
}

//========================================================================

void WavesScalarEqOpOpt::ApplyRHS(double& t)
{

	int nrOuterNodes = sdg.getNoNodes();

	fn_Outer = new real[nrOuterNodes];

	if (USE_FDM)
	{

		if (nsd == 2)
		{

			if (rhs == 0.5)
			{
				initialize2DTimeFunction.applyTimeFunction(fn_Outer, t, dt * dt, elast_pulse1);

			}
		}
		else if (nsd == 3)
		{
			if (rhs == 0.5)
			{
				initialize3DTimeFunction.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs0_5_0_5);
			}
			else if (rhs == 1.0)
			{
				initialize3D_1.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs0_1_0_1);

			}
			else if (rhs == 5.0)
			{

				AssignTimeFunctionOp initialize3D_two(&sdg, D3rhs_two);

				initialize3D_two.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs_two);

			}
			else if (rhs == 6.0)
			{

				AssignTimeFunctionOp initialize3D_test(&sdg, test_exact);

				initialize3D_test.applyTimeFunction(fn_Outer, t, dt * dt, test_exact);

			}
			else if (rhs == 8.0)
			{

				AssignTimeFunctionOp init3D_two(&sdg, D3rhs_onepoint);

				init3D_two.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs_onepoint);

			}

		}
	}
	//======================================================================
	// for FEM/FDM 
	//======================================================================
	if (EXCHANGE)
	{
		if (nsd == 2)
		{

			AssignTimeFunctionOp initialize2D_1(&opt.outerWithHole, elast_pulse1);

			if (rhs == 0.5)
			{
				initialize2D_1.applyTimeFunction(fn_Outer, t, dt * dt, elast_pulse1);

			}
		}
		else if (nsd == 3)
		{

			if (rhs == 0.5)
			{

				AssignTimeFunctionOp initialize3D_1(&opt.outerWithHole, D3rhs0_5_0_5);

				initialize3D_1.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs0_5_0_5);

			}
			else if (rhs == 1.0)
			{

				AssignTimeFunctionOp initialize3D_1(&opt.outerWithHole, D3rhs0_1_0_1);

				initialize3D_1.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs0_1_0_1);

			}
			else if (rhs == 2.0)
			{

				AssignTimeFunctionOp initialize3D_1(&opt.outerWithHole, D3rhs_cube);
				initialize3D_1.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs_cube);

			}
			else if (rhs == 5.0)
			{

				AssignTimeFunctionOp initialize3D_two(&opt.outerWithHole, D3rhs_two);

				initialize3D_two.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs_two);

			}
			else if (rhs == 7.0)
			{

				AssignTimeFunctionOp init_two(&sdg, rhs_two);

				init_two.applyTimeFunction(fn_Outer, t, dt * dt, rhs_two);

			}
			else if (rhs == 8.0)
			{

				AssignTimeFunctionOp init3D_two(&opt.outerWithHole, D3rhs_onepoint);

				init3D_two.applyTimeFunction(fn_Outer, t, dt * dt, D3rhs_onepoint);

			}
			else if (rhs == 9.0)
			{

				AssignTimeFunctionOp init_four(&sdg, rhs_four);

				init_four.applyTimeFunction(fn_Outer, t, dt * dt, rhs_four);

			}
			else if (rhs == 10.0)
			{

				AssignTimeFunctionOp init_four(&sdg, pulse_left_right);

				init_four.applyTimeFunction(fn_Outer, t, dt * dt, pulse_left_right);
			}
		}

	}
}

//=========================================================================

void WavesScalarEqOpOpt::ApplyAdjRHS()
{

	//======================================================================
	// for FDM 
	//======================================================================
	if (EXCHANGE)
	{

		if (nsd == 3)
		{

			AssignAdjTimeFunctionOp initialize3D_1(&opt.outerWithHole, fn_Outer);

			initialize3D_1.applyTimeFunction(fn_Outer, -dt * dt, fn_Outer);

		}
	}
	else
	{
		cout << "other cases don't implemented for adjoint rhs!!!" << endl;
		exit(1);
	}
}

//===========================================================================

void WavesScalarEqOpOpt::AdjScalarFDM(int k, real* adj_f)
{
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double t = k * dt;
	bool marker;

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	if (EXCHANGE)
	{
		if (rank == 0)
//			cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		//cout << "fdm in all domain" << endl;
		timeStepperOuterFDM.apply(v_Outer, u_Outer);
	}

	int i;

	if (USE_RHS)
	{
		fn_Outer = adj_f;
		if (rank == 0)
			//cout << "apply rhs fdm" << endl;
		ApplyAdjRHS();

		if (rank == 0)
			//cout << " after apply rhs" << endl;
		add_sol_and_function.apply(fn_Outer, u_Outer);

	}

	if (USE_DIRICHLET_FDM)
	{
		dirichletBC.apply(u_Outer);
	}
	if (USE_ABSORB)
	{
		if (rank == 0)
			//cout << "with absorbing b.c." << endl;

		cornerBC.apply(u_Outer);

		abs1.apply(v_Outer, u_Outer);
		abs2.apply(v_Outer, u_Outer);
		abs3.apply(v_Outer, u_Outer);
		abs4.apply(v_Outer, u_Outer);
		abs5.apply(v_Outer, u_Outer);
		abs6.apply(v_Outer, u_Outer);

	}

}

//========================================================================

void WavesScalarEqOpOpt::Save_Sol_FEM(int nr_of_it, MV_Vector<double>& Save_u_fem, MV_Vector<int>& Vec_Glob_Num)
{
	int i, j;

	PetscScalar* u21D;
	ierr = VecGetArray(u2,&u21D);
	CHKERRA(ierr);

	int glob_node_numb = Vec_Glob_Num.size();

	for (i = 0; i < glob_node_numb; i++)
	{
		int glob_node = Vec_Glob_Num(i);

		Save_u_fem(i) = u21D[glob_node];
	}

	ierr = VecRestoreArray(u2,&u21D);
	CHKERRA(ierr);

}

void WavesScalarEqOpOpt::Save_Sol_FEM(MV_Vector<double>& Save_u_fem)
{
	int i;

	PetscScalar* u21D;
	ierr = VecGetArray(u2,&u21D);
	CHKERRA(ierr);

	for (i = 0; i < Save_u_fem.size(); i++)
		Save_u_fem(i) = u21D[i];

	ierr = VecRestoreArray(u2,&u21D);
	CHKERRA(ierr);

}

void WavesScalarEqOpOpt::Save_Sol_FEM(int nr_of_it, Mat_real& Save_u_fem)
{
	int i, j;

	PetscScalar* u21D;
	ierr = VecGetArray(u2,&u21D);
	CHKERRA(ierr);

	for (i = 0; i < gg.getNoNodes(); i++)
	{
		Save_u_fem(i, 0) = u21D[i];

	}

	ierr = VecRestoreArray(u2,&u21D);
	CHKERRA(ierr);

}

//=====================================================================

void WavesScalarEqOpOpt::Save_Sol_FDM(int nr_of_it, MV_Vector<double>& Save_u_fdm, MV_Vector<int>& Vec_Glob_Num)
{
	int i, j;

	int glob_node_numb = Vec_Glob_Num.size();

	for (i = 0; i < glob_node_numb; i++)
	{
		int glob_node = Vec_Glob_Num(i);

		Save_u_fdm(i) = u_Outer[glob_node];
	}
}

void WavesScalarEqOpOpt::Save_Sol_FDM(int nr_of_it, MV_Vector<double>& Save_u_fdm)
{
	int i;

	for (i = 0; i < sdg.getNoNodes(); i++)
	{
		Save_u_fdm(i) = u_Outer[i];
	}
}

//=======================================================================

void WavesScalarEqOpOpt::WaveEqSolverFEMforDifMat(double t, int k, int type_of_material, double velocity)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;

	if (USE_RHS)
	{
		//cout << "works rhs for fem" << endl;
		if (nsd == 2)
		{
			if (rhs == 1)
				evalRHSforDifMat(Fn, t, rhs1, lumpedAmass, type_of_material, velocity);
			if (rhs == 2)
				evalRHSforDifMat(Fn, t, rhs2, lumpedAmass, type_of_material, velocity);
			if (rhs == 0.5)
				evalRHSforDifMat(Fn, t, e_3rhs0_5_0_5, lumpedAmass, type_of_material, velocity);

		}
		else
		{
			if (rhs == 0.5)
				evalRHSforDifMat(Fn, t, impl_scalar, lumpedAmass, type_of_material, velocity);
			if (rhs == 5.0)
			{
				evalRHSforDifMat(Fn, t, rhs_amira, lumpedAmass, type_of_material, velocity);

			}
			if (rhs == 3)
			{
				//cout << "rhs = D3rhs0_5_0_5," << " time " << t << endl;
				evalRHSforDifMat(Fn, t, D3rhs0_5_0_5, lumpedAmass, type_of_material, velocity);

			}
			if (rhs == 9.0)
			{
				//cout << "rhs = rhs_four" << endl;
				evalRHSforDifMat(Fn, t, rhs_four, lumpedAmass, type_of_material, velocity);

			}
			if (rhs == 10)
			{
				//cout << "rhs = rhs_four" << endl;
				evalRHSforDifMat(Fn, t, pulse_left_right, lumpedAmass, type_of_material, velocity);
			}
			if (rhs == 11)
			{
				//cout << "rhs = rhs for sphere" << endl;
				evalRHSforDifMat(Fn, t, pulse_for_sphere, lumpedAmass, type_of_material, velocity);
			}

		}
	} // for use_rhs =true

	ierr = MatMult0(A, u1, u2); // u2 = A * u1 - u2

	PetscReal nrm;
	ierr = MatNorm(A, NORM_1, &nrm);
	CHKERRA(ierr);

//	cout << "Norma of matrice A " << nrm << endl;
	if (USE_RHS)
//HRN:	      ierr = VecAXPY(&dtdt,Fn,u2); // u2 = u2 + dt*dt*Fn
		ierr = VecAXPY(u2, dtdt, Fn); // u2 = u2 + dt*dt*Fn

	bool code = 0;
	PetscScalar* uhelp_;
	ierr = VecGetArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	PetscScalar* u2_;
	ierr = VecGetArray(u2,&u2_);
	CHKERRA(ierr);
	PetscScalar* u1_;
	ierr = VecGetArray(u1,&u1_);
	CHKERRA(ierr);

	ierr = VecRestoreArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2,&u2_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1,&u1_);
	CHKERRA(ierr);

	if (USE_DIRICHLET_FEM)
	{
		PetscScalar* u2D;
		ierr = VecGetArray(u2,&u2D);
		CHKERRA(ierr);
		PetscScalar* u1D;
		ierr = VecGetArray(u1,&u1D);
		CHKERRA(ierr);

		for (int i = 0; i < nbn; i++)
		{
			u2D[boundindex[i]] = bound_values[i];
			u1D[boundindex[i]] = bound_values[i];
		}

		ierr = VecRestoreArray(u2,&u2D);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1D);
		CHKERRA(ierr);

	}

}

//=================================================================

void WavesScalarEqOpOpt::WaveEqGLK(double t, int k, MV_Vector<double>& bb)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;
	PetscScalar one = 1.0;

	//cout << "nsd = " << nsd << endl;
	ierr = MatMult(A, u1, uhelp); // uhelp = A*u1

//HRN  ierr = VecAYPX(&minusone,uhelp,u2); // u2 = uhelp - 1*u2
	ierr = VecAYPX(u2, minusone, uhelp); // u2 = uhelp - 1*u2

	PetscReal nrm;
	ierr = MatNorm(A, NORM_1, &nrm);
	CHKERRA(ierr);

//	cout << "Norma of matrice A " << nrm << endl;

	//HRN	      ierr = VecAXPY(&dtdt,Fn,u2); // u2 = u2 + dt*dt*Fn
	if (USE_RHS)
		ierr = VecAXPY(u2, dtdt, Fn); // u2 = u2 + dt*dt*Fn
}
//=======================================================================

void WavesScalarEqOpOpt::PlaneWaveFEMforDifMat(double t, int k, int type_of_material, double velocity, double y_fix)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;

	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;
	//cout << "nsd = " << nsd << endl;

	ierr = MatMult0(A, u1, u2); // u2 = A * u1 - u2

	MV_Vector<int> markNodes;

	applyPlaneWaveFEM(markNodes, t, y_fix);

	PetscReal nrm;
	ierr = MatNorm(A, NORM_1, &nrm);
	CHKERRA(ierr);

//	cout << "Norma of matrice A " << nrm << endl;

	bool code = 0;
	PetscScalar* uhelp_;
	ierr = VecGetArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	PetscScalar* u2_;
	ierr = VecGetArray(u2,&u2_);
	CHKERRA(ierr);
	PetscScalar* u1_;
	ierr = VecGetArray(u1,&u1_);
	CHKERRA(ierr);

	for (int ii = 0; ii < nno; ii++)
	{
		if (markNodes(ii) == 5)
		{
			u1_[ii] = bound_values[ii];
			u2_[ii] = bound_values[ii];

		}
	}
	ierr = VecRestoreArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2,&u2_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1,&u1_);
	CHKERRA(ierr);

	if (USE_DIRICHLET_FEM)
	{
		PetscScalar* u2D;
		ierr = VecGetArray(u2,&u2D);
		CHKERRA(ierr);
		PetscScalar* u1D;
		ierr = VecGetArray(u1,&u1D);
		CHKERRA(ierr);

		for (int i = 0; i < nbn; i++)
		{
			u2D[boundindex[i]] = bound_values[i];
			u1D[boundindex[i]] = bound_values[i];
		}

		ierr = VecRestoreArray(u2,&u2D);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1D);
		CHKERRA(ierr);

	}

}

// ==============  plane wave for lens  ================================

//=======================================================================

void WavesScalarEqOpOpt::PlaneWaveFEM_LENS(double t, int k, MV_Vector<int>& markNodes)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;
	//cout << "nsd = " << nsd << endl;

	ierr = MatMult0(A, u1_x, u2_x);
	ierr = MatMult0(A, u1_y, u2_y);
	ierr = MatMult0(A, u1_z, u2_z);

	ierr = MatMult0(A, u1, u2);

	PetscReal nrm;
	ierr = MatNorm(A, NORM_1, &nrm);
	CHKERRA(ierr);

//	cout << "Norma of matrice A " << nrm << endl;

	bool code = 0;
	PetscScalar* uhelp_;
	ierr = VecGetArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	PetscScalar* u2_;
	ierr = VecGetArray(u2,&u2_);
	CHKERRA(ierr);
	PetscScalar* u1_;
	ierr = VecGetArray(u1,&u1_);
	CHKERRA(ierr);

	PetscScalar* u1_x_;
	ierr = VecGetArray(u1_x,&u1_x_);
	CHKERRA(ierr);
	PetscScalar* u1_y_;
	ierr = VecGetArray(u1_y,&u1_y_);
	CHKERRA(ierr);
	PetscScalar* u1_z_;
	ierr = VecGetArray(u1_z,&u1_z_);
	CHKERRA(ierr);

	PetscScalar* u2_x_;
	ierr = VecGetArray(u2_x,&u2_x_);
	CHKERRA(ierr);
	PetscScalar* u2_y_;
	ierr = VecGetArray(u2_y,&u2_y_);
	CHKERRA(ierr);
	PetscScalar* u2_z_;
	ierr = VecGetArray(u2_z,&u2_z_);
	CHKERRA(ierr);

	double n_x, n_y, n_z, norma;

	if (t < 2 * M_PI / 100.0)

	{
		for (int ii = 0; ii < nno; ii++)
		{

			if (markNodes(ii) == 5)
			{

				n_x = -gg.getCoor(ii, 0);
				n_y = -gg.getCoor(ii, 1);
				n_z = -gg.getCoor(ii, 2);

				norma = sqrt(n_x * n_x + n_y * n_y + n_z * n_z);

				u2_x_[ii] = (n_x * planeWavetest(0, 0, 0, t));
				u2_y_[ii] = (n_y * planeWavetest(0, 0, 0, t));
				u2_z_[ii] = (n_z * planeWavetest(0, 0, 0, t));

				u2_[ii] = norma * planeWavetest(0, 0, 0, t);
			}
		}
	}
	ierr = VecRestoreArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2,&u2_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1,&u1_);
	CHKERRA(ierr);

	ierr = VecRestoreArray(u1_x,&u1_x_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1_y,&u1_y_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1_z,&u1_z_);
	CHKERRA(ierr);

	ierr = VecRestoreArray(u2_x,&u2_x_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2_y,&u2_y_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2_z,&u2_z_);
	CHKERRA(ierr);

}

//*************** plane wave at the earth ******************

//=======================================================================

void WavesScalarEqOpOpt::PlaneWaveFEM_Earth(double t, int k, MV_Vector<int>& markNodes)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;
	//cout << "nsd = " << nsd << endl;

	ierr = MatMult0(A, u1, u2);

	PetscReal nrm;
	ierr = MatNorm(A, NORM_1, &nrm);
	CHKERRA(ierr);

//	cout << "Norma of matrice A " << nrm << endl;

	bool code = 0;
	PetscScalar* uhelp_;
	ierr = VecGetArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	PetscScalar* u2_;
	ierr = VecGetArray(u2,&u2_);
	CHKERRA(ierr);
	PetscScalar* u1_;
	ierr = VecGetArray(u1,&u1_);
	CHKERRA(ierr);

	double norma;

	if (t < 2 * M_PI / 100.0)

	{
		for (int ii = 0; ii < nno; ii++)
		{

			if (markNodes(ii) == 5)
			{

				u2_[ii] = planeWavetest(0, 0, 0, t);
				//cout << " u2 " << u2_[ii] << endl;
			}
		}
	}

	ierr = VecRestoreArray(uhelp,&uhelp_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u2,&u2_);
	CHKERRA(ierr);
	ierr = VecRestoreArray(u1,&u1_);
	CHKERRA(ierr);

}

//**** end of plane wave at the earth ************************

//=======================================================================

void WavesScalarEqOpOpt::WaveEqSolverFEMforDifMat(double t, int k, int type_of_material, MV_Vector<double>& velocity)
{
	int nno = gg.getNoNodes();
//	cout << "fem time is  " << t << endl;
	PetscScalar timest = dt;

	PetscScalar dt2 = 2 * dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar minusone = -1.0;

	if (USE_RHS)
	{
		//cout << "works rhs for fem" << endl;
		if (nsd == 2)
		{
			if (rhs == 1)
				evalRHSforDifMat(Fn, t, rhs1, lumpedAmass, type_of_material, velocity);
			if (rhs == 2)
				evalRHSforDifMat(Fn, t, rhs2, lumpedAmass, type_of_material, velocity);
			if (rhs == 0.5)
				evalRHSforDifMat(Fn, t, e_3rhs0_5_0_5, lumpedAmass, type_of_material, velocity);
		}
		else
		{
			if (rhs == 0.5)
				evalRHSforDifMat(Fn, t, D3rhs0_5_0_5, lumpedAmass, type_of_material, velocity);

			if (rhs == 3)
			{
				//cout << "rhs = D3rhs0_5_0_5" << endl;
				evalRHSforDifMat(Fn, t, D3rhs0_5_0_5, lumpedAmass, type_of_material, velocity);
			}

		}
		if (rhs == 9.0)
		{
			//cout << "rhs = rhs_four" << endl;
			evalRHSforDifMat(Fn, t, rhs_four, lumpedAmass, type_of_material, velocity);

		}
		if (rhs == 10)
		{
			//cout << "rhs = rhs_four" << endl;
			evalRHSforDifMat(Fn, t, pulse_left_right, lumpedAmass, type_of_material, velocity);

		}
		if (rhs == 11)
		{
			//cout << "rhs = rhs for sphere" << endl;
			evalRHSforDifMat(Fn, t, pulse_for_sphere, lumpedAmass, type_of_material, velocity);
		}

	} // for use_rhs =true

	ierr = MatMult0(A, u1, u2); // u2 = A * u1 - u01
	if (USE_RHS)
//HRN	   ierr = VecAXPY(&dtdt,Fn,u2); // u2 = u2 + dt*dt*Fn
		ierr = VecAXPY(u2, dtdt, Fn); // u2 = u2 + dt*dt*Fn

	if (USE_DIRICHLET_FEM)
	{
		PetscScalar* u2D;
		ierr = VecGetArray(u2,&u2D);
		CHKERRA(ierr);
		PetscScalar* u1D;
		ierr = VecGetArray(u1,&u1D);
		CHKERRA(ierr);

		for (int i = 0; i < nbn; i++)
		{
			u2D[boundindex[i]] = bound_values[i];
			u1D[boundindex[i]] = bound_values[i];
		}
		ierr = VecRestoreArray(u2,&u2D);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1D);
		CHKERRA(ierr);
	}

}

//==================================================
void WavesScalarEqOpOpt::SaveSolutionsFEM(int nr_of_it, MV_Vector<double>& Save_fem)
{
	int i, j;
	PetscScalar* u21D;
	ierr = VecGetArray(u2,&u21D);
	CHKERRA(ierr);

	//	cout << " works savesolutions" << endl;

	for (i = 0; i < gg.getNoNodes(); i++)
	{
		Save_fem(i) = u21D[i];
		//	cout << " FEM-solution[" << i << "]=" << 	Save_fem(i) << endl;
	}

	ierr = VecRestoreArray(u2,&u21D);
	CHKERRA(ierr);
}

//==========================================

void WavesScalarEqOpOpt::ApplySwap()
{

	if (EXCHANGE)
	{

		tmpOuter = u_Outer;
		u_Outer = v_Outer;
		v_Outer = tmpOuter;
		ierr = VecSwap(u1, u2);
		CHKERRA(ierr);

	}
	else if (USE_FDM)
	{
		tmpOuter = u_Outer;
		u_Outer = v_Outer;
		v_Outer = tmpOuter;
	//	cout << "after swap" << endl;
	}
	else if (USE_FEM)
	{

		bool code = 0;
		PetscScalar* u2_;
		ierr = VecGetArray(u2,&u2_);
		CHKERRA(ierr);
		PetscScalar* u1_;
		ierr = VecGetArray(u1,&u1_);
		CHKERRA(ierr);

		for (int ii = 0; ii < nno; ii++)
		{
			if (u2_[ii] > 1e+30 || u2_[ii] < -1e30)
			{
				code = 1;
		//		cout << " u2[" << ii << "]=" << u2_[ii] << " u1[" << ii << "]=" << u1_[ii] << endl;
			}
		}

		ierr = VecRestoreArray(u2,&u2_);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1_);
		CHKERRA(ierr);

		ierr = VecSwap(u1, u2);
		CHKERRA(ierr);
	//	cout << "after FEM swap" << endl;
		ierr = VecGetArray(u2,&u2_);
		CHKERRA(ierr);

		ierr = VecGetArray(u1,&u1_);
		CHKERRA(ierr);

		for (int ii = 0; ii < nno; ii++)
		{
			if (u2_[ii] > 1e+30 || u2_[ii] < -1e30 || u1_[ii] > 1e+30 || u1_[ii] < -1e30)
			{
				code = 1;
			//	cout << " u2[" << ii << "]=" << u2_[ii] << " u1[" << ii << "]=" << u1_[ii] << endl;
			}
		}

		ierr = VecRestoreArray(u2,&u2_);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1_);
		CHKERRA(ierr);
	}

}

void WavesScalarEqOpOpt::ApplySwap_LENS()
{

	if (EXCHANGE)
	{

		tmpOuter = u_x_Outer;
		u_x_Outer = v_x_Outer;
		v_x_Outer = tmpOuter;
		tmpOuter = u_y_Outer;
		u_y_Outer = v_y_Outer;
		v_y_Outer = tmpOuter;
		tmpOuter = u_z_Outer;
		u_z_Outer = v_z_Outer;
		v_z_Outer = tmpOuter;

		tmpOuter = u_Outer;
		u_Outer = v_Outer;
		v_Outer = tmpOuter;

		ierr = VecSwap(u1_x, u2_x);
		CHKERRA(ierr);
		ierr = VecSwap(u1_y, u2_y);
		CHKERRA(ierr);
		ierr = VecSwap(u1_z, u2_z);
		CHKERRA(ierr);

		ierr = VecSwap(u1, u2);
		CHKERRA(ierr);

	}

}

//=====================================================================
void WavesScalarEqOpOpt::print_GID_MESH(char* filename)
{

	opt.print_gid_mesh(filename);
}

int WavesScalarEqOpOpt::print_AVS_FEM(int sch)
{

	opt.ApplyAvsOut(sch, u2);

	return ierr;
}

int WavesScalarEqOpOpt::print_AVS_FDM(int sch)
{

	AVSOutputOp avs2dFDM(&sdg, ofs);

	AVSOutputOp avs2d(&opt.outerWithHole, ofs);

	opt.myOpen(ofs);

	if (USE_FDM)
		avs2dFDM.apply(u_Outer);
	else if (EXCHANGE)
		avs2d.apply(u_Outer);

	ofs.close();

	return ierr;
}

int WavesScalarEqOpOpt::print_AVS_ElType(int sch, int EL_TYPE)
{

	opt.ApplyAvsOutElType(sch, u2, EL_TYPE);

	return ierr;
}

void WavesScalarEqOpOpt::print_GID_MESH_FEM(WavesGridB& grid, char* filename)
{
	opt.print_gid_mesh_Amira(grid, filename);
}

void WavesScalarEqOpOpt::print_GID_MESH_FDM(char* filename)
{

	bool pr;

	if (EXCHANGE || USE_FEM)
	{
		GIDOutputMesh gid_out(&opt.outerWithHole, filename);

		pr = gid_out.printMesh();

	}
	if (USE_FDM)
	{
		GIDOutputMesh gid_out(&sdg, filename);

		pr = gid_out.printMesh();

	}

}

void WavesScalarEqOpOpt::print_mesh_common(char* filename, WavesGridB& grid)
{

	bool pr;
	int nno_gg = gg.getNoNodes();
	int nel_gg = gg.getNoElms();

	if (nsd == 2)
	{

		opt.print_gid_nodes_FEM(grid, filename);

		GIDOutputNodes nodes(&opt.outerWithHole, filename, nno_gg);

		pr = nodes.printMesh();

		opt.print_gid_elements_FEM(grid, filename);

		GIDOutputWavesElements elements(&opt.outerWithHole, filename, nel_gg, nno_gg);

		pr = elements.printMesh();

	}

}

void WavesScalarEqOpOpt::print_mesh_common(char* filename)
{

	bool pr;

	int nsd = gg.getNoSpaceDim();

	// old version
	//	GIDOutputMesh  gid_out(&sdg, filename);
	//	pr = gid_out.printCommonMesh();


  if (nsd == 2)
    {
     
    GIDOutputMesh   gid_out(&opt.outerWithHole,
	  		    filename, gg.getNoNodes(), gg.getNoElms());
	   
     pr =  gid_out.printCommonMesh();
	  
  
	
    }
  else if (nsd == 3)
    {
     
    GIDOutputMesh   gid_out(&opt.outerWithHole,
	  		    filename, gg.getNoNodes(), gg.getNoElms());
	   
     pr =  gid_out.printCommonMesh();
	  
  
	
    }


}

void WavesScalarEqOpOpt::InitPlaneWave()
{

	int nrOuterNodes = sdg.getNoNodes();

	// An operator that initializes time 0 on outer grid
	WavesAssignmentOp initializerOuter(&sdg, 0.0);

	u_Outer = new real[nrOuterNodes];
	v_Outer = new real[nrOuterNodes];
}

void WavesScalarEqOpOpt::InitPlaneWave1()
{

	int nrOuterNodes = sdg.getNoNodes();
	//cout << " nr outer nodes " << nrOuterNodes << endl;
	// An operator that initializes time 0 on outer grid
	WavesAssignmentOp initializerOuterFDM(&sdg, 1.0);

	WavesAssignmentOp initializerOuter(&outerWithHole, 1.0);

	//cout << "  outerWithHole.presentWavesSDIndexes in initplanewave" << endl;
	outerWithHole.presentWavesSDIndexes();

	PlaneWaveOpexpl planeWaveLeft(&leftBoundary, planeWave);

	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);
	// arrays for the unknowns on Outer
	u_Outer = new real[nrOuterNodes];
	v_Outer = new real[nrOuterNodes];

	b_Outer = new real[nrOuterNodes];

	if (EXCHANGE)
		initializerOuter.apply(b_Outer); // b is inited 
	else if (USE_FDM)
		initializerOuterFDM.apply(b_Outer);

}

void WavesScalarEqOpOpt::InitPlaneWave_GLK()
{

	int nrOuterNodes = sdg.getNoNodes();
	//cout << " nr outer nodes " << nrOuterNodes << endl;
	// An operator that initializes time 0 on outer grid
	WavesAssignmentOp initializerOuterFDM(&sdg, 1.0);

	WavesAssignmentOp initializerOuter(&opt.outerWithHole, 1.0);

	PlaneWaveOpexpl planeWaveLeft(&leftBoundary, planeWave);
	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);
	// arrays for the unknowns on Outer
	u_Outer = new real[nrOuterNodes];
	v_Outer = new real[nrOuterNodes];

	b_Outer = new real[nrOuterNodes];

	if (EXCHANGE)
		initializerOuter.apply(b_Outer); // b is inited 
	else if (USE_FDM)
		initializerOuterFDM.apply(b_Outer);
}

void WavesScalarEqOpOpt::InitPlaneWave_LENS()
{

	int nrOuterNodes = sdg.getNoNodes();

	// An operator that initializes time 0 on outer grid
	WavesAssignmentOp initializerOuter(&sdg, 0.0);

	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);
	// arrays for the unknowns on Outer

	u_Outer = new real[nrOuterNodes];
	v_Outer = new real[nrOuterNodes];

	u_x_Outer = new real[nrOuterNodes];
	v_x_Outer = new real[nrOuterNodes];

	u_y_Outer = new real[nrOuterNodes];
	v_y_Outer = new real[nrOuterNodes];

	u_z_Outer = new real[nrOuterNodes];
	v_z_Outer = new real[nrOuterNodes];

	initializerOuter.apply(u_Outer); // u is inited
	initializerOuter.apply(v_Outer);

	initializerOuter.apply(u_x_Outer); // u is inited
	initializerOuter.apply(v_x_Outer);

	initializerOuter.apply(u_y_Outer); // u is inited
	initializerOuter.apply(v_y_Outer);

	initializerOuter.apply(u_z_Outer); // u is inited
	initializerOuter.apply(v_z_Outer);

}

//======================================================================

void WavesScalarEqOpOpt::InitFDM()
{
	int nrOuterNodes = sdg.getNoNodes();
	nsd = sdg.getNoSpaceDim();

	// An operator that initializes time 0 on outer grid
	WavesAssignmentOp initializerOuter(&sdg, 0.0);

	// arrays for the unknowns on Outer
	u_Outer = new real[nrOuterNodes];
	v_Outer = new real[nrOuterNodes];

	initializerOuter.apply(u_Outer); // u is inited
	initializerOuter.apply(v_Outer);

}

//==================================================================
void WavesScalarEqOpOpt::PlaneWaveLeftFDM(int k)
{

	double t = k * dt;

	PlaneWaveOpexpl planeWaveLeft(&leftBoundary, planeWave);

	PlaneWaveOpexpl planeWaveLeftvel(&leftBoundary, derplaneWave);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absRight(&sdg, SDiHigh, dt);
	WavesWaveEqBoundary absLeft(&sdg, SDiLow, dt);

	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveLeft.applyPlaneWaveOpexpl(u_Outer, t, planeWave);

	cornerBC.apply(u_Outer);

	absRight.apply(v_Outer, u_Outer); // uOuter is set on boundary
	mirrorTop.apply(v_Outer, u_Outer);
	mirrorLow.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==================================================================
void WavesScalarEqOpOpt::PlaneWaveBackFDM(int k)
{

	double t = k * dt;

	PlaneWaveOpexpl planeWaveBack(&backBoundary, planeWave);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absPer(&sdg, SDkLow, dt);
	WavesWaveEqBoundary absBack(&sdg, SDkHigh, dt);

	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveBack.applyPlaneWaveOpexpl(u_Outer, t, planeWave);

	cornerBC.apply(u_Outer);

	absPer.apply(v_Outer, u_Outer); // uOuter is set on boundary
	mirrorTop.apply(v_Outer, u_Outer);
	mirrorLow.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorLeft.apply(v_Outer, u_Outer);
		mirrorRight.apply(v_Outer, u_Outer);
	}

}


void WavesScalarEqOpOpt::PlaneWaveBackFDM(int k, double omega)
{

  double t = k*dt;
  
  PlaneWaveOpPhotonexpl planeWaveBack(&backBoundary,planeWave);
  
  WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole,dt);
  WavesWaveEqInterior timeStepperOuterFDM(&sdg,dt);
  
  // WaveEqBoundary absRight( &sdg, SDiHigh, dt );
  // WaveEqBoundary absLeft( &sdg, SDiLow, dt );
  
  WavesWaveEqBoundary absPer(&sdg, SDkLow, dt );
  WavesWaveEqBoundary absBack(&sdg, SDkHigh, dt );
  
  WavesMirrorBC mirrorLow(&lowBoundary );
  WavesMirrorBC  mirrorTop(&topBoundary );
  
  WavesMirrorBC mirrorLeft(&leftBoundary);
  WavesMirrorBC mirrorRight(&rightBoundary); 
  
  if (EXCHANGE){
//    cout<<"fdm for exchange"<<endl;
    timeStepperOuter.apply( v_Outer, u_Outer );
    
  } // u is obtained, outer grid
  else {
    // cout<<"fdm in all domain"<<endl;
    timeStepperOuterFDM.apply( v_Outer, u_Outer );
    
}  
//       cout<<"BC on outer domain"<<endl;
       // Apply BCs on outer domain

 if (nsd == 2)
	 planeWaveBack.applyPlaneWaveOpexpl(u_Outer,t, omega,planeWaveOmega);
       else if (nsd == 3)
	 planeWaveBack.applyPlaneWaveOpexpl(u_Outer,t, omega,planeWaveOmega2);
      

 cornerBC.apply(u_Outer);

 absPer.apply(v_Outer,u_Outer);        // uOuter is set on boundary
 mirrorTop.apply(v_Outer,u_Outer);        
 mirrorLow.apply(v_Outer,u_Outer);        
 
  
 if (nsd == 3){
    mirrorLeft.apply(v_Outer,u_Outer);        
    mirrorRight.apply(v_Outer,u_Outer);   
  }

  
}



void WavesScalarEqOpOpt::PlaneWaveBackFDM(int k, double omega, double& param)
{

  double t = k*dt;
  
  PlaneWaveOpPhotonexpl planeWaveBack(&backBoundary,planeWave);
  

  WavesWaveEqInter timeStepperOuter(&opt.outerWithHole,dt,param);
  WavesWaveEqInter timeStepperOuterFDM(&sdg,dt,param);

  // WaveEqBoundary absRight( &sdg, SDiHigh, dt );
  // WaveEqBoundary absLeft( &sdg, SDiLow, dt );
  
  WavesWaveEqBoundary absPer(&sdg, SDkLow, dt );
  WavesWaveEqBoundary absBack(&sdg, SDkHigh, dt );
  
  WavesMirrorBC mirrorLow(&lowBoundary );
  WavesMirrorBC  mirrorTop(&topBoundary );
  
  WavesMirrorBC mirrorLeft(&leftBoundary);
  WavesMirrorBC mirrorRight(&rightBoundary); 
  
  if (EXCHANGE){
//    cout<<"fdm for exchange"<<endl;
    timeStepperOuter.apply( v_Outer, u_Outer );
    
  } // u is obtained, outer grid
  else {
    // cout<<"fdm in all domain"<<endl;
    timeStepperOuterFDM.apply( v_Outer, u_Outer );
    
}  
//       cout<<"BC on outer domain"<<endl;
       // Apply BCs on outer domain

 if (nsd == 2)
	 planeWaveBack.applyPlaneWaveOpexpl(u_Outer,t, omega,planeWaveOmega);
       else if (nsd == 3)
	 planeWaveBack.applyPlaneWaveOpexpl(u_Outer,t, omega,planeWaveOmega2);
      

 cornerBC.apply(u_Outer);

 absPer.apply(v_Outer,u_Outer);        // uOuter is set on boundary
 mirrorTop.apply(v_Outer,u_Outer);        
 mirrorLow.apply(v_Outer,u_Outer);        
 
  
 if (nsd == 3){
    mirrorLeft.apply(v_Outer,u_Outer);        
    mirrorRight.apply(v_Outer,u_Outer);   
  }

  
}

//==================================================================
void WavesScalarEqOpOpt::PlaneWavePerFDM(int k)
{

	double t = k * dt;

	PlaneWaveOpexpl planeWavePer(&perBoundary, planeWave);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absPer(&sdg, SDkLow, dt);
	WavesWaveEqBoundary absBack(&sdg, SDkHigh, dt);

	WavesMirrorBC mirrorLow(&lowBoundary);
	WavesMirrorBC mirrorTop(&topBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWavePer.applyPlaneWaveOpexpl(u_Outer, t, planeWave);

	cornerBC.apply(u_Outer);

	absBack.apply(v_Outer, u_Outer); // uOuter is set on boundary
	mirrorTop.apply(v_Outer, u_Outer);
	mirrorLow.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorLeft.apply(v_Outer, u_Outer);
		mirrorRight.apply(v_Outer, u_Outer);
	}

}

//=========================================================================
//==================================================================
void WavesScalarEqOpOpt::PlaneWaveTopFDM(int k, double omega, double& param)
{

	double t = k * dt;

	PlaneWaveOpPhotonexpl planeWaveTop(&topBoundary, planeWaveOmega);

	WavesWaveEqInter timeStepperOuter(&opt.outerWithHole, dt, param);
	WavesWaveEqInter timeStepperOuterFDM(&sdg, dt, param);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	if (nsd == 2)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, delta2);
	else if (nsd == 3)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega2);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}
//==================================================================
void WavesScalarEqOpOpt::PlaneWaveTopFDM(int k, double omega)
{

	double t = k * dt;

	PlaneWaveOpPhotonexpl planeWaveTop(&topBoundary, planeWaveOmega2);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	if (nsd == 2)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega);
	else if (nsd == 3)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega2);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}
//======================================================================
//==================================================================
void WavesScalarEqOpOpt::PlaneWaveTop1FDM(int k, double omega, double& param)
{

	double t = k * dt;

	PlaneWaveOpPhotonexpl planeWaveTop(&topBoundary, planeWaveOmega);

	WavesWaveEqInter timeStepperOuter(&outerWithHole, dt, param);
	WavesWaveEqInter timeStepperOuterFDM(&sdg, dt, param);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	if (nsd == 2)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega);
	else if (nsd == 3)
		planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega2);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}
//==================================================================
void WavesScalarEqOpOpt::PlaneWaveBotFDM(int k, double omega)
{

	double t = k * dt;

	PlaneWaveOpPhotonexpl planeWaveBot(&lowBoundary, planeWaveOmega);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveBot.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}
//======================================================================

//==================================================================
void WavesScalarEqOpOpt::PlaneWaveBotFDM(int k, double omega, double& param)
{

	double t = k * dt;
	PlaneWaveOpPhotonexpl planeWaveBot(&lowBoundary, planeWaveOmega);

	WavesWaveEqInter timeStepperOuter(&opt.outerWithHole, dt, param);
	WavesWaveEqInter timeStepperOuterFDM(&sdg, dt, param);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveBot.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==================================================================
void WavesScalarEqOpOpt::PlaneWaveLeftFDM(int k, double omega)
{

	double t = k * dt;

	PlaneWaveOpPhotonexpl planeWaveLeft(&leftBoundary, planeWaveOmega);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absRight(&sdg, SDiHigh, dt);
	WavesWaveEqBoundary absLeft(&sdg, SDiLow, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorTop(&topBoundary);
	WavesMirrorBC mirrorBot(&lowBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveLeft.applyPlaneWaveOpexpl(u_Outer, t, omega, planeWaveOmega);

	cornerBC.apply(u_Outer);

	absRight.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / omega;

	mirrorTop.apply(v_Outer, u_Outer);
	mirrorBot.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//======================================================================

void WavesScalarEqOpOpt::apply_reflwave_low(int& k, double omega, MV_Vector<double>& array)
{

	double t_max = 2 * M_PI / omega;

	double t = k * dt;

	WavesReflectedWave refl_wave_low(&lowBoundary, array);
	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);
	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;

	refl_wave_low.applyReflectedWaveLow(u_Outer);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==================================================================
void WavesScalarEqOpOpt::ReflAdjPlaneWaveFDM(int k, real *Difsol)
{

	double t = k * dt;

// to init plane wave at the middle of the domain
// at line SDline
	//cout << " before init plane wave at the middle of domain " << endl;

	int node_left = sdg.coord2node(-3, -3);
	int node_right = sdg.coord2node(3, -3);

	int n_i = sdg.node2i(node_left);
	int n_j = sdg.node2j(node_left);

	int nodeLow = sdg.node(n_i, n_j);

	n_i = sdg.node2i(node_right);
	n_j = sdg.node2j(node_right);
	int nodeHigh = sdg.node(n_i, n_j);

	//cout << "nodeLow " << nodeLow << "  nodeHigh  " << nodeHigh << endl;

	WavesSDIndexLine middle(sdg, nodeLow, nodeHigh);
//	cout << " after sdindexline" << endl;

	PlaneWaveOpexpl3 planeWaveTop(&topBoundary, Difsol);
	PlaneWaveOpexpl3 planeWaveBot(&lowBoundary, Difsol);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==================================================================
void WavesScalarEqOpOpt::AdjPlaneWaveTopFDM(int k, real *Difsol)
{

	double t = k * dt;

// to init plane wave at the middle of the domain
// at line SDline
//	cout << " before init plane wave at the middle of domain " << endl;

	int node_left = sdg.coord2node(-3, -3);
	int node_right = sdg.coord2node(3, -3);

	int n_i = sdg.node2i(node_left);
	int n_j = sdg.node2j(node_left);

	int nodeLow = sdg.node(n_i, n_j);

	n_i = sdg.node2i(node_right);
	n_j = sdg.node2j(node_right);
	int nodeHigh = sdg.node(n_i, n_j);

//	cout << "nodeLow " << nodeLow << "  nodeHigh  " << nodeHigh << endl;

	WavesSDIndexLine middle(sdg, nodeLow, nodeHigh);
//	cout << " after sdindexline" << endl;

	PlaneWaveOpexpl3 planeWaveTop(&topBoundary, Difsol);
	PlaneWaveOpexpl3 planeWaveBot(&lowBoundary, Difsol);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		// cout<<"fdm in all domain"<<endl;
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	planeWaveBot.applyPlaneWaveOpexpl(u_Outer, t);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==========================================================================

void WavesScalarEqOpOpt::AdjPlaneWaveFDM(int k, real *Difsol, WavesSDBoundary& inner)
{

	double t = k * dt;

	PlaneWaveOpexpl3 planeWaveTop(&inner, Difsol);
	PlaneWaveOpexpl3 planeWaveBot(&inner, Difsol);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	planeWaveBot.applyPlaneWaveOpexpl(u_Outer, t);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==================================================================
void WavesScalarEqOpOpt::AdjPlaneWaveLeftFDM(int k, real *Difsol)
{

	double t = k * dt;

// to init plane wave at the middle of the domain
// at line SDline
//	cout << " before init plane wave at the middle of domain " << endl;

	int node_left = sdg.coord2node(-3, -3);
	int node_right = sdg.coord2node(3, -3);

	int n_i = sdg.node2i(node_left);
	int n_j = sdg.node2j(node_left);

	int nodeLow = sdg.node(n_i, n_j);

	n_i = sdg.node2i(node_right);
	n_j = sdg.node2j(node_right);
	int nodeHigh = sdg.node(n_i, n_j);

//	cout << "nodeLow " << nodeLow << "  nodeHigh  " << nodeHigh << endl;

	WavesSDIndexLine middle(sdg, nodeLow, nodeHigh);
//	cout << " after sdindexline" << endl;

	PlaneWaveOpexpl3 planeWaveTop(&topBoundary, Difsol);
	PlaneWaveOpexpl3 planeWaveRight(&rightBoundary, Difsol);

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absRight(&sdg, SDiHigh, dt);
	WavesWaveEqBoundary absLeft(&sdg, SDiLow, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorTop(&topBoundary);
	WavesMirrorBC mirrorBot(&lowBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain

	planeWaveRight.applyPlaneWaveOpexpl(u_Outer, t);

	cornerBC.apply(u_Outer);

	absLeft.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorTop.apply(v_Outer, u_Outer);
	mirrorBot.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//================ seismic start ================================

//==================================================================
void WavesScalarEqOpOpt::PlaneWaveTopFDMSeismic(int k)
{

	double t = k * dt;

	PlaneWaveOpexpl planeWaveTop(&topBoundary, planeWave);

	WavesWaveEqInteriorSeismic timeStepperOuter(&opt.outerWithHole, b_Outer, dt);
	WavesWaveEqInteriorSeismic timeStepperOuterFDM(&sdg, b_Outer, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, planeWave);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//==========seismic end ============================================

//======= ADJOINT SEISMIC START ========================0

//==================================================================
void WavesScalarEqOpOpt::AdjPlaneWaveFDMSeismic(int k, real *Difsol)
{

	double t = k * dt;

	PlaneWaveOpexpl3 planeWaveTop(&topBoundary, Difsol);

	PlaneWaveOpexpl3 planeWaveBot(&lowBoundary, Difsol);

	WavesWaveEqInteriorSeismic timeStepperOuter(&opt.outerWithHole, b_Outer, dt);
	WavesWaveEqInteriorSeismic timeStepperOuterFDM(&sdg, b_Outer, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveBot.applyPlaneWaveOpexpl(u_Outer, t);

	cornerBC.apply(u_Outer);

	absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//========ADJOINT SEISMIC END =================================0000

//==========================================================================
// for exchange in one direction
//==========================================================================
//==================================================================
void WavesScalarEqOpOpt::PlaneWaveTopFDM_(int k)
{

	double t = k * dt;

	PlaneWaveOpexpl planeWaveTop(&topBoundary, planeWave);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	WavesWaveEqBoundary absLow(&sdg, SDjLow, dt);
	WavesWaveEqBoundary absTop(&sdg, SDjHigh, dt);

	WavesMirrorBC mirrorPer(&perBoundary);
	WavesMirrorBC mirrorBack(&backBoundary);

	WavesMirrorBC mirrorLeft(&leftBoundary);
	WavesMirrorBC mirrorRight(&rightBoundary);

	timeStepperOuterFDM.apply(v_Outer, u_Outer);

//	cout << "BC on outer domain" << endl;
	// Apply BCs on outer domain
	planeWaveTop.applyPlaneWaveOpexpl(u_Outer, t, planeWave);

	cornerBC.apply(u_Outer);

	absLow.apply(v_Outer, u_Outer); // uOuter is set on boundary

	real t_max = 2 * M_PI / 5.0;

	if (t > t_max)
		absTop.apply(v_Outer, u_Outer); // uOuter is set on boundary

	mirrorLeft.apply(v_Outer, u_Outer);
	mirrorRight.apply(v_Outer, u_Outer);

	if (nsd == 3)
	{
		mirrorPer.apply(v_Outer, u_Outer);
		mirrorBack.apply(v_Outer, u_Outer);
	}

}

//===========================================================================

void WavesScalarEqOpOpt::ScalarFDM(int k)
{
	double t = k * dt;

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else if (USE_FDM)
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}

	if (USE_RHS)
	{
	//	cout << "apply rhs fdm" << endl;
		ApplyRHS(t);
		add_sol_and_function.apply(fn_Outer, u_Outer);

	}

	if (USE_DIRICHLET_FDM)
	{
	//	cout << "with Dirichlet b.c." << endl;
		dirichletBC.apply(u_Outer);
	//	cout << "after dirichlet fdm" << endl;
	}

	if (USE_ABSORB)
	{
		cornerBC.apply(u_Outer);
		abs1.apply(v_Outer, u_Outer);
		abs2.apply(v_Outer, u_Outer);
		abs3.apply(v_Outer, u_Outer);
		abs4.apply(v_Outer, u_Outer);
		abs5.apply(v_Outer, u_Outer);
		abs6.apply(v_Outer, u_Outer);
	}

}

//===========================================================================

void WavesScalarEqOpOpt::ScalarFDM_LENS(int k)
{
	double t = k * dt;

	WavesWaveEqInterior timeStepperOuter(&opt.outerWithHole, dt);
	WavesWaveEqInterior timeStepperOuterFDM(&sdg, dt);

	if (EXCHANGE)
	{
//		cout << "fdm for exchange" << endl;
		timeStepperOuter.apply(v_x_Outer, u_x_Outer);
		timeStepperOuter.apply(v_y_Outer, u_y_Outer);
		timeStepperOuter.apply(v_z_Outer, u_z_Outer);
		timeStepperOuter.apply(v_Outer, u_Outer);

	} // u is obtained, outer grid
	else if (USE_FDM)
	{
		timeStepperOuterFDM.apply(v_Outer, u_Outer);

	}
	if (USE_ABSORB)
	{
		cornerBC.apply(u_x_Outer);
		cornerBC.apply(u_y_Outer);
		cornerBC.apply(u_z_Outer);
		cornerBC.apply(v_x_Outer);
		cornerBC.apply(v_y_Outer);
		cornerBC.apply(v_z_Outer);

		abs1.apply(v_x_Outer, u_x_Outer);
		abs2.apply(v_x_Outer, u_x_Outer);
		abs3.apply(v_x_Outer, u_x_Outer);
		abs4.apply(v_x_Outer, u_x_Outer);
		abs5.apply(v_x_Outer, u_x_Outer);
		abs6.apply(v_x_Outer, u_x_Outer);

		abs1.apply(v_y_Outer, u_y_Outer);
		abs2.apply(v_y_Outer, u_y_Outer);
		abs3.apply(v_y_Outer, u_y_Outer);
		abs4.apply(v_y_Outer, u_y_Outer);
		abs5.apply(v_y_Outer, u_y_Outer);
		abs6.apply(v_y_Outer, u_y_Outer);

		abs1.apply(v_z_Outer, u_z_Outer);
		abs2.apply(v_z_Outer, u_z_Outer);
		abs3.apply(v_z_Outer, u_z_Outer);
		abs4.apply(v_z_Outer, u_z_Outer);
		abs5.apply(v_z_Outer, u_z_Outer);
		abs6.apply(v_z_Outer, u_z_Outer);

		cornerBC.apply(u_Outer);
		abs1.apply(v_Outer, u_Outer);
		abs2.apply(v_Outer, u_Outer);
		abs3.apply(v_Outer, u_Outer);
		abs4.apply(v_Outer, u_Outer);
		abs5.apply(v_Outer, u_Outer);
		abs6.apply(v_Outer, u_Outer);

	}

}

//========================================================================

void WavesScalarEqOpOpt::AdjScalarFEM(double t, int k, int type_of_material, MV_Vector<double>& Dif_sol, MV_Vector<double>& velocity)
{
	int nno = gg.getNoNodes();

	PetscScalar zero = 0.0;

	t = k * dt;

	PetscScalar timest = dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar one = 1.0;
	PetscScalar minus_dtdt = -1 * dt * dt;

	// in the method  evalAdjRHSforDifMat
	// we get Fn= (-tau�(M^L)^-1 )*S/(-tau^2) = (M^L)^-1 * S
	// here S is assembled rhs 

	if (USE_RHS)
	{
		evalAdjRHSforDifMat(Fn, t, Dif_sol, lumpedAmass, type_of_material, velocity);

	} // for use_rhs = true

	ierr = MatMult0(A, u1, u2); // u2 = A * u1 - u01

	if (USE_RHS)
//HRN ierr = VecAXPY(&minus_dtdt,Fn,u2); 	// u2 = u2 - dt*dt*Fn

		if (USE_RHS)
			ierr = VecAXPY(u2, minus_dtdt, Fn); // u2 = u2 - dt*dt*Fn

	PetscScalar* u2D;
	ierr = VecGetArray(u2,&u2D);
	CHKERRA(ierr);

	for (int i = 0; i < nno; i++)
		if (u2D[i] > 0.0)
			cout << " adj(" << i << ")=" << u2D[i] << " USE_RHS " << USE_RHS << endl;

	ierr = VecRestoreArray(u2,&u2D);
	CHKERRA(ierr);

	if (USE_DIRICHLET_FEM)
	{

		PetscScalar* u2D;
		ierr = VecGetArray(u2,&u2D);
		CHKERRA(ierr);
		PetscScalar* u1D;
		ierr = VecGetArray(u1,&u1D);
		CHKERRA(ierr);

		for (int i = 0; i < nbn; i++)
		{
			u2D[boundindex[i]] = bound_values[i];
			u1D[boundindex[i]] = bound_values[i];
		}

		ierr = VecRestoreArray(u2,&u2D);
		CHKERRA(ierr);
		ierr = VecRestoreArray(u1,&u1D);
		CHKERRA(ierr);
	}

}
//*******************************************************************

void WavesScalarEqOpOpt::AdjFEM_LENS(double t, int k, MV_Vector<double>& Dif_sol, MV_Vector<double>& velocity)
{
	int nno = gg.getNoNodes();

	PetscScalar zero = 0.0;

	PetscScalar timest = dt;
	PetscScalar dtdt = dt * dt;
	PetscScalar one = 1.0;
	PetscScalar minus_dtdt = -1 * dt * dt;

	// in the method  evalAdjRHSforDifMat
	// we get Fn= (-tau�(M^L)^-1 )*S/(-tau^2) = (M^L)^-1 * S
	// here S is assembled rhs 

	if (USE_RHS)
	{
		//   cout<<"works rhs for fem"<<endl; 

		evalAdjRHSforDifMat(Fn, t, Dif_sol, lumpedAmass, type_of_material, velocity);

	} // for use_rhs = true

	ierr = MatMult0(A, u1_x, u2_x);
	ierr = MatMult0(A, u1_y, u2_y);
	ierr = MatMult0(A, u1_z, u2_z);

	ierr = MatMult0(A, u1, u2); // u2 = A * u1 - u01

	if (USE_RHS)
	{
//HRN ierr = VecAXPY(&minus_dtdt,Fn,u2); 	// u2 = u2 - dt*dt*Fn
//    ierr = VecAXPY(&minus_dtdt,Fn,u2_x); 	// u2_x = u2_x - dt*dt*Fn
//    ierr = VecAXPY(&minus_dtdt,Fn,u2_y); 	// u2_y = u2_y - dt*dt*Fn
//    ierr = VecAXPY(&minus_dtdt,Fn,u2_z); 	// u2_z = u2_z - dt*dt*Fn
		ierr = VecAXPY(u2, minus_dtdt, Fn); // u2 = u2 - dt*dt*Fn
		ierr = VecAXPY(u2_x, minus_dtdt, Fn); // u2_x = u2_x - dt*dt*Fn
		ierr = VecAXPY(u2_y, minus_dtdt, Fn); // u2_y = u2_y - dt*dt*Fn
		ierr = VecAXPY(u2_z, minus_dtdt, Fn); // u2_z = u2_z - dt*dt*Fn

	}

}

//===================================================================
int WavesScalarEqOpOpt::print_GID_VEC_FEM(char* filename, int sch)
{

	int code, ierr, ii;

	PetscScalar *u_1_;

	//cout << " nonodes " << nno << endl;
	if (USE_FEM || EXCHANGE)
	{

		ierr = VecGetArray(u1, &u_1_);
		CHKERRQ( ierr);

	}

	if (nsd == 2)
		code = opt.print_2Dgid_result(filename, u_1_, u_1_, u_1_, sch);
	else if (nsd == 3)
		code = opt.print_3Dgid_result(filename, u_1_, u_1_, u_1_, u_1_, sch);

	if (USE_FEM || EXCHANGE)
	{
		ierr = VecRestoreArray(u1, &u_1_);
		CHKERRQ( ierr);
	}

	return ierr;
}
//===============================================================
//===================================================================
int WavesScalarEqOpOpt::print_GID_VEC_COMMON(char* filename, int sch)
{

	int code, ierr, ii;
	int nnogg = gg.getNoNodes();

	PetscScalar *u_1_;

	//cout << " nonodes " << nno << endl;

	if (USE_FEM || EXCHANGE)
	{

		ierr = VecGetArray(u1, &u_1_);
		CHKERRQ( ierr);

	}

	if (nsd == 2)
	{
	  // first, we print out results of the FEM solution on the inner FEM mesh
		code = opt.print_2Dgid_result_common(filename, u_1_, u_1_, u_1_, sch);

		// second, we print out results of FDM solution only at the nodes of the outer FDM mesh.
		GIDOutputOp gid_out2(&opt.outerWithHole, filename, u_Outer, u_Outer, u_Outer, sch, nnogg);

		bool pr = gid_out2.printResultsCommon();

	}
	else if (nsd == 3)
	  {

	    // first, we print out results of the FEM solution on the inner FEM mesh
	    code = opt.print_3Dgid_result_common(filename, u_1_, u_1_, u_1_, u_1_, sch);
	    
	    // second, we print out results of FDM solution only at the nodes of the outer FDM mesh.

	    GIDOutputOp3D   gid_out1(&opt.outerWithHole,
				     filename,u_Outer,
				     u_Outer, u_Outer,
				     u_Outer,
				     sch,nnogg);
	    
	    bool  pr = gid_out1.printResultsCommon();
	  }


	if (USE_FEM || EXCHANGE)
	{
		ierr = VecRestoreArray(u1, &u_1_);
		CHKERRQ( ierr);
	}

	return ierr;
}

//=================================================================
int WavesScalarEqOpOpt::print_GID_VEC_FEM_LENS(char* filename, int sch)
{

	int code, ierr, ii;

	PetscScalar *u_x_;
	PetscScalar *u_y_;
	PetscScalar *u_z_;

	//cout << " nonodes " << nno << endl;
	if (USE_FEM || EXCHANGE)
	{

		ierr = VecGetArray(u1_x, &u_x_);
		CHKERRQ( ierr);
		ierr = VecGetArray(u1_y, &u_y_);
		CHKERRQ( ierr);
		ierr = VecGetArray(u1_z, &u_z_);
		CHKERRQ( ierr);

	}

	if (nsd == 3)
		code = opt.print_3Dgid_LENS(filename, u_x_, u_y_, u_z_, sch);

	if (USE_FEM || EXCHANGE)
	{
		ierr = VecRestoreArray(u1_x, &u_x_);
		CHKERRQ( ierr);
		ierr = VecRestoreArray(u1_y, &u_y_);
		CHKERRQ( ierr);
		ierr = VecRestoreArray(u1_z, &u_z_);
		CHKERRQ( ierr);
	}

	return ierr;
}

//==============================================================
int WavesScalarEqOpOpt::print_GID_VEC_FDM(char* filename, int sch)
{

	int code, ierr, ii;
	bool pr;
	//cout << " nsd == " << nsd << endl;

	// HRN
	ierr = 0;

	if (nsd == 2)
	{

		if (USE_FDM)
		{

			GIDOutputOp gid_out1(&sdg, filename, u_Outer, u_Outer, u_Outer, sch);

			pr = gid_out1.printResults();

		}

		if (EXCHANGE)
		{
			GIDOutputOp gid_out2(&opt.outerWithHole, filename, u_Outer, u_Outer, u_Outer, sch);

			pr = gid_out2.printResults();

		}
	}
	else if (nsd == 3)
	{

		if (USE_FDM)
		{
			GIDOutputOp3D gid_out1(&sdg, filename, u_Outer, u_Outer, u_Outer, u_Outer, sch);

			pr = gid_out1.printResults();

		}

		if (EXCHANGE)
		{
			GIDOutputOp3D gid_out1(&opt.outerWithHole, filename, u_Outer, u_Outer, u_Outer, u_Outer, sch);
			pr = gid_out1.printResults();

		}
	}

	return ierr;
}

