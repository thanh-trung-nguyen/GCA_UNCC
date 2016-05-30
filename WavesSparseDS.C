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

#include "include/wavesSparseDS.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

WavesSparseDS::WavesSparseDS()
{
	nrows = ncolumns = nonzeroes = jcol_length = 0;
	diagonal_computed = false;
}

WavesSparseDS::~WavesSparseDS()
{
}

bool WavesSparseDS::ok() const
{
	return (getNoNonzeroes() == jcol_arr.size() && nrows) ? true : false;
}

bool WavesSparseDS::operator ==(const WavesSparseDS& X)
{
	int i;
	if (this == &X) // the same object ?
		return true;
	else if (nrows != X.nrows)
		return false;
	else if (ncolumns != X.ncolumns)
		return false;
	else if (nonzeroes != X.nonzeroes)
		return false;
	else
	{
		int nrows1 = nrows + 1;
		for (i = 0; i < nrows1; ++i)
			if (irow_arr(i) != X.irow_arr(i))
				return false;

		for (i = 0; i < nonzeroes; ++i)
			if (jcol_arr(i) != X.jcol_arr(i))
				return false;

	}
	return true; // If you arrive here, the patterns are identical
}


bool WavesSparseDS::consistent()

{
	if ((nonzeroes < 0) || (nonzeroes > jcol_length))
	{
		printf("hej1");
		return false;
	}

	if ((nrows < 0) || (ncolumns < 0))
	{
		printf("hej2");
		return false;
	}
	if (nonzeroes != irow_arr(nrows))
	{
		printf("hej3 %d %d %d %d\n", nonzeroes, nrows, irow_arr(nrows), irow_arr(nrows - 1));
		return false;
	}

	int i;

	int prev_row = -1;
	int irowarri; //help variable

	for (i = 0; i < nrows; ++i)
	{
		irowarri = irow_arr(i);
		if ((irowarri < 0) || (irowarri > nonzeroes - 1))
		{
			printf("hej4\n");
			return false;
		}

		if (irowarri <= prev_row)
		{ // Check that irow_arr is in ascending order
			printf("hej5 %d %d %d\n", i, irowarri, prev_row);
			return false;
		}
		prev_row = irowarri;
	}

	int jstart, jstop, j;
	int prev_col;
	int jcolarrj; //help variable

	for (i = 0; i < nrows; ++i)
	{
		jstart = irow_arr(i);
		jstop = irow_arr(i + 1);

		prev_col = -1;
		for (j = jstart; j < jstop; ++j)
		{
			jcolarrj = jcol_arr(j);
			if ((jcolarrj < 0) || (jcolarrj > ncolumns - 1))
			{
				printf("hej6\n");
				return false;
			}

			// Check whether the column numbers are in
			// ascending order within each row?
			if (jcolarrj <= prev_col)
			{
				printf("hej7\n");
				return false;
			}
			prev_col = jcolarrj;
		}
	}
	return true; // If you arrive here, the pattern is consistent
}


bool WavesSparseDS::redim(int nrows_, int ncolumns_, int max_nonzeroes)

{
	bool b = ((this->nrows != nrows_) || (jcol_length != max_nonzeroes)) ? true : false;

	nrows = nrows_;
	ncolumns = ncolumns_;

	irow_arr.newsize(nrows + 1);
	nonzeroes = max_nonzeroes;

	jcol_length = max_nonzeroes;
	jcol_arr.newsize(jcol_length);
	diagonal_computed = false; // no diag_arr is computed yet
	diag_arr.newsize(0); // remove diag_arr if the WavesSparseDS is redim-ed
	return b;
}


bool WavesSparseDS::redimIrow(int nrows_)

{
	// redim only the irow part, let the user compute the jcol part by some
	// list tools and then fill jcol by this list using fillJcol

	bool b = (this->nrows != nrows_) ? true : false;
	nrows = nrows_;
	// here we simply assume that nrows and ncolumns are equal (that is why
	// ncolumns is not an argument to this function, cf. redim(int,int,int)
	ncolumns = nrows_;
	irow_arr.newsize(nrows + 1);
	diagonal_computed = false; // no diag_arr is computed yet
	diag_arr.newsize(0); // remove diag_arr if the WavesSparseDS is redim-ed
	return b;
}


bool WavesSparseDS::redimIrow(int nrows_, int ncolumns_)

{
	// redim only the irow part, let the user compute the jcol part by some
	// list tools and then fill jcol by this list using fillJcol

	bool b = (this->nrows != nrows_) ? true : false;
	nrows = nrows_;
	ncolumns = ncolumns_; // this differs from the redimIrow(int) 
	irow_arr.newsize(nrows + 1);
	diagonal_computed = false; // no diag_arr is computed yet
	diag_arr.newsize(0); // remove diag_arr if the WavesSparseDS is redim-ed
	return b;
}


bool WavesSparseDS::fillJcol(MV_Vector<int>& jcol_vector, const int no_nonzero)

{
	// redimIrow must be called prior to fillJcol!
	// general use:  redimIrow - compute jcol_vector - fillJcol

	bool b = (jcol_arr.size() != no_nonzero) ? true : false;

	if (irow(nrows) != no_nonzero || jcol_vector.size() < no_nonzero)
	{
		cout << "WavesSparseDS::fillJcol irow and jcol are not consistent.\n" << flush;
	}

	jcol_arr.newsize(no_nonzero); //redim and copy

	for (int i = 0; i < no_nonzero; i++)
		jcol_arr(i) = jcol_vector(i);
	nonzeroes = jcol_length = jcol_arr.size(); // exact length
	diagonal_computed = false; //the diagonal is not yet computed

	return b;
}


void WavesSparseDS::endInit()

{
	nonzeroes = irow(nrows) - 1;
}


bool WavesSparseDS::computeDiagonal()

{
// Compute the indices of the diagonal entries in WavesSparseDS.
// If the diagonal is not found for each row this routine will return false.

	bool diagexists = true;
	diagonal_computed = false; // must be off when idx is called
	diag_arr.newsize(nrows);
	int i, idxii;
	//  for( i=1; i <= nrows; i++ ){
	for (i = 0; i < nrows; i++)
	{
		idxii = idx(i, i);
		diag(i) = idxii;
		if (!idxii)
			diagexists = false; // if found a idx(i,i) which is zero
	}
	diagonal_computed = true;
	return diagexists;
}


int WavesSparseDS::idx(int i, int j) const

{
#ifdef DP_DEBUG
	if (!this->ok())
	errorFP("WavesSparseDS::idx",
			"No sparsity pattern not ok, cannot compute index");
#endif
	if ((i == j) && diagonalComputed())
		return diag(i);

// find index by the bisection method, the algorithm assumes that the
// column indices in a row are given in increasing order

	int jl = irow(i);
	int ju = irow(i + 1);

	int jm, jam;
	while ((ju - jl) >= 1)
	{
		jm = (ju + jl) >> 1;
		jam = jcol(jm);
		if (j == jam)
			return jm;
		if (j > jam)
			jl = jm + 1;
		else
			ju = jm;
	}
	// if the loop is terminated the index j  was not found:
	return -1;
}


void WavesSparseDS::printPattern() const

{
	cout << "nrows=" << nrows << "\n";
	cout << "ncolumns=" << ncolumns << "\n";
	cout << "jcol_length=" << jcol_length << "\n";
	cout << "nonzeroes=" << nonzeroes << "\n";

	cout << "irow_arr has " << irow_arr.size() << "elements\n";
	cout << "jcol_arr has " << jcol_arr.size() << "elements\n" << flush;

	int i;
	for (i = 0; i < irow_arr.size(); i++)
		cout << "i=" << i << "  irow_arr(i)=" << irow_arr(i) << "\n";

	for (i = 0; i < jcol_arr.size(); i++)
		cout << "i=" << i << "  jcol_arr(i)=" << jcol_arr(i) << "\n";
	cout << flush;
}

/*
 1 2 3 4 5 6 7 8 9
 -----|    double_horizontal_space
 1| X X              |
 2| X X X            |
 3|   X X X          |
 4|     X X X        |
 5|       X X X      |
 6|         X X X    |
 7|           X X X  |
 8|             X X X|
 9|               X X|




 123456789
 ----
 1|XX       |
 2|XXX      |
 3| XXX     |
 4|  XXX    |
 5|   XXX   |
 6|    XXX  |
 7|     XXX |
 8|      XXX|
 9|       XX|
 ------
 */

/* LOG HISTORY of this file:

 $Log: WavesSparseDS.C,v $
 * Revision 1.1  1996/11/15  12:03:51  job
 * Version 2.4.0
 *
 * Revision 1.4  1996/07/22  13:06:27  job
 * Version 2.2.0.
 *
 */
