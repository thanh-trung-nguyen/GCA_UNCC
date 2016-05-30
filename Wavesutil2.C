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

#include "include/wavesutil2.h"

void copyValues(const real *x, real *y, const int *IX, const int *IY, const int &n)
{

	for (int i = 0; i < n; i++)
		y[IY[i]] = x[IX[i]];

}

void copy2(real *x, real *y, const int *indx, const int *indy, const int &n)
{

	for (int i = 0; i < n; i++)
		y[indy[i]] = (y[indy[i]] + x[indx[i]]) / 2;

}

// ======================================================================

void sort(int *ia, const int &n)
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
