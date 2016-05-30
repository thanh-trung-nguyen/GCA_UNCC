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

#ifndef __WAVESUTIL2_H
#define __WAVESUTIL2_H

typedef double real;

void copyValues(const real *x, real *y, const int *IX, const int *IY, const int &n);

void copy2(real *x, real *y, const int *indx, const int *indy, const int &n);

/*
 Sort values in integer array with n elements
 */
void sort(int *II, const int &n);

#endif
