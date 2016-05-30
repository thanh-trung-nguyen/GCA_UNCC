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

#ifndef __WAVESMULTIDIMARRAYS_H
#define __WAVESMULTIDIMARRAYS_H

#include "wavesMVvtp.h"

template<class TYPE>
class Array2d
{

	public:

		// Constructors & destructor
		inline Array2d() :
				rows(0), cols(0)
		{
		}

		inline Array2d(unsigned int i, unsigned int j) :
				rows(i), cols(j), vec(i * j)
		{
		}

		~Array2d()
		{
		}

		// Access & indices
		inline TYPE& operator()(unsigned int i, unsigned int j)
		{
			return vec(j * rows + i);
		}

		inline const TYPE& operator()(unsigned int i, unsigned int j) const
		{
			return vec(j * rows + i);
		}

		inline
		int size(int dim) const
		{
			if (dim == 0)
				return rows;
			else
				return cols;
		}

		inline Array2d<TYPE>& newsize(unsigned int i, unsigned int j)
		{
			if (rows != i || cols != j)
			{
				rows = i;
				cols = j;
				vec.newsize(i * j);
			}
			return *this;
		}

		// Assignment
		inline Array2d<TYPE>& operator=(const TYPE& alfa)
		{
			vec = alfa;
			return *this;
		}

		inline Array2d<TYPE>& operator=(const Array2d<TYPE>& cpvec)
		{
			if (rows != cpvec.size(0) || cols != cpvec.size(1))
			{
				rows = cpvec.size(0);
				cols = cpvec.size(1);
				vec.newsize(cpvec.size(0) * cpvec.size(1));
			}
			vec = cpvec.vec;
			return *this;
		}

	private:

		int rows;
		int cols;
		MV_Vector<TYPE> vec;

};

template<class TYPE>
class Array3d
{

	public:

		// Constructors & destructor
		inline Array3d() :
				dim1(0), dim2(0), dim3(0)
		{
		}

		inline Array3d(unsigned int i, unsigned int j, unsigned int k) :
				dim1(i), dim2(j), dim3(k), vec(i * j * k)
		{
		}

		~Array3d()
		{
		}

		// Access & indices
		inline TYPE& operator()(unsigned int i, unsigned int j, unsigned int k)
		{
			return vec(k * dim1 * dim2 + j * dim1 + i);
		}

		inline const TYPE& operator()(unsigned int i, unsigned int j, unsigned int k) const
		{
			return vec(k * dim1 * dim2 + j * dim1 + i);
		}

		inline
		int size(int dim) const
		{
			switch (dim)
			{
				case 0:
					return dim1;
				case 1:
					return dim2;
				case 2:
					return dim3;
				default:
					assert(0);
					break;
			}
			return 0;
		}

		inline Array3d<TYPE>& newsize(unsigned int i, unsigned int j, unsigned int k)
		{
			if (dim1 != i || dim2 != j || dim3 != k)
			{
				dim1 = i;
				dim2 = j;
				dim3 = k;
				vec.newsize(i * j * k);
			}
			return *this;
		}

		// Assignment
		inline Array3d<TYPE>& operator=(const TYPE& alfa)
		{
			vec = alfa;
			return *this;
		}

		inline Array3d<TYPE>& operator=(const Array3d<TYPE> cpvec)
		{
			int i = cpvec.size(0);
			int j = cpvec.size(1);
			int k = cpvec.size(2);
			if (dim1 != cpvec.size(0) || dim2 != j || dim3 != k)
			{
				dim1 = i;
				dim2 = j;
				dim3 = k;
				vec.newsize(i * j * k);
			}
			vec = cpvec.vec;
			return *this;
		}

	private:

		int dim1;
		int dim2;
		int dim3;
		MV_Vector<TYPE> vec;

};

template<class TYPE>
class Array4d
{

	public:

		// Constructors & destructor
		inline Array4d() :
				dim1(0), dim2(0), dim3(0), dim4(0)
		{
		}

		Array4d(unsigned int i, unsigned int j, unsigned int k, unsigned int l) :
				dim1(i), dim2(j), dim3(k), dim4(l), vec(i * j * k * l)
		{
		}

		~Array4d()
		{
		}

		// Access & indices
		inline TYPE& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
		{
			return vec(l * dim1 * dim2 * dim3 + k * dim1 * dim2 + j * dim1 + i);
		}

		inline const TYPE& operator()(unsigned int i, unsigned int j, unsigned int k, unsigned int l) const
		{
			return vec(l * dim1 * dim2 * dim3 + k * dim1 * dim2 + j * dim1 + i);
		}

		inline
		int size(int dim) const
		{
			switch (dim)
			{
				case 0:
					return dim1;
				case 1:
					return dim2;
				case 2:
					return dim3;
				case 3:
					return dim4;
			}
		}

		inline Array4d<TYPE>& newsize(unsigned int i, unsigned int j, unsigned int k, unsigned int l)
		{
			if (dim1 != i || dim2 != j || dim3 != k || dim4 != l)
			{
				dim1 = i;
				dim2 = j;
				dim3 = k;
				dim4 = l;
				vec.newsize(i * j * k * l);
			}
			return *this;
		}

		// Assignment
		inline Array4d<TYPE>& operator=(const TYPE& alfa)
		{
			vec = alfa;
			return *this;
		}

	private:

		int dim1;
		int dim2;
		int dim3;
		int dim4;
		MV_Vector<TYPE> vec;

};

#endif
