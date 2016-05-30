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

#ifndef __WAVESSPARSEDS_H
#define __WAVESSPARSEDS_H

#include "wavesMVvtp.h"
#include <math.h>

#ifdef CC_42 /* bool is predefined in gcc, but not in other C compilers */
typedef int bool;
#endif
#define true  1
#define false 0

/*<WavesSparseDS:*/

class WavesSparseDS

{

	private:

		int nrows; // no of rows in matrix
		int ncolumns; // no of columns in matrix

		MV_Vector<int> irow_arr; // representation of irow-function
		MV_Vector<int> jcol_arr; // representation of jcol-function

		int jcol_length; // >= nonzeros, used when allocating jcol_arr
		int nonzeroes;

		MV_Vector<int> diag_arr; // index of diagonal elements in the jcol
		bool diagonal_computed;

	public:

		WavesSparseDS();

		WavesSparseDS(int n, int max_nonzeroes)
		{
			redim(n, max_nonzeroes);
		}

		WavesSparseDS(int nrows_, int ncolumns_, int max_nonzeroes)
		{
			redim(nrows_, ncolumns_, max_nonzeroes);
		}

		~WavesSparseDS();

		bool operator ==(const WavesSparseDS& X);

		bool consistent();
		bool diagonalComputed() const
		{
			return diagonal_computed;
		}
		bool computeDiagonal();
		bool removeDiagonalInfo()
		{
			diagonal_computed = false;
			bool empty = diag_arr.null();
			diag_arr.newsize(0);
			return empty;
		} // To be called when the diagonal indices are
		  // inaccurate or not needed more

		bool ok() const; // tests jcol's length vs getNoNonzeroes
		bool redim(int n, int max_nonzeroes)
		{
			return redim(n, n, max_nonzeroes);
		}
		bool redim(int nrows, int ncolumns, int max_nonzeroes);

		bool redimIrow(int nrows); // redim irow, wait for jcol
		bool redimIrow(int nrows, int ncolumns); // redim irow, set ncolumns
		bool fillJcol(MV_Vector<int>& jcol_vector, const int no_nonzero);
// add jcol vector where no_nonzero is the number of nonzeroes

		int& irow(int i)
		{
			return irow_arr(i);
		}
		int& jcol(int i)
		{
			return jcol_arr(i);
		}
		int& diag(int i)
		{
			return diag_arr(i);
		}

		int irow(int i) const
		{
			return irow_arr(i);
		}
		int jcol(int i) const
		{
			return jcol_arr(i);
		}
		int diag(int i) const
		{
			return diag_arr(i);
		}

		int idx(int i, int j) const;

		int getNoRows() const
		{
			return nrows;
		}
		int getNoColumns() const
		{
			return ncolumns;
		}

		int maxNonZeroes() const
		{
			return jcol_length;
		}
		int getNoNonzeroes() const
		{
			return nonzeroes;
		}
		int& getNoNonzeroes()
		{
			return nonzeroes;
		}
		void endInit(); // update nonzeroes

		void printPattern() const;
};


/*Class:WavesSparseDS

 NAME:  "WavesSparseDS" - data structure for general sparse matrix storage

 SYNTAX:     @WavesSparseDS


 KEYWORDS:

 sparse matrix, compressed sparse row storage


 DESCRIPTION:

 The class defines a data structure and a public interface for representing
 a general sparse matrix. All the nonzero entries of the matrix are stored
 row by row in a vector. The matrix entries are not a part of this
 class. Two indexing functions are used for locating
 a matrix entry. It is these two functions, and related procedures,
 that are offered by class "WavesSparseDS".
 The sparse matrix storage scheme is described in detail below.

 Assume that all the `$`n`_`z`$` nonzeroes in an "nrows" by "ncolumns"
 matrix are stored in an array `$`A`_`s`$`, which will typically be a
 part of class "MatSparse". The information needed to subscript a sparse
 matrix is contained in class "WavesSparseDS". To access the matrix entries
 two address arrays "irow" and "jcol" are needed. Here "irow(k)" gives the
 address in `$`A`_`s`$` of the first nonzero entry in row number
 `$`k`$` in the coefficient matrix, while "jcol(l)" gives the matrix column
 number of the entry `$`A`_`s(l)`$`.

 $The following example illustrates this scheme for $n=5$:
 $\[\left(  \begin{array}{ccccc}
 $ a_{1,1}  &   0   &    0  &  a_{1,4} &   0     \\
  $   0     & a_{2,2}& a_{2,3}&   0     & a_{2,5} \\
  $   0     & a_{3,2}& a_{3,3}&   0     & 0       \\
  $ a_{4,1}  &   0   &    0  &  a_{4,4} & a_{4,5} \\
  $   0     & a_{5,2}&    0  &  a_{5,4} & a_{5,5} \end{array}\right) \]
 $\begin{eqnarray*}
 $A_s  & = &(a_{1,1},a_{1,4},a_{2,2},a_{2,3},a_{2,5},a_{3,2},a_{3,3},a_{4,1},
 $a_{4,4},a_{4,5}, a_{5,2},a_{5,4},a_{5,5}) \\
  $\mbox{\tt irow}  &=& (1,3,6,8,11,14) \\
  $\mbox{\tt jcol}  &=& (1,4,2,3,5,2,3,1,4,5,2,4,5)
 $\end{eqnarray*}
 Notice that "irow" has dimension "nrows"+1 whereas `$`A`_`s`$`
 and "jcol" have dimension
 `$`n`_`z = `\mbox{`"irow"`}`(`\mbox{`"nrows"`}`+1)-1`$`.


 CONSTRUCTORS AND INITIALIZATION:

 There are three constructors. With no arguments, the "redim" function must
 be called later to ensure that memory is allocated for the object.
 There is a constructor that takes the number of rows (and columns) in the
 matrix and the maximum number of nonzeroes as arguments.
 With these constructors or "redim" functions, only the internal data
 structures are allocated, and the user has to use the functions
 "irow", "jcol" and "getNoNonzeroes" to initialize the sparsity pattern
 (this is often done by the "makeSparsityPattern" functions).
 Another way of initializing the object is to use the constructor with no
 arguments and then the "redimIrow" function to allocate storage for the
 ""irow"" part of the data structure (the length of this is known when the
 number of rows in the matrix is known). There after the programmer can
 make a list "List(int,ListItemInst)" and add column (""jcol"") entries
 to this list in a flexible and storage minimizing way, along with filling
 ""irow"" values using the function "WavesSparseDS::irow". Finally, the list
 is submitted to the "WavesSparseDS" object by calling the "fillJcol" function.
 In practice, the usage of "redimIrow" and "fillJcol" is the most flexible
 way of creating a "WavesSparseDS" object. See also the documentation of the
 "redimIrow" function.


 MEMBER FUNCTIONS:

 Most member functions are self-explanatory. Some additional comments are
 given below:

 $\begin{description}

 `\item[`"consistent"`]`  -
 checks the current object to see if the index information is valid.

 `\item[`"redim"`]`  -
 sets the new size for the object. There are three arguments to this
 function: the number of rows ("nrows") and columns ("ncolumns") in
 the matrix and the maximum number of nonzeroes
 ("max_nonzeroes") allowed in the matrix. In case of a
 square matrix, there is an overloaded version of "redim" that
 accepts a single matrix dimension parameter "n" in addition to the
 maximum number of nonzeroes.
 Note that there are two ways of computing the sparsity pattern when
 using these "redim" functions: 1) filling the "irow" and "jcol" 
 functions directly, or 2) compute the sparsity pattern in local data
 structures and then fill the "WavesSparseDS" object by calling "irow"
 and "jcol". In case 1, one must have conservative estimates of the 
 number of nonzeroes (this can often lead to waste of memory). In case
 2, the sparsity pattern is known prior to allocating data in the
 "WavesSparseDS" object such that the amount of storage is exact.
 An way of doing this is to allocate the ""irow"" part
 of the data structures in "WavesSparseDS" and fill this with the "irow"
 function, and compute the column (""jcol"") entries locally using
 a list. See below how this can be done.
 
 `\item[`"redimIrow"`]`  -
 allocates storage for the ""irow"" part of the data structure.
 This can be done before the sparsity pattern is known, and during the
 computation of the sparsity pattern, one can then easily fill in
 ""irow"" values using the "irow" function. Simultaneously, the 
 programmer can use a local list "List(int,ListItemList)" to compute
 the column (""jcol"") numbers. Since the number of nonzeroes is not
 known on beforehand, a local list is a flexible way of storing the
 column number information. After the list is computed, it is
 submitted to the "WavesSparseDS" object by calling "fillJcol" and the
 "WavesSparseDS" is completely initialized. Using "redimIrow" and "fillJcol"
 the "WavesSparseDS" object will allocated the exact amount of memory needed
 to store the sparsity pattern. If the "redim(int)" or "redim(int,int)"
 function is used, the internal arrays in "WavesSparseDS" must be large
 enough to hold the expected sparsity pattern. In an overloaded version
 of "redimIrow" the number of columns is set, the default version 
 assumes that the number of columns equals the number of rows.

 `\item[`"fillJcol"`]`  -
 see the documentation of "redimIrow".

 `\item[`"getNoRows" `{\rm ` and `}` "getNoColumns"`]`  -
 return the number of rows and columns, respectively.

 `\item[`"maxNonZeroes"`]`  -
 returns the maximum number of nonzero elements in the current
 "WavesSparseDS" object.

 `\item[`"getNoNonzeroes"`]`  -
 returns the number of nonzeros actually stored by now. It can also
 be used to set this value.

 `\item[`"endInit"`]`  -
 to be called when the sparsity pattern is fully loaded.

 `\item[`"printPattern"`]`  -
 prints a graphical representation of the sparsity pattern on the
 specified output stream. (Works only for reasonable small matrices).

 `\item[`"print"`]`  -
 prints the sparsity pattern data on the specified output stream.

 `\item[`"scan"`]`  -
 scans the sparsity pattern data from the specified input stream.

 `\item[`"getIrowPtr0"`]`  - 
 returns access to the underlying C array
 for the ""irow"" part of the "WavesSparseDS" structure.
 The "0" in the name indicates that the returned pointer is the
 base of a C array where "[0]" is the first element.
 This is the function to be used when, for example, communicating
 with Fortran programs.

 `\item[`"getIrowPtr1"`]`  - 
 as "getIrowPtr0", but the returned pointer has "[1]" has the first
 entry in the array.

 `\item[`"getJcolPtr0"`]`  - 
 as "getIrowPtr0", but access to the ""jcol"" array is returned.

 `\item[`"getJcolPtr1"`]`  - 
 as "getJcolPtr0", but the returned pointer has "[1]" as the first
 entry.


 $\end{description}

 EXAMPLES:

 SEE ALSO:

 class "MatSparse(Type)", class "prm(Matrix(Type))".

 DEVELOPED BY:

 SINTEF Applied Mathematics, Oslo, Norway, and
 University of Oslo, Dept. of Mathematics, Norway

 AUTHOR:

 Original code by Hans Petter Langtangen, UiO/SINTEF Applied Mathematics.
 Modified and extended by Are Magnus Bruaset, SINTEF Applied Mathematics.
 Further extensions and optimizations by Klas Samuelsson, UiO.

 End:
 */

/*
 #define ClassType WavesSparseDS 
 #include <Handle.h>
 #undef  ClassType
 */

#endif

