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

#ifndef __WAVESSDOPERATOR_H
#define __WAVESSDOPERATOR_H

#include "wavesSDDefs.h"
#include "wavesSDIndexes.h"
#include "wavesSDGeometry.h"

/** 
 * WavesSDOperator is an abstract base class. An WavesSDOperator may be applied on a part of a WavesSDGeometry,
 * described by an object of WavesSDIndexes type. An WavesSDOperator may perform anything from a discretization of a
 * partial differential equation operator to output or measuring a norm or just copying values.
 *
 * The WavesSDOperator contains an association to an WavesSDIndexes object. This must be supplied to the WavesSDOperator. If a subclass allocates it, it shall deallocate it, too.
 *
 * Typically, subclasses are finite difference discretizations that may use the members {\tt weights} and {\tt offsets}. These are freed by the base class.
 */
class WavesSDOperator
{
	protected:
		/// Default for freeIndexes is false, which will not free the memory of the WavesSDIndexes object.
		WavesSDOperator(WavesSDIndexes *sdi_) :
				indexes(sdi_), weights(0), offsets(0), nrOfWeights(0), unary(false)
		{
			assert(indexes != 0);
		}
	public:

		WavesSDOperator();

		///
		virtual ~WavesSDOperator()
		{
			delete[] weights;
			delete[] offsets;
		}

		/// Typically, {\tt y} is {\tt D(x)}, where {\tt D} is an operator. Calls doApply.
		bool apply(const real *x, real *y) const
		{
			return doApply(x, y);
		}

		/// If {\tt unary} is {\tt true}, no {\tt x} is needed to operate on {\tt y}. Calls doApply.
		bool apply(real *y) const
		{
			if (unary)
				return doApply(0, y);
			else
				return false;
		}

		/// Sometimes useful not to present indexes.
		virtual void presentOperator(bool presentIndexes = true);

		/// Needed for Jacobi preconditioning. {\em Something more general later?}
		virtual real getDiagonalEntry() const
		{
			if (nrOfWeights > 0)
				return weights[0];
			return 0.0;
		}

	protected:

		/// Has to be provided.
		virtual bool doApply(const real *x, real *y) const = 0;
		// virtual bool doApplyElastic(real *x1, real *y1,real *x2, real *y2) ;

		/// Association to indexes.
		WavesSDIndexes *indexes;

		/// Aggregation.
		real* weights;
		/// Aggregation.
		int* offsets;
		///
		int nrOfWeights;

		//
		bool unary;

	private:
		/// WavesSDOperator is not intended to be copied
		WavesSDOperator(const WavesSDOperator&);
		///
		WavesSDOperator &operator=(const WavesSDOperator&);

};

//======================================================================

/**
 WavesDplusDminusOp supplies DplusDminus-stencil in the prescribed direction.
 Remark: The Indexes that are  provided must take the operator width (of 1) into account.
 */

class WavesDplusDminusOp : public WavesSDOperator
{
	public:
		// If constructed from WavesSDGeometry, an WavesSDIndexes object covering the interior is created. The direction d is either x, y, or z, prescribed by d=0, 1 or 2.
		WavesDplusDminusOp(WavesSDGeometry *sdg, const int &d);

		/// WavesDplusDminusOp on sdi in direction 0,1 or 2, for x, y, and z directions.
		WavesDplusDminusOp(WavesSDIndexes *sdi, const int &d);

		~WavesDplusDminusOp()
		{
			if (freeIndexes)
				delete indexes;
		}
	protected:
		///
		bool doApply(const real *x, real *y) const;

		// utility method
		void doInitialize(const int &d);

		//
		bool freeIndexes;
};

//======================================================================

/**
 DplusDminusVec2Op supplies DplusDminus-stencil in the prescribed direction.
 It is unary, acting on a vector with two components, such that
 component 2 becomes DplusDminus of component 1. 
 */

class DplusDminusVec2Op : public WavesSDOperator
{
	public:
		//# If constructed from WavesSDGeometry, an WavesSDIndexes object covering the interior is created. The direction d is either x, y, or z, prescribed by d=0, 1 or 2.
		//# DplusDminusVec2Op(WavesSDGeometry *sdg, const int &d);

		/// DplusDminusVec2Op on sdi in direction 0,1 or 2 for x, y, and z directions.
		DplusDminusVec2Op(WavesSDIndexes *sdi, const int &d);

	protected:
		///
		bool doApply(const real *, real *y) const;

		// utility method
		void doInitialize(const int &d);

		/// The number of components in the vector. It is 2 here.
		int nComp;

};

//======================================================================

/**
 LaplacianOp supplies fivepointstencil in 2d and sevenpointstencil in 3d.
 */

class LaplacianOp : public WavesSDOperator
{
	public:

		LaplacianOp();

		// If constructed from WavesSDGeometry, an WavesSDIndexes object covering the interior is created.
		LaplacianOp(WavesSDGeometry *sdg);
		/// Laplaican operator on sdi. Remark: Should take the operator width (of 1) into account.
		LaplacianOp(WavesSDIndexes *sdi);

		// Destructor frees indexes if necessary
		~LaplacianOp()
		{
			if (freeIndexes)
				delete indexes;
		}
	protected:
		///
		bool doApply(const real *x, real *y) const;

		// utility method
		void doInitialize();

		//
		bool freeIndexes;
};

//======================================================================

/**
 DirichletOp supplies Dirichlet BCs. The value is a constant. The operator is unary.
 */
class DirichletOp : public WavesSDOperator
{
	public:

		/// Default is all boundaries. An WavesSDIndexes object is created.
		DirichletOp(WavesSDGeometry *sdg, real value = 0.0);

		/// Index may be provided. May be any domain.
		DirichletOp(WavesSDIndexes *sdi, real value = 0.0);

		/// Destructor frees sdi if constructor 1 used
		~DirichletOp()
		{
			if (freeIndexes)
				delete indexes;
		}
	protected:

		///
		bool doApply(const real *x, real *y) const;

		//
		void doInitialize();

	private:
		///
		real value;
		///
		bool freeIndexes;
};

// ======================================================================

/**
 WavesAssignmentOp assigns a value to {\tt y}. The operator is unary.
 The value may be changed by the methods.
 */

class WavesAssignmentOp : public WavesSDOperator
{
	public:

		///
		WavesAssignmentOp(WavesSDGeometry *sdg, real value = 0.0);
		///
		WavesAssignmentOp(WavesSDIndexes *sdi, real value = 0.0);

		/// Destructor frees sdi if constructor 1 used
		~WavesAssignmentOp()
		{
			if (freeIndexes)
				delete indexes;
		}

		/// Assign value and call doApply
		bool assignValue(real *y, const real &v)
		{
			value = v;
			return doApply(0, y);
		}

		///
		void setValue(const real &v)
		{
			value = v;
		}
		///
		real getValue()
		{
			return value;
		}

	protected:

		///
		bool doApply(const real *, real *y) const;

	private:
		///
		real value;

		///
		bool freeIndexes;
};

//=======================================================================
/**
 WavesAssignmentOp assigns a value to {\tt y}. The operator is unary.
 The value may be changed by the methods.
 */
//=======================================================================
class WavesAssignmentArrayOp : public WavesSDOperator
{
	public:

		///
		WavesAssignmentArrayOp(WavesSDGeometry *sdg, double* array);
		///
		WavesAssignmentArrayOp(WavesSDIndexes *sdi, double* array);

		/// Destructor frees sdi if constructor 1 used
		~WavesAssignmentArrayOp()
		{
			if (freeIndexes)
				delete indexes;
		}

		/// Assign value and call doApply
		bool assignArrayValue(real *y, real* array_)
		{
			array = array_;
			return doApply(0, y);
		}

		bool InnerBoundIndex(MV_Vector<int>& boundindex) const;

	protected:

		///
		bool doApply(const real *, real *y) const;

	private:
		///
		real* array;

		///
		bool freeIndexes;
};

//======================================================================

///
typedef real (*Fcn3d)(real x, real y, real z);

/** WavesAssignFunctionOp assigns values to {\tt y} as a function of the coordinates. 
 *  The operator is unary. The constructor uses a function pointer according to
 \begin{verbatim}
 typedef real (*Fcn3d) (real x, real y, real z);
 \end{verbatim}
 */
class WavesAssignFunctionOp : public WavesSDOperator
{
	public:

		///
		WavesAssignFunctionOp(WavesSDGeometry *sdg, Fcn3d f);

		///
		WavesAssignFunctionOp(WavesSDIndexes *sdi, Fcn3d f);

		/// Destructor frees sdi if constructor 1 used
		~WavesAssignFunctionOp()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		bool freeIndexes;
		///
		Fcn3d f;
};

//======================================================================

///typedef real (*Fcn3dtime) (real x, real y, real z, real t);

/** AssignTimeFunctionOp assigns values to {\tt y} as a function of space and time. 
 *  The operator is unary. The constructor uses a function pointer:
 \begin{verbatim}
 typedef real (*Fcn3dtime) (real x, real y, real z, real t);
 \end{verbatim}
 The {\tt time} may be explicitly set with setTime, before a call to {\tt apply()}.
 The {\tt time} is also set by {applyTimeFunction(real *y,  real t)} and {\tt apply()} is called.
 */
class AssignTimeFunctionOp : public WavesSDOperator
{
	public:

		///
		AssignTimeFunctionOp(WavesSDGeometry *sdg, Fcn3dtime f);

		///
		AssignTimeFunctionOp(WavesSDIndexes *sdi, Fcn3dtime f);

		/// Destructor frees indexes if constructor 1 used
		~AssignTimeFunctionOp()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyTimeFunction(real *y, const real &t, const real &dt, const Fcn3dtime func)
		{
			time = t;
			step = dt;
			f = func;
			return doApply(0, y);
		}

		///
		void setTime(const real &t)
		{
			time = t;
		}
		///
		real getTime()
		{
			return time;
		}

	protected:
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		Fcn3dtime f;
		///
		real time;
		real step;
		///
		bool freeIndexes;
};
//=========================================================================
class AssignAdjTimeFunctionOp : public WavesSDOperator
{
	public:

		///
		AssignAdjTimeFunctionOp(WavesSDGeometry *sdg, real* f);

		///
		AssignAdjTimeFunctionOp(WavesSDIndexes *sdi, real* f);

		/// Destructor frees indexes if constructor 1 used
		~AssignAdjTimeFunctionOp()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyTimeFunction(real *y, const real &dt, real* func)
		{

			step = dt;
			f = func;
			return doApply(0, y);
		}

	protected:
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		real* f;
		///

		real step;
		///
		bool freeIndexes;
};

//======================================================================

///

class ComputeTimeDerivative : public WavesSDOperator
{
	public:
		///
		ComputeTimeDerivative(WavesSDGeometry *sdg);

		///
		ComputeTimeDerivative(WavesSDIndexes *sdi);

		/// Destructor frees indexes if constructor 1 used
		~ComputeTimeDerivative()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyTimeDerivative(real *y, real *x, real *z, real *v, const real &dt)
		{
			v_old_ = x;
			uOuter_ = z;
			vOuter_ = v;
			step = dt;
			return doApply(0, y);
		}

	protected:
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		real *v_old_;
		real *uOuter_;
		real *vOuter_;
		real step;
		///
		bool freeIndexes;
};

// ======================================================================

/**
 WavesAssignmentOp assigns a value to {\tt y}. The operator is unary.
 The value may be changed by the methods.
 */

class ApplyFunctionOp : public WavesSDOperator
{
	public:

		///
		ApplyFunctionOp(WavesSDGeometry *sdg);
		///
		ApplyFunctionOp(WavesSDIndexes *sdi);

		/// Destructor frees sdi if constructor 1 used
		~ApplyFunctionOp()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

	private:
		///
		real value;

		///
		bool freeIndexes;
};

//======================================================================

/**OutputOp just outputs y. The operator is unary. */
class OutputOp : public WavesSDOperator
{
	public:
		/// WavesSDIndexes covering all grid is created.
		OutputOp(WavesSDGeometry *sdg);
		///
		OutputOp(WavesSDIndexes *sdi);

		/// Destructor frees indexes if necessary
		~OutputOp()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

	private:
		///
		bool freeIndexes;

};

/** AVSOutputOp outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class AVSOutputOp : public WavesSDOperator
{
	public:

		///
		AVSOutputOp(WavesSDGeometry *sdg, ostream &ostr_);
		///

		AVSOutputOp(WavesSDIndexes *sdi, ostream &ostr_);

		~AVSOutputOp()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool interiorWithHole;
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		bool freeIndexes;
		///
		ostream &ostr;
};

/** AVSOutputOpOuter outputs the grid and the value of y for outer FDM
 * mesh which is added to the main FDM mesh.  It can be used for an
 * interior with hole or for a whole grid.  An ostream shall be
 * provided.  The operator is unary.
 */
class AVSOutputOpOuter : public WavesSDOperator
{
	public:

		AVSOutputOpOuter(WavesSDIndexes *sdi, ostream &ostr_);

		~AVSOutputOpOuter()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool interiorWithHole;
		///
		bool doApply(const real *, real *y) const;

	private:
		///
		bool freeIndexes;
		///
		ostream &ostr;
};

//=====================================================================

/** GIDOutputOp outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class GIDOutputOp
{
	public:

		///
		GIDOutputOp(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_);
		///

		GIDOutputOp(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_);

		GIDOutputOp(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int nnogg_);

		~GIDOutputOp()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printResults();
		bool printResultsCommon();
	protected:
		///
		bool interiorWithHole;
		char *file;
		real *u_array;
		real *u_1;
		real *u_2;
		int k;
		int nnogg;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};
//=====================================================================

/** GIDOutputOp outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class Out2D
{
	public:

		///
		Out2D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int n_);
		///

		Out2D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int n_);

		~Out2D()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printResults();

	protected:
		///
		bool interiorWithHole;
		char *file;
		real *u_array;
		real *u_1;
		real *u_2;
		int k;
		int n;
		///
		//  bool doApply(const real *, real *y) const ;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};

//=====================================================================

/** Out3D outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class Out3D
{
	public:

		///
		Out3D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_, int n_);
		///

		Out3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_, int n_);

		~Out3D()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printResults();

	protected:
		///
		bool interiorWithHole;
		char *file;
		real *u_array;
		real *u_1;
		real *u_2;
		real *u_3;
		int k;
		int n;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};
//=====================================================================

/** GIDOutputOp3D outputs the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class GIDOutputOp3D
{
	public:

		///
		GIDOutputOp3D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_);
		///

		GIDOutputOp3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_);

		GIDOutputOp3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_,  real* u_3_, real* u_array_, int k_, int nnogg_);

		~GIDOutputOp3D()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printResults();
		bool printResultsCommon();

	protected:
		///
		bool interiorWithHole;
		char *file;
		real *u_array;
		real *u_1;
		real *u_2;
		real *u_3;
		int k;
		int nnogg;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};

//=====================================================================

/** GIDOutputMesh outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class GIDOutputMesh
{
	public:

		///
		GIDOutputMesh(WavesSDGeometry *sdg, char *file_);
		GIDOutputMesh(WavesSDIndexes *sdi, char *file_);

		GIDOutputMesh(WavesSDIndexes *sdi, char *file_, int nnogg_, int nelgg_);

		~GIDOutputMesh()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printMesh();
		bool printCommonMesh();
	protected:
		///
		bool interiorWithHole;
		char *file;
		int nnogg, nelgg;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};
//======================================================================

/** GIDOutputCommon outputs the grid and the value of y.
 * It can be used for an interior with hole or for a whole grid. 
 * An ostream shall be provided. 
 * The operator is unary.
 */
class GIDOutputCommon
{
	public:

		GIDOutputCommon(WavesSDIndexes *sdi, char *file_, real* u_array_, real* u_1_, real* u_2_, int k_);

		~GIDOutputCommon()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printResults();
	protected:
		///

		WavesSDIndexes *indexes;
		bool unary;

		bool interiorWithHole;
		char *file;
		real *u_array;
		real *u_1;
		real *u_2;
		int k;

	private:
		///
		bool freeIndexes;

};

//======================================================================
class GIDOutputNodes
{
	public:
		GIDOutputNodes(WavesSDIndexes *sdi, char *file_, int nno_gg_);

		~GIDOutputNodes()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printMesh();

	protected:
		///
		bool interiorWithHole;
		char *file;
		int nno_gg;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};

//======================================================================
class GIDOutputWavesElements
{
	public:

		GIDOutputWavesElements(WavesSDIndexes *sdi, char *file_, int nel_, int nno_);

		~GIDOutputWavesElements()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool printMesh();

	protected:
		///
		bool interiorWithHole;
		char *file;
		int nel;
		int nno;
		WavesSDIndexes *indexes;
		bool unary;

	private:
		///
		bool freeIndexes;

};

//======================================================================

class PlotMTVOutputOp : public WavesSDOperator
{
	public:

		int writeMTV2D();
		///
		PlotMTVOutputOp(WavesSDGeometry *sdg, ostream &ostr_);
		///

		PlotMTVOutputOp(WavesSDIndexes *sdi, ostream &ostr_);

		~PlotMTVOutputOp()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool interiorWithHole;

		///
		bool doApply(const real *, real *y) const;

	private:
		///
		bool freeIndexes;
		///
		ostream &ostr;
};

//======================================================================

/**
 DifferenceCheckOp is a simple operator that compares x[IX] with y[IY].
 The indexes IX and IY are supplied to the constructor.  If only one index I is supplied, IX = IY = I.
 The call {\tt apply} do the check. If the difference is > tol, warning messages are written.
 */
class DifferenceCheckOp : public WavesSDOperator
{
	public:
		///
		DifferenceCheckOp(WavesSDIndexes *IY, WavesSDIndexes *IX = 0, real tol = 1E-6);
		///
		~DifferenceCheckOp();

	protected:
		///
		bool doApply(const real *x, real *y) const;

		///
		int *iyarray, *ixarray;

		/// oneIndex is true if only one index is supplied
		bool oneIndex;

		///
		real tol;
};

//======================================================================

#endif

