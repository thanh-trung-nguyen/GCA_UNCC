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
// First added:  2000-02-23 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESPLANEWAVEOPEXPL_H
#define __WAVESPLANEWAVEOPEXPL_H

#include "wavesSDDefs.h"
#include "wavesFuncDefs.h"
#include "wavesGridB.h"
#include "wavesSDIndexes.h"
#include "wavesSDGeometry.h"
#include "wavesSDOperator.h"

/** PlaneWaveOpexpl  assigns values to {\tt y} as a function of space and time. 
 *  The operator is unary. The constructor uses a function pointer:
 \begin{verbatim}
 typedef real (*Fcn3dtime) (real x, real y, real z, real t);
 \end{verbatim}
 The {\tt time} may be explicitly set with setTime, before a call to {\tt apply()}.
 The {\tt time} is also set by {applyPlaneWaveOpexpl(real *y,  real t)} and {\tt apply()} is called.
 */
class PlaneWaveOpexpl : public WavesSDOperator
{
	public:
		PlaneWaveOpexpl(WavesSDGeometry *sdg, Fcn3dtime f);

		PlaneWaveOpexpl(WavesSDIndexes *sdi, Fcn3dtime f);

		/// Destructor frees indexes if constructor 1 used
		~PlaneWaveOpexpl()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyPlaneWaveOpexpl(real *y, const real &t, const Fcn3dtime func)
		{
			time = t;
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
//  real step;
		real time;
		///
		bool freeIndexes;
};

//===============================================

class PlaneWaveOpexpl2 : public WavesSDOperator
{
	public:
		///
		PlaneWaveOpexpl2(WavesSDGeometry *sdg, Fcn3dtime f);

		///
		PlaneWaveOpexpl2(WavesSDIndexes *sdi, Fcn3dtime f);

		/// Destructor frees indexes if constructor 1 used
		~PlaneWaveOpexpl2()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyPlaneWaveOpexpl(real *y, const real &t, const real &lambda, const real &mu, const real &ro, const Fcn3dtime func)
		{
			time = t;
			lambda_ = lambda;
			mu_ = mu;
			ro_ = ro;

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
//  real step;
		real time;
		real lambda_;
		real mu_;
		real ro_;
		///
		bool freeIndexes;
};

//==================================================0

class PlaneWaveOpexpl3 : public WavesSDOperator
{
	public:
		///
		PlaneWaveOpexpl3(WavesSDGeometry *sdg, real *Difsol_);

		///
		PlaneWaveOpexpl3(WavesSDIndexes *sdi, real *Difsol_);

		/// Destructor frees indexes if constructor 1 used
		~PlaneWaveOpexpl3()
		{
			if (freeIndexes)
				delete indexes;
		}

		///
		bool applyPlaneWaveOpexpl(real *y, const real &t)
		{
			time = t;

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
		real *Difsol;
		///
		real time;

		///
		bool freeIndexes;
};

//================================================================
class PlaneWaveOpPhotonexpl : public WavesSDOperator
{
	public:
		///
		PlaneWaveOpPhotonexpl(WavesSDGeometry *sdg, Fcn3dtime f);

		///
		PlaneWaveOpPhotonexpl(WavesSDIndexes *sdi, Fcn3dtime f);

		/// Destructor frees indexes if constructor 1 used
		~PlaneWaveOpPhotonexpl()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool applyPlaneWaveOpexpl(real *y, const real &t, const real &omega, const Fcn3dtime func)
		{
			time = t;
			omega_ = omega;

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
		real omega_;

		///
		bool freeIndexes;
};

#endif

