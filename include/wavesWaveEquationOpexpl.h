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
// First added:  2000-01-08 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESWAVEEQUATIONOPEXPL_H
#define __WAVESWAVEEQUATIONOPEXPL_H

#include "wavesSDOperator.h"

class WavesWaveEqInteriorexpl1 : public WavesSDOperator
{
	public:

		/// SDindexes covering the interior is created. The time step dt shall be supplied.
		WavesWaveEqInteriorexpl1(WavesSDGeometry *sdg);
		/// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
		WavesWaveEqInteriorexpl1(WavesSDIndexes *sdi);

		/// Destructor frees memory if needed.
		~WavesWaveEqInteriorexpl1()
		{
			if (freeIndexes)
				delete indexes;
		}
		bool applyWaveEqIntexpl1(double *y1, double &dt, double *x1, double *x2)
		{
			dt_ = dt;
			y1_ = y1;
			x1_ = x1;
			x2_ = x2;
			return doApply(0, y1);
		}

	protected:
		///
		bool doApply(const double *, double *y1) const;

		/// Utility
		void doInitialize();

		///
		double dt_;
		double *y1_;
		double *x1_;
		double *x2_;
		//
		bool freeIndexes;
};

class WavesWaveEqInteriorexpl2 : public WavesSDOperator
{
	public:

		/// SDindexes covering the interior is created. The time step dt shall be supplied.
		WavesWaveEqInteriorexpl2(WavesSDGeometry *sdg);
		/// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
		WavesWaveEqInteriorexpl2(WavesSDIndexes *sdi);

		/// Destructor frees memory if needed.
		~WavesWaveEqInteriorexpl2()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool applyWaveEqIntexpl2(double *y2, double &dt, double *x1, double *x2)
		{

			y2_ = y2;
			dt_ = dt;
			x1_ = x1;
			x2_ = x2;

			return doApply(0, y2);
		}

	protected:
		///
		bool doApply(const double *, double *y) const;

		///

		double dt_;
		double *y2_;
		double *x1_;
		double *x2_;
		//
		bool freeIndexes;
};

/**
 WavesWaveEqBoundexpl represents absorbing boundary conditions of Engquist-Majda type.
 */

class WavesWaveEqBoundexpl : public WavesSDOperator
{
	public:

		/// The {\tt SDBoundaryType sdb} must represent a face.
		WavesWaveEqBoundexpl(WavesSDGeometry *sdg, const SDBoundaryType &sdb, double dt_);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesWaveEqBoundexpl(WavesSDBoundary *sdb, double dt_);

		/// Destructor frees memory if needed.
		~WavesWaveEqBoundexpl()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const double *x, double *y) const;

		/// Utility
		void doInitialize();
		/// Utility
		void initWeights();
		///
		SDBoundaryType sdb;
		///
		double dt;
		///
		bool freeIndexes;
};

#endif

