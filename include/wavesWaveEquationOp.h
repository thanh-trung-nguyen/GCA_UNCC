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

#ifndef __WAVESWAVEEQUATIONOP_H
#define __WAVESWAVEEQUATIONOP_H

#include "wavesSDOperator.h"

/* 
 Here, the operators needed for application 2 is collected.

 */

// In the interior, y_i := 2dt laplace x_i - y_i
/**
 WavesWaveEqInterior is an operator for the scalar wave equation.

 The scalar wave equation is discretized with centred differences. The {\tt SDField y} represnting the unknown on time level {\tt n-1} is overwritten to represent values on time level {\tt n+1}, using {\tt SDField x} representing the unknown on time level {\tt n}, according to 
 \begin{verbatim}
 y_i := 2dt laplace x_i - y_i
 \end{verbatim}
 */

class WavesWaveEqInterior : public WavesSDOperator
{
	public:
		/// SDindexes covering the interior is created. The time step dt shall be supplied.
		WavesWaveEqInterior(WavesSDGeometry *sdg, double dt_);

		/// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
		WavesWaveEqInterior(WavesSDIndexes *sdi, double dt_);

		/// Destructor frees memory if needed.
		~WavesWaveEqInterior()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool doApply(const double *x, double *y) const;

		/// Utility
		void doInitialize();

		///
		double dt;

		//
		bool freeIndexes;
};

/**
 WavesWaveEqBoundary represents absorbing boundary conditions of Engquist-Majda type.
 */

class WavesWaveEqBoundary : public WavesSDOperator
{
	public:
		/// The {\tt SDBoundaryType sdb} must represent a face.
		WavesWaveEqBoundary(WavesSDGeometry *sdg, const SDBoundaryType &sdb, double dt_);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesWaveEqBoundary(WavesSDBoundary *sdb, double dt_);

		/// Destructor frees memory if needed.
		~WavesWaveEqBoundary()
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

//=========================================================================
/**
 WavesWaveEqAdjBoundary represents adjoint boundary conditions for Adjoint equation
 in inverse problem
 */

class WavesWaveEqAdjBoundary : public WavesSDOperator
{
	public:
		/// The {\tt SDBoundaryType sdb} must represent a face.
		WavesWaveEqAdjBoundary(WavesSDGeometry *sdg, const SDBoundaryType &sdb);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesWaveEqAdjBoundary(WavesSDBoundary *sdb);

		/// Destructor frees memory if needed.
		~WavesWaveEqAdjBoundary()
		{
			if (freeIndexes)
				delete indexes;

		}

	protected:

		///
		bool doApply(const double *f, double *y) const;

		/// Utility
		void doInitialize();
		/// Utility
		void initWeights();
		///
		SDBoundaryType sdb;
		///
		bool freeIndexes;
};

//============================================================================
class WavesWaveEqInterior2 : public WavesSDOperator
{
	public:
		/// SDindexes covering the interior is created. The time step dt shall be supplied. It uses nine point stencil.
		WavesWaveEqInterior2(WavesSDGeometry *sdg, double dt_);
		/// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
		WavesWaveEqInterior2(WavesSDIndexes *sdi, double dt_);

		/// Destructor frees memory if needed.
		~WavesWaveEqInterior2()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:
		///
		bool doApply(const double *x, double *y) const;

		/// Utility
		void doInitialize();

		///
		double dt;

		//
		bool freeIndexes;
};

#endif

