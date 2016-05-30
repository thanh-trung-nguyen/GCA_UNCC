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

#ifndef __WAVESWAVEEQOPSEISMIC_H
#define __WAVESWAVEEQOPSEISMIC_H

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

class WavesWaveEqInteriorSeismic : public WavesSDOperator
{
	public:
		/// SDindexes covering the interior is created. The time step dt shall be supplied.
		WavesWaveEqInteriorSeismic(WavesSDGeometry *sdg, double *b_, double dt_);

		/// The supplied WavesSDIndexes shall take the operator width into account. The time step dt shall be supplied.
		WavesWaveEqInteriorSeismic(WavesSDIndexes *sdi, double *b_, double dt_);

		/// Destructor frees memory if needed.
		~WavesWaveEqInteriorSeismic()
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

		double *b;
		//
		bool freeIndexes;
};

#endif

