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
// First added:  2002-01-22 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESREFLECTEDWAVE_H
#define __WAVESREFLECTEDWAVE_H

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "wavesSDDefs.h"
#include "wavesFuncDefs.h"
#include "wavesGridB.h"
#include "wavesSDIndexes.h"
#include "wavesSDGeometry.h"
#include "wavesSDOperator.h"

typedef MV_ColMat<real> Mat_real;

class WavesReflectedWave : public WavesSDOperator
{
	public:

		WavesReflectedWave(WavesSDIndexes *sdi, char* fname_, double* array_);

		WavesReflectedWave(WavesSDIndexes *sdi, char* fname_);
		WavesReflectedWave(WavesSDIndexes *sdi, MV_Vector<double>& array_);

		/// Destructor frees indexes if constructor 1 used
		~WavesReflectedWave()
		{
			if (freeIndexes)
				delete indexes;
		}

		bool applyReflectedWave(real *y, const real &t);
		bool applyReflectedWaveLow(real *y);
		bool applyReflectedWaveLeft(real *y, const real &t);

		void write_reflected_field(double* array);
		void read_reflected_field(double& t);
		void read_reflected_(double& t);

	protected:
		///
		bool doApply(const real *, real *y) const;

	private:
		double* output_ar;
		MV_Vector<double> Data_array;

		char* fname;
		real time;

		///
		bool freeIndexes;
};

#endif
