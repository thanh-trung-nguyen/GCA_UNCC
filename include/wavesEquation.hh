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

#ifndef __WAVESEQUATION_H
#define __WAVESEQUATION_H

#include "waveskraftwerk.hh"
#include "wavesBasisFcn.hh"

class WavesEquation
{

	public:

		WavesEquation(int noeq, int nsd_);
		virtual ~WavesEquation()
		{
		};

		void integrateJacobian(WavesQuadRule& q, // Input:  Quadrature to use 
				WavesElement& e, //         Which element to integrate over
				Array1dReal& result, // Output: Where to put the element matrix
				Array1dInt& rowIdx, //         Row index in the stiffness  matrix
				Array1dInt& colIdx, //         Column index in the stiffness  matrix
				int& noElmsInElMat //         Total size of element matrix
				);

		void integrateResidual(WavesQuadRule& q, // Input:  Quadrature to use  
				WavesElement& e, //         Which element to integrate over
				Array1dReal& result, // Output: Where to put the element matrix
				Array1dInt& rowIdx, //         Row index in the stiffness  matrix 
				int& noElmsInElMat //         Total size of element matrix
				);

		int noEq() const
		{
			return noEqns;
		}

		virtual
		void jacobian() = 0; // PDE definition

		virtual
		void residual() = 0; // PDE definition

	protected:

		void setWavesBasisFcns(WavesElement& e, int eq); // Set WavesEquation's basisfcns to point to WavesElement's

		void setWavesBasisFcns(WavesElement& e); // Set WavesEquation's basisfcns to point to WavesElement's

		int nsd;
		int noEqns;
		MV_Vector<WavesBasisFcn*> testFcns;
		MV_Vector<WavesBasisFcn*> trialFcns;
		Array1dReal work;
		Array1dReal coord;
		WavesElement* myElm;

	private:

		Array4d<real> workElMat;
		Array2d<real> workElVec;

};

#endif
