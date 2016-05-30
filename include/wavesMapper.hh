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

#ifndef __WAVESMAPPER_H
#define __WAVESMAPPER_H

#include "waveskraftwerk.hh"

class WavesMapper
{

	public:

		WavesMapper(int nsd_);
		~WavesMapper()
		{
		};
		void mapIsoparametric(WavesQuadRule& q, WavesElement& e);

	private:

		int nsd;

		Array2d<real> tmpValGPts;
		Array2d<real> tmpCoords;
		Array2d<real> tstval;
		Array2d<real> trialval;
		Array2d<real> jac;
		Array2d<real> cramer;

		Array3d<real> tstgrad;
		Array3d<real> trialgrad;

		Array1dReal det;
		Array1dReal wts;

};

#endif
