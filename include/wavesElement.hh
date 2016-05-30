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

#ifndef __WAVESELEMENT_H
#define __WAVESELEMENT_H

#include "waveskraftwerk.hh"
#include "wavesBasisFcn.hh"
class WavesElement
{

	public:

		WavesElement(int nonode, int nobfcn, int nsd);
		virtual ~WavesElement()=0;

		WavesBasisFcn* ptr2zeroBfcn();
		WavesBasisFcn* ptr2test();
		WavesBasisFcn* ptr2trial();
		int getNoBFcns() const;
		WavesElement& setNoBFcns(int noBf);
		virtual real value(int bfcn, int der, real* coords) const = 0;
		void setIdx(int gpt, int testbf, int trialbf);
		void setIdx(int gpt, int bf);

		real getCoord(int node, int coord) const;
		int loc2glob(int locnode) const;
		int glob2loc(int globnode) const;
		int getNoNodes() const;
		int getElmNo() const;
		int getMaterialType() const;

		WavesElement& setNoGPnts(int nogp);
		int getNoGPnts();
		WavesElement& setGPntCoord(int gpt, int coord, real val);
		real& gPntCoord(int gpt, int coord);
		const real& gPntCoord(int gpt, int coord) const;
		WavesElement& setDetJac(int gpt, real val);
		real& detJac(int gpt);
		const real& detJac(int gpt) const;

		WavesElement& update(int elm, Grid& g);

	protected:

		int nsd;

	private:

		int noBFcns;
		WavesBasisFcn zeroBfcn;
		WavesBasisFcn trial;
		WavesBasisFcn test;

		int elmNo;
		int noNodes;
		int materialType;
		MV_ColMat<real> coords;
		Array1dInt globNodeNo;

		int noGPnts;
		MV_ColMat<real> coordGPnts;
		Array1dReal detJacs;

};

#endif
