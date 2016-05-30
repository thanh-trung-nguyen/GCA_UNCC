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

#ifndef __WAVESBASISFCN_H
#define __WAVESBASISFCN_H

#include "waveskraftwerk.hh"
class WavesMapper;
class WavesElement;

class WavesBasisFcn
{

	public:

		WavesBasisFcn()
		{
		}
		WavesBasisFcn(int nsd_);
		WavesBasisFcn(int gptsize, int bfsize, int nsd_) :
				nsd(nsd_), values(gptsize, bfsize), grad(gptsize, bfsize, nsd_)
		{
		}

		~WavesBasisFcn()
		{
		}

		WavesBasisFcn& newsize(int gptsize, int bfsize, int nsd_);

		inline
		void refill(Array2d<real>& v, Array3d<real>& g)
		{
			values = v;
			grad = g;
		}

		inline
		void set(int gptset, int bfset)
		{
			gpt = gptset;
			bf = bfset;
		}

		inline
		double operator()()
		{
			return values(gpt, bf);
		}

		inline
		double x(const int der)
		{ // chose which derivative with argument i
			return grad(gpt, bf, der);
		}

		inline
		double x()
		{
			return grad(gpt, bf, 0);
		}

		inline
		double y()
		{
			return grad(gpt, bf, 1);
		}

		inline
		double z()
		{
			return grad(gpt, bf, 2);
		}

	private:

		int nsd;
		Array2d<real> values;
		Array3d<real> grad;
		int gpt;
		int bf;

		friend class WavesMapper;

};

#endif
