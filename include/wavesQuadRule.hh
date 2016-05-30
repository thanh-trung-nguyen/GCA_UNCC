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

#ifndef __WAVESQUADRULE_H
#define __WAVESQUADRULE_H

#include "waveskraftwerk.hh"

class WavesQuadRule
{

		virtual ~WavesQuadRule();

	public:

		virtual real getGPntCoord(int coord, int gpt) const = 0;
		virtual real operator()(int coord, int gpt) const = 0;
		virtual real weight(int gpt) const = 0;
		virtual int getNoPts() const = 0;

};

#endif
