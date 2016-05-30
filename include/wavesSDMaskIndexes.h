// Copyright (C) 1999 Larisa Beilina
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
// First added:  1999-10-10 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESSDMASKINDEXES_H
#define __WAVESSDMASKINDEXES_H

#include "wavesSDIndexes.h"
#include <iostream>
#include <assert.h>

enum SDIndexType
{
		code3,
		code4,
		code5,
		code6,
		codeb1,
		codeb2,
		code_inter,
		code10,
		code11,
		code12
};

/*
 WavesSDMaskIndexes produces LoopIndexes depending on the mask.
 */

class WavesSDMaskIndexes : public WavesSDIndexes
{
	public:
		int *mask;
		//cond - adres of the function with type bool

		typedef bool(WavesSDMaskIndexes::*cond)(int i, int j, int k);

		~WavesSDMaskIndexes();

		WavesSDMaskIndexes(const WavesSDGeometry &sdg, int *mask, const SDIndexType &sdm, int w = 1);
		int MaskAddLoop(cond condition, const WavesSDGeometry &sdg);

		// added to make inner loops :can be also rewritten with  w=1 in
		int MaskAddLoopInner(cond condition, const WavesSDGeometry &sdg);

		int size_loop_index(const WavesSDGeometry &sdg, int *mask, const SDIndexType &sdm, int w = 1);

		LoopIndex *mask_loop_index(const WavesSDGeometry &sdg, int *mask, const SDIndexType &sdm, int w = 1);

		bool condition3(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 3;
		}

		bool condition4(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 0 || mask[sdg.node(i, j, k)] == 1 || mask[sdg.node(i, j, k)] == 2 || mask[sdg.node(i, j, k)] == 3;
		}

		bool condition5(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 0 || mask[sdg.node(i, j, k)] == 1 || mask[sdg.node(i, j, k)] == 2;
		}

		bool condition6(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 2 || mask[sdg.node(i, j, k)] == 3 || mask[sdg.node(i, j, k)] == 1;
		}

		bool condition7(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 1;
		}

		bool condition8(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 2;
		}

		bool condition9(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 0;
		}

		bool condition10(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 22;
		}

		bool condition11(int i, int j, int k)
		{
			return mask[sdg.node(i, j, k)] == 1 || mask[sdg.node(i, j, k)] == 22 || mask[sdg.node(i, j, k)] == 2 || mask[sdg.node(i, j, k)] == 3;
		}

};

#endif
