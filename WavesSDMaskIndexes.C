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

#include "include/wavesSDMaskIndexes.h"
#include <iostream>
#include <assert.h>

using namespace std;

WavesSDMaskIndexes::WavesSDMaskIndexes(const WavesSDGeometry &sdg1, int *mask1, const SDIndexType &sdm, int w) :
		WavesSDIndexes(sdg1)
{
	mask = mask1;
	sdg = sdg1;

	assert(w >= 0);

	const int nsd = sdg1.getNoSpaceDim();
	const int n_i = sdg1.getN_i();
	const int n_j = sdg1.getN_j();
	const int n_k = sdg1.getN_k();
	if (nsd == 3)
		nOfLoops = n_j * n_k;
	else if (nsd == 2)
		nOfLoops = n_j;
	else
	{
		cout << "Init of WavesSDIndexes for " << nsd << "space dimensions not implemented." << endl;
		nOfLoops = 0;
		return;
	}

	assert(nOfLoops > 0);
	loopIndexes = new LoopIndex[nOfLoops];
	switch (sdm)
	{
		case code3:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition3, sdg1);
			presentWavesSDIndexes();
			break;
		}

		case code4:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition4, sdg1);
			break;
		}

		case code5:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition5, sdg1);
			break;
		}
		case code6:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition6, sdg1);
			presentWavesSDIndexes();
			break;
		}
		case codeb1:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition7, sdg1);
			presentWavesSDIndexes();
			break;
		}
		case codeb2:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition8, sdg1);
			presentWavesSDIndexes();
			break;
		}
		case code_inter:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition9, sdg1);
			break;
		}

		case code10:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition10, sdg1);
			break;
		}
		case code11:
		{
			MaskAddLoop(&WavesSDMaskIndexes::condition11, sdg1);
			//cout << "outerWithHole sdindexes" << endl;
			presentWavesSDIndexes();
			break;
		}
		case code12:
		{
			MaskAddLoopInner(&WavesSDMaskIndexes::condition6, sdg1);
			//cout << "outerWithHole sdindexes" << endl;
			presentWavesSDIndexes();
			break;
		}

		default:
			cout << " Error in SDIndexMask:unknown options." << endl;
			break;

	}

}

// method: add loop to loopIndexes for mask type 

int WavesSDMaskIndexes::MaskAddLoop(cond condition, const WavesSDGeometry &sdg1)
{
	int mask_loop = 0;
	int counter = 0;
	curLoop = 0;
	const int nsd = sdg1.getNoSpaceDim();
	const int n_i = sdg1.getN_i();
	const int n_j = sdg1.getN_j();
	const int n_k = sdg1.getN_k();
	int i, j, k;

	for (k = 0; k < n_k; k++)
	{
		for (j = 0; j < n_j; j++)
		{
			for (i = 0; i < n_i; i++)
			{
				// condition has an implicit parameter pointing to SDIndexMask

				if ((this->*condition)(i, j, k))
				{

					if (counter == 0)
						mask_loop = i;
					counter++;
				}
				else if (counter > 0)
				{
					doAddLoop(mask_loop, i - 1, j, k);
					counter = 0;
				} // for mask = 0
			} // for i

			//  if we meet only mask-elements 

			if (counter > 0)
			{
				doAddLoop(mask_loop, i - 1, j, k);
				counter = 0;
			}
		} // for j
	} // for k

	//  schrinkLoopIndexArray();
	resize(curLoop);
	return curLoop;

}

//======================================================================

int WavesSDMaskIndexes::MaskAddLoopInner(cond condition, const WavesSDGeometry &sdg1)
{
	int mask_loop = 0;
	int counter = 0;
	curLoop = 0;
	const int nsd = sdg1.getNoSpaceDim();
	const int n_i = sdg1.getN_i();
	const int n_j = sdg1.getN_j();
	const int n_k = sdg1.getN_k();
	int i, j, k;

	if (nsd == 3)
	{
		for (k = 1; k < n_k - 1; k++)
		{
			for (j = 1; j < n_j - 1; j++)
			{
				for (i = 1; i < n_i - 1; i++)
				{
					// condition has an implicit parameter pointing to SDIndexMask

					if ((this->*condition)(i, j, k))
					{

						if (counter == 0)
							mask_loop = i;
						counter++;
					}
					else if (counter > 0)
					{
						doAddLoop(mask_loop, i - 1, j, k);
						counter = 0;
					} // for mask = 0
				} // for i

				//  if we meet only mask-elements 

				if (counter > 0)
				{
					doAddLoop(mask_loop, i - 1, j, k);
					counter = 0;
				}
			} // for j
		} // for k
	} // for nsd == 3
	else if (nsd == 2)
	{

		for (j = 1; j < n_j - 1; j++)
		{
			for (i = 1; i < n_i - 1; i++)
			{
				// condition has an implicit parameter pointing to SDIndexMask

				if ((this->*condition)(i, j, 0))
				{

					if (counter == 0)
						mask_loop = i;
					counter++;
				}
				else if (counter > 0)
				{
					doAddLoop(mask_loop, i - 1, j, 0);
					counter = 0;
				} // for mask = 0
			} // for i

			//  if we meet only mask-elements 

			if (counter > 0)
			{
				doAddLoop(mask_loop, i - 1, j, 0);
				counter = 0;
			}
		} // for j

	}

	//  schrinkLoopIndexArray();
	resize(curLoop);
	return curLoop;

}

//======================================================================

WavesSDMaskIndexes::~WavesSDMaskIndexes()
{
	/*  delete [] loopIndexes;    -- this happens in ~WavesSDIndexes() */
}

