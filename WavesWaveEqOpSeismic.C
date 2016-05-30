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
// First added:  2006-01-07 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#include <math.h>
#include <assert.h>
#include "include/wavesWaveEqOpSeismic.h"

using namespace std;

WavesWaveEqInteriorSeismic::WavesWaveEqInteriorSeismic(WavesSDGeometry *sdg, double *b_, double dt_) :
		WavesSDOperator(new WavesSDInterior(*sdg)), b(b_), dt(dt_), freeIndexes(true)
{
	doInitialize();
}

WavesWaveEqInteriorSeismic::WavesWaveEqInteriorSeismic(WavesSDIndexes *sdi, double *b_, double dt_) :
		WavesSDOperator(sdi), b(b_), dt(dt_), freeIndexes(false)
{
	doInitialize();
}

void WavesWaveEqInteriorSeismic::doInitialize()
{
//	cout << "inside doInitialize" << endl;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	// 1) create Laplace operator
	const int nsd = geometry->getNoSpaceDim();

	assert( nsd>=2 && nsd <= 3);

	nrOfWeights = (nsd == 2) ? 7 : 10;

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx2 = geometry->getDx() * geometry->getDx();
	real dy2 = geometry->getDy() * geometry->getDy();
	real dz2 = geometry->getDz() * geometry->getDz();

	real dx = geometry->getDx();
	real dy = geometry->getDy();
	real dz = geometry->getDz();

	if (nsd == 3)
	{
		weights[0] = -2.0 / dx2 - 2.0 / dy2 - 2.0 / dz2;
		weights[1] = weights[2] = 1.0 / dx2;
		weights[3] = weights[4] = 1.0 / dy2;
		weights[5] = weights[6] = 1.0 / dz2;
		weights[7] = 1.0 / dx;
		weights[8] = 1.0 / dy;
		weights[9] = 1.0 / dz;

		offsets[0] = 0;
		offsets[1] = geometry->node(-1, 0, 0);
		offsets[2] = geometry->node(1, 0, 0);
		offsets[3] = geometry->node(0, -1, 0);
		offsets[4] = geometry->node(0, 1, 0);
		offsets[5] = geometry->node(0, 0, -1);
		offsets[6] = geometry->node(0, 0, 1);
	}
	else
	{
		weights[0] = -2.0 / dx2 - 2.0 / dy2;
		weights[1] = weights[2] = 1.0 / dx2;
		weights[3] = weights[4] = 1.0 / dy2;
		weights[5] = 1.0 / dx;
		weights[6] = 1.0 / dy;

		offsets[0] = 0;
		offsets[1] = geometry->node(-1, 0, 0);
		offsets[2] = geometry->node(1, 0, 0);
		offsets[3] = geometry->node(0, -1, 0);
		offsets[4] = geometry->node(0, 1, 0);
	}

	// 2) y := dt*dt Laplace + 2
	for (int i = 0; i < nrOfWeights; i++)
		weights[i] *= (dt * dt);
	weights[0] += 2;
}

bool WavesWaveEqInteriorSeismic::doApply(const double *x, double *y) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();
	const int nsd = geometry->getNoSpaceDim();
	assert( nsd>=2 && nsd <= 3);

	// We apply y <- (dt^2 Laplace + 2) x - y

	// For faster loops:
	real weights0 = weights[0];
	real weights1 = weights[1];
	real weights2 = weights[2];
	real weights3 = weights[3];
	real weights4 = weights[4];
	real weights5 = weights[5];
	real weights6 = weights[6];

	int offsets1 = offsets[1];
	int offsets2 = offsets[2];
	int offsets3 = offsets[3];
	int offsets4 = offsets[4];
	int offsets5 = (nsd == 3) ? offsets[5] : 0;
	int offsets6 = (nsd == 3) ? offsets[6] : 0;

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		if (nsd == 3)
		{
			for (int n = lstart; n <= lstop; n++)
			{
				y[n] = -y[n] + weights0 * x[n] + weights1 * x[n + offsets1] + weights2 * x[n + offsets2] + weights3 * x[n + offsets3] + weights4 * x[n + offsets4] + weights5 * x[n + offsets5] + weights6 * x[n + offsets6];
			}
		}
		else
		{
			for (int n = lstart; n <= lstop; n++)
			{

				if (geometry->getCoor(n, 1) > 0.46 && geometry->getCoor(n, 1) < 0.504 && geometry->getCoor(n, 0) > 0.42 && geometry->getCoor(n, 0) < 0.508)
				{
					y[n] = -y[n] + weights0 * x[n] + weights1 * x[n + offsets1] + weights2 * x[n + offsets2] + weights3 * x[n + offsets3] + weights4 * x[n + offsets4] - weights1 * (x[n + offsets2] - x[n]) - weights3 * (x[n + offsets4] - x[n]);

				}
				else
				{

					y[n] = -y[n] + weights0 * x[n] + weights1 * x[n + offsets1] + weights2 * x[n + offsets2] + weights3 * x[n + offsets3] + weights4 * x[n + offsets4] - weights1 * (x[n + offsets2] - x[n]) * (b[n + offsets2] - b[n])
							- weights3 * (x[n + offsets4] - x[n]) * (b[n + offsets4] - b[n]);
				}

			}
		}

	}

	return true;
}
