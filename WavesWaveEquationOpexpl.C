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
// First added:  2000-01-18 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#include <math.h>
#include "include/wavesWaveEquationOpexpl.h"

WavesWaveEqInteriorexpl1::WavesWaveEqInteriorexpl1(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg)), freeIndexes(true)
{
	doInitialize();
}

WavesWaveEqInteriorexpl1::WavesWaveEqInteriorexpl1(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	doInitialize();
}

void WavesWaveEqInteriorexpl1::doInitialize()
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	// 1) create Laplace operator
	const int nsd = geometry->getNoSpaceDim();

	nrOfWeights = (nsd == 2) ? 5 : 7;

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx2 = geometry->getDx() * geometry->getDx();
	real dy2 = geometry->getDy() * geometry->getDy();
	real dz2 = geometry->getDz() * geometry->getDz();

	if (nsd == 3)
	{
		weights[0] = -2.0 / dx2 - 2.0 / dy2 - 2.0 / dz2;
		weights[1] = weights[2] = 1.0 / dx2;
		weights[3] = weights[4] = 1.0 / dy2;
		weights[5] = weights[6] = 1.0 / dz2;

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

		offsets[0] = 0;
		offsets[1] = geometry->node(-1, 0, 0);
		offsets[2] = geometry->node(1, 0, 0);
		offsets[3] = geometry->node(0, -1, 0);
		offsets[4] = geometry->node(0, 1, 0);
	}

}
// we compute wave equation as system from two equations
// apply      y1 = tau*Laplace x2 + x1
//            y2 = tau*x1 +x2;
// *x1,*x2 - pointers to the solution (x1,x2) on  the previous iteration

bool WavesWaveEqInteriorexpl1::doApply(const double *, double *y1) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	const int nsd = geometry->getNoSpaceDim();
	assert( nsd>=2 && nsd <= 3);

	// We apply y1 <- dt* Laplace x2 + x1

	// For faster loops:
	real weights0 = weights[0];
	real weights1 = weights[1];
	real weights2 = weights[2];
	real weights3 = weights[3];
	real weights4 = weights[4];
	real weights5 = (nsd == 3) ? weights[5] : 0.;
	real weights6 = (nsd == 3) ? weights[6] : 0.;

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
			for (int n = lstart; n <= lstop; n++)
				y1[n] = x1_[n] + dt_ * (weights0 * x2_[n] + weights1 * x2_[n + offsets1] + weights2 * x2_[n + offsets2] + weights3 * x2_[n + offsets3] + weights4 * x2_[n + offsets4] + weights5 * x2_[n + offsets5] + weights6 * x2_[n + offsets6]);

		else
			for (int n = lstart; n <= lstop; n++)
			{
				y1[n] = x1_[n] + dt_ * (weights0 * x2_[n] + weights1 * x2_[n + offsets1] + weights2 * x2_[n + offsets2] + weights3 * x2_[n + offsets3] + weights4 * x2_[n + offsets4]);
			}
	}
	return true;
}

//
// class WavesWaveEqInteriorexpl2
//


WavesWaveEqInteriorexpl2::WavesWaveEqInteriorexpl2(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg)), freeIndexes(true)
{
	unary = true;
}

WavesWaveEqInteriorexpl2::WavesWaveEqInteriorexpl2(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
}

// we compute wave equation as system from two equations
// and apply      y1 = tau*Laplace x2 + x1
//            y2 = tau*x1 +x2;
// *x1,*x2 - pointers to the solution (x1,x2) on  the previous iteration

bool WavesWaveEqInteriorexpl2::doApply(const double *, double *y2) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();
	const int nsd = geometry->getNoSpaceDim();
	assert( nsd>=2 && nsd <= 3);

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		for (int n = lstart; n <= lstop; n++)
		{
			y2_[n] = x1_[n] * dt_ + x2_[n];
		}
	}
	return true;
}



WavesWaveEqBoundexpl::WavesWaveEqBoundexpl(WavesSDGeometry *sdg, const SDBoundaryType &sdb_, double dt_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), sdb(sdb_), dt(dt_), freeIndexes(true)
{
	cout << "//### WavesWaveEqBoundexpl::WavesWaveEqBoundexpl sdb" << int(sdb) << endl;
	doInitialize();
}

WavesWaveEqBoundexpl::WavesWaveEqBoundexpl(WavesSDBoundary *sdi, double dt_) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), dt(dt_), freeIndexes(false)
{
	cout << "//### WavesWaveEqBoundexpl::WavesWaveEqBoundexpl sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesWaveEqBoundexpl::doInitialize()
{
	switch (sdb)
	{
		case SDallBoundary:
			cout << "SDallBoundary not currently supported!\n" << endl;
			break;
		case SDiLow:
		case SDiHigh:
		case SDjLow:
		case SDjHigh:
		case SDkLow:
		case SDkHigh:
			initWeights();
			break;
		default:
			cout << "WavesWaveEqBoundexpl: ERROR, unknown or not explemented SDBoundaryType " << endl;
			break;
	}
}

void WavesWaveEqBoundexpl::initWeights()
{
	// is called when sdb is known to be one side
	nrOfWeights = 3;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx = geometry->getDx();
	real dy = geometry->getDy();
	real dz = geometry->getDz();

	switch (sdb)
	{
		case SDiLow:
			offsets[0] = 0;
			offsets[1] = geometry->node(1, 0, 0);
			offsets[2] = offsets[1];

			weights[0] = (dx - dt) / (dt + dx);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		case SDjLow:
			offsets[0] = 0;
			offsets[1] = geometry->node(0, 1, 0);
			offsets[2] = offsets[1];

			weights[0] = (dy - dt) / (dt + dy);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		case SDkLow:
			offsets[0] = 0;
			offsets[1] = geometry->node(0, 0, 1);
			offsets[2] = offsets[1];

			weights[0] = (dz - dt) / (dt + dz);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		case SDiHigh:
			offsets[0] = 0;
			offsets[1] = geometry->node(-1, 0, 0);
			offsets[2] = offsets[1];

			weights[0] = (dx - dt) / (dt + dx);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		case SDjHigh:
			offsets[0] = 0;
			offsets[1] = geometry->node(0, -1, 0);
			offsets[2] = offsets[1];

			weights[0] = (dy - dt) / (dt + dy);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		case SDkHigh:
			offsets[0] = 0;
			offsets[1] = geometry->node(0, 0, -1);
			offsets[2] = offsets[1];

			weights[0] = (dz - dt) / (dt + dz);
			weights[1] = 1.;
			weights[2] = -weights[0];
			break;

		default:
			cout << "WavesWaveEqBoundexpl: ERROR, unknown or not explemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesWaveEqBoundexpl::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		// Based on x[0], x[1], and y[1], update y[0].
		for (int n = lstart; n <= lstop; n++)
			y[n] = weights[0] * x[n] + weights[1] * x[n + offsets[1]] + weights[2] * y[n + offsets[2]];
	}

	return true;
}

