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

#include <math.h>
#include <assert.h>
#include "include/wavesWaveEquationOp.h"

using namespace std;

WavesWaveEqInterior::WavesWaveEqInterior(WavesSDGeometry *sdg, double dt_) :
		WavesSDOperator(new WavesSDInterior(*sdg)), dt(dt_), freeIndexes(true)
{
	doInitialize();
}

WavesWaveEqInterior::WavesWaveEqInterior(WavesSDIndexes *sdi, double dt_) :
		WavesSDOperator(sdi), dt(dt_), freeIndexes(false)
{
	doInitialize();
}

void WavesWaveEqInterior::doInitialize()
{
//	cout << "inside doInitialize" << endl;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	// 1) create Laplace operator
	const int nsd = geometry->getNoSpaceDim();

	assert( nsd>=2 && nsd <= 3);

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

	// 2) y := dt*dt Laplace + 2
	for (int i = 0; i < nrOfWeights; i++)
		weights[i] *= (dt * dt);
	weights[0] += 2;
}

bool WavesWaveEqInterior::doApply(const double *x, double *y) const
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
			{
				y[n] = -y[n] + weights0 * x[n] + weights1 * x[n + offsets1] + weights2 * x[n + offsets2] + weights3 * x[n + offsets3] + weights4 * x[n + offsets4] + weights5 * x[n + offsets5] + weights6 * x[n + offsets6];
			}
		else
			for (int n = lstart; n <= lstop; n++)
			{
				y[n] = -y[n] + weights0 * x[n] + weights1 * x[n + offsets1] + weights2 * x[n + offsets2] + weights3 * x[n + offsets3] + weights4 * x[n + offsets4];
			}
	}

	return true;
}



WavesWaveEqBoundary::WavesWaveEqBoundary(WavesSDGeometry *sdg, const SDBoundaryType &sdb_, double dt_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), sdb(sdb_), dt(dt_), freeIndexes(true)
{
//	cout << "//### WavesWaveEqBoundary::WavesWaveEqBoundary sdb" << int(sdb) << endl;
	doInitialize();
}

WavesWaveEqBoundary::WavesWaveEqBoundary(WavesSDBoundary *sdi, double dt_) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), dt(dt_), freeIndexes(false)
{
//	cout << "//### WavesWaveEqBoundary::WavesWaveEqBoundary sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesWaveEqBoundary::doInitialize()
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
			cout << "WavesWaveEqBoundary: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

void WavesWaveEqBoundary::initWeights()
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
			cout << "WavesWaveEqBoundary: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesWaveEqBoundary::doApply(const double *x, double *y) const
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

//======================================================================
// Class WavesWaveEqAdjBoundary compute boundary conditions for adjoint
// equation in Inverse Problem for wave equation
// b.c. in the form lambda grad P = -P grad(lambda)
// For approximation of the derivatives central differances hav been used


WavesWaveEqAdjBoundary::WavesWaveEqAdjBoundary(WavesSDGeometry *sdg, const SDBoundaryType &sdb_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), sdb(sdb_), freeIndexes(true)
{
	cout << "//### WavesWaveEqAdjBoundary::WavesWaveEqAdjBoundary sdb" << int(sdb) << endl;
	doInitialize();
}

WavesWaveEqAdjBoundary::WavesWaveEqAdjBoundary(WavesSDBoundary *sdi) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), freeIndexes(false)
{
	cout << "//### WavesWaveEqAdjBoundary::WavesWaveEqAdjBoundary sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesWaveEqAdjBoundary::doInitialize()
{
	switch (sdb)
	{

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
		cout << "WavesWaveEqAdjBoundary: ERROR, unknown or not implemented SDBoundaryType " << endl;
		break;
	}
}

void WavesWaveEqAdjBoundary::initWeights()
{
	// is called when sdb is known to be one side
	nrOfWeights = 2;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	offsets = new int[nrOfWeights];

	switch (sdb)
	{
		// for left boundary
		case SDiLow:
			offsets[0] = geometry->node(0, 1, 0);
			offsets[1] = geometry->node(0, 2, 0);

			break;
			// for low boundary
		case SDjLow:
			offsets[0] = geometry->node(1, 0, 0);
			offsets[1] = geometry->node(2, 0, 0);

			break;

		case SDkLow:
			cout << "not implemented" << endl;
			break;
			// for right boundary
		case SDiHigh:
			offsets[0] = geometry->node(-1, 0, 0);
			offsets[1] = geometry->node(-2, 0, 0);

			break;

			// for top boundary
		case SDjHigh:
			offsets[0] = geometry->node(0, -1, 0);
			offsets[1] = geometry->node(0, -2, 0);

			break;

		case SDkHigh:
			cout << "not implemented" << endl;
			break;

		default:
			cout << "WavesWaveEqAdjBoundary: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesWaveEqAdjBoundary::doApply(const double *f, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		// 
		for (int n = lstart; n <= lstop; n++)
		{
			double forward_sol = (f[n] - 2 * f[n + offsets[0]] + f[n + offsets[1]]) / f[n + offsets[0]];

			y[n] = -y[n + offsets[0]] * forward_sol + 2 * y[n + offsets[0]] - y[n + offsets[1]];

		}
	}

	return true;
}

//======================================================================


//
// class WavesWaveEqInterior2
//


WavesWaveEqInterior2::WavesWaveEqInterior2(WavesSDGeometry *sdg, double dt_) :
		WavesSDOperator(new WavesSDInterior(*sdg)), dt(dt_), freeIndexes(true)
{
	doInitialize();
}

WavesWaveEqInterior2::WavesWaveEqInterior2(WavesSDIndexes *sdi, double dt_) :
		WavesSDOperator(sdi), dt(dt_), freeIndexes(false)
{
	doInitialize();
}

void WavesWaveEqInterior2::doInitialize()
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	// 1) create 3*Laplace operator
	const int nsd = geometry->getNoSpaceDim();
	assert( nsd>=2 && nsd <= 3);

	nrOfWeights = (nsd == 2) ? 9 : 27;

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx2 = geometry->getDx() * geometry->getDx();
	real dy2 = geometry->getDy() * geometry->getDy();
	real dz2 = geometry->getDz() * geometry->getDz();

	if (nsd == 3)
	{
		cout << "//#### WavesWaveEqInterior2 not implemented for 3d yet" << endl;
		assert( nsd!=3);

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
		weights[0] = -2.0 / dx2 - 2.0 / dy2 - 4.0 / (dx2 + dy2);
		weights[1] = weights[2] = 1.0 / dx2;
		weights[3] = weights[4] = 1.0 / dy2;
		weights[5] = weights[6] = weights[7] = weights[8] = 2.0 / (dx2 + dy2);

		offsets[0] = 0;
		offsets[1] = geometry->node(-1, 0, 0);
		offsets[2] = geometry->node(1, 0, 0);
		offsets[3] = geometry->node(0, -1, 0);
		offsets[4] = geometry->node(0, 1, 0);
		offsets[5] = geometry->node(-1, -1, 0);
		offsets[6] = geometry->node(1, -1, 0);
		offsets[7] = geometry->node(1, 1, 0);
		offsets[8] = geometry->node(-1, 1, 0);
	}

	// 2) y := dt*dt (3*Laplace)/3.0 + 2
	for (int i = 0; i < nrOfWeights; i++)
		weights[i] *= (dt * dt) / 3.0;
	weights[0] += 2;
}

bool WavesWaveEqInterior2::doApply(const double *x, double *y) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	const int nsd = geometry->getNoSpaceDim();
	assert( nsd>=2 && nsd <= 3);

	// We apply y <- (dt^2 Laplace + 2) x - y

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		if (nsd == 3)
			for (int n = lstart; n <= lstop; n++)
				y[n] = -y[n] ;
		else
			for (int n = lstart; n <= lstop; n++)
				y[n] = -y[n] + weights[0] * x[n] + weights[1] * x[n + offsets[1]] + weights[2] * x[n + offsets[2]] + weights[3] * x[n + offsets[3]] + weights[4] * x[n + offsets[4]] + weights[5] * x[n + offsets[5]] + weights[6] * x[n + offsets[6]]
						+ weights[7] * x[n + offsets[7]] + weights[8] * x[n + offsets[8]];
	}
	return true;
}

