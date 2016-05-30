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

#include "include/wavesBCOperators.h"
#include <assert.h>

using namespace std;



WavesMirrorBC::WavesMirrorBC(WavesSDGeometry *sdg, const SDBoundaryType &sdb_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), freeIndexes(true), sdb(sdb_)
{
//	cout << "//### WavesMirrorBC::WavesMirrorBC sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

WavesMirrorBC::WavesMirrorBC(WavesSDBoundary *sdi) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), freeIndexes(false)
{
//	cout << "//### WavesMirrorBC::WavesMirrorBC sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesMirrorBC::doInitialize()
{
//	int i;
//	const WavesSDGeometry *geometry = &indexes->getGeometry();
//
//	assert(geometry->getNoSpaceDim() >= 2);

	switch (sdb)
	{
		case SDallBoundary:
			cout << "WavesMirrorBC:SDallBoundary not currently supported!\n" << endl;
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
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

void WavesMirrorBC::initWeights()
{
	// is called when sdb is known to be one side
	nrOfWeights = 1;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	weights[0] = 1.;

	switch (sdb)
	{
		case SDiLow:
			offsets[0] = geometry->node(2, 0, 0);
			break;

		case SDjLow:
			offsets[0] = geometry->node(0, 2, 0);
			break;

		case SDkLow:
			offsets[0] = geometry->node(0, 0, 2);
			break;

		case SDiHigh:
			offsets[0] = geometry->node(-2, 0, 0);
			break;

		case SDjHigh:
			offsets[0] = geometry->node(0, -2, 0);
			break;

		case SDkHigh:
			offsets[0] = geometry->node(0, 0, -2);
			break;

		default:
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesMirrorBC::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			y[n] = weights[0] * y[n + offsets[0]];
	}

	return true;
}

//======================================================================



WavesMirrorBC_YEE::WavesMirrorBC_YEE(WavesSDGeometry *sdg, const SDBoundaryType &sdb_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), freeIndexes(true), sdb(sdb_)
{
	cout << "//### WavesMirrorBC_YEE::WavesMirrorBC_YEE sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

WavesMirrorBC_YEE::WavesMirrorBC_YEE(WavesSDBoundary *sdi) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), freeIndexes(false)
{
	cout << "//### WavesMirrorBC_YEE::WavesMirrorBC_YEE sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesMirrorBC_YEE::doInitialize()
{
//	int i;
//	const WavesSDGeometry *geometry = &indexes->getGeometry();
//
//	assert(geometry->getNoSpaceDim() >= 2);

	switch (sdb)
	{
		case SDallBoundary:
			cout << "WavesMirrorBC:SDallBoundary not currently supported!\n" << endl;
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
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

void WavesMirrorBC_YEE::initWeights()
{
	// is called when sdb is known to be one side
	nrOfWeights = 1;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	weights[0] = 1.;

	switch (sdb)
	{
		case SDiLow:
			offsets[0] = geometry->node(3, 0, 0);
			break;

		case SDjLow:
			offsets[0] = geometry->node(0, 3, 0);
			break;

		case SDkLow:
			offsets[0] = geometry->node(0, 0, 3);
			break;

		case SDiHigh:
			offsets[0] = geometry->node(-3, 0, 0);
			break;

		case SDjHigh:
			offsets[0] = geometry->node(0, -3, 0);
			break;

		case SDkHigh:
			offsets[0] = geometry->node(0, 0, -3);
			break;

		default:
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesMirrorBC_YEE::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			y[n] = weights[0] * y[n + offsets[0]];
	}

	return true;
}



WavesMirrorBC1::WavesMirrorBC1(WavesSDGeometry *sdg, const SDBoundaryType &sdb_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), freeIndexes(true), sdb(sdb_)
{
	cout << "//### WavesMirrorBC1::WavesMirrorBC1 sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

WavesMirrorBC1::WavesMirrorBC1(WavesSDBoundary *sdi) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), freeIndexes(false)
{
	cout << "//### WavesMirrorBC1::WavesMirrorBC1 sdb" << int(sdb) << endl;
	doInitialize();
	unary = true;
}

void WavesMirrorBC1::doInitialize()
{
//	int i;
//	const WavesSDGeometry *geometry = &indexes->getGeometry();

//	assert(geometry->getNoSpaceDim() >= 2);

	switch (sdb)
	{
		case SDallBoundary:
			cout << "WavesMirrorBC1:SDallBoundary not currently supported!\n" << endl;
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
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

void WavesMirrorBC1::initWeights()
{
	// is called when sdb is known to be one side
	nrOfWeights = 1;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	weights[0] = 1.;

	switch (sdb)
	{
		case SDiLow:
			offsets[0] = geometry->node(1, 0, 0);
			break;

		case SDjLow:
			offsets[0] = geometry->node(0, 1, 0);
			break;

		case SDkLow:
			offsets[0] = geometry->node(0, 0, 1);
			break;

		case SDiHigh:
			offsets[0] = geometry->node(-1, 0, 0);
			break;

		case SDjHigh:
			offsets[0] = geometry->node(0, -1, 0);
			break;

		case SDkHigh:
			offsets[0] = geometry->node(0, 0, -1);
			break;

		default:
			cout << "WavesMirrorBC: ERROR, unknown or not implemented SDBoundaryType " << endl;
			break;
	}
}

bool WavesMirrorBC1::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			y[n] = weights[0] * y[n + offsets[0]];
	}

	return true;
}

//======================================================================
// smoothing for fdm solution
//========================================================================

SmoothingOperator::SmoothingOperator(WavesSDGeometry *sdg, const SDBoundaryType &sdb_) :
		WavesSDOperator(new WavesSDBoundary(*sdg, sdb_)), freeIndexes(true), sdb(sdb_)
{

	unary = true;
}

SmoothingOperator::SmoothingOperator(WavesSDBoundary *sdi) :
		WavesSDOperator(sdi), sdb(sdi->getSDBoundaryType()), freeIndexes(false)
{

	unary = true;
}

SmoothingOperator::SmoothingOperator(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(false)
{

	unary = true;
}

bool SmoothingOperator::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			if (y[n] < 1.0)
			{
				y[n] = 1.0;
				cout << " after smoothing " << y[n] << endl;
			}
	}

	return true;
}
