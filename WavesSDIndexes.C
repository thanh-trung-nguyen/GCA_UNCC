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

#include "include/wavesSDIndexes.h"
#include "include/wavesSDGeometry.h"

#include <string.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <assert.h>

using namespace std;

#define MAX_I_INDEX 300000 /*any number smaller than a typical pointer*/

//
//--  Constructors etc
//
#include "stdio.h"
WavesSDIndexes::WavesSDIndexes(const WavesSDGeometry &sdg_) :
		sdg(sdg_), nOfLoops(0), loopIndexes(0), curLoop(0), nOfNodes(0)
{
//	printf("creating SDindexes %x, %x 1\n", (unsigned int) this, (unsigned int) loopIndexes);
}

WavesSDIndexes::WavesSDIndexes(const WavesSDIndexes &sdi) :
		sdg(sdi.sdg), nOfLoops(sdi.nOfLoops), curLoop(0), nOfNodes(sdi.nOfNodes)
{
	loopIndexes = new LoopIndex[nOfLoops];
	assert(nOfLoops > 0);
	copyLoopIndexes(sdi);
//	printf("creating SDindexes %x, %x 2\n", (unsigned int) this, (unsigned int) loopIndexes);
}

WavesSDIndexes::WavesSDIndexes() :
		nOfLoops(0), loopIndexes(0), curLoop(0), nOfNodes(0)
{
//	printf("creating SDindexes %x, %x 3\n", (unsigned int) this, (unsigned int) loopIndexes);
}

WavesSDIndexes::~WavesSDIndexes()
{
	/*  WavesSDIndexes has two constructors, one of them allocates loopIndexes,
	 and the other initialises it with NULL.
	 The destructor should free the memory only if loopIndex is not NULL.
	 Therefore it is important that loopIndex never contains a freed pointer.
	 */
	if (loopIndexes != 0)
	{
		assert(nOfLoops > 0);

		/* Protection against double free.
		 * Most frees overwrite the first words, so
		 * just checking iStart and iStop catches many problems */
		assert( loopIndexes->iStart >= 0 && loopIndexes->iStart <= loopIndexes->iStop && loopIndexes->iStop < MAX_I_INDEX);

		/* the old protection against double free */

		assert(loopIndexes->k >= 0);
		loopIndexes->k = -1;
		delete[] loopIndexes;
	}
}

WavesSDIndexes &WavesSDIndexes::operator=(const WavesSDIndexes &sdi)
{
	if (this == &sdi)
		return *this;

	printf("assigning sdindexes %x, %x := \n", (unsigned int) this, (unsigned int) loopIndexes, &sdi);

	assert(loopIndexes != 0);
	assert(nOfLoops > 0);
	loopIndexes->k = -2;
	delete[] loopIndexes;

	sdg = sdi.sdg;
	nOfLoops = sdi.nOfLoops;
	loopIndexes = new LoopIndex[nOfLoops];
	copyLoopIndexes(sdi);
	printf("assigned sdindexes %x, %x\n", (unsigned int) this, (unsigned int) loopIndexes);
	nOfNodes = sdi.nOfNodes;
	return *this;
}

void WavesSDIndexes::copyLI(WavesSDIndexes &sdi, int nOfLoops_)
{
	curLoop = 0;
	nOfLoops = nOfLoops_;
	assert(loopIndexes == 0);
	loopIndexes = new LoopIndex[nOfLoops];

	while (curLoop < nOfLoops)
	{

		loopIndexes[curLoop] = sdi.loopIndexes[curLoop];
		curLoop++;
	}

	nOfNodes = sdi.nOfNodes;
}

void WavesSDIndexes::copyinnerLI(WavesSDIndexes &sdi, int nOfLoops_)
{
	curLoop = 1;
	nOfLoops = nOfLoops_;
	assert(loopIndexes == 0);
	loopIndexes = new LoopIndex[nOfLoops];

	while (curLoop < nOfLoops)
	{

		loopIndexes[curLoop] = sdi.loopIndexes[curLoop];
		curLoop++;
	}

	nOfNodes = sdi.nOfNodes;
	sdi.presentWavesSDIndexes();
	exit(1);
}

// Protected methods for adding loops during creation
void WavesSDIndexes::copyLoopIndexes(const WavesSDIndexes &sdi)
{
	curLoop = 0;
	while (curLoop < nOfLoops)
	{
		loopIndexes[curLoop] = sdi.loopIndexes[curLoop];
		curLoop++;
	}
}

// method  for nOfLoops 

void WavesSDIndexes::resize(int n)
{

	tmp = new LoopIndex[n];
	if (!tmp)
	{
		cerr << "********** out of memory *********\n";
		exit(1);
	}
	for (int i = 0; i < nOfLoops && i < n; i++)
	{
		tmp[i] = loopIndexes[i];
	}

	assert(loopIndexes);
	assert(nOfLoops > 0);
	loopIndexes->k = -3;
	delete[] loopIndexes;
	loopIndexes = tmp;
	nOfLoops = n;
}

void WavesSDIndexes::doAddLoop(const int &iStart_, const int &iStop_, const int &j_, const int &k_)
{
	assert(iStart_ <= iStop_);

	if (curLoop >= nOfLoops)
	{
		resize(nOfLoops < 8 ? 8 : 2 * nOfLoops);
	}assert(curLoop < nOfLoops);
	loopIndexes[curLoop].iStart = iStart_;
	loopIndexes[curLoop].iStop = iStop_;
	loopIndexes[curLoop].j = j_;
	loopIndexes[curLoop].k = k_;
	nOfNodes += iStop_ - iStart_ + 1;
	curLoop++;
}

void WavesSDIndexes::insert(const int &nodeLow, const int &nodeHigh)
{

	const int nsd = sdg.getNoSpaceDim();
	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
	const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);

	assert(iStop - iStart >= 0);
	assert(jStop - jStart >= 0);
	assert(kStop - kStart >= 0);

	nOfLoops = (jStop - jStart + 1) * (kStop - kStart + 1);

	assert(nOfLoops > 0);
	assert(loopIndexes == 0);
	loopIndexes = new LoopIndex[nOfLoops];

	insertOnly(nodeLow, nodeHigh);

	assert(curLoop == nOfLoops);

}

void WavesSDIndexes::insertOnly(const int &nodeLow, const int &nodeHigh)
{

	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
	const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);
	assert(iStop - iStart >= 0);
	assert(jStop - jStart >= 0);
	assert(kStop - kStart >= 0);

	assert(loopIndexes);

	for (int k = kStart; k <= kStop; k++)
		for (int j = jStart; j <= jStop; j++)
			doAddLoop(iStart, iStop, j, k);

}

// If LoopIndexArray is allocated for nOfLoops elements, shrinkLoopIndexArray
// resizes the array to get curLoop elements, and sets nOfLoops to curLoop.
void WavesSDIndexes::schrinkLoopIndexArray()
{
	assert(curLoop <= nOfLoops);

	if (curLoop >= nOfLoops + 10) // No need to resize
		return;

	nOfLoops = curLoop;
	LoopIndexArray tmp = loopIndexes;
	loopIndexes = new LoopIndex[nOfLoops];
	assert(loopIndexes);
	assert(tmp);

	for (int i = 0; i < nOfLoops; i++)
		loopIndexes[i] = tmp[i];

	delete[] tmp;
}

//
//-- Methods for loops
//

int WavesSDIndexes::loopStart(const int &loop) const
{
	assert(loop >= 0 && loop < nOfLoops);

	return sdg.node(loopIndexes[loop].iStart, loopIndexes[loop].j, loopIndexes[loop].k);
}

int WavesSDIndexes::loopStop(const int &loop) const
{
	assert(loop >= 0 && loop < nOfLoops);
	return sdg.node(loopIndexes[loop].iStop, loopIndexes[loop].j, loopIndexes[loop].k);
}

// - Misc

ostream &operator<<(ostream &os, const LoopIndex &li)
{
	if (0) // ilya 08.12.02:too much output
		os << "[" << li.iStart << " - " << li.iStop << ", " << li.j << ", " << li.k << "]" << endl;
	return os;
}

void WavesSDIndexes::presentWavesSDIndexes() const
{
//	cout << "----- Present WavesSDIndexes " << endl;
//	cout << "  Loops (" << nOfLoops << "):" << endl;
//	cout << "  curLoop (" << curLoop << "):" << endl;
//	cout << "  Nodes (" << nOfNodes << "):" << endl;
//	for (int loop = 0; loop < nOfLoops; loop++)
//		cout << getLoop(loop);
//	cout << endl;
	//cout << "- " << endl;
}

// Present all global node indices
void WavesSDIndexes::presentGlobalIndexes() const
{
	for (int loop = 0; loop < nOfLoops; loop++)
		for (int i = loopIndexes[loop].iStart; i <= loopIndexes[loop].iStop; i++)
		{
			cout << " i = " << i << "  loop " << loop << endl;
			cout << sdg.node(i, loopIndexes[loop].j, loopIndexes[loop].k) << " ";
		}
}

// Make an array with indices of global node numbers
// indexArray shall be allocated with space for nOfNodeIndex elements 
void WavesSDIndexes::makeIndexArray(int *indexArray) const
{
	int index = 0;
	for (int loop = 0; loop < nOfLoops; loop++)
		for (int i = loopIndexes[loop].iStart; i <= loopIndexes[loop].iStop; i++)
		{
			indexArray[index] = sdg.node(i, loopIndexes[loop].j, loopIndexes[loop].k);
			index++;
		}
}

// Make an array with indices of global node numbers
// indexArray shall be allocated with space for nOfNodeIndex elements; assign these values to array AssignArray 
void WavesSDIndexes::AssignIndexArray(int *indexArray, double* Array, double *AssignArray) const
{
	int index = 0;
	for (int loop = 0; loop < nOfLoops; loop++)
		for (int i = loopIndexes[loop].iStart; i <= loopIndexes[loop].iStop; i++)
		{
			AssignArray[index] = Array[indexArray[index]];
			index++;
		}
}

// write above indxline value in AssignArray to file
void WavesSDIndexes::WriteIndexArray(char *fname, double* Array)
{
	ofstream outp;

	outp.open(fname);

	int index = 0;
	for (int loop = 0; loop < nOfLoops; loop++)
	{
		for (int i = loopIndexes[loop].iStart; i <= loopIndexes[loop].iStop; i++)
		{
			outp << Array[index] << "   ";
			index++;
		}
	}

	outp.close();

	printf("Solution at one point is written to file: %s\n", fname);
}

// Useful for output. In general slow, but can be fast for some cases.
bool WavesSDIndexes::containWavesElement(int e)
{
	if (sdg.getNoSpaceDim() == 2)
		return containNode(sdg.loc2glob(e, 0)) && containNode(sdg.loc2glob(e, 1)) && containNode(sdg.loc2glob(e, 2)) && containNode(sdg.loc2glob(e, 3));
	else if (sdg.getNoSpaceDim() == 3)
		return containNode(sdg.loc2glob(e, 0)) && containNode(sdg.loc2glob(e, 1)) && containNode(sdg.loc2glob(e, 2)) && containNode(sdg.loc2glob(e, 3)) && containNode(sdg.loc2glob(e, 4)) && containNode(sdg.loc2glob(e, 5)) && containNode(sdg.loc2glob(e, 6))
				&& containNode(sdg.loc2glob(e, 7));
	else
		return false;
}
// Useful for output. In general slow, but can be fast for some cases.
bool WavesSDIndexes::containNode(int node)
{
	int n = 0;
	int i, j, k;
	i = sdg.node2i(node);
	j = sdg.node2j(node);
	k = sdg.node2k(node);

	while (n < nOfLoops)
	{
		if (loopIndexes[n].j == j && loopIndexes[n].k == k && loopIndexes[n].iStart <= i && loopIndexes[n].iStop >= i)
			return true;
		n++;
	}
	return false;
}

//======================================================================

/*
 WavesSDInterior produces LoopIndexes for the interior of a grid. 
 The default width=1 means that the boundary is excluded.
 */

WavesSDInterior::WavesSDInterior(const WavesSDGeometry &sdg_, const int &w) :
		WavesSDIndexes(sdg_)
{
	assert(w >= 0);

	const int nsd = sdg.getNoSpaceDim();
	const int n_i = sdg.getN_i();
	const int n_j = sdg.getN_j();
	const int n_k = sdg.getN_k();

	int nodeLow, nodeHigh;

	if (nsd == 2)
	{
		nodeLow = sdg.node(w, w);
		nodeHigh = sdg.node(n_i - 1 - w, n_j - 1 - w);
		insert(nodeLow, nodeHigh);
	}
	else if (nsd == 3)
	{
		nodeLow = sdg.node(w, w, w);
		nodeHigh = sdg.node(n_i - 1 - w, n_j - 1 - w, n_k - 1 - w);
		insert(nodeLow, nodeHigh);
	}
	else
	{
		cout << "Init of WavesSDIndexes for " << nsd << "space dimensions not implemented." << endl;
	}
}

WavesSDInterior::WavesSDInterior(const WavesSDGeometry &sdg_, const int &nodeLow, const int &nodeHigh) :
		WavesSDIndexes(sdg_)
{
	insert(nodeLow, nodeHigh);
}

WavesSDInterior2::WavesSDInterior2(const WavesSDGeometry &sdg_, const int &w) :
		WavesSDIndexes(sdg_)
{
	assert(w >= 0);

	const int nsd = sdg.getNoSpaceDim();
	const int n_i = sdg.getN_i();
	const int n_j = sdg.getN_j();
	const int n_k = sdg.getN_k();

	int nodeLow, nodeHigh;

	if (nsd == 2)
	{
		nodeLow = sdg.node(w, w);
		nodeHigh = sdg.node(n_i - 1 - w, n_j - 1 - w);
		init(nodeLow, nodeHigh);
	}
	else if (nsd == 3)
	{
		nodeLow = sdg.node(w, w, w);
		nodeHigh = sdg.node(n_i - 1 - w, n_j - 1 - w, n_k - 1 - w);
		init(nodeLow, nodeHigh);
	}
	else
	{
		cout << "Init of WavesSDIndexes for " << nsd << "space dimensions not implemented." << endl;
	}
}

WavesSDInterior2::WavesSDInterior2(const WavesSDGeometry &sdg_, const int &nodeLow, const int &nodeHigh) :
		WavesSDIndexes(sdg_)
{
	init(nodeLow, nodeHigh);
}

void WavesSDInterior2::init(const int &nodeLow, const int &nodeHigh)
{
	iStart = sdg.node2i(nodeLow);
	iStop = sdg.node2i(nodeHigh);
	jStart = sdg.node2j(nodeLow);
	jStop = sdg.node2j(nodeHigh);
	kStart = sdg.node2k(nodeLow);
	kStop = sdg.node2k(nodeHigh);
	assert(iStop - iStart >= 0);
	assert(jStop - jStart >= 0);
	assert(kStop - kStart >= 0);

	nOfLoops = (jStop - jStart + 1) * (kStop - kStart + 1);

	assert(nOfLoops > 0);

	nOfNodes = nOfLoops * (iStop - iStart + 1);

	jLength = jStop - jStart + 1;
}

/// Returns "global" node number for where a loop starts: (Overrides base class version)
int WavesSDInterior2::loopStart(const int &loop) const
{
	return sdg.node(iStart, jStart + loop % jLength, kStart + loop / jLength);
}

/// Returns "global" node number for where a loop stops: (Overrides base class version)
int WavesSDInterior2::loopStop(const int &loop) const
{
	return sdg.node(iStop, jStart + loop % jLength, kStart + loop / jLength);

}

LoopIndex WavesSDInterior2::getLoop(const int &loop) const
{
	return LoopIndex(iStart, iStop, jStart + loop % jLength, kStart + loop / jLength);
}

//======================================================================

/*
 WavesSDBoundary produces LoopIndexes for the boundary depending on the enum 
 variable. 


 */

WavesSDBoundary::WavesSDBoundary()
{
}

// By default, this constructor *does* include nodeLow and nodeHigh for a face.

WavesSDBoundary::WavesSDBoundary(const WavesSDGeometry &sdg_, const int &nodeLow, const int &nodeHigh, const SDBoundaryType &sdb_, const bool &doIncludeCorners) :
		WavesSDIndexes(sdg_), sdb(sdb_)
{
	const int nsd = sdg.getNoSpaceDim();
	;

	cout << "  nodeLow " << nodeLow << " ndeHigh" << nodeHigh << endl;

	if (nsd == 2)
		initialize2d(nodeLow, nodeHigh, doIncludeCorners);
	else if (nsd == 3)
		initialize3d(nodeLow, nodeHigh, doIncludeCorners);
	else
	{
		cout << "Init of WavesSDBoundary for " << nsd << "space dimensions not implemented." << endl;
	}

}

// By default, this constructor *does not* include nodeLow and nodeHigh for a face.

WavesSDBoundary::WavesSDBoundary(const WavesSDGeometry &sdg_, const SDBoundaryType &sdb_, const bool &doIncludeCorners) :
		WavesSDIndexes(sdg_), sdb(sdb_)
{
	const int nsd = sdg.getNoSpaceDim();
	const int n_i = sdg.getN_i();
	const int n_j = sdg.getN_j();
	const int n_k = sdg.getN_k();

	if (nsd == 2)
		initialize2d(sdg.node(0, 0), sdg.node(n_i - 1, n_j - 1), doIncludeCorners);
	else if (nsd == 3)
		initialize3d(sdg.node(0, 0, 0), sdg.node(n_i - 1, n_j - 1, n_k - 1), doIncludeCorners);
	else
	{
		cout << "Init of WavesSDBoundary for " << nsd << "space dimensions not implemented." << endl;
	}
}
//========================================

WavesSDBoundary::WavesSDBoundary(const WavesSDGeometry &sdg_, int code, const SDBoundaryType &sdb_, const bool &doIncludeCorners) :
		WavesSDIndexes(sdg_), sdb(sdb_)
{
	const int nsd = sdg.getNoSpaceDim();
	const int n_i = sdg.getN_i();
	const int n_j = sdg.getN_j();
	const int n_k = sdg.getN_k();

	if (nsd == 2)
		initialize2d(sdg.node(0, 1), sdg.node(n_i - 1, n_j - 1), doIncludeCorners);
	else if (nsd == 3)
		initialize3d(sdg.node(0, 1, 0), sdg.node(n_i - 1, n_j - 1, n_k - 1), doIncludeCorners);
	else
	{
		cout << "Init of WavesSDBoundary for " << nsd << "space dimensions not implemented." << endl;
	}
}

//==========================================================

void WavesSDBoundary::initProblem(const WavesSDGeometry& sdg_, const SDBoundaryType &sdb_, const bool &doIncludeCorners)
{
	sdg = sdg_;
	sdb = sdb_;
	cout << "works initialize" << endl;
	const int n_i = sdg.getN_i();
	const int n_j = sdg.getN_j();
	const int n_k = sdg.getN_k();

	int nodeLow = sdg.node(0, 0);
	int nodeHigh = sdg.node(n_i - 1, n_j - 1);

	if (sdg.getNoSpaceDim() == 2)
	{
		const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
		const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
		const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);

		assert(iStop - iStart >= 2);
		assert(jStop - jStart >= 2);
		assert(kStop - kStart >= 0);

		// If corners shall not be included, D=1 adjusts the region.

		nOfLoops = (jStop - jStart + 1) * (kStop - kStart + 1);

		assert(nOfLoops > 0);

		assert(loopIndexes == 0);
		loopIndexes = new LoopIndex[nOfLoops];
		curLoop = 0;
		presentWavesSDIndexes();

		assert(curLoop == nOfLoops);

		assert(loopIndexes);

		for (int k = kStart; k <= kStop; k++)
			for (int j = jStart; j <= jStop; j++)
				doAddLoop(iStart, iStop, j, k);

	}
}
void WavesSDBoundary::initialize2d(const int &nodeLow, const int &nodeHigh, const bool &doIncludeCorners)
{

	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);

	assert(iStop - iStart >= 2);
	assert(jStop - jStart >= 2);
	int D = doIncludeCorners ? 0 : 1;

	// If corners shall not be included, D=1 adjusts the region.

	switch (sdb)
	{
		case SDiLow:
			insert(sdg.node(iStart, jStart + D), sdg.node(iStart, jStop - D));
			break;
		case SDiHigh:
			insert(sdg.node(iStop, jStart + D), sdg.node(iStop, jStop - D));
			break;

		case SDjLow:
			insert(sdg.node(iStart + D, jStart), sdg.node(iStop - D, jStart));
			break;
		case SDjHigh:
			insert(sdg.node(iStart + D, jStop), sdg.node(iStop - D, jStop));
			break;

		case SDCorners:
			nOfLoops = 4;
			loopIndexes = new LoopIndex[nOfLoops];
			curLoop = 0;
			insertOnly(sdg.node(iStart, jStart), sdg.node(iStart, jStart));
			insertOnly(sdg.node(iStart, jStop), sdg.node(iStart, jStop));
			insertOnly(sdg.node(iStop, jStart), sdg.node(iStop, jStart));
			insertOnly(sdg.node(iStop, jStop), sdg.node(iStop, jStop));
			assert(curLoop == nOfLoops);
			break;

		case SDallBoundary:
			nOfLoops = 2 * (jStop - jStart + 1) + 2 * 1; // 2*iface + 2*jface

			loopIndexes = new LoopIndex[nOfLoops];
			insertOnly(sdg.node(iStart, jStart), sdg.node(iStart, jStop));
			insertOnly(sdg.node(iStop, jStart), sdg.node(iStop, jStop));

			insertOnly(sdg.node(iStart + 1, jStart), sdg.node(iStop - 1, jStart));
			insertOnly(sdg.node(iStart + 1, jStop), sdg.node(iStop - 1, jStop));

			assert(curLoop == nOfLoops);
			break;

		default:
			cout << "WavesSDBoundary::init2d, Error: Unkown option. sdb = " << sdb << endl;
			break;
	}

}

void WavesSDBoundary::initialize3d(const int &nodeLow, const int &nodeHigh, const bool &doIncludeCorners)
{
	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
	const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);

	assert(iStop - iStart >= 2);
	assert(jStop - jStart >= 2);
	assert(kStop - kStart >= 2);

	// If corners shall not be included, D=1 adjusts the region.
	int D = doIncludeCorners ? 0 : 1;

	switch (sdb)
	{
		case SDiLow:
			insert(sdg.node(iStart, jStart + D, kStart + D), sdg.node(iStart, jStop - D, kStop - D));
			break;
		case SDiHigh:
			insert(sdg.node(iStop, jStart + D, kStart + D), sdg.node(iStop, jStop - D, kStop - D));
			break;

		case SDjLow:
			insert(sdg.node(iStart + D, jStart, kStart + D), sdg.node(iStop - D, jStart, kStop - D));
			break;
		case SDjHigh:
			insert(sdg.node(iStart + D, jStop, kStart + D), sdg.node(iStop - D, jStop, kStop - D));
			break;

		case SDkLow:
			insert(sdg.node(iStart + D, jStart + D, kStart), sdg.node(iStop - D, jStop - D, kStart));
			break;
		case SDkHigh:
			insert(sdg.node(iStart + D, jStart + D, kStop), sdg.node(iStop - D, jStop - D, kStop));
			break;

		case SDCorners:
			nOfLoops = 8;
			loopIndexes = new LoopIndex[nOfLoops];
			curLoop = 0;
			insertOnly(sdg.node(iStart, jStart, kStart), sdg.node(iStart, jStart, kStart));
			insertOnly(sdg.node(iStart, jStop, kStart), sdg.node(iStart, jStop, kStart));
			insertOnly(sdg.node(iStop, jStart, kStart), sdg.node(iStop, jStart, kStart));
			insertOnly(sdg.node(iStop, jStop, kStart), sdg.node(iStop, jStop, kStart));

			insertOnly(sdg.node(iStart, jStart, kStop), sdg.node(iStart, jStart, kStop));
			insertOnly(sdg.node(iStart, jStop, kStop), sdg.node(iStart, jStop, kStop));
			insertOnly(sdg.node(iStop, jStart, kStop), sdg.node(iStop, jStart, kStop));
			insertOnly(sdg.node(iStop, jStop, kStop), sdg.node(iStop, jStop, kStop));
			assert(curLoop == nOfLoops);
			break;

		case SDallBoundary:
			nOfLoops = 2 * (jStop - jStart + 1) * (kStop - kStart + 1) + 2 * (kStop - kStart + 1) + 2 * (jStop - jStart - 1); // k-faces
			loopIndexes = new LoopIndex[nOfLoops];
			curLoop = 0;
			insertOnly(sdg.node(iStart, jStart, kStart), sdg.node(iStart, jStop, kStop));
			insertOnly(sdg.node(iStop, jStart, kStart), sdg.node(iStop, jStop, kStop));
			insertOnly(sdg.node(iStart + 1, jStart, kStart), sdg.node(iStop - 1, jStart, kStop));
			insertOnly(sdg.node(iStart + 1, jStop, kStart), sdg.node(iStop - 1, jStop, kStop));
			insertOnly(sdg.node(iStart + 1, jStart + 1, kStart), sdg.node(iStop - 1, jStop - 1, kStart));
			insertOnly(sdg.node(iStart + 1, jStart + 1, kStop), sdg.node(iStop - 1, jStop - 1, kStop));
			assert(curLoop == nOfLoops);
			break;
		default:
			cout << "WavesSDBoundary::Error::init3d, Error: Unkown option. sdb = " << sdb << endl;
			break;
	}

}

WavesSDInteriorWithHole::WavesSDInteriorWithHole()
{
}

WavesSDInteriorWithHole::WavesSDInteriorWithHole(const WavesSDGeometry &sdg_, const int &nodeLow, const int &nodeHigh, const int &nodeLowI, const int &nodeHighI) :
		WavesSDIndexes(sdg_), nodeLowInner(nodeLowI), nodeHighInner(nodeHighI)
{
	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
	const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);

	const int iStartI = sdg.node2i(nodeLowI), iStopI = sdg.node2i(nodeHighI);
	const int jStartI = sdg.node2j(nodeLowI), jStopI = sdg.node2j(nodeHighI);
	const int kStartI = sdg.node2k(nodeLowI), kStopI = sdg.node2k(nodeHighI);

	assert(iStop - iStart >= 0);
	assert(jStop - jStart >= 0);
	assert(kStop - kStart >= 0);

	assert(iStopI - iStartI >= 0);
	assert(jStopI - jStartI >= 0);
	assert(kStopI - kStartI >= 0);

	assert(iStart < iStartI);
	assert(iStop > iStopI);
	assert(jStart <= jStartI);
	assert(jStop >= jStopI);
	assert(kStart <= kStartI);
	assert(kStop >= kStopI);

	nOfLoops = (jStop - jStart + 1) * (kStop - kStart + 1) + (jStopI - jStartI + 1) * (kStopI - kStartI + 1);

	assert(nOfLoops > 0);
	loopIndexes = new LoopIndex[nOfLoops];
	curLoop = 0;

	assert(loopIndexes);

	// insert procedure
	for (int k = kStart; k <= kStop; k++)
		for (int j = jStart; j <= jStop; j++)
			if (k < kStartI || k > kStopI || j < jStartI || j > jStopI)
				// No overlap, just add
				doAddLoop(iStart, iStop, j, k);
			else
			{
				// Add two loops
				doAddLoop(iStart, iStartI - 1, j, k);
				doAddLoop(iStopI + 1, iStop, j, k);
			}

	assert(curLoop == nOfLoops);

}

WavesSDInteriorWithHole::WavesSDInteriorWithHole(const WavesSDInteriorWithHole &sdi) :
		WavesSDIndexes(sdi), nodeLowInner(sdi.nodeLowInner), nodeHighInner(sdi.nodeHighInner)
{
}

bool WavesSDInteriorWithHole::containNode(int node)
{
	const int iStartI = sdg.node2i(nodeLowInner), iStopI = sdg.node2i(nodeHighInner);
	const int jStartI = sdg.node2j(nodeLowInner), jStopI = sdg.node2j(nodeHighInner);
	const int kStartI = sdg.node2k(nodeLowInner), kStopI = sdg.node2k(nodeHighInner);

	int i, j, k;
	i = sdg.node2i(node);
	j = sdg.node2j(node);
	k = sdg.node2k(node);

	if (i >= iStartI && i <= iStopI && j >= jStartI && j <= jStopI && k >= kStartI && k <= kStopI)
		return false; // node is in the hole.

	if (node >= 0 && node < sdg.getNoNodes())
		return true;

	return false;
}

//======================================================================

WavesSDIndexLine::WavesSDIndexLine()
{
}

WavesSDIndexLine::WavesSDIndexLine(const WavesSDGeometry &sdg_, const int &nodeLow, const int &nodeHigh) :
		WavesSDIndexes(sdg_)
{
	// Check that nodes parallel to one axis.
	const int nsd = sdg.getNoSpaceDim();
	const int iStart = sdg.node2i(nodeLow), iStop = sdg.node2i(nodeHigh);
	const int jStart = sdg.node2j(nodeLow), jStop = sdg.node2j(nodeHigh);
	const int kStart = sdg.node2k(nodeLow), kStop = sdg.node2k(nodeHigh);
	assert(nsd == 2 || nsd == 3);

	if (nsd == 2)
		assert((iStart == iStop || jStart == jStop) && kStart == kStop);
	else // nsd==3
	if (iStart != iStop)
		assert(jStart == jStop && kStart == kStop);
	else if (jStart != jStop)
		assert(iStart == iStop && kStart == kStop);
	else if (kStart != kStop)
		assert(jStart == jStop && iStart == iStop);

	insert(nodeLow, nodeHigh);
}
