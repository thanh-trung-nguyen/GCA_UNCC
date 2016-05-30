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

#include "include/wavesElement.hh"

using namespace std;

WavesElement::WavesElement(int nonode, int nobfcn, int nsd_)
{
	cout << " inside WavesElement::WavesElement\n" << flush;
	nsd = nsd_;
	noNodes = nonode;
	noBFcns = nobfcn;
	coords.newsize(nonode, nsd_);
	globNodeNo.newsize(nonode);
	zeroBfcn.newsize(1, 1, nsd_);
	trial.newsize(1, 1, nsd_);
	test.newsize(1, 1, nsd_);
}

WavesBasisFcn* WavesElement::ptr2zeroBfcn()
{
	return &zeroBfcn;
}

WavesBasisFcn* WavesElement::ptr2test()
{
	return &test;
}

WavesBasisFcn* WavesElement::ptr2trial()
{
	return &trial;
}

int WavesElement::getNoBFcns() const
{
	return noBFcns;
}

void WavesElement::setIdx(int gpt, int trialIdx, int testIdx)
{
	test.set(gpt, testIdx);
	trial.set(gpt, trialIdx);
	zeroBfcn.set(gpt, testIdx);
}

void WavesElement::setIdx(int gpt, int testIdx)
{
	test.set(gpt, testIdx);
	zeroBfcn.set(gpt, testIdx);
}

real WavesElement::getCoord(int node, int coord) const
{
	return coords(node, coord);
}

int WavesElement::loc2glob(int locnode) const
{
	return globNodeNo(locnode);
}

int WavesElement::glob2loc(int globnode) const
{
	int locnode = -1;
	for (int i = 0; i < noNodes; i++)
		if (globNodeNo(i) == globnode)
			locnode = i;
	return locnode;
}

int WavesElement::getNoNodes() const
{
	return noNodes;
}

int WavesElement::getElmNo() const
{
	return elmNo;
}

int WavesElement::getMaterialType() const
{
	return materialType;
}

WavesElement& WavesElement::setNoGPnts(int nogp)
{
	noGPnts = nogp;
	coordGPnts.newsize(nogp, nsd);
	detJacs.newsize(nogp);
	zeroBfcn.newsize(nogp, noBFcns, nsd);
	test.newsize(nogp, noBFcns, nsd);
	trial.newsize(nogp, noBFcns, nsd);
	return *this;
}

int WavesElement::getNoGPnts()
{
	return noGPnts;
}

WavesElement& WavesElement::setGPntCoord(int gpt, int coor, real val)
{
	coordGPnts(gpt, coor) = val;
	return *this;
}

real& WavesElement::gPntCoord(int gpt, int coor)
{
	return coordGPnts(gpt, coor);
}

const real& WavesElement::gPntCoord(int gpt, int coor) const
{
	return coordGPnts(gpt, coor);
}

WavesElement& WavesElement::setDetJac(int gpt, real val)
{
	detJacs(gpt) = val;
	return *this;
}

real& WavesElement::detJac(int gpt)
{
	return detJacs(gpt);
}

const real& WavesElement::detJac(int gpt) const
{
	return detJacs(gpt);
}

WavesElement& WavesElement::update(int el, Grid& g)
{

	elmNo = el;
	noNodes = g.getNoNodesInElm(el);
	materialType = g.getMaterialType(el);

	globNodeNo.newsize(noNodes);
	coords.newsize(noNodes, nsd);

	for (int i = 0; i < noNodes; i++)
	{
		int tmp = g.loc2glob(elmNo, i);
		globNodeNo(i) = tmp;
		for (int j = 0; j < nsd; j++)
			coords(i, j) = g.getCoor(tmp, j);
	}
	return *this;

}

