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

#include "include/wavesSDOperator.h"
#include "include/wavesOutputs.h"

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <assert.h>

using namespace std;

typedef MV_ColMat<int> Mat_int;

void WavesSDOperator::presentOperator(bool presentIndexes)
{
	cout << "WavesSDOperator::PresentStencil()" << endl;
	cout << "nrOfWeights: " << nrOfWeights << endl;
	cout << "unary" << (unary ? "T" : "F") << endl;

	cout << "Offset";
	cout << "Weight " << endl;
	for (int i = 0; i < nrOfWeights; i++)
	{
		cout << offsets[i] << "    " << weights[i] << endl;
	}
	cout << endl;

	if (presentIndexes)
		indexes->presentWavesSDIndexes();

	cout << " End of WavesSDOperator::PresentStencil() " << endl;
}

//======================================================================
//
// Class LaplacianOp
//
//======================================================================

LaplacianOp::LaplacianOp(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg)), freeIndexes(true)
{
	assert(sdg != 0);
	assert(indexes != 0);

	doInitialize();
}

LaplacianOp::LaplacianOp(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	doInitialize();
}

void LaplacianOp::doInitialize()
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	const int nsd = geometry->getNoSpaceDim();
	assert(nsd >= 2 && nsd <= 3);

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

bool LaplacianOp::doApply(const double *x, double *y) const
{
	const int nsd = indexes->getGeometry().getNoSpaceDim();
	assert(nsd >= 2 && nsd <= 3);

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		if (nsd == 3)
			for (int n = lstart; n <= lstop; n++)
				y[n] = weights[0] * x[n] + weights[1] * x[n + offsets[1]] + weights[2] * x[n + offsets[2]] + weights[3] * x[n + offsets[3]] + weights[4] * x[n + offsets[4]] + weights[5] * x[n + offsets[5]] + weights[6] * x[n + offsets[6]];
		else
			for (int n = lstart; n <= lstop; n++)
			{
				y[n] = weights[0] * x[n] + weights[1] * x[n + offsets[1]] + weights[2] * x[n + offsets[2]] + weights[3] * x[n + offsets[3]] + weights[4] * x[n + offsets[4]];
			}
	}

	return true;
}

//======================================================================
//
// Class WavesDplusDminusOp
//
//======================================================================

WavesDplusDminusOp::WavesDplusDminusOp(WavesSDGeometry *sdg, const int &d) :
		WavesSDOperator(new WavesSDInterior(*sdg)), freeIndexes(false)
{
	doInitialize(d);
}

WavesDplusDminusOp::WavesDplusDminusOp(WavesSDIndexes *sdi, const int &d) :
		WavesSDOperator(sdi), freeIndexes(true)
{
	doInitialize(d);
}

void WavesDplusDminusOp::doInitialize(const int &d)
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	const int nsd = geometry->getNoSpaceDim();
	assert(nsd >= 1 && nsd <= 3);

	nrOfWeights = 3;

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx2 = geometry->getDx() * geometry->getDx();
	real dy2 = geometry->getDy() * geometry->getDy();
	real dz2 = geometry->getDz() * geometry->getDz();

	if (d == 2)
	{
		assert(nsd == 3);
		weights[0] = -2.0 / dz2;
		weights[1] = weights[2] = 1.0 / dz2;

		offsets[0] = 0;
		offsets[1] = geometry->node(0, 0, -1);
		offsets[2] = geometry->node(0, 0, 1);
	}
	else if (d == 1)
	{
		assert(nsd >= 2);
		weights[0] = -2.0 / dy2;
		weights[1] = weights[2] = 1.0 / dy2;

		offsets[0] = 0;
		offsets[1] = geometry->node(0, -1, 0);
		offsets[2] = geometry->node(0, 1, 0);
	}
	else
	{
		assert(d == 0);
		weights[0] = -2.0 / dx2;
		weights[1] = weights[2] = 1.0 / dx2;

		offsets[0] = 0;
		offsets[1] = geometry->node(-1, 0, 0);
		offsets[2] = geometry->node(1, 0, 0);
	}
}

bool WavesDplusDminusOp::doApply(const double *x, double *y) const
{
	const int offsets1 = offsets[1];
	const int offsets2 = offsets[2];
	const real w0 = weights[0];
	const real w1 = weights[1];
	const real w2 = weights[2];

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			y[n] = w0 * x[n] + w1 * x[n + offsets1] + w2 * x[n + offsets2];
	}
	return true;
}

//======================================================================
//
// Class DplusDminusVec2Op
//
//======================================================================

DplusDminusVec2Op::DplusDminusVec2Op(WavesSDIndexes *sdi, const int &d) :
		WavesSDOperator(sdi)
{
	doInitialize(d);
}

void DplusDminusVec2Op::doInitialize(const int &d)
{
	unary = true;

	const WavesSDGeometry *geometry = &indexes->getGeometry();

	const int nsd = geometry->getNoSpaceDim();
	assert(nsd >= 1 && nsd <= 3);

	nrOfWeights = 3;

	weights = new real[nrOfWeights];
	offsets = new int[nrOfWeights];

	real dx2 = geometry->getDx() * geometry->getDx();
	real dy2 = geometry->getDy() * geometry->getDy();
	real dz2 = geometry->getDz() * geometry->getDz();

	// The offsets are multiplied with vector number of components 
	nComp = 2;

	if (d == 2)
	{
		assert(nsd == 3);
		weights[0] = -2.0 / dz2;
		weights[1] = weights[2] = 1.0 / dz2;

		offsets[0] = 0;
		offsets[1] = nComp * geometry->node(0, 0, -1);
		offsets[2] = nComp * geometry->node(0, 0, 1);
	}
	else if (d == 1)
	{
		assert(nsd >= 2);
		weights[0] = -2.0 / dy2;
		weights[1] = weights[2] = 1.0 / dy2;

		offsets[0] = 0;
		offsets[1] = nComp * geometry->node(0, -1, 0);
		offsets[2] = nComp * geometry->node(0, 1, 0);
	}
	else
	{
		assert(d == 0);
		weights[0] = -2.0 / dx2;
		weights[1] = weights[2] = 1.0 / dx2;

		offsets[0] = 0;
		offsets[1] = nComp * geometry->node(-1, 0, 0);
		offsets[2] = nComp * geometry->node(1, 0, 0);
	}
}

bool DplusDminusVec2Op::doApply(const double *, double *y) const
{
	const int offsets1 = offsets[1];
	const int offsets2 = offsets[2];
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		// The actual loop start and stop depends on nComp.
		// Note that component 1 is updated.
		for (int n = nComp * lstart; n <= nComp * lstop; n += nComp)
			y[n + 1] = weights[0] * y[n] + weights[1] * y[n + offsets1] + weights[2] * y[n + offsets2];
	}
	return true;
}

//======================================================================
//
// Class DirichletOp
//
//======================================================================

//## Should inherit from assignment?

DirichletOp::DirichletOp(WavesSDGeometry *sdg, real value_) :
		WavesSDOperator(new WavesSDBoundary(*sdg)), freeIndexes(true), value(value_)
{
	unary = true;
}

DirichletOp::DirichletOp(WavesSDIndexes *sdi, real value_) :
		WavesSDOperator(sdi), freeIndexes(false), value(value_)
{
	unary = true;
}

bool DirichletOp::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			y[n] = value;
	}
	return true;
}

//======================================================================
//
// Class WavesAssignmentOp
//
//======================================================================

WavesAssignmentOp::WavesAssignmentOp(WavesSDGeometry *sdg, real value_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), value(value_)
{
	unary = true;
}

WavesAssignmentOp::WavesAssignmentOp(WavesSDIndexes *sdi, real value_) :
		WavesSDOperator(sdi), freeIndexes(false), value(value_)
{
	unary = true;
}

bool WavesAssignmentOp::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = value;
		}
	}
	return true;
}

//======================================================================
//
// Class WavesAssignmentArrayOp
//
//======================================================================

WavesAssignmentArrayOp::WavesAssignmentArrayOp(WavesSDGeometry *sdg, real* array_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), array(array_)
{
	unary = true;
}

WavesAssignmentArrayOp::WavesAssignmentArrayOp(WavesSDIndexes *sdi, real* error) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
}

bool WavesAssignmentArrayOp::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = array[n];
			cout << "uOuter(" << n << ")=" << y[n] << endl;
		}
	}
	return true;
}

bool WavesAssignmentArrayOp::InnerBoundIndex(MV_Vector<int>& boundindex) const
{
	WavesOutputs out;

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();

		for (int n = lstart; n <= lstop; n++)
		{
			boundindex(n) = n;
		}
	}
	out.WriteToFile((char *) "InnerBoundIndex.m", boundindex);
	return true;

}

//======================================================================
//
// Class WavesAssignFunctionOp
//
//======================================================================

WavesAssignFunctionOp::WavesAssignFunctionOp(WavesSDGeometry *sdg, Fcn3d f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;
}

WavesAssignFunctionOp::WavesAssignFunctionOp(WavesSDIndexes *sdi, Fcn3d f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

bool WavesAssignFunctionOp::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();

		for (int n = lstart; n <= lstop; n++)
			y[n] = f(g->node2x(n), g->node2y(n), g->node2z(n));

	}
	return true;
}

//======================================================================
//
// Class AssignTimeFunctionOp
//
//======================================================================

AssignTimeFunctionOp::AssignTimeFunctionOp(WavesSDGeometry *sdg, Fcn3dtime f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;
}

AssignTimeFunctionOp::AssignTimeFunctionOp(WavesSDIndexes *sdi, Fcn3dtime f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

bool AssignTimeFunctionOp::doApply(const double *, double *y) const
{

	int n;
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();

		for (n = lstart; n <= lstop; n++)
		{
			y[n] = step * f(g->node2x(n), g->node2y(n), g->node2z(n), time);
		}
	}
	return true;
}

//======================================================================
//
// Class AssignAdjTimeFunctionOp
//
//======================================================================

AssignAdjTimeFunctionOp::AssignAdjTimeFunctionOp(WavesSDGeometry *sdg, real* f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;
}

AssignAdjTimeFunctionOp::AssignAdjTimeFunctionOp(WavesSDIndexes *sdi, real* f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

bool AssignAdjTimeFunctionOp::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();
		for (int n = lstart; n <= lstop; n++)
			y[n] = step * f[n];
	}
	return true;
}

//======================================================================
//
// Class ApplyFunctionOp
//
//======================================================================

ApplyFunctionOp::ApplyFunctionOp(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg)), freeIndexes(true)
{
	assert(sdg != 0);
	assert(indexes != 0);
}

ApplyFunctionOp::ApplyFunctionOp(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{

	assert(indexes != 0);
}

bool ApplyFunctionOp::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = x[n] + y[n];
		}
	}
	return true;
}

//======================================================================
//
// Class ComputeTimeDerivative
//
//======================================================================

ComputeTimeDerivative::ComputeTimeDerivative(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true)
{
	unary = true;
}

ComputeTimeDerivative::ComputeTimeDerivative(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
}

bool ComputeTimeDerivative::doApply(const double *, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);
		const WavesSDGeometry *g = &indexes->getGeometry();

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = (v_old_[n] + 3 * uOuter_[n] - 4 * vOuter_[n]) / (2 * step);
			//cout<<"duOuter/dt="<<y[n]<<endl;
		}
	}
	return true;
}

//======================================================================
//
// Class OutputOp
//
//======================================================================

OutputOp::OutputOp(WavesSDGeometry *sdg) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true)
{
	unary = true;
}

OutputOp::OutputOp(WavesSDIndexes *sdi) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
}

bool OutputOp::doApply(const double *x, double *y) const
{
	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			cout << n << " " << y[n] << endl;
	}
	return true;
}

//======================================================================
//
// Class AVSOutputOp
//
//======================================================================

AVSOutputOp::AVSOutputOp(WavesSDGeometry *sdg, ostream &ostr_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), ostr(ostr_)
{
	unary = true;
	interiorWithHole = false;
}

AVSOutputOp::AVSOutputOp(WavesSDIndexes *sdi, ostream &ostr_) :
		WavesSDOperator(sdi), freeIndexes(false), ostr(ostr_)
{
	unary = true;

	interiorWithHole = true;
}

bool AVSOutputOp::doApply(const double *x, double *y) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//*****  Output of geometry    *************

	cout << "### AVSOutputOp::doApply " << endl;
	int n, e;
	// Output of grid
	ostr << geometry->getNoNodes() << " " << geometry->getNoElms() << " " << 1 << " " // data per node
			<< 0 << " " // data per cell
			<< 0 << endl; // model data
	if (geometry->getNoSpaceDim() == 3)
		for (n = 0; n < geometry->getNoNodes(); n++)
			ostr << n << " " << geometry->node2x(n) << " " << geometry->node2y(n) << " " << geometry->node2z(n) << endl;
	else if (geometry->getNoSpaceDim() == 2)
		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			ostr << n << " " << geometry->node2x(n) << " "
					<< geometry->node2y(n) << " " << y[n] << endl;
		}
	else
	{
		cerr << "AVSOutputOp, error: AVS output not defined for 1D." << endl;
		return false;
	}

	//*****  Output of elements   *************
	int m = 1;
	for (e = 0; e < geometry->getNoElms(); e++)
	{
		if (interiorWithHole && indexes->containWavesElement(e))
			m = 1;
		else
			m = 99;
		if (geometry->getNoSpaceDim() == 3)
			ostr << e << " " << m << " hex " << geometry->loc2glob(e, 4) << " " << geometry->loc2glob(e, 5) << " " << geometry->loc2glob(e, 7) << " " << geometry->loc2glob(e, 6) << " " << geometry->loc2glob(e, 0) << " " << geometry->loc2glob(e, 1) << " "
					<< geometry->loc2glob(e, 3) << " " << geometry->loc2glob(e, 2) << endl;
		else
			ostr << e << " " << m << " quad " << geometry->loc2glob(e, 0) << " " << geometry->loc2glob(e, 1) << " " << geometry->loc2glob(e, 3) << " " << geometry->loc2glob(e, 2) << endl;
	}

	//*****   Output of data   *************
	ostr << "1 1" << endl; // 1 data 1 component
	ostr << "data, none" << endl; // label, unit

	/***
	 Even for grid with hole, I think all nodes shall be output.
	 ***/

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
			ostr << n << " " << y[n] << endl;
	}

	return true;
}

//  ====================  AVS output operator outer ============================

AVSOutputOpOuter::AVSOutputOpOuter(WavesSDIndexes *sdi, ostream &ostr_) :
		WavesSDOperator(sdi), freeIndexes(false), ostr(ostr_)
{
	unary = true;

	interiorWithHole = true;
}

bool AVSOutputOpOuter::doApply(const double *x, double *y) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//*****  Output of geometry    *************

	cout << "### AVSOutputOp::doApply " << endl;
	int n, e;
	// Output of grid
	ostr << geometry->getNoNodes() << " " << geometry->getNoElms() << " " << 1 << " " // data per node
			<< 0 << " " // data per cell
			<< 0 << endl; // model data
	if (geometry->getNoSpaceDim() == 3)
		for (n = 0; n < geometry->getNoNodes(); n++)
			ostr << n + 1 << " " << geometry->node2x(n) << " " << geometry->node2y(n) << " " << geometry->node2z(n) << endl;
	else if (geometry->getNoSpaceDim() == 2)
		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			ostr << n + 1 << " " << geometry->node2x(n) << " "
					<< geometry->node2y(n) << " " << "0.0" << endl;
		}
	else
	{
		cerr << "AVSOutputOp, error: AVS output not defined for 1D." << endl;
		return false;
	}

	//*****  Output of elements   *************
	int m = 1;
	for (e = 0; e < geometry->getNoElms(); e++)
	{
		if (interiorWithHole && indexes->containWavesElement(e))
			m = 1;
		else
			m = 99;

		if (geometry->getNoSpaceDim() == 3)
			ostr << e + 1 << " " << m << " hex " << geometry->loc2glob(e, 4) + 1 << " " << geometry->loc2glob(e, 5) + 1 << " " << geometry->loc2glob(e, 7) + 1 << " " << geometry->loc2glob(e, 6) + 1 << " " << geometry->loc2glob(e, 0) + 1 << " "
					<< geometry->loc2glob(e, 1) + 1 << " " << geometry->loc2glob(e, 3) + 1 << " " << geometry->loc2glob(e, 2) + 1 << endl;
		else
			ostr << e + 1 << " " << m << " quad " << geometry->loc2glob(e, 0) + 1 << " " << geometry->loc2glob(e, 1) + 1 << " " << geometry->loc2glob(e, 3) + 1 << " " << geometry->loc2glob(e, 2) + 1 << endl;

	}
	return true;
}

//========================================================================
// Class  GIDOutputOp
// GIDOutputOp  is used to printommon out 2D vectoral FDM solution at the nodes of FDM mesh.
// It can be used to print out  an interior mesh with hole or  a whole grid. 
// This class also allows print out  hybrid FDM solution when FEM solution is already printed out.
// Name of the output file with extension *.res should be provided.
//======================================================================

GIDOutputOp::GIDOutputOp(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_) :
		freeIndexes(true), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), indexes(new WavesSDInterior(*sdg, 0))
{
	cout << " in constructor GIDOutputOp(WavesSDGeometry *sdg,...) " << endl;

	unary = true;
	interiorWithHole = false;
}

GIDOutputOp::GIDOutputOp(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}

GIDOutputOp::GIDOutputOp(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int nnogg_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), nnogg(nnogg_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}

bool GIDOutputOp::printResults()
{

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int n, e;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//=========================================================
	if (geometry->getNoSpaceDim() == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside


		if (k == 1)
			fprintf(fp, "GiD Post Results File 1.0");

		fprintf(fp, "\n");
		fprintf(fp, "Result  \"Displacements\"  \"Wave movement\"  %i    Vector OnNodes   \n", k);
		fprintf(fp, "Values");
		fprintf(fp, "\n");

		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			fprintf(fp, "%i %f %f %f \n ", n + 1, u_1[n], u_2[n], u_array[n]);
		}

		fprintf(fp, "end values");
		fprintf(fp, "\n");

	}

	fclose(fp);

	return true;
}

//========================================================================================================================
// printResultsCommon() prints out hybrid FDM solution in *.res file after FEM solution (this solution should be already
// printed out *.res file).
//==========================================================================================================================
bool GIDOutputOp::printResultsCommon()
{

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int n, e;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//=========================================================
	if (geometry->getNoSpaceDim() == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			fprintf(fp, "%i %f %f %f \n ", nnogg + n + 1, u_1[n], u_2[n], u_array[n]);
		}
		fprintf(fp, "end values");
		fprintf(fp, "\n");
	}
	fclose(fp);

	return true;
}
//========================================================================
// Class GIDOutputOp3D
// GIDOutputOp3D is used to print out 3D vectoral FDM solution at the nodes of FDM mesh.
// It can be used to print out  an interior mesh with hole or  a whole grid. 
// This class also allows print out  hybrid FDM solution when FEM solution is already printed out.
// Name of the output file with extension *.res should be provided.
//======================================================================

GIDOutputOp3D::GIDOutputOp3D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_) :
		freeIndexes(true), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), u_3(u_3_), k(k_), indexes(new WavesSDInterior(*sdg))
{
	cout << " in constructor GIDOutputOp(WavesSDGeometry *sdg,...) " << endl;

	unary = true;
	interiorWithHole = false;
}

GIDOutputOp3D::GIDOutputOp3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), u_3(u_3_), k(k_), indexes(sdi)
{
	unary = true;
	interiorWithHole = true;
}


GIDOutputOp3D::GIDOutputOp3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_, int nnogg_) :
  freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_),  u_3(u_3_), k(k_), nnogg(nnogg_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}


bool GIDOutputOp3D::printResults()
{

 // Print out only FDM solution.
  // This method is used to print out 3D vector results of the FDM solution 
  // using GID output format
  

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int n, e;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//=========================================================
	if (geometry->getNoSpaceDim() == 3)
	{

		// here, we print out displacement of the solution vector in  3D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		fprintf(fp, "DISPLACEMENT   2   %i    2    1    0 \n", k);

		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			fprintf(fp, "%i %f %f %f  \n ", n + 1, u_1[n], u_2[n], u_3[n]);
		}

		fprintf(fp, "\n");

	}

	fclose(fp);

	return true;
}

bool GIDOutputOp3D::printResultsCommon() 
{
  
 FILE *fp;  
 fp = fopen(file, "a+");
  
  if (!fp) 
    {
      cout<< " cannot open file \""<<file<<"\"\n";
      exit(1);
    }
    

  fseek(fp,1,SEEK_END);

  int n;

  const WavesSDGeometry *geometry = &indexes->getGeometry();

  //=========================================================
  if( geometry->getNoSpaceDim() == 3)
    {
    


// Print out  FDM solution after FEM solution when we need 
// print out hybrid method  solution.
  // This method is used to print out 3D vector results of the FDM solution 
  // using GID output format

      for(n = 0; n < geometry->getNoNodes(); n++)
	{
	  /*
	  fprintf(fp, "%i %f %f %f %f  \n ", n + 1,
		  u_1[n] ,  u_2[n],  u_3[n], u_array[n]);
	  */

	  fprintf(fp, "%i %f %f %f  %f \n ", nnogg + n + 1,
		  u_1[n] ,  u_2[n],  u_3[n], u_array[n]);

	  //   cout<<" u_1("<<n<<") = "<<u_1[n]<<" u_2_"<<u_2[n]<<endl;
	}


      fprintf(fp,"end values");
      fprintf(fp,"\n");          
     
    
    }  


    fclose(fp);
 
  return true;
}

//========================================================================

// GIDOutputCommon
//======================================================================

GIDOutputCommon::GIDOutputCommon(WavesSDIndexes *sdi, char *file_, real* u_array_, real* u_1_, real* u_2_, int k_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}

bool GIDOutputCommon::printResults()
{

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int n, e;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//*****  Output of geometry    *************
	//=========================================================
	if (geometry->getNoSpaceDim() == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		cout << " geometry->getNoNodes() " << geometry->getNoNodes() << endl;

		for (n = 0; n < geometry->getNoNodes(); n++)
		{
			fprintf(fp, "%i %f %f %f \n ", n + k + 1, u_array[n], u_1[n], u_2[n]);
		}

		fprintf(fp, "\n");

	}
	fclose(fp);

	return true;
}
//========================================================================
// Out2D
//======================================================================

Out2D::Out2D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int n_) :
		freeIndexes(true), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), n(n_), indexes(new WavesSDInterior(*sdg, 0))
{
	cout << " in constructor Out2D(WavesSDGeometry *sdg,...) " << endl;

	unary = true;
	interiorWithHole = false;
}

Out2D::Out2D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_array_, int k_, int n_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), k(k_), n(n_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}

bool Out2D::printResults()
{

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	int e;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//=========================================================
	if (geometry->getNoSpaceDim() == 2)
	{

		// here, we print out displacement of the elastic vector in  2D
		// Vector_results is results title
		// 2 - type of analysis (displacement analysis)
		// %i - time step
		// 2 - kind of results - vector
		// 1 - position of the data - in nodes
		// 0 - no description inside

		if (k == 1)
			fprintf(fp, "Point Analysis   2   %i    2    1    0 \n", k);

		cout << " geometry->getNoNodes() " << geometry->getNoNodes() << endl;

		fprintf(fp, "%i %f %f %f  \n ", k, u_1[n], u_2[n], u_array[n]);

		fprintf(fp, "\n");

	}

	fclose(fp);

	return true;
}
//========================================================================
// Out3D
//======================================================================

Out3D::Out3D(WavesSDGeometry *sdg, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_, int n_) :
		freeIndexes(true), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), u_3(u_3_), k(k_), n(n_), indexes(new WavesSDInterior(*sdg, 0))
{
	cout << " in constructor GIDOut3D(WavesSDGeometry *sdg,...) " << endl;

	unary = true;
	interiorWithHole = false;
}

Out3D::Out3D(WavesSDIndexes *sdi, char *file_, real* u_1_, real* u_2_, real* u_3_, real* u_array_, int k_, int n_) :
		freeIndexes(false), file(file_), u_array(u_array_), u_1(u_1_), u_2(u_2_), u_3(u_3_), k(k_), n(n_), indexes(sdi)
{
	unary = true;

	interiorWithHole = true;
}

bool Out3D::printResults()
{
	cout << " u_1(" << n << ") = " << u_1[n] << " u_2_" << u_2[n] << " u_array[n] " << u_array[n] << endl;
	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);

	//=========================================================

	// here, we print out displacement of the elastic vector in  2D
	// Vector_results is results title
	// 2 - type of analysis (displacement analysis)
	// %i - time step
	// 2 - kind of results - vector
	// 1 - position of the data - in nodes
	// 0 - no description inside

	if (k == 1)
		fprintf(fp, "Point Analysis   2   %i    2    1    0 \n", k);

	cout << " u_1(" << n << ") = " << u_1[n] << " u_2_" << u_2[n] << endl;

	fprintf(fp, "%i %f %f %f %f  \n ", k, u_1[n], u_2[n], u_3[n], u_array[n]);

	fprintf(fp, "\n");

	fclose(fp);

	return true;
}

//========================================================================
// GIDOutputMesh
//======================================================================

GIDOutputMesh::GIDOutputMesh(WavesSDGeometry *sdg, char *file_) :
		freeIndexes(true), file(file_), indexes(new WavesSDInterior(*sdg, 0))
{
	unary = true;
	interiorWithHole = false;
}

GIDOutputMesh::GIDOutputMesh(WavesSDIndexes *sdi, char *file_) :
		freeIndexes(false), file(file_), indexes(sdi)
{
	unary = true;
	interiorWithHole = true;
}

GIDOutputMesh::GIDOutputMesh(WavesSDIndexes *sdi, char *file_, int nnogg_, int nelgg_) :
		freeIndexes(false), file(file_), indexes(sdi), nnogg(nnogg_), nelgg(nelgg_)
{
	unary = true;
	interiorWithHole = true;
}

bool GIDOutputMesh::printMesh()
{

	FILE *fp;
	fp = fopen(file, "w");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}
	int i, e, el;

	const WavesSDGeometry *geometry = &indexes->getGeometry();
	int nnode = geometry->getNoNodes();
	int nsd = geometry->getNoSpaceDim();

	cout << "  nnode = " << nnode << "  nsd = " << nsd << endl;

	//*****  Output of geometry    *************
	//=========================================================
	if (nsd == 3)
	{

		fprintf(fp, "MESH   dimension  %i  ElemType Hexahedra Nnode  8\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, geometry->getCoor(i, 0), geometry->getCoor(i, 1), geometry->getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		fprintf(fp, "MESH   dimension  %i  ElemType Quadrilateral  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + 1, geometry->getCoor(i, 0), geometry->getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");
	int m = 1;
	if (nsd == 2)
	{
		for (el = 0; el < geometry->getNoElms(); el++)
		{
			if (interiorWithHole && indexes->containWavesElement(el))
				m = 1;
			else
				m = 99;
			fprintf(fp, "%i   %i    %i   %i  %i %i\n", el + 1, geometry->loc2glob(el, 0) + 1, geometry->loc2glob(el, 1) + 1, geometry->loc2glob(el, 3) + 1, geometry->loc2glob(el, 2) + 1, m);

		}
	}
	else if (geometry->getNoSpaceDim() == 3)
	{

		for (el = 0; el < geometry->getNoElms(); el++)
		{
			if (interiorWithHole && indexes->containWavesElement(el))
				m = 1;
			else
				m = 99;
			fprintf(fp, "%i %i %i %i  %i %i %i %i %i %i \n", el + 1, geometry->loc2glob(el, 4) + 1, geometry->loc2glob(el, 5) + 1, geometry->loc2glob(el, 7) + 1, geometry->loc2glob(el, 6) + 1, geometry->loc2glob(el, 0) + 1, geometry->loc2glob(el, 1) + 1,
					geometry->loc2glob(el, 3) + 1, geometry->loc2glob(el, 2) + 1, m);
		}

	}
	fprintf(fp, "end elements\n");
	fclose(fp);

	return true;
}



bool GIDOutputMesh::printCommonMesh() 
{

 FILE *fp;  

 
fp = fopen(file, "a+"); 

  if (!fp)  
    {
      cout<< " cannot open file \""<<file<<"\"\n";
      exit(1);
    }
 fseek(fp,1,SEEK_END);



  int i,e,el;
 
  const WavesSDGeometry *geometry = &indexes->getGeometry();
  int nnode = geometry->getNoNodes();
  int nsd = geometry->getNoSpaceDim();
  int maxno = geometry->getMaxNoNodesInElm();
  Mat_int new_loc2glob_num(nnode,2);
  Mat_int new_loc2glob_el(nnode,maxno);

  cout<<"  nnogg = "<<nnogg<<",  nelgg ="<<nelgg<<",  nsd = "<<nsd<<endl;

  //*****  Output of geometry    *************
  //=========================================================
  if( nsd == 3)
    {
      
      fprintf(fp,"MESH   dimension  %i  ElemType Hexahedra Nnode  8\n", nsd);
      fprintf(fp,"Coordinates\n");
      
      for(i = 0; i < nnode; i++)
	{
	fprintf(fp, "%i %f %f %f\n ", nnogg + i + 1,
		geometry->getCoor(i,0), 
		geometry->getCoor(i,1),
		geometry->getCoor(i,2));

      //new global numeration for nodes in common mesh
	//	new_loc2glob_num(i,0) = i;
	//	new_loc2glob_num(i,1) = nnogg + i;
   } }
  else if( nsd == 2)
    {     
      fprintf(fp,"MESH   dimension  %i  ElemType Quadrilateral  Nnode  4\n", nsd);
      fprintf(fp,"Coordinates\n");
     
   
      for(i = 0; i < nnode; i++)
      {
	fprintf(fp, "%i %f %f %f\n ", nnogg +  i + 1,
		geometry->getCoor(i,0),
		geometry->getCoor(i,1),
		0.0);
	new_loc2glob_num(i,0) = i;
	new_loc2glob_num(i,1) = nnogg + i;
//	cout<<"	new_loc2glob_num(i,0) "<<new_loc2glob_num(i,0)<<"new_loc2glob_num(i,1)"<<new_loc2glob_num(i,1)<<endl;

      }
    }
 
 fprintf(fp,"end coordinates\n");
 fprintf(fp,"\n");
 fprintf(fp,"Elements\n");
 int m = 1;


	// Below we reassign global nodes numbers for some nodes in
	// FDM mesh.  We reassign global nodes numbers for such 
	// nodes in FDM mesh which have the same numeration as
	// the global nodes in FEM mesh.


 if ( nsd == 2)
   {

       for(el = 0; el < geometry->getNoElms(); el++)
       {
	   
	   for(i = 0; i < nnode; i++)
	   {
	       if ( geometry->loc2glob(el,0) == i)
		   new_loc2glob_el(el,0) =   new_loc2glob_num(i,1);
	  else if  ( geometry->loc2glob(el,1) == i)
	      new_loc2glob_el(el,1) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,2) == i)
		   new_loc2glob_el(el,2) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,3) == i)
		   new_loc2glob_el(el,3) =   new_loc2glob_num(i,1);
	       
	   }   
	   
       }

     for(el = 0; el < geometry->getNoElms(); el++)
       {
	 if( interiorWithHole && indexes->containWavesElement(el))
	   m = 98;
	 else
	   m = 99;

/*
	 fprintf(fp, "%i   %i    %i   %i  %i %i\n", nelgg + el + 1, 
		 geometry->loc2glob(el,0)+1,
		 geometry->loc2glob(el,1)+1, 
		 geometry->loc2glob(el,3)+1,
		 geometry->loc2glob(el,2)+1,m);
*/

	 fprintf(fp, "%i   %i    %i   %i  %i %i\n", nelgg + el + 1, 
		 new_loc2glob_el(el,0)+1,
		 new_loc2glob_el(el,1)+1, 
		 new_loc2glob_el(el,3)+1,
		 new_loc2glob_el(el,2)+1,m);
	
//	 cout<<" new_loc2glob_el(el,0)"<< new_loc2glob_el(el,0)<<" new_loc2glob_el(el,1) "<< new_loc2glob_el(el,1)<<"new_loc2glob_el(el,2)"<< new_loc2glob_el(el,2)<<" new_loc2glob_el(el,3)"<< new_loc2glob_el(el,3)<<endl;
 
       }
   }
 else if ( geometry->getNoSpaceDim() == 3)
   {
     
  for(el = 0; el < geometry->getNoElms(); el++)
       {
	 /*
	   for(i = 0; i < nnode; i++)
	     {
	       if ( geometry->loc2glob(el,0) == i)
		 new_loc2glob_el(el,0) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,1) == i)
		 new_loc2glob_el(el,1) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,2) == i)
		 new_loc2glob_el(el,2) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,3) == i)
		 new_loc2glob_el(el,3) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,4) == i)
		 new_loc2glob_el(el,4) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,5) == i)
		 new_loc2glob_el(el,5) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,6) == i)
		 new_loc2glob_el(el,6) =   new_loc2glob_num(i,1);
	       else if  ( geometry->loc2glob(el,7) == i)
		 new_loc2glob_el(el,7) =   new_loc2glob_num(i,1);
	       
	   }   
	 */

  	 new_loc2glob_el(el,0) =  geometry->loc2glob(el,0) + nnogg;
	 new_loc2glob_el(el,1) =  geometry->loc2glob(el,1) + nnogg;
	 new_loc2glob_el(el,2) =  geometry->loc2glob(el,2) + nnogg;
	 new_loc2glob_el(el,3) =  geometry->loc2glob(el,3) + nnogg;
	 new_loc2glob_el(el,4) =  geometry->loc2glob(el,4) + nnogg;
	 new_loc2glob_el(el,5) =  geometry->loc2glob(el,5) + nnogg;
	 new_loc2glob_el(el,6) =  geometry->loc2glob(el,6) + nnogg;
	 new_loc2glob_el(el,7) =  geometry->loc2glob(el,7) + nnogg;


       }

     
      for(el = 0; el < geometry->getNoElms(); el++ )
	{
	  if( interiorWithHole && indexes->containWavesElement(el) ) 
	    m = 98;
	  else
	    m = 99;



	  fprintf(fp, "%i %i %i %i  %i %i %i %i %i %i \n", nelgg + el + 1,
		  	 new_loc2glob_el(el,4)+1,
		  	 new_loc2glob_el(el,5)+1,
		  	 new_loc2glob_el(el,7)+1,
		  	 new_loc2glob_el(el,6)+1,
		  	 new_loc2glob_el(el,0)+1,
		  	 new_loc2glob_el(el,1)+1,
		  	 new_loc2glob_el(el,3)+1,
		         new_loc2glob_el(el,2)+1,m);
	}


   }
 fprintf(fp,"end elements\n"); 
fclose(fp);

return true;
}


//=====================================================================
//  Class for output of the common gg - sdg geometry into gid-file
//=====================================================================

//========================================================================
// GIDOutputNodes
//======================================================================

GIDOutputNodes::GIDOutputNodes(WavesSDIndexes *sdi, char *file_, int nno_gg_) :
		freeIndexes(false), file(file_), indexes(sdi), nno_gg(nno_gg_)
{
	unary = true;
	interiorWithHole = true;
}

bool GIDOutputNodes::printMesh()
{

	FILE *fp;
	// fp = fopen(file, "w");

	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}
	fseek(fp, 1, SEEK_END);

	int i, e, el;

	const WavesSDGeometry *geometry = &indexes->getGeometry();
	int nnode = geometry->getNoNodes();
	int nsd = geometry->getNoSpaceDim();

	cout << "  nnode = " << nnode << "  nsd = " << nsd << " nno_gg = " << nno_gg << endl;

	//*****  Output of geometry    *************
	//=========================================================
	if (nsd == 3)
	{
		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + nno_gg + 1, geometry->getCoor(i, 0), geometry->getCoor(i, 1), geometry->getCoor(i, 2));
	}
	else if (nsd == 2)
	{
		fprintf(fp, "MESH   dimension  %i  ElemType Squares  Nnode  4\n", nsd);
		fprintf(fp, "Coordinates\n");

		for (i = 0; i < nnode; i++)
			fprintf(fp, "%i %f %f %f\n ", i + nno_gg + 1, geometry->getCoor(i, 0), geometry->getCoor(i, 1), 0.0);
	}

	fprintf(fp, "end coordinates\n");
	fprintf(fp, "\n");

	fclose(fp);

	return true;
}

//========================================================================
// GIDOutputWavesElements
//======================================================================

GIDOutputWavesElements::GIDOutputWavesElements(WavesSDIndexes *sdi, char *file_, int nel_, int nno_) :
		freeIndexes(false), file(file_), indexes(sdi), nel(nel_), nno(nno_)
{
	unary = true;
	interiorWithHole = true;
}

bool GIDOutputWavesElements::printMesh()
{

	FILE *fp;
	fp = fopen(file, "a+");

	if (!fp)
	{
		cout << " cannot open file \"" << file << "\"\n";
		exit(1);
	}

	fseek(fp, 1, SEEK_END);
	int i, e, el;

	const WavesSDGeometry *geometry = &indexes->getGeometry();
	int nnode = geometry->getNoNodes();
	int nsd = geometry->getNoSpaceDim();

	cout << "  nnode = " << nnode << "  nsd = " << nsd << endl;

	//*****  Output of geometry    *************
	//=========================================================

	fprintf(fp, "\n");
	fprintf(fp, "Elements\n");
	int m = 1;
	if (nsd == 2)
	{
		for (el = 0; el < geometry->getNoElms(); el++)
		{
			if (interiorWithHole && indexes->containWavesElement(el))
				m = 1;
			else
				m = 7;
			fprintf(fp, "%i   %i    %i   %i  %i \n", el + nel + 1, geometry->loc2glob(el, 0) + nno + 1, geometry->loc2glob(el, 1) + nno + 1, geometry->loc2glob(el, 3) + nno + 1, geometry->loc2glob(el, 2) + nno + 1);

		}
	}

	fprintf(fp, "end elements\n");
	fclose(fp);

	return true;
}

//======================================================================
//
// Class PlotMTVOutputOp
//
//======================================================================

PlotMTVOutputOp::PlotMTVOutputOp(WavesSDGeometry *sdg, ostream &ostr_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), ostr(ostr_)
{
	unary = true;
	interiorWithHole = false;
}

PlotMTVOutputOp::PlotMTVOutputOp(WavesSDIndexes *sdi, ostream &ostr_) :
		WavesSDOperator(sdi), freeIndexes(false), ostr(ostr_)
{
	unary = true;

	interiorWithHole = true;
}

bool PlotMTVOutputOp::doApply(const double *x, double *y) const
{
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	//*****  Output of geometry    *************

	int n, e, m;
	int nx_ = geometry->getN_i();
	int ny_ = geometry->getN_j();
	int nz_ = geometry->getN_k();
	real xmin_ = geometry->getXmin();
	real ymin_ = geometry->getXmin();
	real zmin_ = geometry->getXmin();

	real xmax_ = geometry->getDx() * (nx_ - 1);
	real ymax_ = geometry->getDy() * (ny_ - 1);
	real zmax_ = geometry->getDz() * (nz_ - 1);

	if (geometry->getNoSpaceDim() == 3)
	{
		ostr << "$ DATA=GRID4D \n" << endl;
		ostr << " %toplabel=\"Exact solution \" \n" << endl;
		ostr << " %meshplot=true \n" << endl;
		ostr << " %contstyle=2 \n" << endl;
		ostr << " %nx= " << nx_ << " [xmin = " << xmin_ << "   xmax = " << xmax_ << "]\n" << endl;
		ostr << " %ny= " << ny_ << " [ymin = " << ymin_ << "   ymax = " << ymax_ << "]\n" << endl;
		ostr << " %nz= " << nz_ << " [zmin = " << zmin_ << "   zmax = " << zmax_ << "]\n" << endl;

		for (n = 0; n < geometry->getNoElms(); n++)
		{
			if (interiorWithHole && indexes->containWavesElement(n))
				m = 1;
			else
				m = 99;

			if (m == 1)
			{
				int n1 = geometry->loc2glob(n, 0);
				int n2 = geometry->loc2glob(n, 1);
				int n3 = geometry->loc2glob(n, 3);
				int n4 = geometry->loc2glob(n, 2);

				int n5 = geometry->loc2glob(n, 4);
				int n6 = geometry->loc2glob(n, 5);
				int n7 = geometry->loc2glob(n, 7);
				int n8 = geometry->loc2glob(n, 6);

				ostr << geometry->node2x(n1) << " " << geometry->node2y(n1) << " " << geometry->node2z(n1) << " " << y[n1] << endl;
				ostr << geometry->node2x(n2) << " " << geometry->node2y(n2) << " " << geometry->node2z(n2) << " " << y[n2] << endl;
				ostr << geometry->node2x(n3) << " " << geometry->node2y(n3) << " " << geometry->node2z(n3) << " " << y[n3] << endl;
				ostr << geometry->node2x(n4) << " " << geometry->node2y(n4) << " " << geometry->node2z(n4) << " " << y[n4] << endl;
				ostr << geometry->node2x(n5) << " " << geometry->node2y(n5) << " " << geometry->node2z(n5) << " " << y[n5] << endl;
				ostr << geometry->node2x(n6) << " " << geometry->node2y(n6) << " " << geometry->node2z(n6) << " " << y[n6] << endl;
				ostr << geometry->node2x(n7) << " " << geometry->node2y(n7) << " " << geometry->node2z(n7) << " " << y[n7] << endl;
				ostr << geometry->node2x(n8) << " " << geometry->node2y(n8) << " " << geometry->node2z(n8) << " " << y[n8] << endl;

				ostr << " " << endl;
			}
		}
	}
	else if (geometry->getNoSpaceDim() == 2)
	{
		ostr << "$ DATA=CONTCURVE " << endl;
		ostr << "%toplabel= \"Exact solution\" " << endl;
		ostr << "%meshplot=true " << endl;
		ostr << "%contstyle=2 " << endl;

		for (n = 0; n < geometry->getNoElms(); n++)
		{
			if (interiorWithHole && indexes->containWavesElement(n))
				m = 1;
			else
				m = 99;

			if (m == 1)
			{
				int n1 = geometry->loc2glob(n, 0);
				int n2 = geometry->loc2glob(n, 1);
				int n3 = geometry->loc2glob(n, 3);
				int n4 = geometry->loc2glob(n, 2);

				ostr << geometry->node2x(n1) << " " << geometry->node2y(n1) << " " << y[n1] << "  " << n1 << endl;
				ostr << geometry->node2x(n2) << " " << geometry->node2y(n2) << " " << y[n2] << "  " << n2 << endl;
				ostr << geometry->node2x(n3) << " " << geometry->node2y(n3) << " " << y[n3] << "  " << n3 << endl;
				ostr << geometry->node2x(n4) << " " << geometry->node2y(n4) << " " << y[n4] << "  " << n4 << endl;

				ostr << " " << endl;
			}
		}
	}
	else
	{
		cerr << "PlotMTVOutputOp, error: PLOTMTV output not defined for 1D." << endl;
		return false;
	}

	ostr << "$ END " << endl;

	return true;
}

int PlotMTVOutputOp::writeMTV2D()
{
	cout << "inside function" << endl;
	const WavesSDGeometry *geometry = &indexes->getGeometry();

	int ierr;
	int i, n, m, j;
	int el;

	ostr << "$ DATA=CURVE2D \n" << endl;
	ostr << " %meshplot=true \n" << endl;
	ostr << " %contstyle=2 \n" << endl;

	for (int el = 0; el < geometry->getNoElms(); el++)
	{
		if (interiorWithHole && indexes->containWavesElement(el))
			m = 1;
		else
			m = 99;

		if (m == 1)
		{
			int n1 = geometry->loc2glob(el, 0);
			int n2 = geometry->loc2glob(el, 1);
			int n3 = geometry->loc2glob(el, 3);
			int n4 = geometry->loc2glob(el, 2);

			ostr << geometry->getCoor(n1, 0) << "  " << geometry->getCoor(n1, 1) << endl;
			ostr << geometry->getCoor(n2, 0) << "  " << geometry->getCoor(n2, 1) << endl;
			ostr << geometry->getCoor(n3, 0) << "  " << geometry->getCoor(n3, 1) << endl;
			ostr << geometry->getCoor(n4, 0) << "  " << geometry->getCoor(n4, 1) << endl;
			ostr << geometry->getCoor(n1, 0) << "  " << geometry->getCoor(n1, 1) << endl;

			ostr << "  " << endl;
		}
	}

	ostr << "$ END " << endl;

	return 0;
}

//======================================================================

DifferenceCheckOp::DifferenceCheckOp(WavesSDIndexes *IY, WavesSDIndexes *IX, real tol_) :
		WavesSDOperator(IY), tol(tol_)
{
	// IY is stored in WavesSDOperator::indexes
	iyarray = new int[indexes->nOfNodeIndex()];
	indexes->makeIndexArray(iyarray);

	if (IX)
	{
		assert(indexes->nOfNodeIndex() == IX->nOfNodeIndex());
		ixarray = new int[indexes->nOfNodeIndex()];
		IX->makeIndexArray(ixarray);
		oneIndex = true;
	}
	else
	{
		oneIndex = false;
		ixarray = iyarray;
	}

}

DifferenceCheckOp::~DifferenceCheckOp()
{
	if (!oneIndex)
		delete[] ixarray;
	delete[] iyarray;
}

bool DifferenceCheckOp::doApply(const double *x, double *y) const
{
	real diff, maxDiff = 0;
	bool diffExist = false;
	bool is3d = indexes->getGeometry().getNoSpaceDim() == 3;

	for (int i = 0; i < indexes->nOfNodeIndex(); i++)
	{
		diff = y[iyarray[i]] - x[ixarray[i]];
		if (fabs(diff) > fabs(maxDiff))
			maxDiff = diff;
		if (fabs(diff) > tol)
		{
			if (!diffExist)
			{
				cout << "\n             x             y";
				if (is3d)
					cout << "             z";
				cout << "          diff" << endl;
			}
			cout << indexes->getGeometry().node2x(iyarray[i]);
			cout << indexes->getGeometry().node2y(iyarray[i]);
			if (is3d)
				cout << indexes->getGeometry().node2z(iyarray[i]);
			cout << diff << endl;
			diffExist = true;
		}
	}
	if (!diffExist)
		cout << "No difference > " << tol << " was encountered! " << endl;
	cout << "Max difference was " << maxDiff << endl;
	return !diffExist;
}

//======================================================================
