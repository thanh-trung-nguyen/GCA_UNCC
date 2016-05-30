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

#include "include/wavesEquation.hh"

WavesEquation::WavesEquation(int noeq, int nsd_) :
		nsd(nsd_), noEqns(noeq), work(noeq), testFcns(noeq), trialFcns(noeq), coord(nsd_)
{
}

void WavesEquation::integrateJacobian(WavesQuadRule& q, // Input:  Quadrature to use 
		WavesElement& e, //         Which element to integrate over
		Array1dReal& result, // Output: Where to put the element matrix
		Array1dInt& rowIdx, //         Row index in the stiffness  matrix
		Array1dInt& colIdx, //         Column index in the stiffness  matrix
		int& noElmsInElMat //         Total size of element matrix
		)
{

	int tstbf, trialbf, i_eq, j_eq, gpts; //counters
	int noBf = e.getNoBFcns();
	int el = e.getElmNo();
	int noGpts = q.getNoPts(); // ?? e.noGPtns() ??

	myElm = &e;
	result.newsize(noBf * noBf * noEqns * noEqns);
	rowIdx.newsize(noBf);
	colIdx.newsize(noBf);

	workElMat.newsize(noBf, noBf, noEqns, noEqns);
	workElMat = 0.0;
	for (gpts = 0; gpts < noGpts; gpts++)
	{
		real tmp = e.detJac(gpts) * q.weight(gpts);
		for (tstbf = 0; tstbf < noBf; tstbf++)
		{
			for (trialbf = 0; trialbf < noBf; trialbf++)
			{
				e.setIdx(gpts, trialbf, tstbf);
				for (j_eq = 0; j_eq < noEqns; j_eq++)
				{
					setWavesBasisFcns(e, j_eq);
					for (int i = 0; i < nsd; i++)
						coord(i) = e.gPntCoord(gpts, i);
					jacobian();
					for (i_eq = 0; i_eq < noEqns; i_eq++)
					{
						workElMat(tstbf, trialbf, i_eq, j_eq) += tmp * work[i_eq];
					}
				}
			}
		}
	}

	int cnt1 = 0;
	int cnt2 = 0;
	for (tstbf = 0; tstbf < noBf; tstbf++)
	{
		for (i_eq = 0; i_eq < noEqns; i_eq++)
		{
			for (trialbf = 0; trialbf < noBf; trialbf++)
			{
				for (j_eq = 0; j_eq < noEqns; j_eq++)
				{
					result(cnt1++) = workElMat(tstbf, trialbf, i_eq, j_eq);
				}
			}
		}
		rowIdx(cnt2) = e.loc2glob(tstbf);
		colIdx(cnt2) = e.loc2glob(tstbf);
		cnt2++;
	}
	noElmsInElMat = cnt2;

}

void WavesEquation::integrateResidual(WavesQuadRule& q, // Input:  Quadrature to use  
		WavesElement& e, //         Which element to integrate over
		Array1dReal& result, // Output: Where to put the element matrix
		Array1dInt& rowIdx, //         Row index in the stiffness  matrix 
		int& noElmsInElMat //         Total size of element matrix
		)
{

	int tstbf, i_eq, i, gpts; //counters
	int noBf = e.getNoBFcns();
	int noGpts = q.getNoPts(); // ?? e.noGPtns() ??

	myElm = &e;
	result.newsize(noBf * noEqns);
	rowIdx.newsize(noBf * noEqns); //Icke-blockad variant!!!

	workElVec.newsize(noBf, noEqns);
	workElVec = 0.0;
	setWavesBasisFcns(e);
	for (gpts = 0; gpts < noGpts; gpts++)
	{
		real tmp = e.detJac(gpts) * q.weight(gpts);
		for (tstbf = 0; tstbf < noBf; tstbf++)
		{
			e.setIdx(gpts, tstbf);
			for (int i = 0; i < nsd; i++)
				coord(i) = e.gPntCoord(gpts, i);
			residual();
			for (i_eq = 0; i_eq < noEqns; i_eq++)
			{
				workElVec(tstbf, i_eq) += tmp * work(i_eq);
			}
		}
	}

	int cnt1 = 0;
	int cnt2 = 0;
	for (tstbf = 0; tstbf < noBf; tstbf++)
	{
		for (i_eq = 0; i_eq < noEqns; i_eq++)
		{
			result(cnt1++) = workElVec(tstbf, i_eq);
			rowIdx(cnt2++) = i_eq + noEqns * e.loc2glob(tstbf);
		}
	}
	noElmsInElMat = cnt2;

}

void WavesEquation::setWavesBasisFcns(WavesElement& e, int eq)
{ // Set WavesEquation's basisfcns to point to WavesElement's
	for (int i = 0; i < noEqns; i++)
		trialFcns(i) = e.ptr2zeroBfcn();
	testFcns(0) = e.ptr2test();
	trialFcns(eq) = e.ptr2trial();
}

void WavesEquation::setWavesBasisFcns(WavesElement& e)
{ // Set WavesEquation's basisfcns to point to WavesElement's
	testFcns(0) = e.ptr2test();
}

