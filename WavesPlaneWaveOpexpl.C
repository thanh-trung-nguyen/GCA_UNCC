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
// First added:  2000-02-23 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#include "include/wavesPlaneWaveOpexpl.h"

PlaneWaveOpexpl::PlaneWaveOpexpl(WavesSDIndexes *sdi, Fcn3dtime f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

PlaneWaveOpexpl::PlaneWaveOpexpl(WavesSDGeometry *sdg, Fcn3dtime f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;

}

bool PlaneWaveOpexpl::doApply(const double *, double *y) const
{

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{

		const WavesSDGeometry *g = &indexes->getGeometry();

		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = f(g->node2x(n), g->node2y(n), g->node2z(n), time);
		}
	}
	return true;
}

//======================================================================
//
// Class PlaneWaveOpexpl2
//
//======================================================================

PlaneWaveOpexpl2::PlaneWaveOpexpl2(WavesSDIndexes *sdi, Fcn3dtime f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

PlaneWaveOpexpl2::PlaneWaveOpexpl2(WavesSDGeometry *sdg, Fcn3dtime f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;

}

bool PlaneWaveOpexpl2::doApply(const double *, double *y) const
{

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = f(lambda_, mu_, ro_, time);
		}
	}
	return true;
}

//================================================

PlaneWaveOpexpl3::PlaneWaveOpexpl3(WavesSDIndexes *sdi, real *Difsol_) :
		WavesSDOperator(sdi), freeIndexes(false), Difsol(Difsol_)
{
	unary = true;

}

PlaneWaveOpexpl3::PlaneWaveOpexpl3(WavesSDGeometry *sdg, real *Difsol_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), Difsol(Difsol_)
{
	unary = true;

}

bool PlaneWaveOpexpl3::doApply(const double *, double *y) const
{

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{

		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = Difsol[n];
		}
	}
	return true;
}

//======================================================================
//
// Class PlaneWaveOpPhotonexpl
//
//======================================================================

PlaneWaveOpPhotonexpl::PlaneWaveOpPhotonexpl(WavesSDIndexes *sdi, Fcn3dtime f_) :
		WavesSDOperator(sdi), freeIndexes(false), f(f_)
{
	unary = true;
}

PlaneWaveOpPhotonexpl::PlaneWaveOpPhotonexpl(WavesSDGeometry *sdg, Fcn3dtime f_) :
		WavesSDOperator(new WavesSDInterior(*sdg, 0)), freeIndexes(true), f(f_)
{
	unary = true;

}

bool PlaneWaveOpPhotonexpl::doApply(const double *, double *y) const
{

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = f(omega_, 0.0, 0.0, time);
			//	cout << " initialized gaussian y[" << n << "]=" << y[n] << endl;
		       	

		}
	}
	return true;
}
