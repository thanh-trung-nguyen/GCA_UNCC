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

#include "include/wavesMapper.hh"

using namespace std;

real determinant(Array2d<real> jac, int nsd_)
{

	if (nsd_ == 3)
	{
		real tmpDetJac = jac(0, 0) * jac(1, 1) * jac(2, 2) - jac(0, 0) * jac(2, 1) * jac(1, 2) + jac(1, 0) * jac(2, 1) * jac(0, 2) - jac(1, 0) * jac(0, 1) * jac(2, 2) + jac(2, 0) * jac(0, 1) * jac(1, 2) - jac(2, 0) * jac(1, 1) * jac(0, 2);
		return tmpDetJac;
	}
	else if (nsd_ == 2)
	{
		real tmpDetJac = jac(0, 0) * jac(1, 1) - jac(1, 0) * jac(0, 1);
		return tmpDetJac;
	}
	else
	{
		cout << "Error in WavesMapper: not implemented for this NDIM." << endl;
		exit(1);
	}

}

WavesMapper::WavesMapper(int nsd_) :
		nsd(nsd_), jac(nsd_, nsd_), cramer(nsd_, nsd_)
{
}

void WavesMapper::mapIsoparametric(WavesQuadRule& q, WavesElement& e)
{

	int gpt, bf, derivative, coord, node; //counters

	int noGPts = q.getNoPts();
	int noNodes = e.getNoNodes();
	int noBFcns = e.getNoBFcns();

	if (noNodes != noBFcns)
	{
		cout << "Error:: noNodes ~= noBFcns in WavesMapper.mapIsoparametric" << endl;
		exit(1);
	}
	e.setNoGPnts(noGPts);
	WavesBasisFcn* tst = e.ptr2test(); //Get pointer to access the element's basis functions
	WavesBasisFcn* trial = e.ptr2trial(); //Get pointer to access the element's basis functions

	// Copy element's coordinates

	tmpCoords.newsize(noNodes, nsd);
	for (node = 0; node < noNodes; node++)
		for (coord = 0; coord < nsd; coord++)
			tmpCoords(node, coord) = e.getCoord(node, coord);

	tmpValGPts.newsize(noBFcns, nsd + 1);
	MV_Vector<double> c(nsd);
	for (gpt = 0; gpt < noGPts; gpt++)
	{
		jac = 0.0;
		for (coord = 0; coord < nsd; coord++)
		{ // Copy this gauss point's coordinates
			c(coord) = q.getGPntCoord(gpt, coord);
			e.gPntCoord(gpt, coord) = 0.0;
		}
		for (bf = 0; bf < noBFcns; bf++)
		{ // For each basis function
			for (derivative = 0; derivative < nsd + 1; derivative++) // Compute the basis functions values in gauss point
				tmpValGPts(bf, derivative) = e.value(bf, derivative, &c(0)); // 0<=>no 1<=>dx, 2<=>dy, 3<=>dz
			for (coord = 0; coord < nsd; coord++)
			{
				e.gPntCoord(gpt, coord) += tmpValGPts(bf, 0) * tmpCoords(bf, coord); // Compute mapped gauss points
				for (derivative = 1; derivative < nsd + 1; derivative++)
				{ // Compute jacobian of mapping
					jac(coord, derivative - 1) += tmpValGPts(bf, derivative) * tmpCoords(bf, coord);
				}
			}
		}

		// Compute determinant of the jacobian
		real tmpDetJac = determinant(jac, nsd);

		e.detJac(gpt) = tmpDetJac;

		MV_Vector<double> tmp(nsd);
		for (bf = 0; bf < noBFcns; bf++)
		{
			tst->values(gpt, bf) = trial->values(gpt, bf) = tmpValGPts(bf, 0);
			for (derivative = 0; derivative < nsd; derivative++)
				tmp(derivative) = tmpValGPts(bf, derivative + 1);

			for (derivative = 0; derivative < nsd; derivative++)
			{
				cramer = jac;
				for (coord = 0; coord < nsd; coord++)
					cramer(derivative, coord) = tmp(coord);
				real cramerDet = determinant(cramer, nsd);
				trial->grad(gpt, bf, derivative) = tst->grad(gpt, bf, derivative) = cramerDet / tmpDetJac;
			}
		}
	}
}
