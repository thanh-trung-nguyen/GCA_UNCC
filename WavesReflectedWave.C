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

#include "include/wavesReflectedWave.h"

WavesReflectedWave::WavesReflectedWave(WavesSDIndexes *sdi, char* fname_, double* array_) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
	fname = fname_;
	output_ar = array_;
}

WavesReflectedWave::WavesReflectedWave(WavesSDIndexes *sdi, MV_Vector<double>& array_) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;

	Data_array = array_;
}

WavesReflectedWave::WavesReflectedWave(WavesSDIndexes *sdi, char* fname_) :
		WavesSDOperator(sdi), freeIndexes(false)
{
	unary = true;
	fname = fname_;
}

void WavesReflectedWave::write_reflected_field(double* array)
{

	ofstream outp;
	outp.open(fname);

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			outp << array[n] << "   ";
		}
		outp << "\n";

	}
	outp.close();

	printf("Solution at one point is written to file: %s\n", fname);

}

void WavesReflectedWave::read_reflected_field(double& t)
{
	double t_max = 3.00 / 100;
	ifstream inp;
	inp.open(fname);
	printf("Open the file: %s\n", fname);

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			if (t < t_max)
			{
				inp >> output_ar[n];
			}
			else
				output_ar[n] = 0.0;
		}

	}
	inp.close();

}

void WavesReflectedWave::read_reflected_(double& t)
{

	ifstream inp;
	inp.open(fname);
	printf("Open the file: %s\n", fname);

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			inp >> output_ar[n];
		}
	}
	inp.close();
}

bool WavesReflectedWave::doApply(const double *, double *y) const
{

	for (int loop = 0; loop < indexes->nOfLoopIndex(); loop++)
	{
		int lstart = indexes->loopStart(loop);
		int lstop = indexes->loopStop(loop);

		for (int n = lstart; n <= lstop; n++)
		{
			y[n] = Data_array(n);

			if (y[n] != 0.0)
				cout << "y ( " << n << ") =" << y[n] << " ar (" << n << ") = " << output_ar[n] << endl;
		}
	}

	return true;
}

bool WavesReflectedWave::applyReflectedWave(real *y, const real &t)
{

	time = t;
	read_reflected_field(time);

	return doApply(0, y);
}

bool WavesReflectedWave::applyReflectedWaveLow(real *y)
{

	return doApply(0, y);
}

bool WavesReflectedWave::applyReflectedWaveLeft(real *y, const real &t)
{

	time = t;
	read_reflected_(time);

	return doApply(0, y);
}
