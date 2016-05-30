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
// First added:  2000-01-14 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

/*
 Here are rhs definitions for PDE equations
 */

#include <stdio.h>   
#include <iostream>
#include <string.h>
#include <fstream>
#include "include/wavesFuncDefs.h"

using namespace std;

real rhs1(real x, real y, real z, real t)
{
	real dx = (x - 0.5);
	real dy = (y - 0.5);
	real norm = dx * dx + dy * dy;

	if (norm < 0.0025)
	{
		const real st = sin(40 * M_PI * t) * 1e3;
		const real cc = 0.5 * M_PI / 0.0025;
		return st * cos(cc * norm);
	}
	else
	{
		return 0.0;
	}

}

real rhs2(real x, real y, real z, real t)
{

	real dx = (x - 0.5);
	real dy = (y - 0.5);
	real norm = dx * dx + dy * dy;

	if (norm < 0.0025)
	{
		const real st = sin(40 * M_PI * t) * 1e3;
		return sin(st) * exp(-10 * norm);
	}
	else
	{
		return 0.0;
	}

}

real rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.05 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.0025)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}

real Lens(real x, real y, real z, real t)
{

	return 1e3 * (sin(M_PI * t) * sin(M_PI * t));

}

real dt_rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.05 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.0025)
	{
		return 2 * M_PI * 1e3 * (sin(M_PI * t) * cos(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}

real rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.01 && (x - 0.5) * (x - 0.5) + (y - 0.35) * (y - 0.35) < 0.0025)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}

}

real dt_rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.01 && (x - 0.5) * (x - 0.5) + (y - 0.35) * (y - 0.35) < 0.0025)
	{
		return 2 * M_PI * 1e3 * (sin(M_PI * t) * cos(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}

real t0_1rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.01)
	{
		return (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}

}

real t0_1rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.3) * (y - 0.3) < 0.01)
	{
		return (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}

}
real e_3rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.3) * (y - 0.3) < 0.01)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}
real e_3rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.01)
	{
		cout << " e_3rhs0_5_0_5 = " << 1e3 * (sin(M_PI * t) * sin(M_PI * t)) << "time is " << t << endl;
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//====================================================
// RHS for Maxwell's equations 
//====================================================

real rhsMaxwell(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.8) * (y - 0.8) < 0.01)
	{
		cout << " rhsMaxwell = " << 2 * t * (sin(M_PI * x) * sin(M_PI * y)) << "time is " << t << endl;
		return 2 * t * (sin(M_PI * x) * sin(M_PI * y));
	}
	else
	{
		return 0;
	}
}

real rhsMaxwell1(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return 2 * t * sin(M_PI * y) + (M_PI * M_PI * t * t * t * sin(M_PI * y)) / 3.0;
	}
	else
	{
		return 0;
	}
}

real rhsMaxwell2(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return 2 * t * sin(M_PI * x) + (M_PI * M_PI * t * t * t * sin(M_PI * x)) / 3.0;
	}
	else
	{
		return 0;
	}
}
//=============  rhs for FEM =======================================0

real dt_rhsMaxwell1(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return (2 * sin(M_PI * y) + M_PI * M_PI * t * t * sin(M_PI * y));
	}
	else
	{
		return 0;
	}
}

real dt_rhsMaxwell2(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return (2 * sin(M_PI * x) + M_PI * M_PI * t * t * sin(M_PI * x));
	}
	else
	{
		return 0;
	}
}

//=====================================================================
// exact solution for Maxwell's equations for above RHS (rhsMaxwell1,2)
//======================================================================

real pulseE_2(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return t * t * sin(M_PI * x);
	}
	else
	{
		return 0;
	}
}

real pulseE_1(real x, real y, real z, real t)
{
	if (t < 10.0)
	{
		return t * t * sin(M_PI * y);
	}
	else
	{
		return 0;
	}
}

//=======================================================================

real elast_pulse1(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.01)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real elast_pulse2(real x, real y, real z, real t)
{

	return 0;
}

//=====================================================

real rhs50_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.01)
	{
		return 3 * 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real e_2rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.3) * (y - 0.3) < 0.01)
	{
		return 2.0 * 1e2 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}
real e_2rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.01)
	{
		return 2.0 * 1e2 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0.0;
	}

}
real sejsmrhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.14 && (x - 0.5) * (x - 0.5) + (y - 0.3) * (y - 0.3) < 0.01)
	{
		return (sin(M_PI * t / 0.07) * sin(M_PI * t / 0.07));
	}
	else
	{
		return 0.0;
	}

}

real constplaneWave(real x, real y, real z, real t)
{
	return 1.0;
}

real planeWave(real x, real y, real z, real t)
{
	real t_min = 0.;

	// this is for Maxwell's tests
	//  real t_max = 2*M_PI/5.0;
	//  real t_max = 2*M_PI/10.0;
	real t_max = 2 * M_PI / 25.0;
	//  real t_max = 2*M_PI/150.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	return (0.1 * sin(25 * t - M_PI / 2) + 0.1);
}

real planeWaveOmega(real omega, real y, real z, real t)
{

  //	cout << "time is " << t << " and omega = " << omega << endl;
	real t_min = 0.;
	real t_max = 2 * M_PI / omega;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	return (0.1 * sin(omega * t - M_PI / 2) + 0.1);

}

real planeWaveOmega2(real omega, real y, real z, real t)
{

  //	cout << "time is " << t << " and omega = " << omega << endl;
	real t_min = 0.;
	real t_max = 2 * M_PI / omega;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	if (t < t_max)
		return (sin(omega * t));
	else
		return 0.0;

}

real rhs_maxwell(real x, real y, real z, real t)
{

	// pulse for cube , geometry maxwell_coarse.dat, ... , maxwell_ref5big.dat

	//   if (t <= 0.1 && (x- 2.05)*(x- 2.05)+(y-3.5)*(y-3.5)+(z-1.25)*(z-1.25) < 
	//     0.1)

	// pulse for lens_geometry: lens_cm.dat, ..., lens_ref3type5.dat

	if (t <= 0.1 && x * x + (y - 11.0) * (y - 11.0) + z * z < 3.0)
	{
		return sin(5.0 * t - M_PI / 2);
	}
	else
		return 0;

}
real init_zero(real x, real y, real z, real t)
{

	return 0;

}

real aposterplaneWave(real x, real y, real z, real t)
{
	real t_min = 0.;

	real t_max = 2 * M_PI / 10.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	if (t > t_min && t < t_max)
		return 1.0;
	else
		return 0.0;

}

real planeWavetest(real x, real y, real z, real t)
{
	return (0.1 * sin(100.0 * t - M_PI / 2) + 0.1);
}

// exact solution for plane wave

real exactplaneWave(real x, real y, real z, real t)
{
	real t_min = 0.;
	real t_max = 2 * M_PI / 5.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	return (0.1 * cos(5.0 * t - M_PI / 2)) * (0.1 * cos(5.0 * x - M_PI / 2)) * (0.1 * cos(5.0 * y - M_PI / 2));
}

//========================================================
// plane wave for elastodynamic system
//=========================================================

real ElasticplaneWave1(real x, real y, real z, real t)
{
	real t_min = 0.;
	real t_max = 2 * M_PI / 2.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	return (0.1 * sin(2.0 * t - M_PI / 2) + 0.1);
}

real ElasticplaneWave2(real x, real y, real z, real t)
{

	return 0.0;

}

//=========================================================

real derplaneWave(real x, real y, real z, real t)
{
	real t_min = 0.;
	real t_max = 2 * M_PI / 25.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	return (0.1 * 25.0 * cos(25.0 * t - M_PI / 2));

}

real planeWave2(real x, real y, real z, real t)
{

	if (t <= 0.05 && (x - 0.3) * (x - 0.3) + (y - 0.3) * (y - 0.3) < 0.0025)
	{
		return 1e4 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//========================================================================
real test_exact(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 2.5) * (x - 2.5) + (y - 1.25) * (y - 1.25) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return -M_PI * M_PI * sin(M_PI * t);
	}
	else
	{
		return 0;
	}
}

real exact(real x, real y, real z, real t)
{

	return sin(M_PI * t);

}

//==========================================================================

real D3rhs0_5_0_5(real x, real y, real z, real t)
{

	if (t <= 0.1 && x * x + (y - 5.0) * (y - 5.0) + (z + 1.5) * (z + 1.5) < 0.75)
	{
		return sin(M_PI * t);

	}
	else
	{
		return 0;
	}
}

real delta_func(real x, real y, real z, real t)
{

	if (t <= 0.1 && x * x + (y - 5.0) * (y - 5.0) + z * z < 0.75)
	{
		return sin(M_PI * t);

	}
	else
	{
		return 0;
	}
}

real delta2(real omega, real y, real z, real t)
{

  //	cout << "time is " << t << " and omega = " << omega << endl;
	real t_min = 0.;
	real t_max = 1.0;

	if (t < t_min)
		t = t_min;
	if (t > t_max)
		t = t_max;

	if (t < t_max)
		return 1.0;
	else
		return 0.0;

}

real D3rhs0(real x, real y, real z, real t)
{

	if (t <= 1.0 && x * x + y * y + z * z < 0.5)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real pulse0(real omega, real x, real y, real t)
{
	return (0.1 * sin(omega * t - M_PI / 2));
}

real centered_pulse(real omega, real x, real y, real t)
{
	return (0.05 * sin(omega * t - M_PI / 2));
}

real pulse2(real omega, real x, real y, real t)
{
	return (0.025 * sin(omega * t - M_PI / 2));
}

real pulse3(real omega, real x, real y, real t)
{
	return (0.0125 * sin(omega * t - M_PI / 2));
}

real pulse4(real omega, real x, real y, real t)
{
	return (0.00625 * sin(omega * t - M_PI / 2));
}

real Off_centered_pulse(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 0.0) * (x - 0.0) + (y - 0.0) * (y - 0.0) + (z - 0.0) * (z - 0.0) < 0.2)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real Off_centered_pulse04(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.8) * (x - 0.8) + (y - 0.0) * (y - 0.0) + (z - 0.0) * (z - 0.0) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_two(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 1.25) * (x - 1.25) + (y - 1.25) * (y - 1.25) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 3.75) * (x - 3.75) + (y - 1.25) * (y - 1.25) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//==================================================
real impl_scalar(real x, real y, real z, real t)
{

	if (t <= 0.2 && (x - 3.5) * (x - 3.5) + (y - 1.3) * (y - 1.3) + (z - 1.3) * (z - 1.3) < 0.1)
	{
		cout << " condition works: time is " << t << " and func 3D = " << 1e3 * (sin(M_PI * t) * sin(M_PI * t)) << endl;
		return -1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//==================================================

real rhs_two(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 1.25) * (x - 1.25) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return -1e3 * M_PI * M_PI * sin(M_PI * t);
	}
	else if (t <= 0.1 && (x - 3.75) * (x - 3.75) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return -1e3 * M_PI * M_PI * sin(M_PI * t);
	}
	else
	{
		return 0;
	}
}

//==================================================

real rhs_amira(real x, real y, real z, real t)
{

	if (t <= 1.0 && (x - 12.0) * (x - 12.0) + (y - 9.0) * (y - 9.0) + (z - 10.0) * (z - 10.0) < 1.0)
	{
		return 1e3 * M_PI * M_PI * sin(M_PI * t);
	}
	else
	{
		return 0;
	}

}

//====================================================================

real rhs_four(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 2.05) * (x - 2.05) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real rhs_six(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 0.45) * (x - 0.45) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 1.25) * (x - 1.25) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 2.05) * (x - 2.05) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 2.95) * (x - 2.95) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));

	}
	else if (t <= 0.1 && (x - 3.75) * (x - 3.75) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 4.55) * (x - 4.55) + (y - 2.2) * (y - 2.2) + (z - 1.25) * (z - 1.25) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//====================================================================

real pulse_left_right(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 1.3) * (y - 1.3) + (z - 1.3) * (z - 1.3) < 0.2)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 4.5) * (x - 4.5) + (y - 1.3) * (y - 1.3) + (z - 1.3) * (z - 1.3) < 0.2)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 2.5) * (x - 2.5) + (y - 1.3) * (y - 1.3) + (z - 0.5) * (z - 0.5) < 0.2)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else if (t <= 0.1 && (x - 2.5) * (x - 2.5) + (y - 1.3) * (y - 1.3) + (z - 2.1) * (z - 2.1) < 0.2)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}

}

real pulse_for_sphere(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 15.0) * (x - 15.0) + (y - 20.0) * (y - 20.0) + (z - 5.0) * (z - 5.0) < 1.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real pulse_left_sphere(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 10.0) * (x - 10.0) + (y - 15.0) * (y - 15.0) + (z - 5.0) * (z - 5.0) < 1.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real pulse_front_sphere(real x, real y, real z, real t)
{
	if (t <= 0.1 && (x - 15.0) * (x - 15.0) + (y - 15.0) * (y - 15.0) + (z - 0.0) * (z - 0.0) < 1.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

//=========================================================================

real D3rhs0_1_0_1(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.125) * (x - 0.125) + (y - 0.125) * (y - 0.125) + (z - 0.125) * (z - 0.125) < 0.1)
	{
		return (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_cube(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.14) * (x - 0.14) + (y - 0.09) * (y - 0.09) + (z - 0.1) * (z - 0.1) < 0.05)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_onepoint(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.7) * (y - 0.7) + (z - 0.7) * (z - 0.7) < 0.1)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_for_gid(real x, real y, real z, real t)
{

	if (t <= 2.0 && (x - 13.319100) * (x - 13.319100) + (y - 13.477600) * (y - 13.477600) + (z - 5.251250) * (z - 5.251250) < 10.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_for_gid2(real x, real y, real z, real t)
{

	if (t <= 2.0 && (x - 13.319100) * (x - 13.319100) + (y - 18.0) * (y - 18.0) + (z - 5.251250) * (z - 5.251250) < 10.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs_for_gid3(real x, real y, real z, real t)
{

	if (t <= 2.0 && (x - 13.319100) * (x - 13.319100) + (y - 15.0) * (y - 15.0) + (z - 2.0) * (z - 2.0) < 10.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs0_5_0_5_Poisson(real x, real y, real z)
{

	int a = 400;
	return (4 / M_PI) * a * a * (1 - a * x * x - a * y * y) * exp(-a * (x * x + y * y));

}

real D3rhs_Poisson(real x, real y, real z)
{

	int a = 400;
	double b = 4 / 6;

	return (6 / M_PI) * a * a * (1 - b * a * x * x - b * a * y * y - b * a * z * z) * exp(-a * (x * x + y * y + z * z));

}

real Apel_rhs2(real x, real y, real z)
{

	return 0.0;
}

real exsol_Apel2(real x, real y, real z)
{

	double fi;

	double r = sqrt(x * x + y * y + z * z);

	if (r > 0.0)
		fi = asin(y / r);
	else
		fi = 0.0;

	r = (r * r) * (cbrt(r));

	return r * sin((2.0 / 3.0) * fi);

}

real Apel_rhs(real x, real y, real z)
{

	double r = sqrt(x * x + y * y + z * z);

	double fi = asin(y / r);

	r = (r * r) * (cbrt(r));

	return 0.25 * (1.0 / (z * z * sqrt(fabs(z)))) * r * sin((2.0 / 3.0) * fi);

}

real exsol_Apel(real x, real y, real z)
{

	double fi;

	double r = sqrt(x * x + y * y + z * z);

	if (r > 0.0)
		fi = asin(y / r);
	else
		fi = 0.0;

	r = (r * r) * (cbrt(r));

	cout << " exsol_Apel , fi = " << fi << " r = " << r << "z = " << z << "sin(2/3*fi) " << double(sin((2.0 / 3.0) * fi)) << " sqrt(z)*r*sin(2/3*fi) " << sqrt(z) * r * sin((2.0 / 3.0) * fi) << endl;

	z = fabs(z);
	return sqrt(z) * r * sin((2.0 / 3.0) * fi);

}

real exsol_Poisson(real x, real y, real z)
{

	int a = 400;

	return (a / M_PI) * exp(-a * (x * x + y * y + z * z));

}

real exsol_Poisson2(real x, real y, real z)
{

	double r = sqrt(x * x + y * y + z * z);

	return sqrt(r);

}

real rhs_Poisson2(real x, real y, real z)
{

	double r = sqrt(x * x + y * y + z * z);
	return (-3.0 / 4.0) * (sqrt(r) / (r * r));

}

//====================================================================
// for 7 common cubes 
//===================================================================

real exsol_Poisson3(real x, real y, real z)
{

	double r, r_;
	r = sqrt(x * x + y * y + z * z);

	r_ = sqrt(r);

	if (r > 0.0)
		return 1.0 / sqrt(r_);
	else
	{
		return 0;
	}
}

real rhs_Poisson3(real x, real y, real z)
{

	double r = sqrt(x * x + y * y + z * z);

	double r_ = sqrt(r);

	return (-3.0 / 16.0) * (1.0 / (sqrt(r_) * (r * r)));

}

real u3_x(real x, real y, real z)
{
	double delta_x = 0.01;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt((x + delta_x) * (x + delta_x) + y * y + z * z);

	double r_ = sqrt(r);
	double delta_r_ = sqrt(delta_r);
	double u = 1.0 / sqrt(r_);
	double delta_u = 1.0 / sqrt(delta_r_);
	return (delta_u - u) / delta_x;
}

real u3_y(real x, real y, real z)
{
	double delta_y = 0.01;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + (y + delta_y) * (y + delta_y) + z * z);

	double r_ = sqrt(r);
	double delta_r_ = sqrt(delta_r);
	double u = 1.0 / sqrt(r_);
	double delta_u = 1.0 / sqrt(delta_r_);
	return (delta_u - u) / delta_y;
}

real u3_z(real x, real y, real z)
{
	double delta_z = 0.01;
	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + y * y + (z + delta_z) * (z + delta_z));

	double r_ = sqrt(r);
	double delta_r_ = sqrt(delta_r);
	double u = 1.0 / sqrt(r_);
	double delta_u = 1.0 / sqrt(delta_r_);
	return (delta_u - u) / delta_z;
}

//=====================================================================

real u_x(real x, real y, real z)
{

	double r = sqrt(x * x + y * y + z * z);

	return x / (2.0 * r * sqrt(r));

}

real u_y(real x, real y, real z)
{
	double r = sqrt(x * x + y * y + z * z);

	return y / (2.0 * r * sqrt(r));

}

real u_z(real x, real y, real z)
{
	double r = sqrt(x * x + y * y + z * z);

	return z / (2.0 * r * sqrt(r));

}

//===================================================0
//=====================================================================

real u_x_A(real x, real y, real z)
{
	double delta_x = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt((x + delta_x) * (x + delta_x) + y * y + z * z);

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin(y / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = sqrt(z) * r * sin((2.0 / 3.0) * fi);
	double delta_u = sqrt(z) * delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_x;

}

real u_y_A(real x, real y, real z)
{

	double delta_y = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + (y + delta_y) * (y + delta_y) + z * z);

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin((y + delta_y) / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = sqrt(z) * r * sin((2.0 / 3.0) * fi);
	double delta_u = sqrt(z) * delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_y;

}

real u_z_A(real x, real y, real z)
{

	double delta_z = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + y * y + (z + delta_z) * (z + delta_z));

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin(y / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = sqrt(z) * r * sin((2.0 / 3.0) * fi);
	double delta_u = sqrt(z + delta_z) * delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_z;

}

//==========================

//=====================================================================

real u_x_A2(real x, real y, real z)
{
	double delta_x = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt((x + delta_x) * (x + delta_x) + y * y + z * z);

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin(y / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = r * sin((2.0 / 3.0) * fi);
	double delta_u = delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_x;

}

real u_y_A2(real x, real y, real z)
{

	double delta_y = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + (y + delta_y) * (y + delta_y) + z * z);

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin((y + delta_y) / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = r * sin((2.0 / 3.0) * fi);
	double delta_u = delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_y;

}

real u_z_A2(real x, real y, real z)
{

	double delta_z = 0.01;
	double fi, delta_fi;

	double r = sqrt(x * x + y * y + z * z);
	double delta_r = sqrt(x * x + y * y + (z + delta_z) * (z + delta_z));

	if (r > 0.0)
	{
		fi = asin(y / r);
		delta_fi = asin(y / delta_r);
	}
	else
	{
		fi = 0.0;
		delta_fi = 0.0;
	}

	r = (r * r) * (cbrt(r));
	delta_r = (delta_r * delta_r) * (cbrt(delta_r));

	double u = r * sin((2.0 / 3.0) * fi);
	double delta_u = delta_r * sin((2.0 / 3.0) * delta_fi);

	return (delta_u - u) / delta_z;

}

//==========================

real D3rhs0_5_0_3(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.3) * (y - 0.3) + (z - 0.5) * (z - 0.5) < 0.01)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real D3rhs2(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 2.0) * (x - 2.0) + (y - 2.0) * (y - 2.0) + (z - 2.0) * (z - 2.0) < 1.0)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}
}

real rhs0505(real x, real y, real z, real t)
{

	if (t <= 0.1 && (x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5) < 0.0025)
	{
		return 1e3 * (sin(M_PI * t) * sin(M_PI * t));
	}
	else
	{
		return 0;
	}

}
