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

#ifndef __WAVESFUNCDEFS_H
#define __WAVESFUNCDEFS_H

#include <math.h>
typedef double real;

real Lens(real x, real y, real z, real t);
real ElasticplaneWave1(real x, real y, real z, real t);
real ElasticplaneWave2(real x, real y, real z, real t);
real pulse_left_right(real x, real y, real z, real t);
real pulse_for_sphere(real x, real y, real z, real t);
real pulse_left_sphere(real x, real y, real z, real t);
real pulse_front_sphere(real x, real y, real z, real t);

real rhs_amira(real x, real y, real z, real t);
real rhs_maxwell(real x, real y, real z, real t);
real init_zero(real x, real y, real z, real t);

real elast_pulse1(real x, real y, real z, real t);
real elast_pulse2(real x, real y, real z, real t);
real exact(real x, real y, real z, real t);
real impl_scalar(real x, real y, real z, real t);

real D3rhs0_5_0_5_Poisson(real x, real y, real z);
real delta_func(real x, real y, real z);

real delta2(real omega, real y, real z, real t);

real Apel_rhs(real x, real y, real z);
real exsol_Apel(real x, real y, real z);
real Apel_rhs2(real x, real y, real z);
real exsol_Apel2(real x, real y, real z);
real D3rhs_for_gid(real x, real y, real z, real t);
real test_exact(real x, real y, real z, real t);

real D3rhs_two(real x, real y, real z, real t);
real D3rhs_onepoint(real x, real y, real z, real t);
real rhs_two(real x, real y, real z, real t);
real rhs_four(real x, real y, real z, real t);
real rhs_six(real x, real y, real z, real t);

real D3rhs_cube(real x, real y, real z, real t);
real D3rhs0_1_0_1(real x, real y, real z, real t);

real D3rhs_for_gid1(real x, real y, real z, real t);
real D3rhs_for_gid2(real x, real y, real z, real t);
real D3rhs_for_gid3(real x, real y, real z, real t);

real exsol_Poisson2(real x, real y, real z);
real rhs_Poisson2(real x, real y, real z);
real u_x(real x, real y, real z);
real u_y(real x, real y, real z);
real u_z(real x, real y, real z);
//================================================
real u_x_A(real x, real y, real z);
real u_y_A(real x, real y, real z);
real u_z_A(real x, real y, real z);
real u_x_A2(real x, real y, real z);
real u_y_A2(real x, real y, real z);
real u_z_A2(real x, real y, real z);
//=================================================
real exsol_Poisson3(real x, real y, real z);
real rhs_Poisson3(real x, real y, real z);
real u3_x(real x, real y, real z);
real u3_y(real x, real y, real z);
real u3_z(real x, real y, real z);

//=================================================

real exsol_Poisson(real x, real y, real z);

real D3rhs_Poisson(real x, real y, real z);

real rhs1(real x, real y, real z, real t);

real rhs2(real x, real y, real z, real t);

real rhs0_5_0_5(real x, real y, real z, real t);

real dt_rhs0_5_0_5(real x, real y, real z, real t);

real rhs0_5_0_3(real x, real y, real z, real t);

real dt_rhs0_5_0_3(real x, real y, real z, real t);

real t0_1rhs0_5_0_5(real x, real y, real z, real t);

real t0_1rhs0_5_0_3(real x, real y, real z, real t);

real e_3rhs0_5_0_3(real x, real y, real z, real t);

real e_3rhs0_5_0_5(real x, real y, real z, real t);

real rhsMaxwell(real x, real y, real z, real t);
real rhsMaxwell1(real x, real y, real z, real t);
real rhsMaxwell2(real x, real y, real z, real t);
real dt_rhsMaxwell1(real x, real y, real z, real t);
real dt_rhsMaxwell2(real x, real y, real z, real t);

real pulseE_1(real x, real y, real z, real t);
real pulseE_2(real x, real y, real z, real t);

real rhs50_5_0_5(real x, real y, real z, real t);

real e_2rhs0_5_0_3(real x, real y, real z, real t);

real e_2rhs0_5_0_5(real x, real y, real z, real t);

real sejsmrhs0_5_0_3(real x, real y, real z, real t);

real planeWave(real x, real y, real z, real t);
real planeWaveOmega(real omega, real y, real z, real t);
real planeWaveOmega2(real omega, real y, real z, real t);

real aposterplaneWave(real x, real y, real z, real t);

real planeWavetest(real x, real y, real z, real t);

real exactplaneWave(real x, real y, real z, real t);

real constplaneWave(real x, real y, real z, real t);

real derplaneWave(real x, real y, real z, real t);

real planeWave2(real x, real y, real z, real t);

real D3rhs0(real x, real y, real z, real t);

real centered_pulse(real omega, real x, real y, real t);
real pulse0(real omega, real x, real y, real t);
real pulse2(real omega, real x, real y, real t);
real pulse3(real omega, real x, real y, real t);
real pulse4(real omega, real x, real y, real t);

real Off_centered_pulse(real x, real y, real z, real t);

real Off_centered_pulse04(real x, real y, real z, real t);

real D3rhs0_5_0_5(real x, real y, real z, real t);

real D3rhs0_5_0_3(real x, real y, real z, real t);

real D3rhs2(real x, real y, real z, real t);

real rhs0505(real x, real y, real z, real t);

#endif
