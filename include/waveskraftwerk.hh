//KRAFTWERK

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

#ifndef __WAVESKRAFTWERK_H
#define __WAVESKRAFTWERK_H

#include "wavesMVvtp.h"
#include "wavesMVvind.h"
#include "wavesMVmtp.h"

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>

#include "petscvec.h"
#include "petscmat.h"
#include "petscoptions.h"
#include "petscsys.h"

#include "wavesMultiDimArrays.hh"
#include "wavesGridB.h"
#include "wavesNeighborFE.h"

typedef WavesGridB Grid;
typedef double real;
typedef MV_Vector<real> Array1dReal;
typedef MV_Vector<int> Array1dInt;

#include "wavesQuadRule.hh"
#include "wavesBasisFcn.hh"
#include "wavesElement.hh"
#include "wavesEquation.hh"
#include "wavesMapper.hh"
#endif
