//KRISTER

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

#ifndef __WAVESSDDEFS_H
#define __WAVESSDDEFS_H

/**@name StandardDefinitions
 The ABCD classes have following standard definitions.

 */

#define true  1
#define false 0

/** real is double.

 Note that real* usually is a scalar field over a grid where the components are ordered in natural ordering. WavesSDOperators can act on SDFields, over an WavesSDIndexes domain on an WavesSDGeometry. */
typedef double real;

typedef real (*Fcn3dtime)(real x, real y, real z, real t);

#endif
