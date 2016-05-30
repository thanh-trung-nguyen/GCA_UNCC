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

#ifndef __ELEMENTTYPE_H
#define __ELEMENTTYPE_H

enum ElementType
{
		NO_ELEMENT = 0, // default dummy for not implemented element
		ELMTRI1,
		ELMQUAD1,
		ELMTET1,
		ELMPYR1,
		ELMPRISM1,
		ELMHEX1,
		ELMLINE2,
		ELMLINE1,
		ELMTRI2,
		ELMTET2
};

#endif
