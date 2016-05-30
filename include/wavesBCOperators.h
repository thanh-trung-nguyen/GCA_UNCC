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

#ifndef __WAVESBCOPERATORS_H
#define __WAVESBCOPERATORS_H

#include "wavesSDOperator.h"

/*
 Some BC operators are declared here.
 */

/**
 WavesMirrorBC is used as a mirror BC. In 3d, if Face 0 is the boundary and Face 1 is the interior neighbour face, and so on,
 {\tt WavesMirrorBC::apply} results in that Face 0 obtains the values of Face 2.
 The orientation is determined by the {\tt SDBoundaryType}. The options {\tt SDAllBoundary} and
 {\tt SDCorners} are not supported. The operator is unary.
 */

class WavesMirrorBC : public WavesSDOperator
{
	public:
		/// An WavesSDBoundary which does not contain the edges of the face is constructed.
		WavesMirrorBC(WavesSDGeometry *sdg, const SDBoundaryType &sdb);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesMirrorBC(WavesSDBoundary *sdb);

		/// Destructor frees memory if needed.
		~WavesMirrorBC()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

		//
		void doInitialize();
		//
		void initWeights();
		//
		SDBoundaryType sdb;

		//
		bool freeIndexes;
};

class WavesMirrorBC_YEE : public WavesSDOperator
{
	public:
		/// An WavesSDBoundary which does not contain the edges of the face is constructed.
		WavesMirrorBC_YEE(WavesSDGeometry *sdg, const SDBoundaryType &sdb);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesMirrorBC_YEE(WavesSDBoundary *sdb);

		/// Destructor frees memory if needed.
		~WavesMirrorBC_YEE()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

		//
		void doInitialize();
		//
		void initWeights();
		//
		SDBoundaryType sdb;

		//
		bool freeIndexes;
};

class WavesMirrorBC1 : public WavesSDOperator
{
	public:
		/// An WavesSDBoundary which does not contain the edges of the face is constructed.
		WavesMirrorBC1(WavesSDGeometry *sdg, const SDBoundaryType &sdb);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		WavesMirrorBC1(WavesSDBoundary *sdb);

		/// Destructor frees memory if needed.
		~WavesMirrorBC1()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

		//
		void doInitialize();
		//
		void initWeights();
		//
		SDBoundaryType sdb;

		//
		bool freeIndexes;
};

class SmoothingOperator : public WavesSDOperator
{
	public:
		/// An WavesSDBoundary which does not contain the edges of the face is constructed.
		SmoothingOperator(WavesSDGeometry *sdg, const SDBoundaryType &sdb);

		/// Here, {\tt sdb} shall have a face as {\tt SDBoundaryType}.
		SmoothingOperator(WavesSDBoundary *sdb);
		SmoothingOperator(WavesSDGeometry *sdg);
		/// Destructor frees memory if needed.
		~ SmoothingOperator()
		{
			if (freeIndexes)
				delete indexes;
		}

	protected:

		///
		bool doApply(const real *x, real *y) const;

		//

		SDBoundaryType sdb;

		//
		bool freeIndexes;
};

#endif
