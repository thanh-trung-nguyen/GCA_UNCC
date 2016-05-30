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

/*
 WavesSDIndexes represents indexes for a subdomain in a WavesSDGeometry

 */

#ifndef __WAVESSDINDEXES_H
#define __WAVESSDINDEXES_H

#include <iostream>
#include <assert.h>

#include "wavesSDDefs.h"
#include "wavesSDGeometry.h"

using namespace std;

/**@name LoopIndex
 
 LoopIndex is an abstraction of indices for a loop in i-direction.
 The loop goes from iStart to iStop at j and k 

 */
struct LoopIndex
{
		///
		LoopIndex()
		{
			iStart = iStop = j = k = 0;
		}
		///
		LoopIndex(const int& iStart_, const int& iStop_, const int& j_, const int& k_) :
				iStart(iStart_), iStop(iStop_), j(j_), k(k_)
		{
		}
		/// 
		int iStart;
		///
		int iStop;
		///
		int j;
		///
		int k;
};

/// Simple output of LoopIndex
ostream &operator<<(ostream &os, const LoopIndex &li);

/// A typedef for an array of LoopIndexes
typedef LoopIndex* LoopIndexArray;

/** @name SDBoundaryType
 *  SDBoundaryType is used to mark which side of a cube an WavesSDBoundary 
 *  represents. It is also possible to create indices for all boundary
 *  and for the corners. In 3d, "corners" also include edges. Here are the options:
 */
enum SDBoundaryType
{
		SDallBoundary,
		SDiLow,
		SDiHigh,
		SDjLow,
		SDjHigh,
		SDkLow,
		SDkHigh,
		SDCorners
};

/** 
 * WavesSDIndexes represent any index domain in a structured grid.
 */
class WavesSDIndexes
{
	public:
		/// Copy constructor
		WavesSDIndexes(const WavesSDIndexes &sdi);

		WavesSDIndexes();
		/// Destructor
		virtual ~WavesSDIndexes();
		/// Assignment 
		WavesSDIndexes &operator=(const WavesSDIndexes &sdi);

		/// Returns number of loops
		int nOfLoopIndex() const
		{
			return nOfLoops;
		}

		/// Returns "global" node number for where a loop starts:
		virtual int loopStart(const int &loop) const;

		/// Returns "global" node number for where a loop stops:
		virtual int loopStop(const int &loop) const;

		/// Returns a particular loop
		virtual LoopIndex getLoop(const int &loop) const
		{
			assert(loop >= 0 && loop < nOfLoops);
			return loopIndexes[loop];
		}

		/// Present misc. info
		void presentWavesSDIndexes() const;

		/// Present all global node indices
		void presentGlobalIndexes() const;

		/// Returns number of nodes
		int nOfNodeIndex() const
		{
			return nOfNodes;
		}

		/** Make an array with indices of global node numbers. The array
		 indexArray shall be allocated with space for nOfNodeIndex elements */
		void makeIndexArray(int *indexArray) const;
		void AssignIndexArray(int *indexArray, double* Array, double *AssignArray) const;
		void WriteIndexArray(char *fname, double* Array);

		///
		const WavesSDGeometry &getGeometry() const
		{
			return sdg;
		}

		/// Returns true if all nodes of an element is contained in the index domain.
		virtual bool containWavesElement(int e);
		///  Returns true if a node is contained in the index domain.
		// This general routine performs line search. Inheritors can often speed it up.
		virtual bool containNode(int n);

		void copyLI(WavesSDIndexes& sdi, int nOfLoops_);
		// the same function,as copyLI, but only for inner of the domain,
		// don't include boundaries
		void copyinnerLI(WavesSDIndexes& sdi, int nOfLoops_);

		/// This constructor produces an empty index domain connected to a SDgeometry
		WavesSDIndexes(const WavesSDGeometry &sdg);

		/// The geometry is contained by value in the SDindex object.
		WavesSDGeometry sdg;

	protected:
		// This array may hold indexes for all nodes that shall be updated
		LoopIndexArray loopIndexes;
		// Total number of loops
		int nOfLoops;
		// Current number of loops, used when creating.
		int curLoop;
		// Total number of nodes.
		int nOfNodes;

		// Protected methods
		void doAddLoop(const int &iStart_, const int &iStop_, const int &j_, const int &k_);

		// 
		void copyLoopIndexes(const WavesSDIndexes &sdi);

		// Allocate loops and insert
		void insert(const int &nodeLow, const int &nodeHigh);

		// Only insert. Memory must be allocated beforehand.
		void insertOnly(const int &nodeLow, const int &nodeHigh);

		// If LoopIndexArray is allocated for nOfLoops elements, shrinkLoopIndexArray
		// resizes the array to get curLoop elements, and sets nOfLoops to curLoop.
		void schrinkLoopIndexArray();

		// If LoopIndexArray is allocated for nOfLoops elements, resize(int n)
		// resizes the array from n elements to get an array for n  elements 
		LoopIndexArray tmp;
		void resize(int n);
};

/**WavesSDInterior covers the interior of a domain.

 WavesSDInterior produces LoopIndexes for the interior of a domain, 
 described either as WavesSDGeometry or as two nodes on an WavesSDGeometry.
 */
class WavesSDInterior : public WavesSDIndexes
{
	public:
		/** Produce LoopIndexes covering all of sdg except width.
		 The default width=1 means that the boundary is excluded.
		 If the boundary is to be included, use width 0. */
		WavesSDInterior(const WavesSDGeometry &sdg, const int &width = 1);
		WavesSDInterior();
		/** Produce LoopIndexes covering all the interior of a cube specified by two
		 global node indexes. These nodes are included. Restrictions:
		 \begin{verbatim}
		 0 <= i(nodeHigh)-i(nodeLow), 
		 0 <= j(nodeHigh)-j(nodeLow), 
		 0 <= k(nodeHigh)-k(nodeLow); 
		 \end{verbatim}
		 */
		WavesSDInterior(const WavesSDGeometry &sdg, const int &nodeLow, const int &nodeHigh);
};

/**WavesSDInterior2 covers the interior of a domain.

 WavesSDInterior2 implements WavesSDInterior with a different interface.
 */

class WavesSDInterior2 : public WavesSDIndexes
{
	public:
		/** Produce LoopIndexes covering all of sdg except width.
		 The default width=1 means that the boundary is excluded.
		 If the boundary is to be included, use width 0. */
		WavesSDInterior2(const WavesSDGeometry &sdg, const int &width = 1);
		WavesSDInterior2();
		/** Produce LoopIndexes covering all the interior of a cube specified by two
		 global node indexes. These nodes are included. Restrictions:
		 \begin{verbatim}
		 0 <= i(nodeHigh)-i(nodeLow), 
		 0 <= j(nodeHigh)-j(nodeLow), 
		 0 <= k(nodeHigh)-k(nodeLow); 
		 \end{verbatim}
		 */
		WavesSDInterior2(const WavesSDGeometry &sdg, const int &nodeLow, const int &nodeHigh);

		/// Returns "global" node number for where a loop starts: (Overrides base class version)
		int loopStart(const int &loop) const;

		/// Returns "global" node number for where a loop stops: (Overrides base class version)
		int loopStop(const int &loop) const;

		// Returns a particular loop. (Overrides base class version)
		LoopIndex getLoop(const int &loop) const;

	private:
		/// Utility
		void init(const int &nodeLow, const int &nodeHigh);
		/// Members needed for calculating loopStart and loopStop
		int iStart, iStop, jStart, jStop, kStart, kStop, jLength;
};

/** WavesSDBoundary covers the boundary depending on the enum variable.   

 For a face, the parameter doIncludeCorners decides whether the boundary shall be included or not. For the first constructor, the default option is not to include corners and edges. For the second constructor, the default choice is to include them. 
 @see SDBoundaryType.
 */
class WavesSDBoundary : public WavesSDIndexes
{
	public:

		WavesSDBoundary();
		/// By default, this constructor *does not* include nodeLow and nodeHigh for a face.
		WavesSDBoundary(const WavesSDGeometry &sdg, const SDBoundaryType &sdb = SDallBoundary, const bool &doIncludeCorners = false);

		WavesSDBoundary(const WavesSDGeometry &sdg, int code, const SDBoundaryType &sdb = SDallBoundary, const bool &doIncludeCorners = false);

		/// By default, this constructor *does* include nodeLow and nodeHigh for a face.
		WavesSDBoundary(const WavesSDGeometry &sdg, const int &nodeLow, const int &nodeHigh, const SDBoundaryType &sdb = SDallBoundary, const bool &doIncludeCorners = true);

		/// Return the SDBoundaryType.
		SDBoundaryType getSDBoundaryType()
		{

			return sdb;

		}
		void initProblem(const WavesSDGeometry &sdg, const SDBoundaryType &sdb = SDallBoundary, const bool &doIncludeCorners = false);

		/// Utility method
		void initialize2d(const int &nodeLow, const int &nodeHigh, const bool &doIncludeCorners);
		/// Utility method
		void initialize3d(const int &nodeLow, const int &nodeHigh, const bool &doIncludeCorners);

		///
		SDBoundaryType sdb;
};

/**
 WavesSDInteriorWithHole produces LoopIndexes for the interior of a grid, 
 where we have a hole.

 The domain is specified by four indexes. The result can be described approximately as
 \begin{verbatim}
 WavesSDInterior(sdg,nodeLowOuter,nodeHighOuter) - 
 WavesSDInterior(sdg,nodeLowInner,nodeHighInner).
 \end{verbatim}
 Thus, nodeLowOuter and nodeLowOuter are included in the domain, whereas
 nodeLowInner and nodeHighInner are not.

 The inner domain must be strictly contained in the outer domain.

 */
class WavesSDInteriorWithHole : public WavesSDIndexes
{
	public:
		///
		WavesSDInteriorWithHole(const WavesSDGeometry &sdg, const int &nodeLowOuter, const int &nodeHighOuter, const int &nodeLowInner, const int &nodeHighInner);

		WavesSDInteriorWithHole();
		///
		int nodeInnerLow()
		{
			return nodeLowInner;
		}
		///
		int nodeInnerHigh()
		{
			return nodeHighInner;
		}

		/// containNode is more efficiently implemented than in base class
		bool containNode(int n);

	private:
		/// Not implemented
		WavesSDInteriorWithHole &operator=(const WavesSDInteriorWithHole &);
		///
		WavesSDInteriorWithHole(const WavesSDInteriorWithHole &sdi);
		/// 
		int nodeLowInner, nodeHighInner;
};

/**
 WavesSDIndexLine produces indexes for a line in a 2d or 3d grid.
 The line must be parallel to one axis. Both end points are included.
 The line may be only one point.

 */

class WavesSDIndexLine : public WavesSDIndexes
{
	public:
		/// 
		WavesSDIndexLine();

		WavesSDIndexLine(const WavesSDGeometry &sdg, const int &nodeLow, const int &nodeHigh);
};

//--class WavesSDIndexes ------

/*

 Example: In 2d, we represent the following schematic 6x5grid:

 j\i: 012345
 0: xxxxxx
 1: xx  xx
 2: x  xxx
 3:   xx
 4: xxx

 with loopIndices (k == 0)
 [0:5;0]
 [0:1;1]
 [1:1;2]
 [3:5;2]
 [2:3;3]
 [1:2;4]

 These indexes are the nodes where we use FD to apply a stencil on a vector. If the vector has one unknown, the FD algorithm is as follows:

 for all loopIndices:
 for(n = loopIndex.start(); n <= loopIndex.end(); n++)
 v_new[n] = Stencil applied on v_old[n] 
 end
 end

 If they are used to assemble a matrix, the algorithm is similar:

 for all loopIndices:
 for(n = loopIndex.start(); n <= loopIndex.end(); n++)
 A[n,:] filled according to Stencil 
 end
 end

 Remark: If we rotate the grid, the loop indexes change. One task is to find the optimal orientation.


 */

#endif
