

#ifndef _WAVESFINDBOUNDARYNODES_
#define _WAVESFINDBOUNDARYNODES_
#include <math.h>
#include "waveskraftwerk.hh"


void findBoundaryNodes(WavesGridB& grid, MV_Vector<int>& boundaryMarkers);
void findBoundaryNodes(WavesGridB& grid, MV_Vector<int>& boundaryMarkers,
		       MV_ColMat<int>& elmneigh);
 
#endif 
