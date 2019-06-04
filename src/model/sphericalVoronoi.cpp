#include "sphericalVoronoi.h"
/*! \file sphericalVoronoi.cpp" */

sphericalVoronoi::sphericalVoronoi(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
    {
    };

void sphericalVoronoi::convexHull()
    {
    ArrayHandle<dVec> pos(positions);
    convexHuller.sphericalConvexHull(pos.data,N,neighbors,numberOfNeighbors,neighborIndex);
    };

