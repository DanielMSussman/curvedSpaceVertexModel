#include "sphericalVoronoi.h"
/*! \file sphericalVoronoi.cpp" */

sphericalVoronoi::sphericalVoronoi(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
    {
    };

void sphericalVoronoi::convexHull()
    {
    ArrayHandle<dVec> pos(positions);
    convexHuller.sphericalConvexHull(pos.data,N,allNeighs,numNeighs);
    /*
    for (int ii = 0; ii < N; ++ii)
        {
        printf("particle %i (%f,%f,%f) neighbors: ",ii, pos.data[ii][0],pos.data[ii][1],pos.data[ii][2]);
        for (int jj = 0; jj < convexHuller.numNeighs[ii]; ++jj)
            printf("%i ",convexHuller.allNeighs[ii][jj]);
        printf("\n");
        }
    */
    };

