#ifndef convexHullCGAL_H
#define convexHullCGAL_H

#include <cmath>
#include <algorithm>
#include <memory>
#include <ctype.h>
#include <random>
#include <stdio.h>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <sys/time.h>
#include <string.h>
#include <stdexcept>
#include <cassert>

#include "vector_types.h"
#include <vector>
#include "gpuarray.h"
#include "indexer.h"

#include "dDimensionalVectorTypes.h"


class convexHullCGALInterface
    {
    public:
        //! compute the convex hull of a set of points on a sphere
        void sphericalConvexHull(dVec *points, int n, GPUArray<int> &allNeighs, GPUArray<unsigned int> &numNeighs, Index2D &nidx);

        //!meant to be used as the vertex model initializer... computes DT, then infers vertex positions and neighbor relations
        void sphericalConvexHullForVertexModel(dVec *cellpoints, int n, GPUArray<int> &allNeighs,
                    GPUArray<unsigned int> &numNeighs, Index2D &nidx, GPUArray<dVec> &vertexPositions, GPUArray<int> &vertexNeighs, 
                    GPUArray<int> &vertexCellNeighs, GPUArray<unsigned int> &numVertexNeighs, Index2D &vnidx);
    };
#endif
