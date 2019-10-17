#ifndef sphericalVertexModel_H
#define sphericalVertexModel_H

#include "sphericalModel.h"
#include "convexHullCGAL.h"

//!set up a spherical vertex model. N is the number of CELLS, not vertices
class sphericalVertexModel : public sphericalModel
    {
    public:
        sphericalVertexModel(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);

        virtual void getNeighbors(){};
        convexHullCGALInterface convexHuller;
        Index2D cellNeighborIndex;
        GPUArray<unsigned int> cellNumberOfNeighbors;
        GPUArray<int> cellNeighbors;
    };
#endif
