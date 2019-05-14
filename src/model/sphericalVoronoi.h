#ifndef sphericalVoronoi_H
#define sphericalVoronoi_H

#include "convexHullCGAL.h"

#include "simpleModel.h"
#include "noiseSource.h"
#include "indexer.h"


class sphericalVoronoi : public simpleModel
    {
    public:
        sphericalVoronoi(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);
        //!move the degrees of freedom
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);

        virtual void setParticlePositionsRandomly(noiseSource &noise);
        //!update the lists of neighbors
        void convexHull();

        Index2D neighborIndex;
        GPUArray<int> numberOfVoronoiNeighbors;
        GPUArray<int> voronoiNeighbors;

        noiseSource noise;
        convexHullCGALInterface convexHuller;
    };

#endif
