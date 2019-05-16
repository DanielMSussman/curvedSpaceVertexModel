#ifndef sphericalVoronoi_H
#define sphericalVoronoi_H

#include "convexHullCGAL.h"

#include "sphericalModel.h"
#include "noiseSource.h"
#include "indexer.h"
#include "sphericalDomain.h"


class sphericalVoronoi : public sphericalModel
    {
    public:
        sphericalVoronoi(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);

        virtual void getNeighbors(){convexHull();};
        //!update the lists of neighbors
        void convexHull();

        convexHullCGALInterface convexHuller;
    };

#endif

