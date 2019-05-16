#ifndef sphericalModel_H
#define sphericalModel_H

#include "simpleModel.h"
#include "noiseSource.h"
#include "indexer.h"
#include "sphericalDomain.h"


class sphericalModel : public simpleModel
    {
    public:
        sphericalModel(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);
        //!move the degrees of freedom
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);

        void setRadius(scalar _r);

        virtual void setParticlePositionsRandomly(noiseSource &noise);

        virtual void getNeighbors(){};
        //!update the lists of neighbors
        void convexHull();
        //!return a reference to the GPUArray of positions
        virtual GPUArray<dVec> & returnDirectors(){return directors;};

        Index2D neighborIndex;
        GPUArray<int> numberOfNeighbors;
        GPUArray<int> neighbors;
        GPUArray<dVec> directors;
        std::vector< std::vector<int> > allNeighs; //!<The list of neighbors of every point in the convex hull
        std::vector<int> numNeighs;

        noiseSource noise;
        sphericalDomain sphere;
    };

#endif
