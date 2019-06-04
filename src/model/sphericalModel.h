#ifndef sphericalModel_H
#define sphericalModel_H

#include "simpleModel.h"
#include "noiseSource.h"
#include "sphericalDomain.h"

class sphericalModel : public simpleModel
    {
    public:
        sphericalModel(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);
        //!move the degrees of freedom
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);

        void setRadius(scalar _r);

        virtual void setParticlePositionsRandomly(noiseSource &noise);
        virtual void setParticlePositionsBandedRandomly(noiseSource &noise,scalar angularExtent);


        virtual void computeForces(bool zeroOutForces=false);

        virtual void setSoftRepulsion(scalar range = 1.0, scalar stiffness = 1.0);


        noiseSource noise;
        sphericalDomain sphere;

        scalar repulsionRange;
        scalar repulsionStiffness;
    };

#endif
