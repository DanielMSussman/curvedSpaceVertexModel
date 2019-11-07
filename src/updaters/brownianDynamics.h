#ifndef brownianDynamics_H
#define brownianDynamics_H

#include "equationOfMotion.h"
#include "simpleModel.h"
#include "noiseSource.h"

class brownianDynamics : public equationOfMotion
    {
    public:
        brownianDynamics(bool _reproducible = true);
        virtual void integrateEOMGPU();
        virtual void integrateEOMCPU();

        void setT(scalar _t);

        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            initializeFromModel();
            noise.setReproducible(reproducible);
            noise.initialize(Ndof);
            if(model->useGPU)
                noise.initializeGPURNGs();
            };

        scalar temperature;
        scalar mu;
        noiseSource noise;
    };
#endif
