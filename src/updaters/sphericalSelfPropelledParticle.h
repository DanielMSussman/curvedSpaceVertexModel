#ifndef sphericalSPP_H
#define sphericalSPP_H

#include "vectorialVicsek.h" //use vector vicsek for noise and parameters contained in that class
#include "sphericalModel.h"
/*! \file sphericalSelfPropelledParticle.h */

/*!
mu = inverse friction
v0 = self-propulsion
tau = D_r^{-1}
*/
class sphericalSelfPropelledParticle : public vectorialVicsek
    {
    public:
        sphericalSelfPropelledParticle(){useGPU = false; mu = 1.0; tau = 1.0;v0=0.01;};
        virtual void integrateEOMGPU(){};
        virtual void integrateEOMCPU();

        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model);
        shared_ptr<sphericalModel> sphereModel;
    };
#endif

