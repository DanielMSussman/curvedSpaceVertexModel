#ifndef vectorialVicsek_H
#define vectorialVicsek_H

#include "equationOfMotion.h"
#include "sphericalModel.h"
#include "noiseSource.h"
/*! \file vectorialVicsek.h */

class vectorialVicsek : public equationOfMotion
    {
    public:
        vectorialVicsek(){useGPU = false; mu = 1.0; Eta = 1.0; tau = 1.0;v0=0.01;};
        virtual void integrateEOMGPU(){};
        virtual void integrateEOMCPU();

        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            voronoiModel = dynamic_pointer_cast<sphericalModel>(model);
            initializeFromModel();
            newVelocityDirector.resize(Ndof);
            };
        shared_ptr<sphericalModel> voronoiModel;
        //!Set the number of degrees of freedom of the equation of motion
        void setMu(scalar _mu){mu=_mu;};
        void setEta(scalar _Eta){Eta=_Eta;};
        void setTau(scalar _tau){tau=_tau;};
        void setV0(scalar _v0){v0=_v0;};
        //!The value of the alignment coupling -- tau = J^-1
        scalar tau;
        //!The value of the inverse friction constant
        scalar mu;
        //!The value of the strength of vectorial noise
        scalar Eta;
        //intrinsic self-propulsion
        scalar v0;

        noiseSource noise;
        vector<dVec> newVelocityDirector;
    };
#endif

