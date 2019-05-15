#ifndef voronoiVicsek_H
#define voronoiVicsek_H

#include "equationOfMotion.h"
#include "sphericalVoronoi.h"
#include "noiseSource.h"
/*! \file voronoiVicsek.h */

class voronoiVicsek : public equationOfMotion
    {
    public:
        voronoiVicsek(){useGPU = false; mu = 1.0; Eta = 1.0; tau = 1.0;v0=0.01;};
        virtual void integrateEOMGPU(){};
        virtual void integrateEOMCPU();

        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            voronoiModel = dynamic_pointer_cast<sphericalVoronoi>(model);
            initializeFromModel();
            newVelocityDirector.resize(Ndof);
            };
        shared_ptr<sphericalVoronoi> voronoiModel;
        //!Set the number of degrees of freedom of the equation of motion
        void setMu(scalar _mu){mu=_mu;};
        void setEta(scalar _Eta){Eta=_Eta;};
        void setTau(scalar _tau){tau=_tau;};
        void setV0(scalar _v0){v0=_v0;};
        //!The value of the alignment coupling
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

