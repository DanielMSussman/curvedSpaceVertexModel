#ifndef sphericalVectorialVicsek_H
#define sphericalVectorialVicsek_H

#include "vectorialVicsek.h"
#include "sphericalModel.h"
/*! \file sphericalVectorialVicsek.h */

class sphericalVectorialVicsek : public vectorialVicsek
    {
    public:
        sphericalVectorialVicsek(){useGPU = false; mu = 1.0; Eta = 1.0; tau = 1.0;v0=0.01;};
        virtual void integrateEOMGPU(){};
        virtual void integrateEOMCPU();

        //! virtual function to allow the model to be a derived class
        virtual void setModel(shared_ptr<simpleModel> _model)
            {
            model=_model;
            voronoiModel = dynamic_pointer_cast<sphericalModel>(model);
            initializeFromModel();
            if(model->neverGPU)
                newVelocityDirector.noGPU = true;
            newVelocityDirector.resize(Ndof);
            };
        shared_ptr<sphericalModel> voronoiModel;
        //!Set the number of degrees of freedom of the equation of motion
    };
#endif

