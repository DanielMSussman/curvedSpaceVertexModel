#ifndef EQUATIONOFMOTION_H
#define EQUATIONOFMOTION_H
#include "baseUpdater.h"
/*! \file equationOfMotion.h */
//! child class just need to implement the CPU and GPU functions
class equationOfMotion : public updater
    {
    public:

        virtual void performUpdate(){integrateEquationOfMotion();};

        virtual void integrateEquationOfMotion()
            {
            if (model->getNumberOfParticles() != Ndof)
                initializeFromModel();
            if (useGPU)
                integrateEOMGPU();
            else
                integrateEOMCPU();
            };

        virtual void integrateEOMGPU(){};
        virtual void integrateEOMCPU(){};

        virtual void initializeFromModel()
            {
            printf("setting model size for updater\n");
            Ndof = model->getNumberOfParticles();
            if(model->neverGPU)
                displacement.noGPU = true;
            displacement.resize(Ndof);
            };
        //!an array of displacements
        GPUArray<dVec> displacement;
    };
#endif
