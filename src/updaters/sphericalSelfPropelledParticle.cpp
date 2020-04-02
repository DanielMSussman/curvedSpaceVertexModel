#include "sphericalSelfPropelledParticle.h"

/*! \file sphericalSelfPropelledParticle.cpp */

void sphericalSelfPropelledParticle::setModel(shared_ptr<simpleModel> _model)
    {
    model=_model;
    sphereModel = dynamic_pointer_cast<sphericalModel>(model);
    initializeFromModel();
    if(model->neverGPU)
            newVelocityDirector.noGPU = true;
    newVelocityDirector.resize(Ndof);
    };


void sphericalSelfPropelledParticle::integrateEOMCPU()
    {
    sim->computeForces();

    {
    ArrayHandle<dVec> n(sphereModel->returnDirectors());
    ArrayHandle<dVec> f(sphereModel->returnForces());
    ArrayHandle<dVec> disp(displacement);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        disp.data[ii] = deltaT*(v0*n.data[ii]+mu*f.data[ii]);
        }

    }

    sim->moveParticles(displacement);
    //update directors
    {//ArrayHandle scope
    ArrayHandle<dVec> p(sphereModel->returnPositions());
    ArrayHandle<dVec> n(sphereModel->returnDirectors());
    ArrayHandle<dVec> nDisp(newVelocityDirector);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        scalar randomNumber = noise.getRealUniform(-.5,.5);
        scalar theta = randomNumber*sqrt(2.0*deltaT*tau);
        sphereModel->sphere->projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
        nDisp.data[ii] = n.data[ii];
        rodriguesRotation(nDisp.data[ii],p.data[ii],theta);
        n.data[ii] = nDisp.data[ii];
        };
    }//arrayhandle scope
    };
