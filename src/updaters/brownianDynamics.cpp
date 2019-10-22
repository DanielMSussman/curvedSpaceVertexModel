#include "brownianDynamics.h"

brownianDynamics::brownianDynamics(bool _reproducible)
    {
    useGPU = false;
    temperature = 0.0;
    deltaT = 0.01;
    mu = 1.0;
    reproducible = _reproducible;
    };

void brownianDynamics::integrateEOMGPU()
    {
    UNWRITTENCODE("no brownian gpu yet");
    };

void brownianDynamics::integrateEOMCPU()
    {
    //compute the forces
    sim->computeForces();

    {//scope for array handles
    scalar forcePrefactor = deltaT*mu;
    scalar noisePrefactor = sqrt(2.0*forcePrefactor*temperature);
    ArrayHandle<dVec> f(model->returnForces());
    ArrayHandle<dVec> disp(displacement);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        for(int dd = 0; dd < DIMENSION; ++dd)
            disp.data[ii][dd] = noise.getRealNormal()*noisePrefactor + f.data[ii][dd]*forcePrefactor;
        }
    };//end array handle scope

    sim->moveParticles(displacement);
    }

void brownianDynamics::setT(scalar _t)
    {
    temperature = _t;
    }

