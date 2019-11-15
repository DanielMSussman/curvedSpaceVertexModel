#include "brownianDynamics.h"
#include "brownianDynamics.cuh"

brownianDynamics::brownianDynamics(bool _reproducible)
    {
    useGPU = false;
    temperature = 0.0;
    deltaT = 0.01;
    mu = 1.0;
    reproducible = _reproducible;
    updateProfiler.setName("brownian dynamics profiler");
    };

void brownianDynamics::integrateEOMGPU()
    {
    sim->computeForces();
    updateProfiler.start();
    {
    ArrayHandle<dVec> d_f(model->returnForces(),access_location::device,access_mode::read);
    ArrayHandle<dVec> d_disp(displacement,access_location::device,access_mode::overwrite);
    ArrayHandle<curandState> d_RNG(noise.RNGs,access_location::device,access_mode::readwrite);
    gpu_brownian_eom_integration(d_f.data,d_disp.data,d_RNG.data,
                               Ndof,deltaT,mu,temperature);
    }
    updateProfiler.end();
    sim->moveParticles(displacement);
    };

void brownianDynamics::integrateEOMCPU()
    {
    //compute the forces
    sim->computeForces();

    updateProfiler.start();
    {//scope for array handles
    scalar forcePrefactor = deltaT*mu;
    scalar noisePrefactor = sqrt(2.0*forcePrefactor*temperature);
    ArrayHandle<dVec> f(model->returnForces());
    ArrayHandle<dVec> disp(displacement);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        for(int dd = 0; dd < DIMENSION; ++dd)
            disp.data[ii][dd] = noise.getRealNormalStandardized()*noisePrefactor + f.data[ii][dd]*forcePrefactor;
        }
    };//end array handle scope
    updateProfiler.end();
    sim->moveParticles(displacement);
    }

void brownianDynamics::setT(scalar _t)
    {
    temperature = _t;
    }

