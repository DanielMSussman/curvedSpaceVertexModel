#include "vectorialVicsek.h"
#include "vectorialVicsek.cuh"

/*! \file vectorialVicsek.cpp */

void vectorialVicsek::integrateEOMGPU()
    {
    sim->computeForces();
    {
    ArrayHandle<dVec> n(model->returnDirectors(),access_location::device,access_mode::read);
    ArrayHandle<dVec> f(model->returnForces(),access_location::device,access_mode::read);
    ArrayHandle<dVec> disp(displacement,access_location::device,access_mode::overwrite);
    gpu_selfPropelledParticleDisplacement(f.data,n.data,disp.data,v0,mu,deltaT,Ndof);
    }

    sim->moveParticles(displacement);

    //find new target directors
    {//ArrayHandle scope
    ArrayHandle<dVec> n(model->returnDirectors(),access_location::device,access_mode::read);
    ArrayHandle<dVec> nDisp(newVelocityDirector,access_location::device,access_mode::overwrite);
    ArrayHandle<unsigned int> nNeighs(model->numberOfNeighbors,access_location::device,access_mode::read);
    ArrayHandle<int> neighs(model->neighbors,access_location::device,access_mode::read);
    ArrayHandle<curandState> rngs(noise.RNGs,access_location::device,access_mode::readwrite);

    gpu_vectorVicsek_target_directors(n.data,nDisp.data,nNeighs.data,neighs.data,rngs.data,
                                      model->neighborIndex,Eta,Ndof);
    }

    //update actual directors
    {
    ArrayHandle<dVec> n(model->returnDirectors(),access_location::device,access_mode::readwrite);
    ArrayHandle<dVec> nDisp(newVelocityDirector,access_location::device,access_mode::read);
    gpu_vectorVicsek_update_directors(n.data,nDisp.data,tau,deltaT,Ndof);
    }

    };

void vectorialVicsek::integrateEOMCPU()
    {
    sim->computeForces();

    {
    ArrayHandle<dVec> n(model->returnDirectors());
    ArrayHandle<dVec> f(model->returnForces());
    ArrayHandle<dVec> disp(displacement);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        disp.data[ii] = deltaT*(v0*n.data[ii]+mu*f.data[ii]);
        }
    }

    sim->moveParticles(displacement);
    //update directors
    {//ArrayHandle scope
    ArrayHandle<dVec> p(model->returnPositions());
    ArrayHandle<dVec> n(model->returnDirectors());
    ArrayHandle<dVec> nDisp(newVelocityDirector);
    ArrayHandle<unsigned int> nNeighs(model->numberOfNeighbors);
    ArrayHandle<int> neighs(model->neighbors);
    dVec spherePoint;
    for(int ii = 0; ii < Ndof; ++ii)
        {
        //get a random point on the sphere
#if DIMENSION == 3
        scalar u = noise.getRealUniform();
        scalar w = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*w-1);
        spherePoint.x[0] = 1.0*sin(theta)*cos(phi);
        spherePoint.x[1] = 1.0*sin(theta)*sin(phi);
        spherePoint.x[2] = 1.0*cos(theta);
#else
        scalar u = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        spherePoint.x[0] = cos(phi);
        spherePoint.x[1] = sin(phi);
#endif
        //average direction of neighbors?
        int m = nNeighs.data[ii];

        nDisp.data[ii] = n.data[ii];
        for (int jj = 0; jj < m; ++jj)
            {
            nDisp.data[ii] += n.data[neighs.data[model->neighborIndex(jj,ii)]];
            }
        m +=1;//account for self-alignment
        scalar mi = 1.0/m;
        nDisp.data[ii] = nDisp.data[ii] * mi + spherePoint*Eta;
        }
    for (int ii = 0; ii < Ndof; ++ii)
        {
//        n.data[ii] += (deltaT/tau)*(nDisp.data[ii] - n.data[ii]);
//        n.data[ii] = n.data[ii]*(1.0/norm(n.data[ii]));
        n.data[ii] = nDisp.data[ii];
        };

    }//arrayhandle scope
    };
