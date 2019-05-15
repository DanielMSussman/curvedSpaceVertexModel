#include "voronoiVicsek.h"

/*! \file voronoiVicsek.cpp */

void voronoiVicsek::integrateEOMCPU()
    {
    sim->computeForces();

    sim->moveParticles(voronoiModel->returnVelocities(),deltaT*v0);
    {//ArrayHandle scope
    ArrayHandle<dVec> p(voronoiModel->returnPositions());
    ArrayHandle<dVec> v(voronoiModel->returnVelocities());
    dVec spherePoint;
    for(int ii = 0; ii < Ndof; ++ii)
        {
        //get a random point on the sphere
        scalar u = noise.getRealUniform();
        scalar w = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*w-1);
        spherePoint.x[0] = 1.0*sin(theta)*cos(phi);
        spherePoint.x[1] = 1.0*sin(theta)*sin(phi);
        spherePoint.x[2] = 1.0*cos(theta);
        //project it onto the tangent plane
        voronoiModel->sphere.projectToTangentPlane(spherePoint,p.data[ii]);
        spherePoint = spherePoint*(1.0/norm(spherePoint));
        //average direction of neighbors?
        int m = voronoiModel->convexHuller.numNeighs[ii];
        newVelocityDirector[ii] = make_dVec(0.0);
        for (int jj = 0; jj < m; ++jj)
            {
            newVelocityDirector[ii] += v.data[voronoiModel->convexHuller.allNeighs[ii][jj]];
            }
        newVelocityDirector[ii] = newVelocityDirector[ii] * (1.0/m) + spherePoint*Eta;
        newVelocityDirector[ii] = newVelocityDirector[ii]*(1.0/norm(newVelocityDirector[ii]));
        }
    for (int ii = 0; ii < Ndof; ++ii)
        v.data[ii] = newVelocityDirector[ii];

    }//arrayhandle scope
    };
