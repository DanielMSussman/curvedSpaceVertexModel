#include "sphericalVectorialVicsek.h"

/*! \file sphericalVectorialVicsek.cpp */

void sphericalVectorialVicsek::integrateEOMCPU()
    {
    sim->computeForces();

    {
    ArrayHandle<dVec> n(voronoiModel->returnDirectors());
    ArrayHandle<dVec> f(voronoiModel->returnForces());
    ArrayHandle<dVec> disp(displacement);
    for(int ii = 0; ii < Ndof; ++ii)
        {
        disp.data[ii] = deltaT*(v0*n.data[ii]+mu*f.data[ii]);
        }

    }

    sim->moveParticles(displacement);
    //update directors
    {//ArrayHandle scope
    ArrayHandle<dVec> p(voronoiModel->returnPositions());
    ArrayHandle<dVec> n(voronoiModel->returnDirectors());
    ArrayHandle<dVec> nDisp(newVelocityDirector);
    ArrayHandle<unsigned int> nNeighs(voronoiModel->numberOfNeighbors);
    ArrayHandle<int> neighs(voronoiModel->neighbors);
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
        voronoiModel->sphere->projectToTangentPlane(spherePoint,p.data[ii]);
        spherePoint = spherePoint*(1.0/norm(spherePoint));
        //average direction of neighbors?
        int m = nNeighs.data[ii];

        nDisp.data[ii] = n.data[ii];
        for (int jj = 0; jj < m; ++jj)
            {
            nDisp.data[ii] += n.data[neighs.data[voronoiModel->neighborIndex(jj,ii)]];
            }
        m +=1;//account for self-alignment
        scalar mi = 1.0/m;
        nDisp.data[ii] = nDisp.data[ii] * mi + spherePoint*Eta;

        voronoiModel->sphere->projectToTangentPlaneAndNormalize(nDisp.data[ii],p.data[ii]);

        }
    for (int ii = 0; ii < Ndof; ++ii)
        {
        n.data[ii] += (deltaT/tau)*(nDisp.data[ii] - n.data[ii]);
        voronoiModel->sphere->projectToTangentPlane(n.data[ii],p.data[ii]);
        n.data[ii] = n.data[ii]*(1.0/norm(n.data[ii]));
//        n.data[ii] = nDisp.data[ii];
        };

    }//arrayhandle scope
    };
