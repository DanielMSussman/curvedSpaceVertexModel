#include "sphericalModel.h"
/*! \file sphericalModel.cpp" */

sphericalModel::sphericalModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : simpleModel(n,_useGPU, _neverGPU)
    {
    if(neverGPU)
        {
        numberOfNeighbors.noGPU = true;
        neighbors.noGPU = true;
        }
    numberOfNeighbors.resize(N);

    //set random positions on the sphere of radius 1
    setParticlePositionsRandomly(_noise);

    getNeighbors();
    }

void sphericalModel::setRadius(scalar _r)
    {
    sphere.radius = _r;
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(velocities);
    printf(" %f %f %f\t%f %f %f\n",p.data[0][0],p.data[0][1],p.data[0][2],V.data[0][0],V.data[0][1],V.data[0][2]);
    for (int ii = 0; ii < N; ++ii)
        {
        //sphere.putInBoxReal(p.data[ii]);
        sphere.move(p.data[ii],V.data[ii],0.00);
        }
    printf(" %f %f %f\t%f %f %f\n",p.data[0][0],p.data[0][1],p.data[0][2],V.data[0][0],V.data[0][1],V.data[0][2]);
    };


void sphericalModel::setParticlePositionsRandomly(noiseSource &noise)
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(velocities);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u = noise.getRealUniform();
        scalar v = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*v-1);
        p.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        p.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        p.data[ii].x[2] = 1.0*cos(theta);
        }
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u2 = noise.getRealUniform();
        scalar v2 = noise.getRealUniform();
        scalar phi = 2.0*PI*u2;
        scalar theta = acos(2.0*v2-1);
        V.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        V.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        V.data[ii].x[2] = 1.0*cos(theta);
        //project the velocity onto the tangent plane
        sphere.move(p.data[ii],V.data[ii],0.00);
        }
    printf(" %f %f %f\t%f %f %f\n",p.data[0][0],p.data[0][1],p.data[0][2],V.data[0][0],V.data[0][1],V.data[0][2]);
    }

void sphericalModel::moveParticles(GPUArray<dVec> &displacements, scalar scale)
    {
        {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(displacements);
    for(int ii = 0; ii < N; ++ii)
        sphere.move(p.data[ii],V.data[ii],scale);
        }
    getNeighbors();
    };
