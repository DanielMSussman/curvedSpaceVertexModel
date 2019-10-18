#include "sphericalModel.h"
/*! \file sphericalModel.cpp" */

sphericalModel::sphericalModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : simpleModel(n,_useGPU, _neverGPU)
    {
    //set random positions on the sphere of radius 1
    setParticlePositionsRandomly(_noise);
    }

void sphericalModel::setRadius(scalar _r)
    {
    sphere.changeRadius(_r);
    inverseRadius = 1.0/sphere.radius;
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(velocities);
    for (int ii = 0; ii < N; ++ii)
        {
        sphere.putInBoxReal(p.data[ii]);
        sphere.putInBoxVirtual(V.data[ii]);
        }

    metricNeighbors.setBasics(1.0,sphere.radius+2.5);
    getNeighbors();
    };


void sphericalModel::setParticlePositionsRandomly(noiseSource &noise)
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> n(directors);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u = noise.getRealUniform();
        scalar v = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*v-1);
        p.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        p.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        p.data[ii].x[2] = 1.0*cos(theta);
        sphere.putInBoxReal(p.data[ii]);
        }
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u2 = noise.getRealUniform();
        scalar v2 = noise.getRealUniform();
        scalar phi = 2.0*PI*u2;
        scalar theta = acos(2.0*v2-1);
        n.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        n.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        n.data[ii].x[2] = 1.0*cos(theta);
        //project the velocity onto the tangent plane
        sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
        }
    //printf("%f %f %f\n", n.data[0][0],n.data[0][1],n.data[0][2]);
    }
void sphericalModel::setParticlePositionsBandedRandomly(noiseSource &noise, scalar angularExtent)
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> n(directors);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u = noise.getRealUniform();
        scalar v = noise.getRealUniform(0.5-angularExtent,0.5+angularExtent);
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*v-1);
        p.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        p.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        p.data[ii].x[2] = 1.0*cos(theta);
        sphere.putInBoxReal(p.data[ii]);
        }
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u2 = noise.getRealUniform();
        scalar v2 = noise.getRealUniform(0.5-angularExtent,0.5+angularExtent);
        scalar phi = 2.0*PI*u2;
        scalar theta = acos(2.0*v2-1);
        n.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        n.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        n.data[ii].x[2] = 1.0*cos(theta);
        //project the velocity onto the tangent plane
        sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
        }
    //printf("%f %f %f\n", n.data[0][0],n.data[0][1],n.data[0][2]);
    }

void sphericalModel::moveParticles(GPUArray<dVec> &displacements, scalar scale)
    {
    if(scale == 1.)
        {
        ArrayHandle<dVec> p(positions);
        ArrayHandle<dVec> V(displacements);
        ArrayHandle<dVec> n(directors);
        for(int ii = 0; ii < N; ++ii)
            {
            sphere.move(p.data[ii],V.data[ii]);
            sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
            }
        }
    else
        {
        ArrayHandle<dVec> p(positions);
        ArrayHandle<dVec> V(displacements);
        ArrayHandle<dVec> n(directors);
        for(int ii = 0; ii < N; ++ii)
            {
            sphere.move(p.data[ii],scale*V.data[ii]);
            sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
            }
        }
    getNeighbors();
    };
