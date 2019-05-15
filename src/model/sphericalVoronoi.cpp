#include "sphericalVoronoi.h"
/*! \file sphericalVoronoi.cpp" */

sphericalVoronoi::sphericalVoronoi(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : simpleModel(n,_useGPU, _neverGPU)
    {
    if(neverGPU)
        {
        numberOfVoronoiNeighbors.noGPU = true;
        voronoiNeighbors.noGPU = true;
        }
    numberOfVoronoiNeighbors.resize(N);

    //set random positions on the sphere of radius 1
    setParticlePositionsRandomly(_noise);

    convexHull();
    }

void sphericalVoronoi::convexHull()
    {
    ArrayHandle<dVec> pos(positions);
    convexHuller.sphericalConvexHull(pos.data,N);
    /*
    for (int ii = 0; ii < N; ++ii)
        {
        printf("particle %i (%f,%f,%f) neighbors: ",ii, pos.data[ii][0],pos.data[ii][1],pos.data[ii][2]);
        for (int jj = 0; jj < convexHuller.numNeighs[ii]; ++jj)
            printf("%i ",convexHuller.allNeighs[ii][jj]);
        printf("\n");
        }
    */
    };

void sphericalVoronoi::setParticlePositionsRandomly(noiseSource &noise)
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
        sphere.move(p.data[ii],V.data[ii],0.001);
        }
    }

void sphericalVoronoi::moveParticles(GPUArray<dVec> &displacements, scalar scale)
    {
        {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(displacements);
    for(int ii = 0; ii < N; ++ii)
        sphere.move(p.data[ii],V.data[ii],scale);
        }
    convexHull();
    };
