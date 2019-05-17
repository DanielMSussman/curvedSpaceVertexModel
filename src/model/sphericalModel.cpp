#include "sphericalModel.h"
/*! \file sphericalModel.cpp" */

sphericalModel::sphericalModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : simpleModel(n,_useGPU, _neverGPU)
    {
    if(neverGPU)
        {
        numberOfNeighbors.noGPU = true;
        neighbors.noGPU = true;
        directors.noGPU=true;
        }
    numberOfNeighbors.resize(N);
    directors.resize(N);

    //set random positions on the sphere of radius 1
    setParticlePositionsRandomly(_noise);

    getNeighbors();
    }

void sphericalModel::setRadius(scalar _r)
    {
    sphere.radius = _r;
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> V(velocities);
    for (int ii = 0; ii < N; ++ii)
        {
        sphere.putInBoxReal(p.data[ii]);
        sphere.putInBoxVirtual(V.data[ii]);
        }
    metricNeighbors.setBasics(1.0,sphere.radius);
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
void sphericalModel::setParticlePositionsBandedRandomly(noiseSource &noise)
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> n(directors);
    for (int ii = 0; ii < N; ++ii)
        {
        scalar u = noise.getRealUniform();
        scalar v = noise.getRealUniform(0.3,0.7);
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

void sphericalModel::setSoftRepulsion(scalar range, scalar stiffness)
    {
    cout << "soft repulsion set" << endl;
    selfForceCompute = true;
    repulsionRange = range;
    repulsionStiffness = stiffness;
    }

void sphericalModel::computeForces(bool zeroOutForces)
    {
    if(!useGPU)
        {
        metricNeighbors.computeNeighborLists(positions);
        ArrayHandle<unsigned int> nNeighs(metricNeighbors.neighborsPerParticle);
        ArrayHandle<int> neighs(metricNeighbors.particleIndices);
        ArrayHandle<dVec> h_nv(metricNeighbors.neighborVectors);
        ArrayHandle<scalar> h_nd(metricNeighbors.neighborDistances);
        ArrayHandle<dVec> h_f(forces);
        ArrayHandle<dVec> p(positions);
        dVec dArrayZero(0.0);
        dVec disp;
        scalar dist;
        dVec newForce;
        for(int ii = 0; ii <N;++ii)
            {
            if(zeroOutForces)
                h_f.data[ii] = dArrayZero;
            int num = nNeighs.data[ii];
            for (int jj = 0; jj < num; ++jj)
                {
                int nIdx = metricNeighbors.neighborIndexer(jj,ii);
                int otherIdx = neighs.data[nIdx];
                if(otherIdx < ii)//newton force summation etc.
                    continue;
                disp = h_nv.data[nIdx];
                dist = h_nd.data[nIdx];
                if(dist < repulsionRange)
                    {
                    scalar delta = (1.0-dist/repulsionRange);
                    newForce = repulsionStiffness*(1.0/repulsionRange)*delta*(1.0/dist)*disp;
                    h_f.data[ii] = h_f.data[ii]+newForce;
                    h_f.data[otherIdx] = h_f.data[otherIdx]-newForce;
                    };
                }
            };

        }
    };
