#include "simpleModel.h"
#include "utilities.cuh"
#include "simpleModel.cuh"
/*! \file simpleModel.cpp" */

/*!
 * Set the size of basic data structures...
*/
simpleModel::simpleModel(int n, bool _useGPU, bool _neverGPU) :
    N(n), useGPU(_useGPU), neverGPU(_neverGPU)
    {
    cout << "initializing a model with "<< N << " particles" << endl;
    initializeSimpleModel(n);
    Box = make_shared<periodicBoundaryConditions>(pow(N,1.0/DIMENSION));
    };

simpleModel::simpleModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : simpleModel(n,_useGPU, _neverGPU)
    {
    //set random positions on the sphere of radius 1
    setParticlePositionsRandomly(_noise);
    }
/*!
 * actually set the array sizes. positions, velocities, forces are zero
 * masses are set to unity
*/
void simpleModel::initializeSimpleModel(int n)
    {
    N=n;
    selfForceCompute = false;
    if(neverGPU)
        {
        positions.noGPU =true;
        velocities.noGPU = true;
        forces.noGPU = true;
        masses.noGPU = true;
        radii.noGPU = true;
        types.noGPU = true;
        numberOfNeighbors.noGPU = true;
        neighbors.noGPU = true;
        directors.noGPU=true;
        }
    numberOfNeighbors.resize(N);
    directors.resize(N);

    positions.resize(n);
    velocities.resize(n);
    forces.resize(n);
    masses.resize(n);
    radii.resize(n);
    types.resize(n);
    vector<dVec> zeroes(N,make_dVec(0.0));
    vector<scalar> ones(N,1.0);
    vector<scalar> halves(N,.5);
    vector<int> units(N,0);
    fillGPUArrayWithVector(units,types);
    fillGPUArrayWithVector(zeroes,positions);
    fillGPUArrayWithVector(zeroes,velocities);
    fillGPUArrayWithVector(zeroes,forces);
    fillGPUArrayWithVector(ones,masses);
    fillGPUArrayWithVector(halves,radii);
    };

scalar simpleModel::computeKineticEnergy()
    {
    ArrayHandle<scalar> h_m(masses,access_location::host,access_mode::read);
    ArrayHandle<dVec> h_v(velocities);
    scalar en = 0.0;
    for (int ii = 0; ii < N; ++ii)
        {
        en += 0.5*h_m.data[ii]*dot(h_v.data[ii],h_v.data[ii]);
        }
    return en;
    };

void simpleModel::getMeanDirection(dVec &meanDir)
    {
    for (int dd = 0; dd < DIMENSION; ++dd)
        meanDir[dd] = 0.0;
    ArrayHandle<dVec> n(directors);
    for (int ii = 0; ii < N; ++ii)
        meanDir += n.data[ii];

    meanDir = meanDir*(1./N);
    };

scalar simpleModel::computeInstantaneousTemperature(bool fixedMomentum)
    {
    ArrayHandle<scalar> h_m(masses,access_location::host,access_mode::read);
    ArrayHandle<dVec> h_v(velocities);
    scalar en = 0.0;
    for (int ii = 0; ii < N; ++ii)
        {
        en += 1.0*h_m.data[ii]*dot(h_v.data[ii],h_v.data[ii]);
        }
    if(fixedMomentum)
        return en /((N-DIMENSION)*DIMENSION);
    else
        return en /(N*DIMENSION);
    };

void simpleModel::setParticlePositions(vector<dVec> &newPositions)
    {
    if(N !=newPositions.size())
        initializeSimpleModel(newPositions.size());
    ArrayHandle<dVec> p(positions);
    for (int pp = 0;pp < N; ++pp)
        {
        p.data[pp] = newPositions[pp];
        Box->putInBoxReal(p.data[pp]);
        };
    };
void simpleModel::setParticlePositions(GPUArray<dVec> &newPositions)
    {
    if(N !=newPositions.getNumElements())
        initializeSimpleModel(newPositions.getNumElements());
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> np(newPositions);
    for (int pp = 0;pp < N; ++pp)
        {
        p.data[pp] = np.data[pp];
        Box->putInBoxReal(p.data[pp]);
        };
    };
/*!
 */
void simpleModel::setParticlePositionsRandomly(noiseSource &noise)
    {
    dVec bDims;
    Box->getBoxDims(bDims);
    ArrayHandle<dVec> n(directors);
    ArrayHandle<dVec> pos(positions);
    printf("setting in a box of total side length %f\n", bDims.x[0]);
    for(int pp = 0; pp < N; ++pp)
        for (int dd = 0; dd <DIMENSION; ++dd)
            {
            pos.data[pp].x[dd] = noise.getRealUniform(-0.5*bDims.x[dd],0.5*bDims.x[dd]);
            Box->putInBoxReal(pos.data[pp]);
            };
    for (int ii = 0; ii < N; ++ii)
        {
#if DIMENSION == 3
        scalar u = noise.getRealUniform();
        scalar w = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        scalar theta = acos(2.0*w-1);
        n.data[ii].x[0] = 1.0*sin(theta)*cos(phi);
        n.data[ii].x[1] = 1.0*sin(theta)*sin(phi);
        n.data[ii].x[2] = 1.0*cos(theta);
#else
        scalar u = noise.getRealUniform();
        scalar phi = 2.0*PI*u;
        n.data[ii].x[0] = cos(phi);
        n.data[ii].x[1] = sin(phi);
#endif
        //project the velocity onto the tangent plane
        }
    };

scalar simpleModel::setVelocitiesMaxwellBoltzmann(scalar T,noiseSource &noise)
    {
    ArrayHandle<scalar> h_m(masses,access_location::host,access_mode::read);
    ArrayHandle<dVec> h_v(velocities);
    scalar KE = 0.0;
    dVec P(0.0);
    for (int ii = 0; ii < N; ++ii)
        {
        for (int dd = 0; dd <DIMENSION; ++dd)
            h_v.data[ii].x[dd] = noise.getRealNormal(0.0,sqrt(T/h_m.data[ii]));
        P += h_m.data[ii]*h_v.data[ii];
        KE += 0.5*h_m.data[ii]*dot(h_v.data[ii],h_v.data[ii]);
        }
    //remove excess momentum, calculate the ke
    KE = 0.0;
    for (int ii = 0; ii < N; ++ii)
        {
        h_v.data[ii] += (-1.0/(N*h_m.data[ii]))*P;
        KE += 0.5*h_m.data[ii]*dot(h_v.data[ii],h_v.data[ii]);
        };
    return KE;
    };

/*!
 * move particles on either CPU or gpu
*/
void simpleModel::moveParticles(GPUArray<dVec> &displacement, scalar scale)
    {
    if(!useGPU)
        {//cpu branch
        ArrayHandle<dVec> h_disp(displacement, access_location::host,access_mode::read);
        ArrayHandle<dVec> h_pos(positions);
        #include "ompParallelLoopDirective.h"
        for(int pp = 0; pp < N; ++pp)
            {
            h_pos.data[pp] += scale*h_disp.data[pp];
            Box->putInBoxReal(h_pos.data[pp]);
            }
        }
    else
        {//gpu branch
        ArrayHandle<dVec> d_disp(displacement,access_location::device,access_mode::readwrite);
        ArrayHandle<dVec> d_pos(positions,access_location::device,access_mode::readwrite);
        gpu_move_particles(d_pos.data,d_disp.data,*(Box),scale,N);
        };
    forcesComputed = false;
    getNeighbors();
    };

void simpleModel::setRadius(scalar _r)
    {
    Box = make_shared<periodicBoundaryConditions>(2.0*_r);
    ArrayHandle<dVec> p(positions);
    for (int ii = 0; ii < N; ++ii)
        {
        Box->putInBoxReal(p.data[ii]);
        }
    metricNeighbors.setBasics(1.0,2.0*_r);
    //metricNeighbors.setGPU(useGPU);
    getNeighbors();
    };
/*!
 *
*/

void simpleModel::computeHarmonicRepulsions(bool zeroOutForces)
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

void simpleModel::computeForces(bool zeroOutForces)
    {
    if(selfForceCompute)
        computeHarmonicRepulsions(zeroOutForces);
    if(zeroOutForces && !selfForceCompute)
        {
        if(!useGPU)
            {//cpu branch
            ArrayHandle<dVec> h_f(forces);
            dVec dArrayZero(0.0);
            for(int pp = 0; pp <N;++pp)
                h_f.data[pp] = dArrayZero;
            }
        else
            {
                ArrayHandle<dVec> d_f(forces);
                gpu_zero_array(d_f.data,N);
            };
        };
    };

void simpleModel::setSoftRepulsion(scalar range, scalar stiffness)
    {
    cout << "soft repulsion set" << endl;
    selfForceCompute = true;
    repulsionRange = range;
    repulsionStiffness = stiffness;
    }



void simpleModel::getNeighbors()
    {
    metricNeighbors.computeNeighborLists(positions);
    numberOfNeighbors = metricNeighbors.neighborsPerParticle;
    neighbors = metricNeighbors.particleIndices;
    neighborIndex=metricNeighbors.neighborIndexer;
    };

