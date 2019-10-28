#ifndef SIMPLEMODEL_H
#define SIMPLEMODEL_H

#include "std_include.h"
#include "gpuarray.h"
#include "periodicBoundaryConditions.h"
#include "functions.h"
#include "noiseSource.h"
#include "indexer.h"
#include "neighborList.h"

/*! \file simpleModel.h
 * \brief defines an interface for models that compute forces
 */

//! A base interfacing class that defines common operations
/*!
This provides an interface, guaranteeing that SimpleModel S will provide access to
S.setGPU();
S.getNumberOfParticles();
S.computeForces();
S.moveParticles();
S.returnForces();
S.returnPositions();
S.returnDirectors();
S.returnVelocities();
S.returnRadii();
S.returnMasses();
S.returnTypes();
S.spatialSorting();
S.returnAdditionalData();
*/
class simpleModel
    {
    public:
        //!The base constructor requires the number of particles
        simpleModel(int n, bool _useGPU = false, bool _neverGPU=true);
        //!The base constructor requires the number of particles
        simpleModel(int n, noiseSource &_noise,bool _useGPU = false, bool _neverGPU=true);
        //!a blank default constructor
        simpleModel(){};
        //!initialize the size of the basic data structure arrays
        void initializeSimpleModel(int n);

        virtual void setRadius(scalar _r);

        //!Enforce GPU operation
        virtual void setGPU(bool _useGPU=true){useGPU = _useGPU;};
        //!get the number of degrees of freedom, defaulting to the number of cells
        virtual int getNumberOfParticles(){return N;};
        //!move the degrees of freedom
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);
        //!do everything unusual to compute additional forces... by default, sets forces to zero
        virtual void computeForces(bool zeroOutForces=false);

        virtual void computeHarmonicRepulsions(bool zeroOutForces);

        void setParticlePositions(GPUArray<dVec> &newPositions);
        void setParticlePositions(vector<dVec> &newPositions);
        virtual void setParticlePositionsRandomly(noiseSource &noise);
        virtual void setParticlePositionsBandedRandomly(noiseSource &noise,scalar angularExtent){setParticlePositionsRandomly(noise);};

        //!Set velocities via a temperature. The return value is the total kinetic energy
        scalar setVelocitiesMaxwellBoltzmann(scalar T,noiseSource &noise);

        //!compute the current KE
        virtual scalar computeKineticEnergy();
        //!compute the dimension-dependent instantaneous temperature
        virtual scalar computeInstantaneousTemperature(bool fixedMomentum=true);

        //!do everything necessary to perform a Hilbert sort
        virtual void spatialSorting(){};

        //!return a reference to the GPUArray of positions
        virtual GPUArray<dVec> & returnPositions(){return positions;};
        //!return a reference to the GPUArray of the current forces
        virtual GPUArray<dVec> & returnForces(){return forces;};
        //!return a reference to the GPUArray of the particle radii
        virtual GPUArray<scalar> & returnRadii(){return radii;};
        //!return a reference to the GPUArray of the integer types
        virtual GPUArray<int> & returnTypes(){return types;};
        //!return a reference to the GPUArray of the masses
        virtual GPUArray<scalar> & returnMasses(){return masses;};
        //!return a reference to the GPUArray of the current velocities
        virtual GPUArray<dVec> & returnVelocities(){return velocities;};

        //!Does this model have a special force it needs to compute itself?
        bool selfForceCompute;

        //!The space in which the particles live
        BoxPtr Box;
        //!Are the forces current? set to false after every call to moveParticles. set to true after the SIMULATION calls computeForces
        bool forcesComputed;

        virtual void setSoftRepulsion(scalar range = 1.0, scalar stiffness = 1.0);

        //!allow for setting multiple threads
        virtual void setNThreads(int n){nThreads = n;};

        virtual void getNeighbors();
        //!return a reference to the GPUArray of positions
        virtual GPUArray<dVec> & returnDirectors(){return directors;};
        neighborList metricNeighbors;

        virtual void getMeanDirection(dVec &meanDir);
        virtual void getMeanForce(dVec &meanForce){};

        Index2D neighborIndex;
        GPUArray<unsigned int> numberOfNeighbors;
        GPUArray<int> neighbors;
        GPUArray<dVec> directors;

        noiseSource noise;
        scalar repulsionRange;
        scalar repulsionStiffness;
    protected:
        //!The number of particles
        int N;
        //!number of threads to use
        int nThreads=1;
        //!particle  positions
        GPUArray<dVec> positions;
        //!particle velocities
        GPUArray<dVec> velocities;
        //!Forces on particles
        GPUArray<dVec> forces;
        //!particle radii
        GPUArray<scalar> radii;
        //!particle masses
        GPUArray<scalar> masses;
        //!particle types
        GPUArray<int> types;

    public:
        virtual void setPreferredParameters(scalar _a0, scalar _p0){};
        virtual void setScalarModelParameter(scalar _param){};
        //!Whether the GPU should be used to compute anything
        bool useGPU;

        bool neverGPU;

    };
typedef shared_ptr<simpleModel> ConfigPtr;
typedef weak_ptr<simpleModel> WeakConfigPtr;
#endif
