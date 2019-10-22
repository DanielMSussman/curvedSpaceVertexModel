#ifndef sphericalVertexModel_H
#define sphericalVertexModel_H

#include "sphericalModel.h"
#include "convexHullCGAL.h"

//!set up a spherical vertex model. N is the number of CELLS, not vertices
class sphericalVertexModel : public sphericalModel
    {
    public:
        sphericalVertexModel(int n, noiseSource &_noise, scalar _area = 1.0, scalar _perimeter = 3.8, bool _useGPU=false, bool _neverGPU = true);

        virtual void getNeighbors(){};

        virtual void computeGeometry()
            {
            if(useGPU)
                computeGeometryGPU();
            else
                computeGeometryCPU();
            };
        virtual void computeForces(bool zeroOutForces)
            {
            if(useGPU)
                {
                computeForceGPU();
                }
            else
                {
                computeForceCPU();
                }
            };
        virtual void moveParticles(GPUArray<dVec> &displacements,scalar scale = 1.);

        virtual void enforceTopology(){};
        virtual void computeForceCPU();
        virtual void computeForceGPU();

        virtual void computeGeometryCPU();
        virtual void computeGeometryGPU();
        convexHullCGALInterface convexHuller;
        Index2D cellNeighborIndex;
        GPUArray<unsigned int> cellNumberOfNeighbors;
        GPUArray<int> cellNeighbors;
        GPUArray<int> vertexCellNeighbors;
        GPUArray<scalar2> areaPerimeter;
        GPUArray<scalar2> areaPerimeterPreference;
        GPUArray<dVec> cellPositions;

        GPUArray<dVec> currentVertexAroundCell;
        GPUArray<dVec> lastVertexAroundCell;
        GPUArray<dVec> nextVertexAroundCell;
        int nCells;
        int maxVNeighs;

        virtual void setPreferredParameters(scalar _a0, scalar _p0);
    };
#endif
