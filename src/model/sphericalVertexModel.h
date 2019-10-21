#ifndef sphericalVertexModel_H
#define sphericalVertexModel_H

#include "sphericalModel.h"
#include "convexHullCGAL.h"

//!set up a spherical vertex model. N is the number of CELLS, not vertices
class sphericalVertexModel : public sphericalModel
    {
    public:
        sphericalVertexModel(int n, noiseSource &_noise, bool _useGPU=false, bool _neverGPU = true);

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
        GPUArray<dVec> cellPositions;

        GPUArray<dVec> currentVertexAroundCell;
        GPUArray<dVec> lastVertexAroundCell;
        GPUArray<dVec> nextVertexAroundCell;
        int nCells;

    };
#endif
