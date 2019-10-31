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
        //!compute the current PE
        virtual scalar computeEnergy();

        virtual void enforceTopology();
        virtual void computeForceCPU();
        virtual void computeForceGPU();

        //! computes areas via included angles
        virtual void computeGeometryCPU();
        //! if included angles are bad due to edge-crossings, approximate via triangles
        virtual void recomputeAreasCPU();
        virtual void computeGeometryGPU();
        convexHullCGALInterface convexHuller;
        Index2D cellNeighborIndex;
        //!number of vertices around each cell
        GPUArray<unsigned int> cellNumberOfNeighbors;
        //!vertex indices composing each cell
        GPUArray<int> cellNeighbors;
        //!cells around each vertex (vertex neighbor number in "numberOfNeighbors")
        GPUArray<int> vertexCellNeighbors;
        GPUArray<scalar2> areaPerimeter;
        GPUArray<scalar2> areaPerimeterPreference;
        GPUArray<dVec> cellPositions;

        GPUArray<dVec> currentVertexAroundCell;
        GPUArray<dVec> lastVertexAroundCell;
        GPUArray<dVec> nextVertexAroundCell;
        int nCells;
        int maxVNeighs;
        int maximumVerticesPerCell;

        scalar Kr;
        scalar t1Threshold;
        scalar energy;
        virtual void setPreferredParameters(scalar _a0, scalar _p0);
        virtual void setScalarModelParameter(scalar _param){Kr=_param;};
        //!Set the length threshold for T1 transitions
        virtual void setT1Threshold(scalar t1t){t1Threshold = t1t;};

    protected:
        void preserveOrientatedFaces();
        //!Simple test for T1 transitions (edge length less than threshold) on the CPU
        void testAndPerformT1TransitionsCPU();
        //!Simple test for T1 transitions (edge length less than threshold) on the GPU...calls the following functions
        void testAndPerformT1TransitionsGPU();
        
        //!For finding T1s on the CPU; find the set of vertices and cells involved in the transition
        void getCellVertexSetForT1(int v1, int v2, int4 &cellSet, int4 &vertexSet, bool &growList);

        //!if the maximum number of vertices per cell increases, grow the cellVertices list
        void growCellVerticesList(int newVertexMax);

        //!Initialize the data structures for edge flipping...should also be called if Nvertices changes
        void initializeEdgeFlipLists();
        //! data structure to help with cell-vertex list
        GPUArray<int> growCellVertexListAssist;

        //!test the edges for a T1 event, and grow the cell-vertex list if necessary
        void testEdgesForT1GPU();
        //!perform the edge flips found in the previous step
        void flipEdgesGPU();
        //! data structure to help with not simultaneously trying to flip nearby edges
        GPUArray<int> finishedFlippingEdges;

        //! data structure per cell for not simulataneously flipping nearby edges
        GPUArray<int> cellEdgeFlips;
        //! data structure per cell for not simulataneously flipping nearby edges
        GPUArray<int4> cellSets;
    };
#endif
