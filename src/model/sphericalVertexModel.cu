#include "sphericalVertexModel.cuh"
/*!
    \addtogroup modelKernels
    @{
*/
__global__ void gpu_move_particles_on_sphere_kernel(dVec *pos,
                      dVec *disp,
                      sphericalDomain sphere,
                      scalar scale,
                      int N
                      )
    {
    // read in the particle that belongs to this thread
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= N)
        return;
    sphere.move(pos[idx],scale * disp[idx]);
    };

bool gpu_move_particles_on_sphere(dVec *pos,
                                  dVec *disp,
                                  sphericalDomain &sphere,
                                  scalar scale,
                                  int N
                                  )
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = N/block_size + 1;
    gpu_move_particles_on_sphere_kernel<<<nblocks,block_size>>>(pos,disp,sphere,scale,N);

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    }

__global__ void gpu_spherical_vertex_model_geometry_kernel(dVec *vertexPos,
                                         dVec *cellPos,
                                         int *cellNeighbors,
                                         int *vertexCellNeighbors,
                                         unsigned int *vertexCellNumberOfNeighbors,
                                         dVec *currentVertexAroundCell,
                                         dVec *lastVertexAroundCell,
                                         dVec *nextVertexAroundCell,
                                         unsigned int *cellNumberOfNeighbors,
                                         scalar2 *areaPerimeter,
                                         Index2D cellNeighborIndex,
                                         Index2D neighborIndex,
                                         sphericalDomain sphere,
                                         int nCells
                                         )
    {
    // read in the particle that belongs to this thread
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= nCells)
        return;

    int neighs = cellNumberOfNeighbors[idx];
    dVec cPos(0.0);
    for (int nn = 0; nn < neighs;++nn)
        cPos = cPos + vertexPos[cellNeighbors[cellNeighborIndex(nn,idx)]];
    sphere.putInBoxReal(cPos);
    cellPos[idx] = cPos;
        
    int lastVertexIdx = cellNeighbors[cellNeighborIndex(neighs-2,idx)];
    int curVertexIdx = cellNeighbors[cellNeighborIndex(neighs-1,idx)];
    int nextVertexIdx;
    dVec lastVertexPos = vertexPos[lastVertexIdx];
    dVec curVertexPos = vertexPos[curVertexIdx];
    dVec nextVertexPos;
    scalar perimeter = 0.; 
    scalar area = 0.;
    scalar tempVal;
    for (int nn = 0; nn < neighs; ++nn)
        {
        int cni = cellNeighborIndex(nn,idx);
        int vNeighs = vertexCellNumberOfNeighbors[curVertexIdx];
        int forceSetIdx = -1;
        for (int vn = 0; vn < vNeighs; ++vn)
            {
            int newIdx = neighborIndex(vn,curVertexIdx);
            if(vertexCellNeighbors[newIdx] == idx)
                forceSetIdx = newIdx;
            }

        nextVertexIdx = cellNeighbors[cni];
        nextVertexPos = vertexPos[nextVertexIdx];

        sphere.geodesicDistance(lastVertexPos,curVertexPos,tempVal);
        perimeter += tempVal;
        sphere.includedAngle(lastVertexPos,curVertexPos,nextVertexPos,tempVal);
        area += tempVal;

        lastVertexAroundCell[forceSetIdx] = lastVertexPos;
        currentVertexAroundCell[forceSetIdx] = curVertexPos;
        nextVertexAroundCell[forceSetIdx] = nextVertexPos;
        
        lastVertexPos = curVertexPos;
        curVertexIdx = nextVertexIdx;
        curVertexPos = nextVertexPos;
        }
    area = (area-(neighs-2)*PI);
    int extraAngularArea = floor(area/(1.0*PI));
    if(extraAngularArea > 0)
        area -= extraAngularArea*PI;
    area *= (sphere.radius*sphere.radius);

    areaPerimeter[idx].x = area;
    areaPerimeter[idx].y = perimeter;
    };
        
bool gpu_spherical_vertex_model_geometry(dVec *vertexPos,
                                         dVec *cellPos,
                                         int *cellNeighbors,
                                         int *vertexCellNeighbors,
                                         unsigned int *vertexCellNumberOfNeighbors,
                                         dVec *currentVertexAroundCell,
                                         dVec *lastVertexAroundCell,
                                         dVec *nextVertexAroundCell,
                                         unsigned int *cellNumberOfNeighbors,
                                         scalar2 *areaPerimeter,
                                         Index2D cellNeighborIndex,
                                         Index2D neighborIndex,
                                         sphericalDomain &sphere,
                                         int nCells
                                         )
    {
    unsigned int block_size = 512;
    if (nCells < 512) block_size = 32;
    unsigned int nblocks  = nCells/block_size + 1;
    gpu_spherical_vertex_model_geometry_kernel<<<nblocks,block_size>>>
        (vertexPos,cellPos,cellNeighbors,vertexCellNeighbors,vertexCellNumberOfNeighbors,
            currentVertexAroundCell,lastVertexAroundCell,nextVertexAroundCell,
            cellNumberOfNeighbors,areaPerimeter,cellNeighborIndex,neighborIndex,sphere,
            nCells);

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

__global__ void  gpu_quadratic_spherical_cellular_force_kernel(dVec *cellPos,
                                            dVec *vertexPos,
                                            dVec *forces,
                                            int *vertexCellNeighbors,
                                            unsigned int *vertexCellNeighborNumber,
                                            dVec *currentVertexAroundCell,
                                            dVec *lastVertexAroundCell,
                                            dVec *nextVertexAroundCell,
                                            unsigned int *cellNumberOfNeighbors,
                                            scalar2 *areaPerimeter,
                                            scalar2 *areaPerimeterPreference,
                                            Index2D neighborIndex,
                                            scalar Kr,
                                            sphericalDomain sphere,
                                            int N)
    {
    // read in the particle that belongs to this thread
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= N)
        return;
    dVec vLast,vCur,vNext,cPos,tempVar;
    dVec f(0.0);
    int vNeighs = vertexCellNeighborNumber[idx];
    for (int cc = 0; cc < vNeighs; ++cc)
        {
        dVec fSet(0.0);
        int vni = neighborIndex(cc,idx);
        int cellIndex = vertexCellNeighbors[vni];
        cPos = cellPos[cellIndex];
        vLast = lastVertexAroundCell[vni];
        vCur = currentVertexAroundCell[vni];
        vNext =nextVertexAroundCell[vni];

        scalar areaDifference = areaPerimeter[cellIndex].x - areaPerimeterPreference[cellIndex].x;
        scalar perimeterDifference = areaPerimeter[cellIndex].y - areaPerimeterPreference[cellIndex].y;
            
        dVec thetaHat, phiHat;
        scalar r0, t0, p0;
        sphere.getAngularCoordinates(vCur,r0,t0,p0);
        sphere.cartesianSphericalBasisChange(t0,p0,thetaHat,phiHat);

        sphere.gradientGeodesicDistance(vCur,vLast,tempVar,thetaHat,phiHat);
        fSet -= 2.0*Kr*perimeterDifference*tempVar;
        sphere.gradientGeodesicDistance(vCur,vNext,tempVar,thetaHat,phiHat);
        fSet -= 2.0*Kr*perimeterDifference*tempVar;

        sphere.gradientTriangleArea(vCur,vLast,cPos,tempVar,thetaHat,phiHat);
        fSet -= 2.0*areaDifference*tempVar;
        sphere.gradientTriangleArea(vCur,cPos,vNext,tempVar,thetaHat,phiHat);
        fSet -= 2.0*areaDifference*tempVar;
            
        if(!isnan(fSet[0]))
            f += fSet;
        };

    forces[idx] = f;
    };


bool gpu_quadratic_spherical_cellular_force(dVec *cellPos,
                                            dVec *vertexPos,
                                            dVec *forces,
                                            int *vertexCellNeighbors,
                                            unsigned int *vertexCellNeighborNumber,
                                            dVec *currentVertexAroundCell,
                                            dVec *lastVertexAroundCell,
                                            dVec *nextVertexAroundCell,
                                            unsigned int *cellNumberOfNeighbors,
                                            scalar2 *areaPerimeter,
                                            scalar2 *areaPerimeterPreference,
                                            Index2D neighborIndex,
                                            scalar Kr,
                                            sphericalDomain &sphere,
                                            int N)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = N/block_size + 1;
    gpu_quadratic_spherical_cellular_force_kernel<<<nblocks,block_size>>>(cellPos,vertexPos,forces,
                vertexCellNeighbors,vertexCellNeighborNumber,currentVertexAroundCell,lastVertexAroundCell,nextVertexAroundCell,
                cellNumberOfNeighbors,areaPerimeter,areaPerimeterPreference,neighborIndex,Kr,sphere,N);

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

__global__ void vm_simple_T1_test_kernel(dVec *d_vertexPositions,
                int *d_vertexNeighbors,
                int *d_vertexEdgeFlips,
                int      *d_vertexCellNeighbors,
                unsigned int      *d_cellVertexNum,
                int      *d_cellVertices,
                sphericalDomain &sphere,
                scalar  T1THRESHOLD,
                int      NvTimes3,
                int      vertexMax,
                int      *d_grow,
                Index2D  &cellNeighborIndex)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= NvTimes3)
        return;
    int vertex1 = idx/3;
    int vertex2 = d_vertexNeighbors[idx];
    scalar arcLength;
    if(vertex1 < vertex2)
        {
        sphere.geodesicDistance(d_vertexPositions[vertex1],d_vertexPositions[vertex2],arcLength);
        if(norm(edge) < T1THRESHOLD)
            {
            d_vertexEdgeFlips[idx]=1;

            //test the number of neighbors of the cells connected to v1 and v2 to see if the
            //cell list should grow. This is kind of slow, and I wish I could optimize it away,
            //or at least not test for it during every time step. The latter seems pretty doable.
            //But this is boring, so we'll revisit if optimizations require it
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex1]] == vertexMax)
                d_grow[0] = 1;
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex1+1]] == vertexMax)
                d_grow[0] = 1;
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex1+2]] == vertexMax)
                d_grow[0] = 1;
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex2]] == vertexMax)
                d_grow[0] = 1;
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex2+1]] == vertexMax)
                d_grow[0] = 1;
            if(d_cellVertexNum[d_vertexCellNeighbors[3*vertex2+2]] == vertexMax)
                d_grow[0] = 1;
            }
        else
            d_vertexEdgeFlips[idx]=0;
        }
    else
        d_vertexEdgeFlips[idx] = 0;
    };

//!Test every edge for a potential T1 event; see if vertexMax needs to increase
bool gpu_vm_test_edges_for_T1(dVec *d_vertexPositions,
                int *d_vertexNeighbors,
                int *d_vertexEdgeFlips,
                int      *d_vertexCellNeighbors,
                unsigned int      *d_cellVertexNum,
                int      *d_cellVertices,
                sphericalDomain &sphere,
                scalar  T1THRESHOLD,
                int      Nvertices,
                int      vertexMax,
                int      *d_grow,
                Index2D  &cellNeighborIndex)
    {
    unsigned int blockSize = 512;
    int nV3 = Nvertices*3;
    if (nV3 < 512) blockSize = 32;
    unsigned int nBlocks = nV3/blockSize + 1;

    vm_simple_T1_test_kernel<<<nBlocks,blockSize>>>(d_vertexPositions,d_vertexNeighbors,
                                                      d_vertexEdgeFlips,d_vertexCellNeighbors,
                                                      d_cellVertexNum,
                                                      Box,T1THRESHOLD,
                                                      nV3,vertexMax,d_grow,cellNeighborIndex);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

bool gpu_vm_parse_multiple_flips(
                    int      *d_vertexEdgeFlips,
                    int      *d_vertexEdgeFlipsCurrent,
                    int      *d_vertexNeighbors,
                    int      *d_vertexCellNeighbors,
                    unsigned int      *d_cellVertexNum,
                    int      *d_cellVertices,
                    int      *d_finishedFlippingEdges,
                    int      *d_edgeFlips,
                    int4     *d_cellSets,
                    Index2D  &cellNeighborIndex,
                    int      Ncells)
    {
    }

bool gpu_vm_flip_edges(
                    int      *d_vertexEdgeFlipsCurrent,
                    dVec *d_vertexPositions,
                    int      *d_vertexNeighbors,
                    int      *d_vertexCellNeighbors,
                    unsigned int      *d_cellVertexNum,
                    int      *d_cellVertices,
                    int      *d_edgeFlips,
                    int4     *d_cellSets,
                    sphericalDomain   &sphere,
                    Index2D  &cellNeighborIndex,
                    int      Nvertices,
                    int      Ncells)
    {
    };
/** @} */ //end of group declaration
