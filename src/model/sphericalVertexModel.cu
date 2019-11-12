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
        sphere.cartesianSphericalBasisChange(vCur,thetaHat,phiHat);

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
                sphericalDomain sphere,
                scalar  T1THRESHOLD,
                int      NvTimes3,
                int      vertexMax,
                int      *d_grow,
                Index2D  cellNeighborIndex)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx >= NvTimes3)
        return;
    int vertex1 = idx/3;
    int vertex2 = d_vertexNeighbors[idx];
    scalar arcLength;
    if(vertex1 < vertex2)
        {
        sphere.geodesicDistance( d_vertexPositions[vertex1], d_vertexPositions[vertex2],arcLength);
        if(arcLength < T1THRESHOLD)
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
                                                      d_cellVertexNum,d_cellVertices,
                                                      sphere,T1THRESHOLD,
                                                      nV3,vertexMax,d_grow,cellNeighborIndex);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

/*!
  There will be severe topology mismatches if a cell is involved in more than one T1 transition
  simultaneously (due to incoherent updates of the cellVertices structure). So, go through the
  current list of edges that are marked to take part in a T1 transition and select one edge per
  cell to be flipped on this trip through the functions.
  */
__global__ void vm_one_T1_per_cell_per_vertex_kernel(
                                        int* __restrict__ d_vertexEdgeFlips,
                                        int* __restrict__ d_vertexEdgeFlipsCurrent,
                                        const int* __restrict__ d_vertexNeighbors,
                                        const int* __restrict__ d_vertexCellNeighbors,
                                        const unsigned int* __restrict__ d_cellVertexNum,
                                        const int * __restrict__ d_cellVertices,
                                        int *d_finishedFlippingEdges,
                                        int *d_cellEdgeFlips,
                                        int4 *d_cellSets,
                                        Index2D cellNeighborIndex,
                                        int Ncells)
    {
    unsigned int cell = blockDim.x * blockIdx.x + threadIdx.x;
    if (cell >= Ncells)
        return;

    //look through every vertex of the cell
    int cneigh = d_cellVertexNum[cell];
    int vertex;
    bool flipFound = false;
    bool moreFlipsFound = false;
    for (int cc = 0; cc < cneigh; ++cc)
        {
        vertex = d_cellVertices[cellNeighborIndex(cc,cell)];
        //what are the other cells attached to this vertex? For correctness, only one cell should
        //own each vertex here. For simplicity, only the lowest-indexed cell gets to do any work.
        int c1,c2,c3,c4;
        c1 = d_vertexCellNeighbors[3*vertex];
        c2 = d_vertexCellNeighbors[3*vertex+1];
        c3 = d_vertexCellNeighbors[3*vertex+2];
        if(c1 < cell || c2 < cell || c3 < cell)
            continue;

        for (int idx = 3*vertex; idx < 3*vertex+3; ++idx)
            {
            if(d_vertexEdgeFlips[idx] == 1)
                {
                int vertex2 = d_vertexNeighbors[idx];
                int ctest;
                for (int ff = 0; ff < 3; ++ff)
                    {
                    ctest = d_vertexCellNeighbors[3*vertex2+ff];
                    if(ctest != c1 && ctest != c2 && ctest != c3)
                    c4=ctest;
                    }

                if (flipFound)
                    {
                    moreFlipsFound = true;
                    break;
                    }
                //check if the cells have been reserved; if not reserve them
                int cc1 = atomicExch(&(d_cellEdgeFlips[c1]),1);
                int cc2 = atomicExch(&(d_cellEdgeFlips[c2]),1);
                int cc3 = atomicExch(&(d_cellEdgeFlips[c3]),1);
                int cc4 = atomicExch(&(d_cellEdgeFlips[c4]),1);
                flipFound = true;
                if(cc1 ==0 && cc2 ==0 &&cc3==0&&cc4==0)
                    {
//                printf("(%i,%i,%i,%i)\t(%i,%i)\n",c1,c2,c3,c4,vertex,vertex2);
                    atomicExch(&d_vertexEdgeFlipsCurrent[idx],1);
                    atomicExch(&d_vertexEdgeFlips[idx],0);
                    int4 cs;cs.x=c1;cs.y=c2;cs.z=c3;cs.w=c4;
                    d_cellSets[idx] = cs;
                    };
                }
            };
        };
    if (flipFound)
        {
        d_finishedFlippingEdges[0] = 1;
        if(moreFlipsFound)
            d_finishedFlippingEdges[1] = 1;
        };
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
    unsigned int block_size = 512;

    /*The issue is that if a cell is involved in two edge flips done by different threads, the resulting
    data structure for what vertices belong to cells and what cells border which vertex will be
    inconsistently updated.

    The strategy will be to take the d_vertexEdgeFlips list, put at most one T1 per cell per vertex into the
    d_vertexEdgeFlipsCurrent list (erasing it from the d_vertexEdgeFlips list), and swap the edges specified
    by the "current" list. If d_vertexEdgeFlips is empty, we will set d_finishedFlippingEdges[0] to 1,
     and if any cell has multiple edges to flip, we set d_finishedFlippingEdges[1] to 1. As long
    as the zeroth entry is 1, the flip edges kernel is called; as long as the first entry is 1 the cpp code will continue calling this gpu_avm_flip_edges function.
    */

    //first select a few edges to flip...
    if(Ncells <512) block_size = 32;
    unsigned int nblocks = Ncells/block_size + 1;
    vm_one_T1_per_cell_per_vertex_kernel<<<nblocks,block_size>>>(
                                                                d_vertexEdgeFlips,
                                                                d_vertexEdgeFlipsCurrent,
                                                                d_vertexNeighbors,
                                                                d_vertexCellNeighbors,
                                                                d_cellVertexNum,
                                                                d_cellVertices,
                                                                d_finishedFlippingEdges,
                                                                d_edgeFlips,
                                                                d_cellSets,
                                                                cellNeighborIndex,
                                                                Ncells);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    }

/*!
  Flip any edge labeled for re-wiring in the vertexEdgeFlipsCurrent list
  */
__global__ void vm_flip_edges_kernel(int *d_vertexEdgeFlipsCurrent,
                    dVec *d_vertexPositions,
                    int      *d_vertexNeighbors,
                    int      *d_vertexCellNeighbors,
                    unsigned int      *d_cellVertexNum,
                    int      *d_cellVertices,
                    int      *d_cellEdgeFlips,
                    int4     *d_cellSets,
                    sphericalDomain   sphere,
                    Index2D  cellNeighborIndex,
                    int      NvTimes3
                    )
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    //return if the index is out of bounds or if the edge isn't marked for flipping
    if (idx >= NvTimes3 || d_vertexEdgeFlipsCurrent[idx] == 0)
        return;
    //identify the vertices and reset the flag
    int vertex1 = idx/3;
    int vertex2 = d_vertexNeighbors[idx];
    d_vertexEdgeFlipsCurrent[idx] = 0;

    //first, identify the cell and vertex set involved...
    int4 cellSet;cellSet.x=-1;cellSet.y=-1;cellSet.z=-1;cellSet.w=-1;
    //int4 vertexSet;
    int2 vertexSet;//vertexSet.x = "b", vertexSet.y = "a"
    /*
    The following is fairly terrible GPU code, and should be considered for refactoring
    */
    int4 cells = d_cellSets[idx];
    int cell1,cell2,cell3;
    int vlast, vcur, vnext, cneigh;
    cell1 = cells.x;
    cell2 = cells.y;
    cell3 = cells.z;
    cellSet.w = cells.w;

    //classify cell1
    cneigh = d_cellVertexNum[cell1];
    vlast = d_cellVertices[ cellNeighborIndex(cneigh-2,cell1) ];
    vcur = d_cellVertices[ cellNeighborIndex(cneigh-1,cell1) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell1)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell1;
    else if(vnext == vertex2)
        cellSet.z = cell1;
    else
        {
        cellSet.y = cell1;
        };

    //classify cell2
    cneigh = d_cellVertexNum[cell2];
    vlast = d_cellVertices[ cellNeighborIndex(cneigh-2,cell2) ];
    vcur = d_cellVertices[ cellNeighborIndex(cneigh-1,cell2) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell2)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell2;
    else if(vnext == vertex2)
        cellSet.z = cell2;
    else
        {
        cellSet.y = cell2;
        };
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell1)];
        }
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell2)];
        }
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell3)];
        }

    //classify cell3
    cneigh = d_cellVertexNum[cell3];
    vlast = d_cellVertices[ cellNeighborIndex(cneigh-2,cell3) ];
    vcur = d_cellVertices[ cellNeighborIndex(cneigh-1,cell3) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cell3)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    if(vlast == vertex2)
        cellSet.x = cell3;
    else if(vnext == vertex2)
        cellSet.z = cell3;
    else
        {
        cellSet.y = cell3;
        };

    //get the vertexSet by examining cells j and l
    cneigh = d_cellVertexNum[cellSet.y];
    vlast = d_cellVertices[ cellNeighborIndex(cneigh-2,cellSet.y) ];
    vcur = d_cellVertices[ cellNeighborIndex(cneigh-1,cellSet.y) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cellSet.y)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.x=vnext;
    cneigh = d_cellVertexNum[cellSet.w];
    vlast = d_cellVertices[ cellNeighborIndex(cneigh-2,cellSet.w) ];
    vcur = d_cellVertices[ cellNeighborIndex(cneigh-1,cellSet.w) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = d_cellVertices[cellNeighborIndex(cn,cellSet.w)];
        if(vcur == vertex2) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.y=vnext;

    d_cellEdgeFlips[cells.x] = 0;
    d_cellEdgeFlips[cells.y] = 0;
    d_cellEdgeFlips[cells.z] = 0;
    d_cellEdgeFlips[cells.w] = 0;

    //forbid a T1 transition that would shrink a triangular cell
    if (d_cellVertexNum[cellSet.x] ==3 || d_cellVertexNum[cellSet.z] ==3)
        return;
    if(cellSet.x <0 || cellSet.y < 0 || cellSet.z <0 || cellSet.w <0)
        return;

    //okay, we're ready to go. First, rotate the vertices in the edge and set them at twice their original distance
    dVec v1 = d_vertexPositions[vertex1];
    dVec v2 = d_vertexPositions[vertex2];
    dVec midpoint = 0.5*(v1+v2);
    sphere.putInBoxVirtual(midpoint);
    //chose the angle of rotation based on whether the edges are currently crossed...
    dVec vC = d_vertexPositions[vertexSet.y];//vSet.y is vSet.z

    scalar determinant = vC[0]*(v1[1]*v2[2]-v1[2]*v2[1])
                        +vC[1]*(v1[2]*v2[0]-v1[0]*v2[2])
                        +vC[2]*(v1[0]*v2[1]-v1[1]*v2[0]);
    determinant = determinant > 0 ? 1. : -1. ;

    rodriguesRotation(v1,midpoint,-0.5*determinant*PI);
    rodriguesRotation(v2,midpoint,-0.5*determinant*PI);

    dVec diff = 0.5*(v1-v2);
    v1 = v1 + diff;
    v2 = v2 - diff;
    sphere.putInBoxReal(v1);
    sphere.putInBoxReal(v2);

    d_vertexPositions[vertex1] = v1;
    d_vertexPositions[vertex2] = v2;

    //now, re-wire the cells and vertices
    //start with the vertex-vertex and vertex-cell  neighbors
    for (int vert = 0; vert < 3; ++vert)
        {
        //vertex-cell neighbors
        if(d_vertexCellNeighbors[3*vertex1+vert] == cellSet.z)
            d_vertexCellNeighbors[3*vertex1+vert] = cellSet.w;
        if(d_vertexCellNeighbors[3*vertex2+vert] == cellSet.x)
            d_vertexCellNeighbors[3*vertex2+vert] = cellSet.y;

        //vertex-vertex neighbors
        if(d_vertexNeighbors[3*vertex1+vert] == vertexSet.x)
            d_vertexNeighbors[3*vertex1+vert] = vertexSet.y;
        if(d_vertexNeighbors[3*vertex2+vert] == vertexSet.y)
            d_vertexNeighbors[3*vertex2+vert] = vertexSet.x;

        if(d_vertexNeighbors[3*vertexSet.x+vert] == vertex1)
            d_vertexNeighbors[3*vertexSet.x+vert] = vertex2;
        if(d_vertexNeighbors[3*vertexSet.y+vert] == vertex2)
            d_vertexNeighbors[3*vertexSet.y+vert] = vertex1;
        };
    //now rewire the cells...
    //cell i loses v2 as a neighbor
    cneigh = d_cellVertexNum[cellSet.x];
    int cidx = 0;
    for (int cc = 0; cc < cneigh-1; ++cc)
        {
        if(d_cellVertices[cellNeighborIndex(cc,cellSet.x)] == vertex2)
            cidx +=1;
        d_cellVertices[cellNeighborIndex(cc,cellSet.x)] = d_cellVertices[cellNeighborIndex(cidx,cellSet.x)];
        cidx +=1;
        };
    d_cellVertexNum[cellSet.x] -= 1;

    //cell j gains v2 in between v1 and b, so step through list backwards and insert
    cneigh = d_cellVertexNum[cellSet.y];
    bool found0 = false;
    for (int cc = cneigh-1; cc >= 0; --cc)
        {
        int cellIndex = d_cellVertices[cellNeighborIndex(cc,cellSet.y)];
        if(!found0)
            d_cellVertices[cellNeighborIndex(cc+1,cellSet.y)] = cellIndex;
        if(cellIndex == vertexSet.x)
            {
            found0 = true;
            d_cellVertices[cellNeighborIndex(cc,cellSet.y)] = vertex2;
            }
        }
    d_cellVertexNum[cellSet.y] += 1;

    //cell k loses v1 as a neighbor
    cneigh = d_cellVertexNum[cellSet.z];
    cidx = 0;
    for (int cc = 0; cc < cneigh-1; ++cc)
        {
        if(d_cellVertices[cellNeighborIndex(cc,cellSet.z)] == vertex1)
            cidx +=1;
        d_cellVertices[cellNeighborIndex(cc,cellSet.z)] = d_cellVertices[cellNeighborIndex(cidx,cellSet.z)];
        cidx +=1;
        };
    d_cellVertexNum[cellSet.z] -= 1;

    //cell l gains v1 in between v2 and a...copy the logic of cell j
    cneigh = d_cellVertexNum[cellSet.w];
    bool found = false;
    for (int cc = cneigh-1; cc >= 0; --cc)
        {
        int cellIndex = d_cellVertices[cellNeighborIndex(cc,cellSet.w)];
        if(!found)
            d_cellVertices[cellNeighborIndex(cc+1,cellSet.w)] = cellIndex;
        if(cellIndex == vertexSet.y)
            {
            found = true;
            d_cellVertices[cellNeighborIndex(cc,cellSet.w)] = vertex1;
            }
        }
    d_cellVertexNum[cellSet.w] += 1;
    }

bool gpu_vm_flip_edges(
                    int      *d_vertexEdgeFlipsCurrent,
                    dVec *d_vertexPositions,
                    int      *d_vertexNeighbors,
                    int      *d_vertexCellNeighbors,
                    unsigned int      *d_cellVertexNum,
                    int      *d_cellVertices,
                    int      *d_cellEdgeFlips,
                    int4     *d_cellSets,
                    sphericalDomain   &sphere,
                    Index2D  &cellNeighborIndex,
                    int      Nvertices,
                    int      Ncells)
    {
    unsigned int block_size = 128;
    int NvTimes3 = Nvertices*3;
    if (NvTimes3 < 128) block_size = 32;
    unsigned int nblocks  = NvTimes3/block_size + 1;

    vm_flip_edges_kernel<<<nblocks,block_size>>>(
                                                  d_vertexEdgeFlipsCurrent,d_vertexPositions,d_vertexNeighbors,
                                                  d_vertexCellNeighbors,d_cellVertexNum,d_cellVertices,d_cellEdgeFlips,d_cellSets,
                                                  sphere,
                                                  cellNeighborIndex,NvTimes3);

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };
/** @} */ //end of group declaration
