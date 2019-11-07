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
                                            int N)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = N/block_size + 1;
    //gpu_move_particles_on_sphere_kernel<<<nblocks,block_size>>>(pos,disp,sphere,scale,N);
    UNWRITTENCODE("Force ON GPU");

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
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
                                         int nCells
                                         )
    {
    unsigned int block_size = 512;
    if (nCells < 512) block_size = 32;
    unsigned int nblocks  = nCells/block_size + 1;
    //gpu_move_particles_on_sphere_kernel<<<nblocks,block_size>>>(pos,disp,sphere,scale,N);
    UNWRITTENCODE("GEOMETRY ON GPU");

    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };
/** @} */ //end of group declaration
