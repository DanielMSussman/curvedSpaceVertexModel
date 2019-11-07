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

/** @} */ //end of group declaration
