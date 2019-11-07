#include "brownianDynamics.cuh"

/** \file brownianDynamics.cu
    * Defines kernel callers and kernels for GPU calculations of overdamped brownian dynamics
*/

/*!
    \addtogroup updaterKernels
    @{
*/

/*!
Each thread calculates the displacement of an individual cell
*/
__global__ void brownian_eom_integration_kernel(dVec *forces,
                                           dVec *displacements,
                                           curandState *RNGs,
                                           int N,
                                           scalar noisePrefactor,
                                           scalar forcePrefactor)
    {
    unsigned int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if (idx >=N)
        return;
    curandState_t randState;

    randState=RNGs[idx];
    for (int dd = 0; dd < DIMENSION; ++dd)
        displacements[idx][dd] = noisePrefactor*cur_norm(&randState) + forcePrefactor*forces[idx][dd];

    RNGs[idx] = randState;
    return;
    };

//!get the current timesteps vector of displacements into the displacement vector
bool gpu_brownian_eom_integration(
                    dVec *forces,
                    dVec *displacements,
                    curandState *RNGs,
                    int N,
                    scalar deltaT,
                    scalar mu,
                    scalar T)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = N/block_size + 1;

    scalar forcePrefactor = deltaT*mu;
    scalar noisePrefactor = sqrt(2.0*forcePrefactor*T);

    brownian_eom_integration_kernel<<<nblocks,block_size>>>(
                                forces,displacements,
                                RNGs,
                                N,noisePrefactor,forcePrefactor);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };

/** @} */ //end of group declaration

