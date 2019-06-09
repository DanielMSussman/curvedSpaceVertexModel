#include "vectorialVicsek.cuh"
/*! \file vectorialVicsek.cu
\addtogroup updaterKernels
@{
*/

__global__ void gpu_spp_displacement_kernel(dVec *force,
                                           dVec *directors,
                                           dVec *displacements,
                                           scalar v0,
                                           scalar mu,
                                           scalar deltaT,
                                           int N)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx >= N)
        return;
    displacements[idx] = deltaT*(v0*directors[idx]+mu*force[idx]);
    };

bool gpu_selfPropelledParticleDisplacement(dVec *force,
                                           dVec *directors,
                                           dVec *displacements,
                                           scalar v0,
                                           scalar mu,
                                           scalar deltaT,
                                           int N)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = (N)/block_size + 1;
    gpu_spp_displacement_kernel<<<nblocks,block_size>>>(force,directors,displacements,v0,mu,deltaT,N);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    }

__global__ void gpu_vectorVicsek_update_directors_kernel(dVec *director,
                                           dVec *disp,
                                           scalar tau,
                                           scalar deltaT,
                                           int N)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx >= N)
        return;
    director[idx] = disp[idx];
    }
bool gpu_vectorVicsek_update_directors(dVec *director,
                                       dVec *disp,
                                       scalar tau,
                                       scalar deltaT,
                                       int N)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = (N)/block_size + 1;
    gpu_vectorVicsek_update_directors_kernel<<<nblocks,block_size>>>(director,disp,tau,deltaT,N);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    }

__global__ void gpu_vectorVicsek_target_directors_kernel(dVec *director,
                                           dVec *disp,
                                           unsigned int *nNeighs,
                                           int *neighs,
                                           curandState *rngs,
                                           Index2D neighborIndex,
                                           scalar eta,
                                           int N)
    {
    unsigned int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if(idx >= N)
        return;
    curandState_t randState;
    randState = rngs[idx];
    dVec spherePoint;
#if DIMENSION == 3
    scalar u = curand_uniform_double(&randState);
    scalar w = curand_uniform_double(&randState);
    scalar phi = 2.0*PI*u;
    scalar theta = acos(2.0*w-1);
    spherePoint.x[0] = 1.0*sin(theta)*cos(phi);
    spherePoint.x[1] = 1.0*sin(theta)*sin(phi);
    spherePoint.x[2] = 1.0*cos(theta);
#else
    scalar u = curand_uniform_double(&randState);
    scalar phi = 2.0*PI*u;
    spherePoint.x[0] = cos(phi);
    spherePoint.x[1] = sin(phi);
#endif

    int m = nNeighs[idx];
    disp[idx] = director[idx];
    for (int jj = 0; jj < m; ++jj)
        disp[idx] += director[neighs[neighborIndex(jj,idx)]];
    scalar mi = 1./(m+1);
    disp[idx] = disp[idx]*(mi) + spherePoint*eta;

    rngs[idx] = randState;
    }

bool gpu_vectorVicsek_target_directors(dVec *director,
                                       dVec *disp,
                                       unsigned int *nNeighs,
                                       int *neighs,
                                       curandState *rngs,
                                       Index2D &neighborIndex,
                                       scalar eta,
                                       int N)
    {
    unsigned int block_size = 512;
    if (N < 512) block_size = 32;
    unsigned int nblocks  = (N)/block_size + 1;
    gpu_vectorVicsek_target_directors_kernel<<<nblocks,block_size>>>(director,disp,nNeighs,neighs,rngs,neighborIndex,eta,N);
    HANDLE_ERROR(cudaGetLastError());
    return cudaSuccess;
    };
/** @} */ //end of group declaration
