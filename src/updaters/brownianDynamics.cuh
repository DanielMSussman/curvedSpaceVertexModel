#ifndef brownianDynamics_CUH
#define brownianDynamics_CUH

#include "std_include.h"
#include <cuda_runtime.h>
#include "curand_kernel.h"
/*! \file brownianDynamics.cuh */

/** @addtogroup updaterKernels updater Kernels
 * @{
 * \brief CUDA kernels and callers
 */

//!set the vector of displacements from forces and noise
bool gpu_brownian_eom_integration(
                    dVec *forces,
                    dVec *displacements,
                    curandState *RNGs,
                    int N,
                    scalar deltaT,
                    scalar mu,
                    scalar T);

/** @} */ //end of group declaration

#endif
