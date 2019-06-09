#ifndef vectorialVicsek_CUH
#define vectorialVicsek_CUH

#include "std_include.h"
#include "indexer.h"
#include "noiseSource.h"
#include <cuda_runtime.h>
/*! \file vectorialVicsek.cuh */

/** @addtogroup updaterKernels updater Kernels
 * @{
 * \brief CUDA kernels and callers
 */

bool gpu_selfPropelledParticleDisplacement(dVec *force,
                                           dVec *directors,
                                           dVec *displacements,
                                           scalar v0,
                                           scalar mu,
                                           scalar deltaT,
                                           int N);

bool gpu_vectorVicsek_target_directors(dVec *director,
                                       dVec *disp,
                                       unsigned int *nNeighs,
                                       int *neighs,
                                       curandState *rngs,
                                       Index2D &neighborIndex,
                                       scalar eta,
                                       int N);

bool gpu_vectorVicsek_update_directors(dVec *director,
                                       dVec *disp,
                                       scalar tau,
                                       scalar deltaT,
                                       int N);

/** @} */ //end of group declaration
#endif
