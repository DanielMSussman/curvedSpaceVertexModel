#ifndef sphericalVertexModel_CUH
#define sphericalVertexModel_CUH

#include "std_include.h"
#include "sphericalDomain.h"
#include <cuda_runtime.h>

/** @addtogroup modelKernels model Kernels
 * @{
 * \brief CUDA kernels and callers
 */

bool gpu_move_particles_on_sphere(dVec *pos,
                                  dVec *disp,
                                  sphericalDomain &sphere,
                                  scalar scale,
                                  int N
                                  );

/** @} */ //end of group declaration
#endif
