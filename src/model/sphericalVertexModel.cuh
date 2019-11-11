#ifndef sphericalVertexModel_CUH
#define sphericalVertexModel_CUH

#include "std_include.h"
#include "sphericalDomain.h"
#include "indexer.h"
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
                                         );

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
                                            int N);

bool gpu_vm_test_edges_for_T1(
                    dVec *d_vertexPositions,
                    int      *d_vertexNeighbors,
                    int      *d_vertexEdgeFlips,
                    int      *d_vertexCellNeighbors,
                    unsigned int      *d_cellVertexNum,
                    int      *d_cellVertices,
                    sphericalDomain &sphere,
                    scalar  T1THRESHOLD,
                    int      Nvertices,
                    int      vertexMax,
                    int      *d_grow,
                    Index2D  &cellNeighborIndex);

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
                    int      Ncells);

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
                    int      Ncells);
/** @} */ //end of group declaration
#endif
