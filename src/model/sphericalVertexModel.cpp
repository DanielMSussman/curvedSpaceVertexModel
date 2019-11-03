#include "sphericalVertexModel.h"

sphericalVertexModel::sphericalVertexModel(int n, noiseSource &_noise, scalar _area, scalar _perimeter, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
    {
    selfForceCompute =true;
    //the initialization sets n cell positions randomly
    //use the convex huller to determine an initial delaunay triangulation on the sphere
    cout << "extracting vertex positions from convex huller" << endl;
    nCells = n;
    setRadius(sqrt(nCells/(4.0*PI)));
    t1Threshold = 0.1;
    cellPositions = positions;
        {
        ArrayHandle<dVec> cellPos(cellPositions);
        convexHuller.sphericalConvexHullForVertexModel(cellPos.data,N,cellNeighbors,cellNumberOfNeighbors,cellNeighborIndex,positions,neighbors,vertexCellNeighbors,numberOfNeighbors,neighborIndex);
        };
    maximumVerticesPerCell = cellNeighborIndex.getW();

    int cnnArraySize = vertexCellNeighbors.getNumElements();
    currentVertexAroundCell.resize(cnnArraySize);
    lastVertexAroundCell.resize(cnnArraySize);
    nextVertexAroundCell.resize(cnnArraySize);
    int nVertices = positions.getNumElements();
    maxVNeighs = cnnArraySize / nVertices;

    printf("initialized a system with %i cells and %i vertices\n",nCells, nVertices);


    N=nVertices;
    //here is the set of data structures to be resized
    velocities.resize(nVertices);
    directors.resize(nVertices);
    forces.resize(nVertices);
    masses.resize(nVertices);
    radii.resize(nVertices);
    types.resize(nVertices);
    vector<dVec> zeroes(nVertices,make_dVec(0.0));
    vector<dVec> dirs(nVertices,make_dVec(0.577));
    vector<scalar> ones(nVertices,1.0);
    //vector<scalar> halves(nVertices,.5);
    vector<int> units(N,0);
    fillGPUArrayWithVector(units,types);
    fillGPUArrayWithVector(zeroes,velocities);
    fillGPUArrayWithVector(zeroes,forces);
    fillGPUArrayWithVector(ones,masses);
    fillGPUArrayWithVector(zeroes,directors);
    //fillGPUArrayWithVector(halves,radii);
    
    areaPerimeter.resize(nCells);
    areaPerimeterPreference.resize(nCells);

    printf("(a0,p0)=%f\t%f\n",_area,_perimeter);
    setPreferredParameters(_area,_perimeter);
    setScalarModelParameter(1.0);
    computeGeometry();
    preserveOrientatedFaces();
    };

void sphericalVertexModel::preserveOrientatedFaces()
    {
    /*
    ArrayHandle<int> cvn(cellNeighbors);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    scalar totalArea = 0;
    scalar totalPerimeter = 0.;

    for (int cc = 0; cc < nCells; ++cc)
        {
        int neighs = cnn.data[cc];
        dVec CP = cp.data[cc];
        for (int nn = 0; nn < neighs; ++nn)
            {
            int cni = cellNeighborIndex(nn,cc);
            dVec pt1 = lastVert.data[cni];
            dVec pt2 = curVert.data[cni];
            dVec pt3 = nextVert.data[cni];
            scalar determinant = pt1[0]*(pt2[1]*pt3[2]-pt2[2]*pt3[1])
                        +pt1[1]*(pt2[2]*pt3[0]-pt2[0]*pt3[2])
                        +pt1[2]*(pt2[0]*pt3[1]-pt2[1]*pt3[0]);
            pt1 = CP;
            scalar determinant2 = pt1[0]*(pt2[1]*pt3[2]-pt2[2]*pt3[1])
                        +pt1[1]*(pt2[2]*pt3[0]-pt2[0]*pt3[2])
                        +pt1[2]*(pt2[0]*pt3[1]-pt2[1]*pt3[0]);
           // printf("%i\t%i\t%i\t%i\n",cc,nn,determinant > 0 ? 1:-1, determinant2 > 0 ? 1:-1);
            };
        };
    */
    }
    
void sphericalVertexModel::setPreferredParameters(scalar _a0, scalar _p0)
    {
    ArrayHandle<scalar2> app(areaPerimeterPreference);
    scalar2 prefs; prefs.x = _a0; prefs.y = _p0;
    for (int cc = 0; cc < nCells; ++cc)
        app.data[cc] = prefs;
    }

void sphericalVertexModel::moveParticles(GPUArray<dVec> &displacements, scalar scale)
    {
    if(scale == 1.)
        {
        ArrayHandle<dVec> p(positions);
        ArrayHandle<dVec> V(displacements);
        for(int ii = 0; ii < N; ++ii)
            {
            sphere.move(p.data[ii],V.data[ii]);
        //    sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
            }
        }
    else
        {
        ArrayHandle<dVec> p(positions);
        ArrayHandle<dVec> V(displacements);
        for(int ii = 0; ii < N; ++ii)
            {
            sphere.move(p.data[ii],scale*V.data[ii]);
        //    sphere.projectToTangentPlaneAndNormalize(n.data[ii],p.data[ii]);
            }
        }
    if(!restrictedMotion)
        enforceTopology();
    };

scalar sphericalVertexModel::computeEnergy()
    {
    ArrayHandle<scalar2> ap(areaPerimeter);
    ArrayHandle<scalar2> app(areaPerimeterPreference);
    energy = 0.0;
    for (int i = 0; i < nCells; ++i)
        {
        energy += (ap.data[i].x-app.data[i].x)*(ap.data[i].x-app.data[i].x) + Kr*(ap.data[i].y-app.data[i].y)*(ap.data[i].y-app.data[i].y);
        }
    //printf("current energy = %g\n",energy);
    return energy;
    }

void sphericalVertexModel::recomputeAreasCPU()
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<int> cvn(cellNeighbors);
    ArrayHandle<int> vcn(vertexCellNeighbors);
    ArrayHandle<unsigned int> vcnn(numberOfNeighbors);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<scalar2> ap(areaPerimeter);
    scalar totalArea = 0;
    scalar totalPerimeter = 0.;
    for (int cc = 0; cc < nCells; ++cc)
        {
        int neighs = cnn.data[cc];
        dVec cellPos = cp.data[cc];;

        int lastVertexIdx = cvn.data[cellNeighborIndex(neighs-2,cc)];
        int curVertexIdx = cvn.data[cellNeighborIndex(neighs-1,cc)];
        dVec lastVertexPos = p.data[lastVertexIdx];
        dVec curVertexPos = p.data[curVertexIdx];
        int nextVertexIdx;
        dVec nextVertexPos;
        scalar perimeter = 0;
        scalar area = 0;
        scalar tempVal;

        for (int nn = 0; nn < neighs; ++nn)
            {
            int cni = cellNeighborIndex(nn,cc);
            nextVertexIdx = cvn.data[cni];
            nextVertexPos = p.data[nextVertexIdx];
            sphere.sphericalTriangleArea(cellPos,lastVertexPos,curVertexPos,tempVal);
            area +=tempVal;

            lastVertexPos = curVertexPos;
            curVertexIdx = nextVertexIdx;
            curVertexPos = nextVertexPos;
            }
        ap.data[cc].x = area;
        totalArea += area;
        }
 //       printf("recomputed area = %f \n",totalArea);
    }

void sphericalVertexModel::computeGeometryCPU()
    {
    scalar excessArea;
    scalar totalArea = 0;
    scalar totalPerimeter = 0.;
    {//arrayHandle scope
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<int> cvn(cellNeighbors);
    ArrayHandle<int> vcn(vertexCellNeighbors);
    ArrayHandle<unsigned int> vcnn(numberOfNeighbors);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<scalar2> ap(areaPerimeter);
    for (int cc = 0; cc < nCells; ++cc)
        {
        int neighs = cnn.data[cc];
        dVec cellPos(0.0);
        for (int nn = 0; nn < neighs; ++nn)
            {
            cellPos = cellPos + p.data[cvn.data[cellNeighborIndex(nn,cc)]];
            }
        sphere.putInBoxReal(cellPos);
        cp.data[cc]=cellPos;

        int lastVertexIdx = cvn.data[cellNeighborIndex(neighs-2,cc)];
        int curVertexIdx = cvn.data[cellNeighborIndex(neighs-1,cc)];
        dVec lastVertexPos = p.data[lastVertexIdx];
        dVec curVertexPos = p.data[curVertexIdx];
        int nextVertexIdx;
        dVec nextVertexPos;
        scalar perimeter = 0;
        scalar area = 0;
        scalar tempVal;

        for (int nn = 0; nn < neighs; ++nn)
            {
            int cni = cellNeighborIndex(nn,cc);
            //what index is the cell in the vertex-cell neighbor list?
            int vNeighs = vcnn.data[curVertexIdx];
            int forceSetIdx = -1;
            for (int vn = 0; vn < vNeighs; ++vn)
                {
                int newIdx = neighborIndex(vn,curVertexIdx);
                if(vcn.data[newIdx] == cc)
                    forceSetIdx = newIdx;
                }
            nextVertexIdx = cvn.data[cni];
            nextVertexPos = p.data[nextVertexIdx];
            sphere.geodesicDistance(lastVertexPos,curVertexPos,tempVal);
            perimeter += tempVal;
            sphere.includedAngle(lastVertexPos,curVertexPos,nextVertexPos,tempVal);
            area +=tempVal;

            curVert.data[forceSetIdx] = curVertexPos;
            lastVert.data[forceSetIdx] = lastVertexPos;
            nextVert.data[forceSetIdx] = nextVertexPos;

            lastVertexPos = curVertexPos;
            curVertexIdx = nextVertexIdx;
            curVertexPos = nextVertexPos;
            }
        area = (area-(neighs-2)*PI);
        int extraAngularArea = floor(area/(1.0*PI));
        if(extraAngularArea > 0)
            area -= extraAngularArea*PI;
        area = area * sphere.radius*sphere.radius; 
        ap.data[cc].x = area;
        ap.data[cc].y = perimeter;
        totalArea += area;
        totalPerimeter += perimeter;
        }
        excessArea = totalArea - 4.0*PI*sphere.radius*sphere.radius;
        }//arrayHandle scope
        if(fabs(excessArea)> 1e-6)
            {
 //           printf("excess area = %g\t total peri = %f ..recomputing area via triangles: ",excessArea,totalPerimeter);
            recomputeAreasCPU();
            }
    }

void sphericalVertexModel::computeGeometryGPU()
    {
    }
void sphericalVertexModel::computeForceGPU()
    {
    }
void sphericalVertexModel::computeForceCPU()
    {
    //printf("computing forces\n");
    computeGeometry();
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> force(forces);
    ArrayHandle<int> vcn(vertexCellNeighbors);
    ArrayHandle<unsigned int> vcnn(numberOfNeighbors);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<scalar2> ap(areaPerimeter);
    ArrayHandle<scalar2> app(areaPerimeterPreference);
    
    scalar forceNorm = 0.0;
    dVec vLast,vCur,vNext,cPos,tempVar;
    dVec meanForce(0.0);
    for (int vertexIndex = 0; vertexIndex < N; ++vertexIndex)
        {
        dVec f(0.0);
        int vNeighs = vcnn.data[vertexIndex];
        for (int cc = 0; cc < vNeighs; ++cc)
            {
            int cellIndex = vcn.data[neighborIndex(cc,vertexIndex)];
            cPos = cp.data[cellIndex];
            vLast = lastVert.data[neighborIndex(cc,vertexIndex)];
            vCur = curVert.data[neighborIndex(cc,vertexIndex)];
            vNext = nextVert.data[neighborIndex(cc,vertexIndex)];
            scalar areaDifference = ap.data[cellIndex].x - app.data[cellIndex].x;
            scalar perimeterDifference = ap.data[cellIndex].y - app.data[cellIndex].y;

            sphere.gradientGeodesicDistance(vCur,vLast,tempVar);
            f -= 2.0*Kr*perimeterDifference*tempVar;
            sphere.gradientGeodesicDistance(vCur,vNext,tempVar);
            f -= 2.0*Kr*perimeterDifference*tempVar;
            sphere.gradientTriangleArea(vCur,vLast,cPos,tempVar);
            f -= 2.0*areaDifference*tempVar;
            sphere.gradientTriangleArea(vCur,cPos,vNext,tempVar);
            f -= 2.0*areaDifference*tempVar;
//if(isnan(tempVar[0])) {printf("area last nan %f\t (%f,%f,%f), (%f,%f,%f), (%f,%f,%f) \n",areaDifference, vCur[0],vCur[1],vCur[2],vLast[0],vLast[1],vLast[2],cPos[0],cPos[1],cPos[2]);}

            //sphere.gradientTriangleArea(vCur,cPos,vNext,tempVar);
            //f -= 2.0*areaDifference*tempVar;
//if(isnan(tempVar[0])) {printf("area next nan %f\t (%f,%f,%f), (%f,%f,%f), (%f,%f,%f) \n",areaDifference, vCur[0],vCur[1],vCur[2],vNext[0],vNext[1],vNext[2],cPos[0],cPos[1],cPos[2]);}

            };
        //printf("%f,%f,%f\n",f[0],f[1],f[2]);
        if(isnan(f[0]))
            {
            f = dVec(0);
            printf("forceNan on vidx %i\n",vertexIndex);
            }
        force.data[vertexIndex] = f;
        meanForce = meanForce + f;
        forceNorm += dot(f,f);
        //printf("vertex %i, force (%f,%f,%f)\n",vertexIndex, f[0],f[1],f[2]);
        };

//    printf("total force norm =  %g, mean force = (%f,%f,%f)\n",forceNorm/N,meanForce[0]/N,meanForce[1]/N,meanForce[2]/N);
//    getMeanForce(meanForce);
//    printf("projected mean force = (%f,%f,%f)\n",meanForce[0],meanForce[1],meanForce[2]);
    }

void sphericalVertexModel::enforceTopology()
    {
        //see if vertex motion leads to T1 transitions
    if(useGPU)
        {
        testAndPerformT1TransitionsGPU();
        }
    else
        {
        testAndPerformT1TransitionsCPU();
        };
    }

/*!
Test whether a T1 needs to be performed on any edge by simply checking if the edge length is beneath a threshold.
This function also performs the transition and maintains the auxiliary data structures
 */
void sphericalVertexModel::testAndPerformT1TransitionsCPU()
    {
    ArrayHandle<dVec> h_v(positions,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_vn(neighbors,access_location::host,access_mode::readwrite);
    ArrayHandle<unsigned int> h_cvn(cellNumberOfNeighbors,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_cv(cellNeighbors,access_location::host,access_mode::readwrite);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::readwrite);

    dVec edge;
    //first, scan through the list for any T1 transitions...
    int vertex2;
    //keep track of whether vertexMax needs to be increased
    int vMax = maximumVerticesPerCell;
    /*
     The following is the convention:
     cell i: contains both vertex 1 and vertex 2, in CW order
     cell j: contains only vertex 1
     cell k: contains both vertex 1 and vertex 2, in CCW order
     cell l: contains only vertex 2
     */
    int4 cellSet;
    /*
    vertexSet (a,b,c,d) have those indices in which before the transition
    cell i has CCW vertices: ..., c, v2, v1, a, ...
    and
    cell k has CCW vertices: ..., b,v1,v2,d, ...
    */
    int4 vertexSet;
    dVec v1,v2;
    scalar arcLength;
    for (int vertex1 = 0; vertex1 < N; ++vertex1)
        {
        v1 = h_v.data[vertex1];
        //look at vertexNeighbors list for neighbors of vertex, compute edge length
        for (int vv = 0; vv < 3; ++vv)
            {
            vertex2 = h_vn.data[3*vertex1+vv];
            //only look at each pair once
            if(vertex1 < vertex2)
                {
                v2 = h_v.data[vertex2];
                sphere.geodesicDistance(v1,v2,arcLength);
                if(arcLength < t1Threshold)
                    {
                    bool growCellVertexList = false;
                    getCellVertexSetForT1(vertex1,vertex2,cellSet,vertexSet,growCellVertexList);
                    //forbid a T1 transition that would shrink a triangular cell
                    if( h_cvn.data[cellSet.x] == 3 || h_cvn.data[cellSet.z] == 3)
                        continue;
                    //Does the cell-vertex-neighbor data structure need to be bigger?
                    if(growCellVertexList)
                        {
                        vMax +=1;
                        growCellVerticesList(vMax);
                        h_cv = ArrayHandle<int>(cellNeighbors,access_location::host,access_mode::readwrite);
                        };

                    //Rotate the vertices in the edge and set them at twice their original distance
                    dVec midpoint = 0.5*(v1+v2);
                    sphere.putInBoxVirtual(midpoint);
                    rodriguesRotation(v1,midpoint,-0.5*PI);
                    rodriguesRotation(v2,midpoint,-0.5*PI);

                    dVec diff = 0.5*(v1-v2);
                    v1 = v1 + diff;
                    v2 = v2 - diff;
                    sphere.putInBoxReal(v1);
                    sphere.putInBoxReal(v2);
                    h_v.data[vertex1] = v1;
                    h_v.data[vertex2] = v2;

                    //re-wire the cells and vertices
                    //start with the vertex-vertex and vertex-cell  neighbors
                    for (int vert = 0; vert < 3; ++vert)
                        {
                        //vertex-cell neighbors
                        if(h_vcn.data[3*vertex1+vert] == cellSet.z)
                            h_vcn.data[3*vertex1+vert] = cellSet.w;
                        if(h_vcn.data[3*vertex2+vert] == cellSet.x)
                            h_vcn.data[3*vertex2+vert] = cellSet.y;
                        //vertex-vertex neighbors
                        if(h_vn.data[3*vertexSet.y+vert] == vertex1)
                            h_vn.data[3*vertexSet.y+vert] = vertex2;
                        if(h_vn.data[3*vertexSet.z+vert] == vertex2)
                            h_vn.data[3*vertexSet.z+vert] = vertex1;
                        if(h_vn.data[3*vertex1+vert] == vertexSet.y)
                            h_vn.data[3*vertex1+vert] = vertexSet.z;
                        if(h_vn.data[3*vertex2+vert] == vertexSet.z)
                            h_vn.data[3*vertex2+vert] = vertexSet.y;
                        };
                    //now rewire the cells
                    //cell i loses v2 as a neighbor
                    int cneigh = h_cvn.data[cellSet.x];
                    int cidx = 0;
                    for (int cc = 0; cc < cneigh-1; ++cc)
                        {
                        if(h_cv.data[cellNeighborIndex(cc,cellSet.x)] == vertex2)
                            cidx +=1;
                        h_cv.data[cellNeighborIndex(cc,cellSet.x)] = h_cv.data[cellNeighborIndex(cidx,cellSet.x)];
                        cidx +=1;
                        };
                    h_cvn.data[cellSet.x] -= 1;

                    //cell j gains v2 in between v1 and b
                    cneigh = h_cvn.data[cellSet.y];
                    vector<int> cvcopy1(cneigh+1);
                    cidx = 0;
                    for (int cc = 0; cc < cneigh; ++cc)
                        {
                        int cellIndex = h_cv.data[cellNeighborIndex(cc,cellSet.y)];
                        cvcopy1[cidx] = cellIndex;
                        cidx +=1;
                        if(cellIndex == vertex1)
                            {
                            cvcopy1[cidx] = vertex2;
                            cidx +=1;
                            };
                        };
                    for (int cc = 0; cc < cneigh+1; ++cc)
                        h_cv.data[cellNeighborIndex(cc,cellSet.y)] = cvcopy1[cc];
                    h_cvn.data[cellSet.y] += 1;

                    //cell k loses v1 as a neighbor
                    cneigh = h_cvn.data[cellSet.z];
                    cidx = 0;
                    for (int cc = 0; cc < cneigh-1; ++cc)
                        {
                        if(h_cv.data[cellNeighborIndex(cc,cellSet.z)] == vertex1)
                            cidx +=1;
                        h_cv.data[cellNeighborIndex(cc,cellSet.z)] = h_cv.data[cellNeighborIndex(cidx,cellSet.z)];
                        cidx +=1;
                        };
                    h_cvn.data[cellSet.z] -= 1;

                    //cell l gains v1 in between v2 and a
                    cneigh = h_cvn.data[cellSet.w];
                    vector<int> cvcopy2(cneigh+1);
                    cidx = 0;
                    for (int cc = 0; cc < cneigh; ++cc)
                        {
                        int cellIndex = h_cv.data[cellNeighborIndex(cc,cellSet.w)];
                        cvcopy2[cidx] = cellIndex;
                        cidx +=1;
                        if(cellIndex == vertex2)
                            {
                            cvcopy2[cidx] = vertex1;
                            cidx +=1;
                            };
                        };
                    for (int cc = 0; cc < cneigh+1; ++cc)
                        h_cv.data[cellNeighborIndex(cc,cellSet.w)] = cvcopy2[cc];
                    h_cvn.data[cellSet.w] = cneigh + 1;

                    };//end condition that a T1 transition should occur
                };
            };//end loop over vertex2
        };//end loop over vertices
    };

/*!
A utility function for the CPU T1 transition routine. Given two vertex indices representing an edge that will undergo
a T1 transition, return in the pass-by-reference variables a helpful representation of the cells in the T1
and the vertices to be re-wired...see the comments in "testAndPerformT1TransitionsCPU" for what that representation is
*/
void sphericalVertexModel::getCellVertexSetForT1(int vertex1, int vertex2, int4 &cellSet, int4 &vertexSet, bool &growList)
    {
    int cell1,cell2,cell3,ctest;
    int vlast, vcur, vnext, cneigh;
    ArrayHandle<int> h_cv(cellNeighbors,access_location::host, access_mode::read);
    ArrayHandle<unsigned int> h_cvn(cellNumberOfNeighbors,access_location::host,access_mode::read);
    ArrayHandle<int> h_vcn(vertexCellNeighbors,access_location::host,access_mode::read);
    ArrayHandle<unsigned int> vcnn(numberOfNeighbors,access_location::host,access_mode::read); 

    cell1 = h_vcn.data[3*vertex1];
    cell2 = h_vcn.data[3*vertex1+1];
    cell3 = h_vcn.data[3*vertex1+2];
    //cell_l doesn't contain vertex 1, so it is the cell neighbor of vertex 2 we haven't found yet
    for (int ff = 0; ff < 3; ++ff)
        {
        ctest = h_vcn.data[3*vertex2+ff];
        if(ctest != cell1 && ctest != cell2 && ctest != cell3)
            cellSet.w=ctest;
        };
    //find vertices "c" and "d"
    cneigh = h_cvn.data[cellSet.w];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cellSet.w) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cellSet.w) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cell1)];
        if(vcur == vertex2) break;
        vlast = vcur;
        vcur = vnext;
        };

    //classify cell1
    cneigh = h_cvn.data[cell1];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cell1) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cell1) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cell1)];
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
    cneigh = h_cvn.data[cell2];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cell2) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cell2) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cell2)];
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

    //classify cell3
    cneigh = h_cvn.data[cell3];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cell3) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cell3) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cell3)];
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
    cneigh = h_cvn.data[cellSet.y];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cellSet.y) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cellSet.y) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cellSet.y)];
        if(vcur == vertex1) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.x=vlast;
    vertexSet.y=vnext;
    cneigh = h_cvn.data[cellSet.w];
    vlast = h_cv.data[ cellNeighborIndex(cneigh-2,cellSet.w) ];
    vcur = h_cv.data[ cellNeighborIndex(cneigh-1,cellSet.w) ];
    for (int cn = 0; cn < cneigh; ++cn)
        {
        vnext = h_cv.data[cellNeighborIndex(cn,cellSet.w)];
        if(vcur == vertex2) break;
        vlast = vcur;
        vcur = vnext;
        };
    vertexSet.w=vlast;
    vertexSet.z=vnext;

    //Does the cell-vertex-neighbor data structure need to be bigger...for safety check all cell-vertex numbers, even if it won't be incremented?
    if(h_cvn.data[cellSet.x] == maximumVerticesPerCell || h_cvn.data[cellSet.y] == maximumVerticesPerCell || h_cvn.data[cellSet.z] == maximumVerticesPerCell || h_cvn.data[cellSet.w] == maximumVerticesPerCell)
        growList = true;
    }

void sphericalVertexModel::growCellVerticesList(int newVertexMax)
    {
    cout << "maximum number of vertices per cell grew from " <<maximumVerticesPerCell << " to " << newVertexMax << endl;
    maximumVerticesPerCell = newVertexMax+1;
    Index2D old_idx = cellNeighborIndex;
    cellNeighborIndex = Index2D(maximumVerticesPerCell,nCells);

    GPUArray<int> newCellNeighbors;
    newCellNeighbors.resize(maximumVerticesPerCell*nCells);
    {//scope for array handles
    ArrayHandle<unsigned int> h_nn(cellNumberOfNeighbors,access_location::host,access_mode::read);
    ArrayHandle<int> h_n_old(cellNeighbors,access_location::host,access_mode::read);
    ArrayHandle<int> h_n(newCellNeighbors,access_location::host,access_mode::readwrite);

    for(int cell = 0; cell < nCells; ++cell)
        {
        int neighs = h_nn.data[cell];
        for (int n = 0; n < neighs; ++n)
            {
            h_n.data[cellNeighborIndex(n,cell)] = h_n_old.data[old_idx(n,cell)];
            };
        };
    };//scope for array handles
    cellNeighbors.resize(maximumVerticesPerCell*nCells);
    cellNeighbors.swap(newCellNeighbors);
    }

void sphericalVertexModel::testAndPerformT1TransitionsGPU()
    {
    UNWRITTENCODE("t1 unfinished");
    /*
        GPUArray<int> growCellVertexListAssist;

        GPUArray<int> finishedFlippingEdges;

        //! data structure per cell for not simulataneously flipping nearby edges
        GPUArray<int> cellEdgeFlips;
        //! data structure per cell for not simulataneously flipping nearby edges
        GPUArray<int4> cellSets;
        */
    }

void sphericalVertexModel::initializeEdgeFlipLists()
    {
    UNWRITTENCODE("t1 unfinished");
    }
void sphericalVertexModel::testEdgesForT1GPU()
    {
    UNWRITTENCODE("t1 unfinished");
    }
void sphericalVertexModel::flipEdgesGPU()
    {
    UNWRITTENCODE("t1 unfinished");
    }

