#include "sphericalVertexModel.h"

sphericalVertexModel::sphericalVertexModel(int n, noiseSource &_noise, scalar _area, scalar _perimeter, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
    {
    selfForceCompute =true;
    //the initialization sets n cell positions randomly
    //use the convex huller to determine an initial delaunay triangulation on the sphere
    cout << "extracting vertex positions from convex huller" << endl;
    nCells = n;
    setRadius(sqrt(nCells/(4.0*PI)));
    cellPositions = positions;
        {
        ArrayHandle<dVec> cellPos(cellPositions);
        convexHuller.sphericalConvexHullForVertexModel(cellPos.data,N,cellNeighbors,cellNumberOfNeighbors,cellNeighborIndex,positions,neighbors,vertexCellNeighbors,numberOfNeighbors,neighborIndex);
        };
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

    computeGeometry();
    };

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
    enforceTopology();
    };



void sphericalVertexModel::computeGeometryCPU()
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
            sphere.sphericalTriangleArea(cellPos,lastVertexPos,curVertexPos,tempVal);
            area +=tempVal;

            curVert.data[forceSetIdx] = curVertexPos;
            lastVert.data[forceSetIdx] = lastVertexPos;
            nextVert.data[forceSetIdx] = nextVertexPos;

            lastVertexPos = curVertexPos;
            curVertexIdx = nextVertexIdx;
            curVertexPos = nextVertexPos;
            }

        ap.data[cc].x = area;
        ap.data[cc].y = perimeter;
        totalArea += area;
        totalPerimeter += perimeter;
//        printf("%i, n=%i: %f\t%f\n",cc,neighs,area,perimeter);
        }
        scalar excessArea = totalArea - 4.0*PI*sphere.radius*sphere.radius;
        printf("excess area = %g\t total peri = %f \n",excessArea,totalPerimeter);
    }

void sphericalVertexModel::computeGeometryGPU()
    {
    }
void sphericalVertexModel::computeForceGPU()
    {
    }
void sphericalVertexModel::computeForceCPU()
    {
    printf("computing forces\n");
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
            sphere.dGeodesicDistanceDVertex(vCur,vLast,tempVar);
            f -= 2.0*perimeterDifference*tempVar;
            sphere.dGeodesicDistanceDVertex(vCur,vNext,tempVar);
            f -= 2.0*perimeterDifference*tempVar;

            sphere.dSphericalTriangleAreaDVertex(vCur,vLast,vNext,tempVar);
//            printf("vertex %i, cell %i, area force (%f,%f,%f)\n",vertexIndex, cc,tempVar[0],tempVar[1],tempVar[2]);

            f -= 2.0*areaDifference*tempVar;
            };
        //only allow forces in the tangent plane?
        sphere.projectToTangentPlane(f,vCur);
        force.data[vertexIndex] = f;
        forceNorm += dot(f,f);
//        printf("vertex %i, force (%f,%f,%f)\n",vertexIndex, f[0],f[1],f[2]);
        };

    printf("total force norm =  %g\n",forceNorm/N);
    }
