#include "sphericalVertexModel.h"

sphericalVertexModel::sphericalVertexModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
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
    int cnnArraySize = cellNeighbors.getNumElements();
    currentVertexAroundCell.resize(cnnArraySize);
    lastVertexAroundCell.resize(cnnArraySize);
    nextVertexAroundCell.resize(cnnArraySize);
    int nVertices = positions.getNumElements();

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
    computeGeometry();
    };

void sphericalVertexModel::computeGeometryCPU()
    {
    ArrayHandle<dVec> p(positions);
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<int> cvn(cellNeighbors);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<scalar2> ap(areaPerimeter);
    scalar totalArea = 0;
    for (int cc = 0; cc < nCells; ++cc)
        {
        int neighs = cnn.data[cc];
        dVec cellPos(0.0);
        for (int nn = 0; nn < neighs; ++nn)
            {
            cellPos = cellPos + p.data[cvn.data[cellNeighborIndex(nn,cc)]];
            /// move curLastNext vertex stuff here...
            }
        sphere.putInBoxReal(cellPos);
        cp.data[cc]=cellPos;

        int lastVertexIdx = cvn.data[cellNeighborIndex(neighs-1,cc)];
        dVec lastVertexPos = p.data[lastVertexIdx];
        dVec curVertexPos;
        int curVertexIdx;
        scalar perimeter = 0;
        scalar area = 0;
        scalar tempVal;

        for (int nn = 0; nn < neighs; ++nn)
            {
            int cni = cellNeighborIndex(nn,cc);
            curVertexIdx = cvn.data[cellNeighborIndex(nn,cc)];
            curVertexPos = p.data[curVertexIdx];
            sphere.geodesicDistance(lastVertexPos,curVertexPos,tempVal);
            perimeter += tempVal;
            sphere.sphericalTriangleArea(cellPos,lastVertexPos,curVertexPos,tempVal);
            area +=tempVal;

            curVert.data[cni] = curVertexPos;
            lastVert.data[cni] = lastVertexPos;

            lastVertexIdx = curVertexIdx;
            lastVertexPos = curVertexPos;
            }
        for (int nn = 0; nn < neighs - 1; ++nn)
            nextVert.data[cellNeighborIndex(nn,cc)] = curVert.data[cellNeighborIndex(nn+1,cc)];
        nextVert.data[cellNeighborIndex(neighs-1,cc)] = curVert.data[cellNeighborIndex(0,cc)];

        ap.data[cc].x = area;
        ap.data[cc].y = perimeter;
        totalArea += area;
        }
        printf("total area = %f\n",totalArea);
    }

void sphericalVertexModel::computeGeometryGPU()
    {
    }
void sphericalVertexModel::computeForceGPU()
    {
    }
void sphericalVertexModel::computeForceCPU()
    {
    ArrayHandle<dVec> cp(cellPositions);
    ArrayHandle<int> cvn(cellNeighbors);
    ArrayHandle<dVec> curVert(currentVertexAroundCell);
    ArrayHandle<dVec> lastVert(lastVertexAroundCell);
    ArrayHandle<dVec> nextVert(nextVertexAroundCell);
    ArrayHandle<unsigned int> cnn(cellNumberOfNeighbors);
    ArrayHandle<scalar2> ap(areaPerimeter);
    }
