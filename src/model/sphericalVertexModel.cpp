#include "sphericalVertexModel.h"

sphericalVertexModel::sphericalVertexModel(int n, noiseSource &_noise, bool _useGPU, bool _neverGPU) : sphericalModel(n,_noise,_useGPU,_neverGPU)
    {
    //the initialization sets n cell positions randomly
    //use the convex huller to determine an initial delaunay triangulation on the sphere
    cout << "extracting vertex positions from convex huller" << endl;
        {
        GPUArray<dVec> cellPositions(positions);

        ArrayHandle<dVec> cellPos(cellPositions);
        convexHuller.sphericalConvexHullForVertexModel(cellPos.data,N,cellNeighbors,cellNumberOfNeighbors,cellNeighborIndex,positions,neighbors,vertexCellNeighbors,numberOfNeighbors,neighborIndex);
        };
    int nVertices = positions.getNumElements();
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
    };
