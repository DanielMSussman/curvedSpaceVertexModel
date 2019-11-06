#define ENABLE_CUDA

#include "structuralFeatures.h"
/*! \file structuralFeatures.cpp */

/*!
A brute-force, O(N^2) computation of the radial distribution function for the point pattern. The
answer is stored in the GofR vector.
*/
void structuralFeatures::computeRadialDistributionFunction(vector<dVec> &points,vector<scalar2> &GofR, scalar binWidth)
    {
    int N = points.size();
    scalar L,b2,b3,b4;
    Box->getBoxDims(L,b2,b3,b4);

    //Initialize the answer vector
    int totalBins = floor(0.5*L/binWidth);
    GofR.resize(totalBins);
    for (int bb = 0; bb < totalBins; ++bb)
        GofR[bb] = make_scalar2((bb+0.5)*binWidth,0.0);

    //loop through points
    dVec dist;
    for (int ii = 0; ii < N-1; ++ii)
        {
        for (int jj = ii+1; jj < N; ++jj)
            {
            Box->minDist(points[ii],points[jj],dist);
            scalar d=norm(dist);
            int ibin = floor(d/binWidth);
            if (ibin < totalBins)
                GofR[ibin].y += 1.0;
            };
        };
    //finally, normalize the function appropriately
    for (int bb = 0; bb < totalBins; ++bb)
        {
        scalar annulusArea = PI*(((bb+1)*binWidth)*((bb+1)*binWidth)-(bb*binWidth)*(bb*binWidth));
        scalar yVal = (2.0*GofR[bb].y/N) / annulusArea;
        GofR[bb].y=yVal;
        };
    };
