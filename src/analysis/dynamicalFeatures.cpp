#define ENABLE_CUDA

#include "dynamicalFeatures.h"
#include "functions.h"
/*! \file dynamicalFeatures.cpp */

dynamicalFeatures::dynamicalFeatures(GPUArray<dVec> &initialPos, shared_ptr<sphericalDomain> _sphere, scalar fractionAnalyzed)
    {
    sphere = _sphere;
    copyGPUArrayData(initialPos,iPos);
    N = iPos.size();
    if(fractionAnalyzed < 1)
        N = floor(N*fractionAnalyzed);
    cout << "dynamical analysis package pointing at sphere of radius " << sphere->radius << " and with " << N << "points to analyze" << endl; cout.flush();
    };

scalar dynamicalFeatures::computeMSD(GPUArray<dVec> &currentPos)
    {
    scalar msd = 0.0;
    ArrayHandle<dVec> fPos(currentPos,access_location::host,access_mode::read);
    dVec cur,init;
    scalar disp;
    for (int ii = 0; ii < N; ++ii)
        {
        cur = fPos.data[ii];
        init = iPos[ii];
        sphere->geodesicDistance(init,cur,disp);
        if(!isnan(disp))
            msd += disp*disp;
        else
            printf("%i %g %g %g\n",ii ,init[0]-cur[0],init[1]-cur[1],init[2]-cur[2]);
        };
    msd = msd / N;
    return msd;
    };

scalar dynamicalFeatures::computeOverlapFunction(GPUArray<dVec> &currentPos, scalar cutoff)
    {
    scalar overlap = 0.0;
    ArrayHandle<dVec> fPos(currentPos,access_location::host,access_mode::read);
    dVec cur,init;
    scalar disp;
    for (int ii = 0; ii < N; ++ii)
        {
        cur = fPos.data[ii];
        init = iPos[ii];
        sphere->geodesicDistance(init,cur,disp);
        if(disp < cutoff)
            overlap += 1;
        };
    overlap = overlap / N;
    return overlap;
    }
