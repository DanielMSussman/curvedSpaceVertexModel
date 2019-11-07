#ifndef dynamicalFeatures_H
#define dynamicalFeatures_H

#include "std_include.h"
#include "sphericalDomain.h"

/*! \file dynamicalFeatures.h */

//! A class that calculates various dynamical features for 2D systems
class dynamicalFeatures
    {
    public:
        //!The constructor takes in a defining set of boundary conditions
        dynamicalFeatures(GPUArray<dVec> &initialPos, shared_ptr<sphericalDomain> _sphere, scalar fractionAnalyzed = 1.0);

        //!Compute the mean squared displacement of the passed vector from the initial positions
        scalar computeMSD(GPUArray<dVec> &currentPos);

        //!compute the overlap function
        scalar computeOverlapFunction(GPUArray<dVec> &currentPos, scalar cutoff = 0.5);
    protected:
        //!the box defining the periodic domain
        shared_ptr<sphericalDomain> sphere;
        //!the initial positions
        vector<dVec> iPos;
        //!the number of dVecs
        int N;
    };
#endif
