#ifndef structuralFeatures_H
#define structuralFeatures_H

#include "std_include.h"
#include "functions.h"
#include "sphericalDomain.h"

/*! \file structuralFeatures.h */

//! A class that calculates various structural features of 2D point patterns
class structuralFeatures
    {
    public:
        //!The constructor takes in a defining set of boundary conditions
        structuralFeatures(BoxPtr _bx){Box = _bx;};

        //!Compute the (isotropic) radial distribution function of the point pattern
        void computeRadialDistributionFunction(vector<dVec> &points,vector<scalar2> &GofR, scalar binWidth = 0.1);

    protected:
        //!the box defining the periodic domain
        sphericalDomain Box;
    };
#endif
