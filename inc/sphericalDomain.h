#ifndef sphericalDomain_H
#define sphericalDomain_H

#include "std_include.h"

#ifdef NVCC
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#endif

class sphericalDomain
    {
    public:
        sphericalDomain(scalar _radius=1.0){radius = _radius;};

        HOSTDEVICE void putInBoxReal(dVec &p);
        //returns euclidean distance, not geodesic
        HOSTDEVICE void minDist(dVec &p1, dVec &p2, dVec &pans);
        HOSTDEVICE void move(dVec &p1, dVec &velocityDirection, scalar magnitude);
        HOSTDEVICE void projectToTangentPlane(dVec &vec, dVec &normal);

        scalar radius=1.0;
        dVec tangentPlaneProjection;
        dVec pt;
    };

void sphericalDomain::projectToTangentPlane(dVec &vec, dVec &normal)
    {
    vec = vec - dot(vec,normal)*normal;
    }

void sphericalDomain::putInBoxReal(dVec &p)
    {
    pt = p*(1.0/norm(p));
    p = radius*pt;
    };

void sphericalDomain::minDist(dVec &p1, dVec &p2, dVec &pans)
    {
    pans = p1-p2;
    };

void sphericalDomain::move(dVec &p1, dVec &velocityDirection, scalar magnitude)
    {
    projectToTangentPlane(velocityDirection,p1);
    velocityDirection = velocityDirection*(1.0/ norm(velocityDirection));
    p1 = p1+magnitude*velocityDirection;
    putInBoxReal(p1);
    projectToTangentPlane(velocityDirection,p1);
    velocityDirection = velocityDirection*(1.0/ norm(velocityDirection));
    }

#undef HOSTDEVICE
#endif
