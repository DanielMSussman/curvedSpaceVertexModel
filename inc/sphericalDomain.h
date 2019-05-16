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

        HOSTDEVICE void putInBoxVirtual(dVec &p);
        HOSTDEVICE void putInBoxReal(dVec &p);
        //returns euclidean distance, not geodesic
        HOSTDEVICE void minDist(dVec &p1, dVec &p2, dVec &pans);
        HOSTDEVICE void move(dVec &p1, dVec &velocityDirection, scalar magnitude);
        HOSTDEVICE void move(dVec &p1, const dVec &velocityDirection);
        HOSTDEVICE void projectToTangentPlane(dVec &vec, const dVec &normal);
        HOSTDEVICE void projectToTangentPlaneAndNormalize(dVec &vec, const dVec &normal);

        scalar radius=1.0;
        dVec tangentPlaneProjection;
        dVec pt;
        dVec disp;
    };

void sphericalDomain::projectToTangentPlane(dVec &vec, const dVec &normal)
    {
    pt = normal*(1.0/norm(normal));
    vec = vec - dot(vec,pt)*pt;
    }
void sphericalDomain::projectToTangentPlaneAndNormalize(dVec &vec, const dVec &normal)
    {
    pt = normal*(1.0/norm(normal));
    vec = vec - dot(vec,pt)*pt;
    vec = vec*(1.0/norm(vec));
    }

void sphericalDomain::putInBoxVirtual(dVec &p)
    {
    p = p*(1.0/norm(p));
    };
void sphericalDomain::putInBoxReal(dVec &p)
    {
    pt = p*(1.0/norm(p))*radius;
    p = pt;
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
    }

void sphericalDomain::move(dVec &p1, const dVec &velocityDirection)
    {
    disp = velocityDirection;
    projectToTangentPlane(disp,p1);
    p1 = p1+disp;
    putInBoxReal(p1);
    }

#undef HOSTDEVICE
#endif
