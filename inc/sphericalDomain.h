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
        sphericalDomain(scalar _radius=1.0){radius = _radius; inverseRadius = 1.0/_radius;};

        HOSTDEVICE void putInBoxVirtual(dVec &p);
        HOSTDEVICE void putInBoxReal(dVec &p);
        //returns euclidean distance, not geodesic
        HOSTDEVICE void minDist(dVec &p1, dVec &p2, dVec &pans);
        HOSTDEVICE void move(dVec &p1, dVec &velocityDirection, scalar magnitude);
        HOSTDEVICE void move(dVec &p1, const dVec &velocityDirection);
        HOSTDEVICE void projectToTangentPlane(dVec &vec, const dVec &normal);
        HOSTDEVICE void projectToTangentPlaneAndNormalize(dVec &vec, const dVec &normal);

        HOSTDEVICE void changeRadius(scalar _r){radius = _r; inverseRadius = 1.0/_r;};

        HOSTDEVICE void geodesicDistance(dVec &p1, dVec &p2, scalar &dist);
        HOSTDEVICE void sphericalTriangleArea(dVec &p1, dVec &p2, dVec &p3, scalar &area);

        scalar radius=1.0;
        scalar inverseRadius = 1.0;
        dVec tangentPlaneProjection;
        dVec pt;
        dVec pt1,pt2,pt3;
        dVec disp;
    };

void sphericalDomain::geodesicDistance(dVec &p1, dVec &p2, scalar &dist)
    {
    pt1 = p1;
    pt2 = p2;
    putInBoxVirtual(pt1);
    putInBoxVirtual(pt2);
    dist = radius*acos(dot(pt1,pt2));
    }

void sphericalDomain::sphericalTriangleArea(dVec &p1, dVec &p2, dVec &p3, scalar &area)
    {
    pt1 = p1;
    pt2 = p2;
    pt3= p3;
    putInBoxVirtual(pt1);
    putInBoxVirtual(pt2);
    putInBoxVirtual(pt3);

    scalar p1Dotp2 = dot(pt1,pt2);
    scalar p1Dotp3 = dot(pt1,pt3);
    scalar p2Dotp3 = dot(pt2,pt3);
    area = -PI;
    area += acos((p2Dotp3-p1Dotp2*p1Dotp3) / (sqrt(1-p1Dotp2*p1Dotp2)*sqrt(1-p1Dotp3*p1Dotp3)) );
    area += acos((p1Dotp2-p1Dotp3*p2Dotp3) / (sqrt(1-p1Dotp3*p1Dotp3)*sqrt(1-p2Dotp3*p2Dotp3)) );
    area += acos((p1Dotp3-p1Dotp2*p2Dotp3) / (sqrt(1-p1Dotp2*p1Dotp2)*sqrt(1-p2Dotp3*p2Dotp3)) );
    area *= (radius*radius);
    }

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
