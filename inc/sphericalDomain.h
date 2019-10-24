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
        HOSTDEVICE void dGeodesicDistanceDVertex(dVec &p, dVec &other, dVec &derivative);
        HOSTDEVICE void sphericalTriangleArea(dVec &p1, dVec &p2, dVec &p3, scalar &area);
        HOSTDEVICE void dSphericalTriangleAreaDVertex(dVec &p1, dVec &p2, dVec &p3, dVec &derivative);

        HOSTDEVICE void dGeodesicDistanceDVertex2(dVec &p, dVec &other, dVec &derivative);
        HOSTDEVICE void dSphericalTriangleAreaDVertex2(dVec &p1, dVec &p2, dVec &p3, dVec &derivative);
        HOSTDEVICE scalar normCross(dVec &p1, dVec &p2);
        scalar radius=1.0;
        scalar inverseRadius = 1.0;
        dVec tangentPlaneProjection;
        dVec pt;
        dVec pt1,pt2,pt3;
        dVec disp;
    };

scalar sphericalDomain::normCross(dVec &p1, dVec &p2)
    {
    scalar term1 = (p1[1]*p2[0]-p1[0]*p2[1]);
    scalar term2 = (p1[2]*p2[0]-p1[0]*p2[2]);
    scalar term3 = (p1[2]*p2[1]-p1[1]*p2[2]);

    return sqrt(term1*term1+term2*term2+term3*term3);
    }

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
    pt3 = p3;
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
/*
    scalar a,b,c,alpha,beta,gamma;
    //a=asin(normCross(pt2,pt3));
    //b=asin(normCross(pt3,pt1));
    //c=asin(normCross(pt1,pt2));
    a=atan(normCross(pt2,pt3)/dot(pt2,pt3));
    b=atan(normCross(pt3,pt1)/dot(pt3,pt1));
    c=atan(normCross(pt1,pt2)/dot(pt1,pt2));
    alpha = acos((cos(a)-cos(b)*cos(c))/(sin(b)*sin(c)));
    beta = acos((cos(b)-cos(a)*cos(c))/(sin(a)*sin(c)));
    gamma = acos((cos(c)-cos(a)*cos(b))/(sin(a)*sin(b)));
    area = radius*radius*(alpha+beta+gamma-PI);

*/
    }

//!Change vec so that only the component in the tangent plane to normal remains
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

//!Project velocit to tangent plane and displace
void sphericalDomain::move(dVec &p1, const dVec &velocityDirection)
    {
    disp = velocityDirection;
    projectToTangentPlane(disp,p1);
    p1 = p1+disp;
    putInBoxReal(p1);
    }

void sphericalDomain::dSphericalTriangleAreaDVertex(dVec &p1, dVec &p2, dVec &p3, dVec &derivative)
    {
    pt1 = p1;
    pt2 = p2;
    pt3 = p3;
    putInBoxVirtual(pt1);
    putInBoxVirtual(pt2);
    putInBoxVirtual(pt3);
    scalar p1Dotp2 = dot(pt1,pt2);
    scalar p1Dotp3 = dot(pt1,pt3);
    scalar p2Dotp3 = dot(pt2,pt3);

    derivative[0] = -(((pt3[0]*(p1Dotp3)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(sqrt(1 - pow(p1Dotp2,2))*pow(1 - pow(p1Dotp3,2),1.5)) + (-((p1Dotp2)*pt3[0]) - pt2[0]*(p1Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p1Dotp3,2))) + (pt2[0]*(p1Dotp2)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p1Dotp3,2))))/sqrt(1 - pow(p2Dotp3 - (p1Dotp2)*(p1Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p1Dotp3,2))))) - ((pt3[0] - pt2[0]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt2[0]*(p1Dotp2)*(p1Dotp3 - (p1Dotp2)*(p2Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp3 - (p1Dotp2)*(p2Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p2Dotp3,2)))) - ((pt2[0] - pt3[0]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp3,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt3[0]*(p1Dotp3)*(p1Dotp2 - (p1Dotp3)*(p2Dotp3)))/(pow(1 - pow(p1Dotp3,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp2 - (p1Dotp3)*(p2Dotp3),2)/((1 - pow(p1Dotp3,2))*(1 - pow(p2Dotp3,2))));
    derivative[1]=-(((pt3[1]*(p1Dotp3)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(sqrt(1 - pow(p1Dotp2,2))*pow(1 - pow(p1Dotp3,2),1.5)) + (-((p1Dotp2)*pt3[1]) - pt2[1]*(p1Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p1Dotp3,2))) + (pt2[1]*(p1Dotp2)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p1Dotp3,2))))/sqrt(1 - pow(p2Dotp3 - (p1Dotp2)*(p1Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p1Dotp3,2))))) - ((pt3[1] - pt2[1]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt2[1]*(p1Dotp2)*(p1Dotp3 - (p1Dotp2)*(p2Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp3 - (p1Dotp2)*(p2Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p2Dotp3,2)))) - ((pt2[1] - pt3[1]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp3,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt3[1]*(p1Dotp3)*(p1Dotp2 - (p1Dotp3)*(p2Dotp3)))/(pow(1 - pow(p1Dotp3,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp2 - (p1Dotp3)*(p2Dotp3),2)/((1 - pow(p1Dotp3,2))*(1 - pow(p2Dotp3,2))));
    derivative[2]=-(((pt3[2]*(p1Dotp3)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(sqrt(1 - pow(p1Dotp2,2))*pow(1 - pow(p1Dotp3,2),1.5)) + (-((p1Dotp2)*pt3[2]) - pt2[2]*(p1Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p1Dotp3,2))) + (pt2[2]*(p1Dotp2)*(p2Dotp3 - (p1Dotp2)*(p1Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p1Dotp3,2))))/sqrt(1 - pow(p2Dotp3 - (p1Dotp2)*(p1Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p1Dotp3,2))))) - ((pt3[2] - pt2[2]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp2,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt2[2]*(p1Dotp2)*(p1Dotp3 - (p1Dotp2)*(p2Dotp3)))/(pow(1 - pow(p1Dotp2,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp3 - (p1Dotp2)*(p2Dotp3),2)/((1 - pow(p1Dotp2,2))*(1 - pow(p2Dotp3,2)))) - ((pt2[2] - pt3[2]*(p2Dotp3))/(sqrt(1 - pow(p1Dotp3,2))*sqrt(1 - pow(p2Dotp3,2))) + (pt3[2]*(p1Dotp3)*(p1Dotp2 - (p1Dotp3)*(p2Dotp3)))/(pow(1 - pow(p1Dotp3,2),1.5)*sqrt(1 - pow(p2Dotp3,2))))/sqrt(1 - pow(p1Dotp2 - (p1Dotp3)*(p2Dotp3),2)/((1 - pow(p1Dotp3,2))*(1 - pow(p2Dotp3,2))));

    derivative = radius*radius*derivative;
    }

void sphericalDomain::dGeodesicDistanceDVertex(dVec &p, dVec &other, dVec &derivative)
    {
    pt1 = p;
    pt2 = other;
    putInBoxVirtual(pt1);
    putInBoxVirtual(pt2);
    scalar denomInverse = -1.*radius/sqrt(1-dot(pt1,pt2));
    derivative = pt2*denomInverse;
    }

void sphericalDomain::dSphericalTriangleAreaDVertex2(dVec &p1, dVec &p2, dVec &p3, dVec &derivative)
    {
    pt1 = p1;
    pt2 = p2;
    pt3 = p3;
    derivative.x = (2*(pt2[1] - pt3[1])*(-(pt2[1]*pt3[0]) + pt1[1]*(-pt2[0] + pt3[0]) + pt1[0]*(pt2[1] - pt3[1]) + pt2[0]*pt3[1]) + 2*(pt2[2] - pt3[2])*(-(pt2[2]*pt3[0]) + pt1[2]*(-pt2[0] + pt3[0]) + pt1[0]*(pt2[2] - pt3[2]) + pt2[0]*pt3[2]))/(4.*sqrt(pow(pt1[1]*(pt2[0] - pt3[0]) + pt2[1]*pt3[0] - pt2[0]*pt3[1] + pt1[0]*(-pt2[1] + pt3[1]),2) + pow(pt1[2]*(pt2[0] - pt3[0]) + pt2[2]*pt3[0] - pt2[0]*pt3[2] + pt1[0]*(-pt2[2] + pt3[2]),2) + pow(pt1[2]*(pt2[1] - pt3[1]) + pt2[2]*pt3[1] - pt2[1]*pt3[2] + pt1[1]*(-pt2[2] + pt3[2]),2)));
    derivative.y = (2*(pt2[0] - pt3[0])*(pt1[1]*(pt2[0] - pt3[0]) + pt2[1]*pt3[0] - pt2[0]*pt3[1] + pt1[0]*(-pt2[1] + pt3[1])) + 2*(pt2[2] - pt3[2])*(-(pt2[2]*pt3[1]) + pt1[2]*(-pt2[1] + pt3[1]) + pt1[1]*(pt2[2] - pt3[2]) + pt2[1]*pt3[2]))/(4.*sqrt(pow(pt1[1]*(pt2[0] - pt3[0]) + pt2[1]*pt3[0] - pt2[0]*pt3[1] + pt1[0]*(-pt2[1] + pt3[1]),2) + pow(pt1[2]*(pt2[0] - pt3[0]) + pt2[2]*pt3[0] - pt2[0]*pt3[2] + pt1[0]*(-pt2[2] + pt3[2]),2) + pow(pt1[2]*(pt2[1] - pt3[1]) + pt2[2]*pt3[1] - pt2[1]*pt3[2] + pt1[1]*(-pt2[2] + pt3[2]),2)));
    derivative.z = (2*(pt2[0] - pt3[0])*(pt1[2]*(pt2[0] - pt3[0]) + pt2[2]*pt3[0] - pt2[0]*pt3[2] + pt1[0]*(-pt2[2] + pt3[2])) + 2*(pt2[1] - pt3[1])*(pt1[2]*(pt2[1] - pt3[1]) + pt2[2]*pt3[1] - pt2[1]*pt3[2] + pt1[1]*(-pt2[2] + pt3[2])))/(4.*sqrt(pow(pt1[1]*(pt2[0] - pt3[0]) + pt2[1]*pt3[0] - pt2[0]*pt3[1] + pt1[0]*(-pt2[1] + pt3[1]),2) + pow(pt1[2]*(pt2[0] - pt3[0]) + pt2[2]*pt3[0] - pt2[0]*pt3[2] + pt1[0]*(-pt2[2] + pt3[2]),2) + pow(pt1[2]*(pt2[1] - pt3[1]) + pt2[2]*pt3[1] - pt2[1]*pt3[2] + pt1[1]*(-pt2[2] + pt3[2]),2)));
    }

void sphericalDomain::dGeodesicDistanceDVertex2(dVec &p, dVec &other, dVec &derivative)
    {
    pt1 = p;
    pt2 = other;
    scalar denomInverse = norm(pt1-pt2);
    derivative=denomInverse*(pt1-pt2);
    }
#undef HOSTDEVICE
#endif
