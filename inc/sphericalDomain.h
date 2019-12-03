#ifndef sphericalDomain_H
#define sphericalDomain_H

#include "std_include.h"
#include "functions.h"

#ifdef NVCC
#define HOSTDEVICE __host__ __device__ inline
#else
#define HOSTDEVICE inline __attribute__((always_inline))
#endif

class sphericalDomain
    {
    public:
        sphericalDomain(scalar _radius=1.0){radius = _radius; inverseRadius = 1.0/_radius;};

        HOSTDEVICE void changeRadius(scalar _r){radius = _r; inverseRadius = 1.0/_r;};

        HOSTDEVICE void putInBoxVirtual(dVec &p);
        HOSTDEVICE void putInBoxReal(dVec &p);
        
        //returns euclidean distance, not geodesic
        HOSTDEVICE void minDist(dVec &p1, dVec &p2, dVec &pans);
        HOSTDEVICE void move(dVec &p1, dVec &velocityDirection, scalar magnitude);
        HOSTDEVICE void move(dVec &p1, const dVec &velocityDirection);
        
        HOSTDEVICE void projectToTangentPlane(dVec &vec, const dVec &normal);
        HOSTDEVICE void projectToTangentPlaneAndNormalize(dVec &vec, const dVec &normal);

        HOSTDEVICE void getAngularCoordinates(dVec &pos, scalar &radius, scalar &theta, scalar &phi);
        HOSTDEVICE void cartesianSphericalBasisChange(scalar t, scalar p, dVec &thetaHat, dVec &phiHat); 
        HOSTDEVICE void cartesianSphericalBasisChange(dVec &cartesianPosition, dVec &thetaHat, dVec &phiHat); 

        HOSTDEVICE void geodesicDistance(dVec &p1, dVec &p2, scalar &dist);
        //!version well-conditioned for all angles
        HOSTDEVICE void geodesicDistanceTan(dVec &p1, dVec &p2, scalar &dist);
        
        //!take the gradient in spherical coordinates, even though p1 and p3 are 3-vectors, then project back. Grad is definitionally in the tangent plane
        HOSTDEVICE void gradientGeodesicDistance(dVec &p, dVec &other, dVec &derivative);
        HOSTDEVICE void gradientTriangleArea(dVec &v1, dVec &v2, dVec &v3, dVec &derivative);
        HOSTDEVICE void gradientGeodesicDistance(dVec &p, dVec &other, dVec &derivative, dVec &thetaHat, dVec &phiHat);
        HOSTDEVICE void gradientTriangleArea(dVec &v1, dVec &v2, dVec &v3, dVec &derivative, dVec &thetaHat, dVec &phiHat);
        //!take the gradient relative to v1 of changing the included angles it is part of
        HOSTDEVICE void gradientIncludedAngleSet(dVec &v1, quadAngularPosition &angleSet, dVec &derivative);

        //!given an ordered set of vertices (p1,p2,p3), what is the included angle?
        HOSTDEVICE void includedAngle(dVec &p1, dVec &p2, dVec &p3, scalar &angle);
        HOSTDEVICE void sphericalTriangleArea(dVec &p1, dVec &p2, dVec &p3, scalar &area);

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

void sphericalDomain::geodesicDistanceTan(dVec &p1, dVec &p2, scalar &dist)
    {
    //pt1 = p1;
    //pt2 = p2;
    if(p1==p2)
        {
        dist = 0.;
        return;
        }
    putInBoxVirtual(pt1);
    putInBoxVirtual(pt2);
    //acos formulation
    //dist = radius*acos(dot(pt1,pt2));
    //scalar numerator = p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
    //scalar denominator = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
    //dist = radius*acos(numerator/denominator);
    //
    //asin formulation
    dVec crossP = cross(pt1,pt2);
    //dist = radius*asin(norm(crossP));

    //atan formulation... best conditioned for all angles
    dist = radius*atan2(norm(crossP),dot(pt1,pt2));
    }

void sphericalDomain::geodesicDistance(dVec &p1, dVec &p2, scalar &dist)
    {
    //pt1 = p1;
    //pt2 = p2;
    //putInBoxVirtual(pt1);
    //putInBoxVirtual(pt2);
    //acos formulation
    //dist = radius*acos(dot(pt1,pt2));
    scalar numerator = p1[0]*p2[0]+p1[1]*p2[1]+p1[2]*p2[2];
    scalar denominator = sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])*sqrt(p2[0]*p2[0]+p2[1]*p2[1]+p2[2]*p2[2]);
    dist = radius*acos(numerator/denominator);
    //
    //asin formulation
    //dVec crossP = cross(p1,p2);
    //dist = asin(norm(crossP));

    //atan formulation... best conditioned for all angles
    //dist = atan2(norm(crossP),dot(p1,p2));
    }

void sphericalDomain::includedAngle(dVec &p1, dVec &p2, dVec &p3, scalar &angle)
    {
    pt1 = p1;
    pt2 = p2;
    pt3 = p3;

    dVec crossI = cross(pt2,pt3);
    dVec crossIm1 = cross(pt1,pt2);

    scalar determinant = pt1[0]*(pt2[1]*pt3[2]-pt2[2]*pt3[1])
                        +pt1[1]*(pt2[2]*pt3[0]-pt2[0]*pt3[2])
                        +pt1[2]*(pt2[0]*pt3[1]-pt2[1]*pt3[0]);


    angle = acos(-dot(crossI,crossIm1)/(norm(crossI)*norm(crossIm1)));
    if(determinant > 0)
        angle *= -1;
    if(angle < 0)
        angle += 2.*PI;
    };

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
    vec = vec - (dot(vec,normal)/dot(normal,normal))*normal;
    }
void sphericalDomain::projectToTangentPlaneAndNormalize(dVec &vec, const dVec &normal)
    {
    vec = vec - (dot(vec,normal)/dot(normal,normal))*normal;
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

void sphericalDomain::getAngularCoordinates(dVec &pos, scalar &radius, scalar &theta, scalar &phi)
    {
    radius = sqrt(dot(pos,pos));
    theta = acos(pos[2]/radius);
    phi = atan2(pos[1],pos[0]);
    }

void sphericalDomain::cartesianSphericalBasisChange(dVec &cartesianPosition, dVec &thetaHat, dVec &phiHat)
    {
    scalar radius = sqrt(cartesianPosition[0]*cartesianPosition[0]+cartesianPosition[1]*cartesianPosition[1]+cartesianPosition[2]*cartesianPosition[2]);
    scalar cosT = cartesianPosition[2]/radius;
    scalar sinT = sqrt(1.-cosT*cosT);
    scalar denom = sqrt(cartesianPosition[0]*cartesianPosition[0] + cartesianPosition[1]*cartesianPosition[1]);
    scalar cosP = cartesianPosition[0] / denom;
    scalar sinP = cartesianPosition[1] / denom;
    thetaHat[0] = cosT*cosP;
    thetaHat[1] = cosT*sinP;
    thetaHat[2] = -sinT;
    phiHat[0] = -sinP;
    phiHat[1] = cosP;
    phiHat[2] = 0;
    }

void sphericalDomain::cartesianSphericalBasisChange(scalar t, scalar p, dVec &thetaHat, dVec &phiHat)
    {
    thetaHat[0] = cos(t)*cos(p);
    thetaHat[1] = cos(t)*sin(p);
    thetaHat[2] = -sin(t);
    phiHat[0] = -sin(p);
    phiHat[1] = cos(p);
    phiHat[2] = 0;
    }

/*!
given a geodesic with endpoints v1 and v2, what is the gradient (spherical coords) with respect to v1?
This version assumes that the \hat{theta} and \hat{phi} directions at point v1 are already known
*/
void sphericalDomain::gradientGeodesicDistance(dVec &v1, dVec &v2, dVec &derivative, dVec &thetaHat, dVec &phiHat)
    {
    scalar r1,r2;
    r1 = sqrt(dot(v1,v1));
    r2 = sqrt(dot(v2,v2));

    scalar cosT1 = v1[2]/r1;
    scalar cosT2 = v2[2]/r2;
    scalar sinT1 = sqrt(1-cosT1*cosT1);
    scalar sinT2 = sqrt(1-cosT2*cosT2);
    scalar cosP1 = v1[0]/sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
    scalar cosP2 = v2[0]/sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
    scalar sinP1 = v1[1]/sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
    scalar sinP2 = v2[1]/sqrt(v2[0]*v2[0]+v2[1]*v2[1]);

    scalar sinP1MinusP2 = sinP1*cosP2-cosP1*sinP2;
    scalar cosP1MinusP2 = cosP1*cosP2+sinP1*sinP2;

    scalar denomPart = cosT1*cosT2 + cosP1MinusP2*sinT1*sinT2;
    scalar denomInverse = 1./sqrt(1-denomPart*denomPart);

    scalar gradTheta = 1.0*(cosT2*sinT1-cosT1*cosP1MinusP2*sinT2)*denomInverse;
    scalar gradPhi = sinT2*sinP1MinusP2 * denomInverse;

    derivative = gradTheta*thetaHat + gradPhi*phiHat;
    }

/*!
given a geodesic with endpoints v1 and v2, what is the gradient (spherical coords) with respect to v1?
This version computes the local basis at v1, then calls the other function
*/
void sphericalDomain::gradientGeodesicDistance(dVec &p, dVec &other, dVec &derivative)
    {
    scalar r1,t1,ph1;
    getAngularCoordinates(pt1 ,r1,t1,ph1);
    dVec thetaHat, phiHat;
    cartesianSphericalBasisChange(t1,ph1,thetaHat,phiHat);
    gradientGeodesicDistance(p,other,derivative,thetaHat,phiHat);
    }

void sphericalDomain::gradientIncludedAngleSet(dVec &v1, quadAngularPosition &angleSet, dVec &derivative)
    {
    pt1 = v1;
    scalar r0,t0,p0,tn1,pn1,tn2,pn2,t1,p1,t2,p2;
    scalar gradTheta, gradPhi;
    getAngularCoordinates(pt1 ,r0,t0,p0);
    tn2 = angleSet[0];
    pn2 = angleSet[1];
    tn1 = angleSet[2];
    pn1 = angleSet[3];
    t1 = angleSet[4];
    p1 = angleSet[5];
    t2 = angleSet[6];
    p2 = angleSet[7];


    scalar s01 = cos(t0)*cos(t1)+cos(p0-p1)*sin(t0)*sin(t1);
    scalar s02 = cos(t0)*cos(t2)+cos(p0-p2)*sin(t0)*sin(t2);
    scalar s12 = cos(t1)*cos(t2)+cos(p1-p2)*sin(t1)*sin(t2);
    scalar s0n1 = cos(t0)*cos(tn1)+cos(p0-pn1)*sin(t0)*sin(tn1);
    scalar s0n2 = cos(t0)*cos(tn2)+cos(p0-pn2)*sin(t0)*sin(tn2);
    scalar sn1n2 = cos(tn1)*cos(tn2)+cos(pn1-pn2)*sin(tn1)*sin(tn2);
    scalar s1n1 = cos(t1)*cos(tn1)+cos(p1-pn1)*sin(t1)*sin(tn1);
    scalar alt01 = cos(p0-p1)*cos(t0)*sin(t1) - cos(t1)*sin(t0);
    scalar alt02 = cos(p0-p2)*cos(t0)*sin(t2) - cos(t2)*sin(t0);
    scalar alt12 = cos(p1-p2)*cos(t1)*sin(t2) - cos(t2)*sin(t1);
    scalar alt0n1 = cos(p0-pn1)*cos(t0)*sin(tn1) - cos(tn1)*sin(t0);
    scalar alt0n2 = cos(p0-pn2)*cos(t0)*sin(tn2) - cos(tn2)*sin(t0);
    scalar altn1n2 = cos(pn1-pn2)*cos(tn1)*sin(tn2) - cos(tn2)*sin(tn1);
    scalar d01 = sin(p0-p1)*sin(t1);
    scalar d02 = sin(p0-p2)*sin(t2);
    scalar d0n1 = sin(p0-pn1)*sin(tn1);
    scalar d0n2 = sin(p0-pn2)*sin(tn2);

    scalar denom1 = 1.0/(pow(1-s01*s01,1.5)*sqrt(1.0-s12*s12)
        *sqrt((1-s01*s01-s02*s02+2.*s01*s02*s12 - s12*s12)/((s01*s01-1.)*(s12*s12-1.))));
    scalar denom2 = 1.0/(pow(1-s01*s01,1.5)*pow(1.0-s0n1*s0n1,1.5)
        *sqrt((1-s01*s01-s0n1*s0n1+2.*s01*s0n1*s1n1 - s1n1*s1n1)/((s01*s01-1.)*(s0n1*s0n1-1.))));
    scalar denom3 = 1.0/(pow(1-s0n1*s0n1,1.5)*sqrt(1.0-sn1n2*sn1n2)
        *sqrt((1-s0n1*s0n1-s0n2*s0n2+2.*s0n1*s0n2*sn1n2 - sn1n2*sn1n2)/((s0n1*s0n1-1.)*(sn1n2*sn1n2-1.))));


    gradTheta = denom1*(alt02*(s01*s01-1.)+alt01*(s12-s01*s02))
               +denom2*(alt01*(1.-s0n1*s0n1)*(s0n1-s01*s1n1) - alt0n1*(s01*s01-1.)*(s01-s0n1*s1n1))
               +denom3*(alt0n2*(s0n1*s0n1-1.)+alt0n1*(sn1n2-s0n1*s0n2));
    gradTheta *= radius;
    
    gradPhi = -denom1*(d02*(s01*s02-1)+d01*(s12-s01*s02))
              +denom2*(d01*(s0n1*s0n1-1.)*(s0n1-s01*s1n1)+d0n1*(s01*s01-1.)*(s01-s0n1*s1n1))
              -denom3*(d0n2*(s0n1*s0n1-1.)+d0n1*(sn1n2-s0n1*s0n2));
    gradPhi *= radius;

    dVec thetaHat, phiHat;
    cartesianSphericalBasisChange(t0,p0,thetaHat,phiHat);

    derivative = gradTheta*thetaHat + gradPhi*phiHat;

    }

void sphericalDomain::gradientTriangleArea(dVec &v1, dVec &v2, dVec &v3, dVec &derivative, dVec &thetaHat, dVec &phiHat)
    {
    scalar r1,r2,r3;
    r1 = sqrt(dot(v1,v1));
    r2 = sqrt(dot(v2,v2));
    r3 = sqrt(dot(v3,v3));

    scalar cosT1 = v1[2]/r1;
    scalar cosT2 = v2[2]/r2;
    scalar cosT3 = v3[2]/r3;
    scalar sinT1 = sqrt(1-cosT1*cosT1);
    scalar sinT2 = sqrt(1-cosT2*cosT2);
    scalar sinT3 = sqrt(1-cosT3*cosT3);
    scalar cosP1 = v1[0]/sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
    scalar cosP2 = v2[0]/sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
    scalar cosP3 = v3[0]/sqrt(v3[0]*v3[0]+v3[1]*v3[1]);
    scalar sinP1 = v1[1]/sqrt(v1[0]*v1[0]+v1[1]*v1[1]);
    scalar sinP2 = v2[1]/sqrt(v2[0]*v2[0]+v2[1]*v2[1]);
    scalar sinP3 = v3[1]/sqrt(v3[0]*v3[0]+v3[1]*v3[1]);

    scalar sinP1MinusP2 = sinP1*cosP2-cosP1*sinP2;
    scalar sinP1MinusP3 = sinP1*cosP3-cosP1*sinP3;
    scalar cosP1MinusP2 = cosP1*cosP2+sinP1*sinP2;
    scalar cosP1MinusP3 = cosP1*cosP3+sinP1*sinP3;
    scalar cosP2MinusP3 = cosP2*cosP3+sinP2*sinP3;

    double s12,s13,s23,d12,d13,denom1,denom2,denom3;
    //double tempNum;
    s12 = cosT1*cosT2+cosP1MinusP2*sinT1*sinT2;
    s13 = cosT1*cosT3+cosP1MinusP3*sinT1*sinT3;
    s23 = cosT2*cosT3+cosP2MinusP3*sinT2*sinT3;
    d12 = cosT1*cosP1MinusP2*sinT2 - cosT2*sinT1;
    d13 = cosT1*cosP1MinusP3*sinT3 - cosT3*sinT1;
    denom1 = sqrt((-1.+s12*s12)*(-1.+s12*s12)*(1-s12*s12-s13*s13-s23*s23+2.0*s12*s13*s23));
    denom2 = sqrt((-1.+s13*s13)*(-1.+s13*s13)*(1-s12*s12-s13*s13-s23*s23+2.0*s12*s13*s23));
    denom3 = sqrt((s12*s12-1.)*(s13*s13-1.)*(s12*s12-1.)*(s13*s13-1.)*(1-s12*s12-s13*s13-s23*s23+2.0*s12*s13*s23));

    scalar gradTheta = (d13*(s12*s12-1.0)+d12*(s23-s12*s13))/denom1
                      +(d12*(s13*s13-1.0)+d13*(s23-s12*s13))/denom2
                      -(d12*(s13*s13-1.)*(s13-s12*s23)+d13*(s12*s12-1.0)*(s12-s13*s23))/denom3;
    gradTheta *= radius;

    scalar gradPhi =((s12*s13 - s23)*sinP1MinusP2*sinT2 - (-1 + s12*s12)*sinP1MinusP3*sinT3)/(denom1)
                     - ((-1 + s13*s13)*sinP1MinusP2*sinT2 + (-(s12*s13) + s23)*sinP1MinusP3*sinT3)/(denom2) 
                    +((-1 + s13*s13)*(s13 - s12*s23)*sinP1MinusP2*sinT2 + (-1 + s12*s12)*(s12 - s13*s23)*
                    sinP1MinusP3*sinT3)/(denom3);
    gradPhi *=radius;

    scalar determinant = v1[0]*(v2[1]*v3[2]-v2[2]*v3[1])
                        +v1[1]*(v2[2]*v3[0]-v2[0]*v3[2])
                        +v1[2]*(v2[0]*v3[1]-v2[1]*v3[0]);

    derivative = gradTheta*thetaHat + gradPhi*phiHat;

    if(determinant < 0)
        derivative = -1.*derivative;
/*
    if(isnan(derivative[0]))
        {
        printf("%f %f %f\t %f %f %f\t %f %f %f\n",v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],v3[0],v3[1],v3[2]);
        printf("%f %f %f\n",denom1, denom2, denom3);
        printf("%.8g %.8g %.8g\n",s12,s13,s23);
        }
*/
  }

void sphericalDomain::gradientTriangleArea(dVec &v1, dVec &v2, dVec &v3, dVec &derivative)
    {
    scalar r1,t1,p1;
    getAngularCoordinates(v1 ,r1,t1,p1);
    dVec thetaHat, phiHat;
    cartesianSphericalBasisChange(t1,p1,thetaHat,phiHat);
    gradientTriangleArea(v1,v2,v3,derivative,thetaHat,phiHat);
    }
#undef HOSTDEVICE
#endif
