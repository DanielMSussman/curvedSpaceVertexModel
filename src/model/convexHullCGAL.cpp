#include "convexHullCGAL.h"
/*! \file convexHullCGAL.cpp" */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point_3;
typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

void convexHullCGALInterface::sphericalConvexHull(dVec *points, int n, GPUArray<int> &allNeighs, GPUArray<unsigned int> &numNeighs, Index2D &nidx)
    {
    if(numNeighs.getNumElements() < n)
        numNeighs.resize(n);

    //convex hull part
    std::vector<Point_3 > p(n);
    std::map<Point_3,int> pointToIndex;
    for(int ii = 0; ii < p.size();++ii)
        {
        p[ii] = Point_3(points[ii].x[0],points[ii].x[1],points[ii].x[2]);
        pointToIndex[p[ii]] = ii;
        }
    Polyhedron_3 poly;
    CGAL::convex_hull_3(p.begin(),p.end(),poly);
    int maximumDegree = 0;

    ArrayHandle<unsigned int> nNeigh(numNeighs);
    for(Polyhedron_3::Vertex_iterator v = poly.vertices_begin(); v!= poly.vertices_end(); ++v)
        {
        int idx = pointToIndex[v->point()];
        int curDegree = v->vertex_degree();
        if(curDegree > maximumDegree)
            maximumDegree = curDegree;
        nNeigh.data[idx]= curDegree;
        }

    if(allNeighs.getNumElements() < n*maximumDegree)
        {
        allNeighs.resize(n*maximumDegree);
        nidx = Index2D(maximumDegree,n);
        }

    ArrayHandle<int> neighs(allNeighs);
    for(Polyhedron_3::Vertex_iterator v = poly.vertices_begin(); v!= poly.vertices_end(); ++v)
        {
  //      std::cout << v->point() << std::endl;
        int idx = pointToIndex[v->point()];
//        printf("\t\t (%f,%f,%f)\n",points[idx].x[0], points[idx].x[1], points[idx].x[2]);

        std::vector<int> nbs;nbs.reserve(6);
        Polyhedron_3::Halfedge_around_vertex_circulator h = v->vertex_begin();
        int ii = 0;
        do
            {
            neighs.data[nidx(ii,idx)] = pointToIndex[h->opposite()->vertex()->point()];
            ++ii;
            ++h;
            }while (h!= v->vertex_begin());
        };
    };
