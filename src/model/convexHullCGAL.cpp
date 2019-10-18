#include "convexHullCGAL.h"
/*! \file convexHullCGAL.cpp" */

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point_3;
typedef Polyhedron_3::Facet_handle Facet;
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

void getFacetCenter(Polyhedron_3 &P, Polyhedron_3::Facet_iterator f, dVec &center, int &neighs)
    {
    neighs = 0;
    center[0]=0.0;center[1]=0.0;center[2]=0.0;
    Polyhedron_3::Halfedge_around_facet_circulator h = f->facet_begin();
    do {
        Point_3 pp = h->vertex()->point();
        center[0] += pp.x();
        center[1] += pp.y();
        center[2] += pp.z();
        neighs += 1;
        } while(++h != f->facet_begin());
    center = (1.0/neighs)*center;
    }

void convexHullCGALInterface::sphericalConvexHullForVertexModel(dVec *points, int n, GPUArray<int> &allNeighs, GPUArray<unsigned int> &numNeighs, Index2D &nidx, GPUArray<dVec> &vertexPositions, GPUArray<int> &vertexNeighs, GPUArray<int> &vertexCellNeighs, GPUArray<unsigned int> &numVertexNeighs, Index2D &vnidx)
    {
    //first, sadly copyPaste the above routine to get the convex hull and cellNeighbor positions
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

    std::map<Facet,int> facetToIndex;
    int faceIdx= 0;
    for(Polyhedron_3::Facet_iterator f = poly.facets_begin(); f != poly.facets_end();++f)
        {
        facetToIndex[f] = faceIdx;
        faceIdx+=1;
        }
    int nFaces = faceIdx + 1;

    vertexPositions.resize(nFaces);
    numVertexNeighs.resize(nFaces);
    int maxVns = 0;
    {//get vertex positions as the facet centers
    int vn;
    ArrayHandle<dVec> vp(vertexPositions);
    ArrayHandle<unsigned int> vns(numVertexNeighs);
    int idx = 0;
    for(Polyhedron_3::Facet_iterator f = poly.facets_begin(); f != poly.facets_end();++f)
        {
        getFacetCenter(poly, f, vp.data[idx], vn);
        vns.data[idx] = vn;
        if(vn > maxVns)
            maxVns = vn;
        idx +=1;
        }
    }

    printf("max vns = %i\n",maxVns);
    vertexNeighs.resize(nFaces*maxVns);
    vertexCellNeighs.resize(nFaces*maxVns);
    vnidx = Index2D(maxVns,nFaces);
    {//get vertex neighbors
    ArrayHandle<int> vnset(vertexNeighs);
    ArrayHandle<int> vcset(vertexCellNeighs);
    for(Polyhedron_3::Facet_iterator f = poly.facets_begin(); f != poly.facets_end();++f)
        {
        int vIdx = facetToIndex[f];
        Polyhedron_3::Halfedge_around_facet_circulator h = f->facet_begin();
        int neighNum = 0;
        do {
            int otherVidx = facetToIndex[h->opposite()->facet()];
            int cellIndex = pointToIndex[h->vertex()->point()];
            vnset.data[vnidx(neighNum,vIdx)] = otherVidx;
            vcset.data[vnidx(neighNum,vIdx)] = cellIndex;
            neighNum +=1;
            } while(++h != f->facet_begin());
        }
    }

    //determine the vertices that compose each cell...first the cell sizes
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

    //...and next the actual neighbors
    ArrayHandle<int> neighs(allNeighs);
    for(Polyhedron_3::Vertex_iterator v = poly.vertices_begin(); v!= poly.vertices_end(); ++v)
        {
        int idx = pointToIndex[v->point()];

        std::vector<int> nbs;nbs.reserve(6);
        Polyhedron_3::Halfedge_around_vertex_circulator h = v->vertex_begin();
        int ii = 0;
        do
            {
            int vertexIdx = facetToIndex[h->facet()];
            neighs.data[nidx(ii,idx)] = vertexIdx;
            ++ii;
            ++h;
            }while (h!= v->vertex_begin());
        };
    };
