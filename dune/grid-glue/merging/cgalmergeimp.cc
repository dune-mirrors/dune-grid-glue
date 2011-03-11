// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CGAL_EXTERN
#include "config.h"
#endif

#include <dune/common/bitsetvector.hh>

#include <dune/grid-glue/merging/cgalmergeimp.hh>



#ifdef HAVE_CGAL  // without CGAL we can still handle 1d problems
// 2d
//#include "bso_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>


#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#endif

#ifdef HAVE_CGAL

#ifdef CGAL_USE_GMP

// GMP is installed. Use the GMP rational number-type.
  #include <CGAL/Gmpq.h>

typedef CGAL::Gmpq Number_type;

#else

// GMP is not installed. Use CGAL's exact rational number-type.
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif

// for debugging: print a CGAL polygon to the screen
template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator vit;

  std::cout << "[ " << P.size() << " vertices:";
  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit)
    std::cout << " (" << *vit << ')';
  std::cout << " ]" << std::endl;
}

// for debugging: print a CGAL polygon with holes to the screen
template<class Kernel, class Container>
void print_polygon_with_holes(const CGAL::Polygon_with_holes_2<Kernel, Container> & pwh)
{
  if (! pwh.is_unbounded()) {
    std::cout << "{ Outer boundary = ";
    print_polygon (pwh.outer_boundary());
  }
  else
    std::cout << "{ Unbounded polygon." << std::endl;

  typename CGAL::Polygon_with_holes_2<Kernel,Container>::Hole_const_iterator hit;
  unsigned int k = 1;

  std::cout << "  " << pwh.number_of_holes() << " holes:" << std::endl;
  for (hit = pwh.holes_begin(); hit != pwh.holes_end(); ++hit, ++k) {
    std::cout << "    Hole #" << k << " = ";
    print_polygon (*hit);
  }
  std::cout << " }" << std::endl;
}

/** See CGAL manual Section 25.3.7 to see how this works
 * \tparam T The type used by Dune for coordinates
 * \tparam dim The element dimension.  It must be 3.
 *         It is a template parameter only to make the code compile
 */
template <int dim, class T>
void CGALMergeImp<dim,T>::makeHexahedron(Polyhedron_3& P,
                                         const std::vector<Dune::FieldVector<T,dim> >& c)
{
  // appends a cube of size [0,1]ˆ3 to the polyhedron P.
  CGAL_precondition( P.is_valid());
  typedef typename Polyhedron_3::Point_3 Point;
  typedef typename Polyhedron_3::Halfedge_handle Halfedge_handle;
  Halfedge_handle h = P.make_tetrahedron( Point( c[1][0], c[1][1], c[1][2] ),       // vertex 1
                                          Point( c[4][0], c[4][1], c[4][2] ),       // vertex 4
                                          Point( c[0][0], c[0][1], c[0][2] ),       // vertex 0
                                          Point( c[2][0], c[2][1], c[2][2] ));      // vertex 2
  Halfedge_handle g = h->next()->opposite()->next();

  P.split_edge( h->next());
  P.split_edge( g->next());
  P.split_edge( g);

  h->next()->vertex()->point()     = Point( c[5][0], c[5][1], c[5][2] );            // vertex 5
  g->next()->vertex()->point()     = Point( c[6][0], c[6][1], c[6][2] );            // vertex 6
  g->opposite()->vertex()->point() = Point( c[3][0], c[3][1], c[3][2] );            // vertex 3

  Halfedge_handle f = P.split_facet( g->next(),
                                     g->next()->next()->next());

  Halfedge_handle e = P.split_edge( f);
  e->vertex()->point() = Point( c[7][0], c[7][1], c[7][2] );                        // vertex 7

  P.split_facet( e, f->next()->next());

  CGAL_postcondition( P.is_valid());
}
#endif

//////////////////////////////////////////////////////////////////////
//  Compute for each face of each element whether it will also intersect
//////////////////////////////////////////////////////////////////////

/** \todo This method should really be a private class method, but then we'd get in trouble
 * because we don't want cgal stuff in cgalmerge.hh at all (because of compile time).
 * I guess the proper solution is a separate implementation class.
 */
template<int dim, typename T>
void CGALMergeImp<dim,T>::computeNeighborIntersections(const Dune::GeometryType& elementType,
                                                       const std::vector<Dune::FieldVector<T,dim> >& elementCorners,
                                                       const Nef_Polyhedron_3& NQ,
                                                       const Nef_Polyhedron_3& intersection,
                                                       std::bitset<(1<<dim)>& neighborIntersects
                                                       )
{
  const Dune::GenericReferenceElement<T,dim>& refElement = Dune::GenericReferenceElements<T,dim>::general(elementType);

  // for each vertex: is it also a vertex of the intersection?
  Dune::BitSetVector<1> isVertexInIntersection(elementCorners.size(),false);
  Dune::BitSetVector<1> isContainedInOtherElement(elementCorners.size(),false);

  for (int i=0; i!=elementCorners.size(); ++i) {

    Point_3 v(elementCorners[i][0], elementCorners[i][1], elementCorners[i][2]);

    typename Nef_Polyhedron_3::Volume_const_handle volumeHandle;
    typename Nef_Polyhedron_3::Halffacet_const_handle facetHandle;
    typename Nef_Polyhedron_3::Halfedge_const_handle edgeHandle;
    typename Nef_Polyhedron_3::Vertex_const_handle vertexHandle;
    isVertexInIntersection[i] = assign(vertexHandle, intersection.locate(v));

    if (isVertexInIntersection[i][0]) {
      // if the vertex is a vertex of the intersection it must necessarily
      // be contained in the other element
      isContainedInOtherElement[i] = true;
    } else {
      typename Nef_Polyhedron_3::Object_handle q_handle = NQ.locate(v);

      /*            std::cout << "volume: " << (assign(volumeHandle, q_handle) && volumeHandle->mark()) << std::endl;
                  std::cout << "facet: " << assign(facetHandle, q_handle) << std::endl;
                  std::cout << "edge: " << assign(edgeHandle, q_handle) << std::endl;
                  std::cout << "vertex: " << assign(vertexHandle, q_handle) << std::endl;*/

      isContainedInOtherElement[i] = (assign(volumeHandle, q_handle) && volumeHandle->mark())
                                     or assign(facetHandle, q_handle)
                                     or assign(edgeHandle, q_handle)
                                     or assign(vertexHandle, q_handle);
    }
  }

  /*    std::cout << "isVertexInIntersection1: " << isVertexInIntersection << std::endl;
      std::cout << "isContainedInOtherElement1: " << isContainedInOtherElement << std::endl;*/

  //
  for (size_t i=0; i<refElement.size(1); i++) {

    // how many vertices of this face are a vertex of the intersection?
    int count = 0;
    for (size_t j=0; j<refElement.size(i,1,dim); j++)
      count += isVertexInIntersection[refElement.subEntity(i,1,j,dim)][0];

    if (count != 0 and count != refElement.size(i,1,dim)) {
      neighborIntersects[i] = true;
    } else {

      int countInOther = 0;
      for (size_t j=0; j<refElement.size(i,1,dim); j++)
        countInOther += isContainedInOtherElement[refElement.subEntity(i,1,j,dim)][0];

      neighborIntersects[i] = (countInOther == refElement.size(i,1,dim));

    }

  }

  //std::cout << "neighborIntersects: " << neighborIntersects << std::endl;

}


template<int dim, typename T>
void CGALMergeImp<dim,T>::compute1dIntersection(const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<T,dim,dim> >& grid1Geometry,
                                                const std::vector<Dune::FieldVector<T,dim> >& grid1ElementCorners,
                                                unsigned int grid1Index,
                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<T,dim,dim> >& grid2Geometry,
                                                const std::vector<Dune::FieldVector<T,dim> >& grid2ElementCorners,
                                                unsigned int grid2Index,
                                                std::vector<RemoteSimplicialIntersection<T,dim,dim,dim> >& intersections
                                                )
{
  // Check consistent orientation
  // \todo Reverse the orientation if this check fails
  assert(grid1ElementCorners[0][0] <= grid1ElementCorners[1][0]);
  assert(grid2ElementCorners[0][0] <= grid2ElementCorners[1][0]);

  T lowerBound = std::max(grid1ElementCorners[0][0], grid2ElementCorners[0][0]);
  T upperBound = std::min(grid1ElementCorners[1][0], grid2ElementCorners[1][0]);

  if (lowerBound <= upperBound) {    // Intersection is non-empty

    intersections.push_back(RemoteSimplicialIntersection<T,dim,dim,dim>());

    // Compute local coordinates in the grid1 element
    intersections.back().grid1Local_[0] = grid1Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
    intersections.back().grid1Local_[1] = grid1Geometry.local(Dune::FieldVector<T,dim>(upperBound));

    // Compute local coordinates in the grid2 element
    intersections.back().grid2Local_[0] = grid2Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
    intersections.back().grid2Local_[1] = grid2Geometry.local(Dune::FieldVector<T,dim>(upperBound));

    // Set indices
    intersections.back().grid1Entity_ = grid1Index;
    intersections.back().grid2Entity_ = grid2Index;

    //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

  }

}


// Explicit instantiation
#ifdef CGAL_EXTERN
#define DECL extern
#else
#define DECL
#endif
template class CGALMergeImp<1,double>;
template class CGALMergeImp<2,double>;
template class CGALMergeImp<3,double>;
#undef DECL
