// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CGAL_EXTERN
#include "config.h"
#endif

#include <dune/common/bitsetvector.hh>

#include <dune/grid-glue/merging/cgalmergeimp.hh>



#if HAVE_CGAL  // without CGAL we can still handle 1d problems
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

#ifdef CGAL_USE_GMP

// GMP is installed. Use the GMP rational number-type.
  #include <CGAL/Gmpq.h>

typedef CGAL::Gmpq Exact_number_type;

#else

// GMP is not installed. Use CGAL's exact rational number-type.
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Exact_number_type;

#endif // CGAL_USE_GMP

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
 * \tparam Dune_number_type The type used by Dune for coordinates
 * \tparam CGAL_number_type The type used by CGAL for coordinates
 * \tparam dim The element dimension.  It must be 3.
 *         It is a template parameter only to make the code compile
 */
template <int dim, class Dune_number_type, class CGAL_number_type>
void CGALMergeImp<dim,Dune_number_type,CGAL_number_type>::makeHexahedron(Polyhedron_3& P,
                                                                         const std::vector<Dune::FieldVector<Dune_number_type,dim> >& c)
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

//////////////////////////////////////////////////////////////////////
//  Compute for each face of each element whether it will also intersect
//////////////////////////////////////////////////////////////////////

/** \todo This method should really be a private class method, but then we'd get in trouble
 * because we don't want cgal stuff in cgalmerge.hh at all (because of compile time).
 * I guess the proper solution is a separate implementation class.
 */
template<int dim, class Dune_number_type, class CGAL_number_type>
void CGALMergeImp<dim,Dune_number_type,CGAL_number_type>::computeNeighborIntersections(const Dune::GeometryType& elementType,
                                                                                       const std::vector<Dune::FieldVector<Dune_number_type,dim> >& elementCorners,
                                                                                       const Nef_Polyhedron_3& NQ,
                                                                                       const Nef_Polyhedron_3& intersection,
                                                                                       std::bitset<(1<<dim)>& neighborIntersects
                                                                                       )
{
  const Dune::GenericReferenceElement<Dune_number_type,dim>& refElement = Dune::GenericReferenceElements<Dune_number_type,dim>::general(elementType);

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

#endif // HAVE_CGAL

template<int dim, class Dune_number_type, class CGAL_number_type>
void CGALMergeImp<dim,Dune_number_type,CGAL_number_type>::compute1dIntersection(const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                                                                unsigned int grid1Index,
                                                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                                                                unsigned int grid2Index,
                                                                                std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                                                                )
{
  // Check consistent orientation
  // \todo Reverse the orientation if this check fails
  assert(grid1ElementCorners[0][0] <= grid1ElementCorners[1][0]);
  assert(grid2ElementCorners[0][0] <= grid2ElementCorners[1][0]);

  Dune_number_type lowerBound = std::max(grid1ElementCorners[0][0], grid2ElementCorners[0][0]);
  Dune_number_type upperBound = std::min(grid1ElementCorners[1][0], grid2ElementCorners[1][0]);

  if (lowerBound <= upperBound) {    // Intersection is non-empty

    intersections.push_back(RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim>());

    // Compute local coordinates in the grid1 element
    intersections.back().grid1Local_[0] = grid1Geometry.local(Dune::FieldVector<Dune_number_type,dim>(lowerBound));
    intersections.back().grid1Local_[1] = grid1Geometry.local(Dune::FieldVector<Dune_number_type,dim>(upperBound));

    // Compute local coordinates in the grid2 element
    intersections.back().grid2Local_[0] = grid2Geometry.local(Dune::FieldVector<Dune_number_type,dim>(lowerBound));
    intersections.back().grid2Local_[1] = grid2Geometry.local(Dune::FieldVector<Dune_number_type,dim>(upperBound));

    // Set indices
    intersections.back().grid1Entity_ = grid1Index;
    intersections.back().grid2Entity_ = grid2Index;

    //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

  }

}

#if HAVE_CGAL

template<int dim, class Dune_number_type, class CGAL_number_type>
void CGALMergeImp<dim,Dune_number_type,CGAL_number_type>::compute2dIntersection(const Dune::GeometryType& grid1ElementType,
                                                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                                                                unsigned int grid1Index,
                                                                                const Dune::GeometryType& grid2ElementType,
                                                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                                                                unsigned int grid2Index,
                                                                                std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                                                                )
{
  typedef std::list<Polygon_with_holes_2>            Pwh_list_2;

  // Construct the two input polygons.
  Polygon_2 P;
  if (grid1ElementType.isQuadrilateral()) {
    // Vertex renumbering Dune --> CGAL
    P.push_back( Point_2(grid1ElementCorners[0][0], grid1ElementCorners[0][1]));
    P.push_back( Point_2(grid1ElementCorners[1][0], grid1ElementCorners[1][1]));
    P.push_back( Point_2(grid1ElementCorners[3][0], grid1ElementCorners[3][1]));
    P.push_back( Point_2(grid1ElementCorners[2][0], grid1ElementCorners[2][1]));

  } else
    for (std::size_t i=0; i<grid1ElementCorners.size(); i++)
      P.push_back (Point_2 (grid1ElementCorners[i][0], grid1ElementCorners[i][1]));

  //std::cout << "P = "; print_polygon (P);

  Polygon_2 Q;
  if (grid2ElementType.isQuadrilateral()) {
    // Vertex renumbering Dune --> CGAL
    Q.push_back( Point_2(grid2ElementCorners[0][0], grid2ElementCorners[0][1]));
    Q.push_back( Point_2(grid2ElementCorners[1][0], grid2ElementCorners[1][1]));
    Q.push_back( Point_2(grid2ElementCorners[3][0], grid2ElementCorners[3][1]));
    Q.push_back( Point_2(grid2ElementCorners[2][0], grid2ElementCorners[2][1]));

  } else
    for (std::size_t i=0; i<grid2ElementCorners.size(); i++)
      Q.push_back (Point_2 (grid2ElementCorners[i][0], grid2ElementCorners[i][1]));

  //std::cout << "Q = "; print_polygon (Q);

  // Compute the intersection of P and Q.
  Pwh_list_2 intR;
  typename Pwh_list_2::const_iterator it;

  CGAL::intersection (P, Q, std::back_inserter(intR));

  //std::cout << "The intersection:" << std::endl;
  for (it = intR.begin(); it != intR.end(); ++it) {
    //             std::cout << "--> ";
    //             print_polygon_with_holes (*it);

    // Intersections of bounded polygons must be bounded
    assert(!it->is_unbounded());

    // We cannot have holes
    assert (it->number_of_holes() == 0);

    // Today I am too lazy to program anything more general than convex intersections
    if (not it->outer_boundary().is_convex())
      DUNE_THROW(Dune::NotImplemented, "CGAL merging with non-convex intersections");

    // Less than 3 vertices would be a bug
    assert(it->outer_boundary().size() >= 3);

    // Now we split the intersection into triangles.  We use the simplest algorithm
    // conceivable, because triangle quality is not important here.

    typename Polygon_2::Vertex_const_iterator anchor = it->outer_boundary().vertices_begin();

    typename Polygon_2::Vertex_const_iterator next = anchor;
    ++next;

    typename Polygon_2::Vertex_const_iterator nextNext = next;
    ++nextNext;

    do {

      // Make Dune types from CGAL types

      // The following check is there, because we are using CGAL::to_double().
      // If necessary the code can be generalized somewhat.
      dune_static_assert((Dune::is_same<Dune_number_type,double>::value), "Dune_number_type must be 'double'");


      Dune::FieldVector<Dune_number_type,dim> anchorFieldVector;
      anchorFieldVector[0] = CGAL::to_double(anchor->x());
      anchorFieldVector[1] = CGAL::to_double(anchor->y());

      Dune::FieldVector<Dune_number_type,dim> nextFieldVector;
      nextFieldVector[0] = CGAL::to_double(next->x());
      nextFieldVector[1] = CGAL::to_double(next->y());

      Dune::FieldVector<Dune_number_type,dim> nextNextFieldVector;
      nextNextFieldVector[0] = CGAL::to_double(nextNext->x());
      nextNextFieldVector[1] = CGAL::to_double(nextNext->y());

      // ///////////////////////////////////////////////////
      // Output the triangle (anchor, next, nextNext)
      // ///////////////////////////////////////////////////

      intersections.push_back(RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim>());

      // Compute local coordinates in the grid1 element
      intersections.back().grid1Local_[0] = grid1Geometry.local(anchorFieldVector);
      intersections.back().grid1Local_[1] = grid1Geometry.local(nextFieldVector);
      intersections.back().grid1Local_[2] = grid1Geometry.local(nextNextFieldVector);

      // Compute local coordinates in the grid1 element
      intersections.back().grid2Local_[0] = grid2Geometry.local(anchorFieldVector);
      intersections.back().grid2Local_[1] = grid2Geometry.local(nextFieldVector);
      intersections.back().grid2Local_[2] = grid2Geometry.local(nextNextFieldVector);

      // Set indices
      intersections.back().grid1Entity_ = grid1Index;
      intersections.back().grid2Entity_ = grid2Index;

      //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

      // move to the next triangle
      ++next;
      ++nextNext;

    } while (nextNext != it->outer_boundary().vertices_end());

  }

}


template<int dim, class Dune_number_type, class CGAL_number_type>
void CGALMergeImp<dim,Dune_number_type,CGAL_number_type>::compute3dIntersection(const Dune::GeometryType& grid1ElementType,
                                                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                                                                unsigned int grid1Index,
                                                                                std::bitset<(1<<dim)>& neighborIntersects1,
                                                                                const Dune::GeometryType& grid2ElementType,
                                                                                const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                                                                const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                                                                unsigned int grid2Index,
                                                                                std::bitset<(1<<dim)>& neighborIntersects2,
                                                                                std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                                                                )
{
  // Construct the two input polyhedra.
  Polyhedron_3 P;

  if (grid1ElementType.isSimplex()) {
    P.make_tetrahedron(Point_3 (grid1ElementCorners[0][0], grid1ElementCorners[0][1], grid1ElementCorners[0][2]),
                       Point_3 (grid1ElementCorners[1][0], grid1ElementCorners[1][1], grid1ElementCorners[1][2]),
                       Point_3 (grid1ElementCorners[2][0], grid1ElementCorners[2][1], grid1ElementCorners[2][2]),
                       Point_3 (grid1ElementCorners[3][0], grid1ElementCorners[3][1], grid1ElementCorners[3][2]));
  } if (grid1ElementType.isCube()) {

    makeHexahedron(P, grid1ElementCorners);

  } else
    DUNE_THROW(Dune::GridError, "Type " << grid1ElementType << " not supported by CGALMerge yet");
  //std::cout << "P = " << P << std::endl;

  // Turn polyhedron into a Nef polyhedron
  Nef_Polyhedron_3 NP(P);
  Polyhedron_3 Q;

  if (grid2ElementType.isSimplex()) {
    Q.make_tetrahedron(Point_3 (grid2ElementCorners[0][0], grid2ElementCorners[0][1], grid2ElementCorners[0][2]),
                       Point_3 (grid2ElementCorners[1][0], grid2ElementCorners[1][1], grid2ElementCorners[1][2]),
                       Point_3 (grid2ElementCorners[2][0], grid2ElementCorners[2][1], grid2ElementCorners[2][2]),
                       Point_3 (grid2ElementCorners[3][0], grid2ElementCorners[3][1], grid2ElementCorners[3][2]));
  } if (grid2ElementType.isCube()) {

    makeHexahedron(Q, grid2ElementCorners);

  } else
    DUNE_THROW(Dune::GridError, "Type " << grid2ElementType << " not supported by CGALMerge yet");

  Nef_Polyhedron_3 NQ(Q);
  //std::cout << "Q = " << Q << std::endl;

  //////////////////////////////////////////////////////////
  // Compute the intersection of P and Q.
  //////////////////////////////////////////////////////////
  Nef_Polyhedron_3 intersection = NP * NQ;

  if (intersection.is_empty())
    return;

  Polyhedron_3 intersectionP;
  if(intersection.is_simple()) {
    intersection.convert_to_polyhedron(intersectionP);
    //std::cout << intersectionP;
  } else
    std::cerr << "N1 is not a 2-manifold." << std::endl;

  //std::cout << "Intersection has " << intersection.number_of_vertices() << " vertices\n";

  //////////////////////////////////////////////////////////////////////
  //  Compute for each face of each element whether it will also intersect
  //////////////////////////////////////////////////////////////////////

  computeNeighborIntersections(grid1ElementType, grid1ElementCorners, NQ, intersection, neighborIntersects1);

  computeNeighborIntersections(grid2ElementType, grid2ElementCorners, NP, intersection, neighborIntersects2);

  //////////////////////////////////////////////////////////////////////
  //  Triangulate the intersection polyhedron.
  //  For this we first compute the centroid to have a point that
  //  is certainly within the intersection.  (this of course presupposes
  //  that the intersection is convex).  Then we connect the centroid
  //  to all facets.
  //////////////////////////////////////////////////////////////////////

  assert(intersectionP.is_closed());

  // If the intersection is a tetrahedron we take a shortcut
  if (intersectionP.is_tetrahedron(intersectionP.halfedges_begin())) {

    // Make Dune types from CGAL types

    // The following check is there, because we are using CGAL::to_double().
    // If necessary the code can be generalized somewhat.
    dune_static_assert((Dune::is_same<Dune_number_type,double>::value), "Dune_number_type must be 'double'");

    Dune::array<Dune::FieldVector<Dune_number_type,dim>, 4> duneVertexPos;

    typename Polyhedron_3::Point_iterator vIt = intersectionP.points_begin();

    for (int i=0; i<4; ++i, ++vIt)
      for (int j=0; j<3; j++)
        duneVertexPos[i][j] = CGAL::to_double((*vIt)[j]);

    assert(vIt==intersectionP.points_end());

    // ///////////////////////////////////////////////////
    // Output the tetrahedron
    // ///////////////////////////////////////////////////

    intersections.push_back(RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim>());

    // Compute local coordinates in the grid1 and grid2 elements
    for (int i=0; i<4; i++) {
      intersections.back().grid1Local_[i] = grid1Geometry.local(duneVertexPos[i]);
      intersections.back().grid2Local_[i] = grid2Geometry.local(duneVertexPos[i]);
    }

    // Set indices
    intersections.back().grid1Entity_ = grid1Index;
    intersections.back().grid2Entity_ = grid2Index;

    //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

  } else {

    //////////////////////////////////////////////////////////////////
    //   Compute the centroid
    //////////////////////////////////////////////////////////////////

    Dune::FieldVector<Dune_number_type,dim> centroid(0);

    for (typename Polyhedron_3::Point_iterator vIt = intersectionP.points_begin();
         vIt != intersectionP.points_end();
         ++vIt) {
      for (int i=0; i<dim; i++)
        centroid[i] += CGAL::to_double((*vIt)[i]);
    }

    centroid /= intersectionP.size_of_vertices();

    //////////////////////////////////////////////////////////////////////////
    //  Loop over all facets, triangulate the facet and create
    //  an intersection from each such triangle together with the centroid
    //////////////////////////////////////////////////////////////////////////

    for (typename Polyhedron_3::Facet_iterator fIt = intersectionP.facets_begin();
         fIt != intersectionP.facets_end();
         ++fIt) {

      ////////////////////////////////////////////////////////////
      // Get the vertices of this facet in circular order
      ////////////////////////////////////////////////////////////
      unsigned int nVertices = fIt->facet_degree();

      std::vector<Dune::FieldVector<Dune_number_type,dim> > facetVertices(nVertices);

      typename Polyhedron_3::Facet::Halfedge_around_facet_circulator fVCirc = fIt->facet_begin();
      assert (fVCirc != 0);        // an empty circulator would be an error
      int i=0;         // array iterator

      do {

        for (int j=0; j<dim; j++)
          facetVertices[i][j] = CGAL::to_double(fVCirc->vertex()->point()[j]);

        ++i;

      } while (++fVCirc != fIt->facet_begin());

      /////////////////////////////////////////////////////////////////////////
      //  Triangulate the facet and enter an intersection for each triangle
      /////////////////////////////////////////////////////////////////////////

      typename Polyhedron_3::Facet::Halfedge_around_facet_circulator anchor = fIt->facet_begin();

      typename Polyhedron_3::Facet::Halfedge_around_facet_circulator next = anchor;
      ++next;

      typename Polyhedron_3::Facet::Halfedge_around_facet_circulator nextNext = next;
      ++nextNext;

      do {

        // Make Dune types from CGAL types

        // The following check is there, because we are using CGAL::to_double().
        // If necessary the code can be generalized somewhat.
        dune_static_assert((Dune::is_same<Dune_number_type,double>::value), "Dune_number_type must be 'double'");

        Dune::FieldVector<Dune_number_type,dim> anchorFieldVector;
        Dune::FieldVector<Dune_number_type,dim> nextFieldVector;
        Dune::FieldVector<Dune_number_type,dim> nextNextFieldVector;

        for (int i=0; i<dim; i++) {
          anchorFieldVector[i]   = CGAL::to_double(anchor->vertex()->point()[i]);
          nextFieldVector[i]     = CGAL::to_double(next->vertex()->point()[i]);
          nextNextFieldVector[i] = CGAL::to_double(nextNext->vertex()->point()[i]);
        }

        // ////////////////////////////////////////////////////////////
        // Output the tetrahedron (anchor, next, nextNext, centroid)
        // ////////////////////////////////////////////////////////////

        intersections.push_back(RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim>());

        // Compute local coordinates in the grid1 element
        intersections.back().grid1Local_[0] = grid1Geometry.local(anchorFieldVector);
        intersections.back().grid1Local_[1] = grid1Geometry.local(nextFieldVector);
        intersections.back().grid1Local_[2] = grid1Geometry.local(nextNextFieldVector);
        intersections.back().grid1Local_[3] = grid1Geometry.local(centroid);

        // Compute local coordinates in the grid1 element
        intersections.back().grid2Local_[0] = grid2Geometry.local(anchorFieldVector);
        intersections.back().grid2Local_[1] = grid2Geometry.local(nextFieldVector);
        intersections.back().grid2Local_[2] = grid2Geometry.local(nextNextFieldVector);
        intersections.back().grid2Local_[3] = grid2Geometry.local(centroid);

        // Set indices
        intersections.back().grid1Entity_ = grid1Index;
        intersections.back().grid2Entity_ = grid2Index;

        //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

        // move to the next triangle
        ++next;
        ++nextNext;

      } while (nextNext != fIt->facet_begin());

    }

  }

}

#endif // HAVE_CGAL

// Explicit instantiation
#if HAVE_CGAL
#ifndef CGAL_EXTERN
/*template class CGALMergeImp<1,double,double>;
   template class CGALMergeImp<2,double,double>;
   template class CGALMergeImp<3,double,double>;*/
template class CGALMergeImp<1,double,Exact_number_type>;
template class CGALMergeImp<2,double,Exact_number_type>;
template class CGALMergeImp<3,double,Exact_number_type>;
#endif
#else
template class CGALMergeImp<1,double,double>;
#endif // HAVE_CGAL
