// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CGAL_EXTERN
#include "config.h"
#endif

#include <dune/grid-glue/merging/cgalmerge.hh>

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
template <class Polyhedron, class T, int dim>
void makeHexahedron(Polyhedron& P,
                    const std::vector<Dune::FieldVector<T,dim> >& c)
{
  // appends a cube of size [0,1]ˆ3 to the polyhedron P.
  CGAL_precondition( P.is_valid());
  typedef typename Polyhedron::Point_3 Point;
  typedef typename Polyhedron::Halfedge_handle Halfedge_handle;
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


/* IMPLEMENTATION */

template<int dim, typename T>
void CGALMerge<dim, T>::
computeIntersection(const Dune::GeometryType& grid1ElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& grid1ElementCorners,
                    unsigned int grid1Index,
                    const Dune::GeometryType& grid2ElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& grid2ElementCorners,
                    unsigned int grid2Index)
{

  // A few consistency checks
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid1ElementType).size(dim)) == grid1ElementCorners.size());
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid2ElementType).size(dim)) == grid2ElementCorners.size());

  // Make generic geometries representing the grid1- and grid2 element.
  // this eases computation of local coordinates.
  typedef Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<T,dim,dim> > Geometry;

  Geometry grid1Geometry(grid1ElementType, grid1ElementCorners);
  Geometry grid2Geometry(grid2ElementType, grid2ElementCorners);

  // /////////////////////////////////////////////////////////////////////////////////////
  //   Compute the intersection between the two elements.  The 1d case is implemented
  //   by hand; 2d and 3d use CGAL.
  // /////////////////////////////////////////////////////////////////////////////////////

  switch (dim) {
  case 1 : {

    // Check consistent orientation
    // \todo Reverse the orientation if this check fails
    assert(grid1ElementCorners[0][0] <= grid1ElementCorners[1][0]);
    assert(grid2ElementCorners[0][0] <= grid2ElementCorners[1][0]);

    T lowerBound = std::max(grid1ElementCorners[0][0], grid2ElementCorners[0][0]);
    T upperBound = std::min(grid1ElementCorners[1][0], grid2ElementCorners[1][0]);

    if (lowerBound <= upperBound) {      // Intersection is non-empty

      this->intersections_.push_back(RemoteSimplicialIntersection());

      // Compute local coordinates in the grid1 element
      this->intersections_.back().grid1Local_[0] = grid1Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().grid1Local_[1] = grid1Geometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Compute local coordinates in the grid2 element
      this->intersections_.back().grid2Local_[0] = grid2Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().grid2Local_[1] = grid2Geometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Set indices
      this->intersections_.back().grid1Entity_ = grid1Index;
      this->intersections_.back().grid2Entity_ = grid2Index;

      //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

    }

    break;
  }
#ifdef HAVE_CGAL
  case 2 : {

    //typedef T Number_type;

    typedef CGAL::Cartesian<Number_type>               Kernel;
    typedef typename Kernel::Point_2 Point_2;
    typedef CGAL::Polygon_2<Kernel>                    Polygon_2;
    typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;
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
        dune_static_assert((Dune::is_same<T,double>::value), "T must be 'double'");


        Dune::FieldVector<T,dim> anchorFieldVector;
        anchorFieldVector[0] = CGAL::to_double(anchor->x());
        anchorFieldVector[1] = CGAL::to_double(anchor->y());

        Dune::FieldVector<T,dim> nextFieldVector;
        nextFieldVector[0] = CGAL::to_double(next->x());
        nextFieldVector[1] = CGAL::to_double(next->y());

        Dune::FieldVector<T,dim> nextNextFieldVector;
        nextNextFieldVector[0] = CGAL::to_double(nextNext->x());
        nextNextFieldVector[1] = CGAL::to_double(nextNext->y());

        // ///////////////////////////////////////////////////
        // Output the triangle (anchor, next, nextNext)
        // ///////////////////////////////////////////////////

        this->intersections_.push_back(RemoteSimplicialIntersection());

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0] = grid1Geometry.local(anchorFieldVector);
        this->intersections_.back().grid1Local_[1] = grid1Geometry.local(nextFieldVector);
        this->intersections_.back().grid1Local_[2] = grid1Geometry.local(nextNextFieldVector);

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid2Local_[0] = grid2Geometry.local(anchorFieldVector);
        this->intersections_.back().grid2Local_[1] = grid2Geometry.local(nextFieldVector);
        this->intersections_.back().grid2Local_[2] = grid2Geometry.local(nextNextFieldVector);

        // Set indices
        this->intersections_.back().grid1Entity_ = grid1Index;
        this->intersections_.back().grid2Entity_ = grid2Index;

        //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

        // move to the next triangle
        ++next;
        ++nextNext;

      } while (nextNext != it->outer_boundary().vertices_end());

    }

    break;
  }
  case 3 : {

    typedef CGAL::Cartesian<Number_type>   Kernel;
    typedef typename Kernel::Point_3 Point_3;
    typedef CGAL::Polyhedron_3<Kernel>     Polyhedron_3;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_Polyhedron_3;

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
    //std::cout << "P = "; print_polygon (P);

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
    //std::cout << "Q = "; print_polygon (Q);

    //////////////////////////////////////////////////////////
    // Compute the intersection of P and Q.
    //////////////////////////////////////////////////////////
    Nef_Polyhedron_3 intersection = NP * NQ;

    Polyhedron_3 intersectionP;
    if(intersection.is_simple()) {
      intersection.convert_to_polyhedron(intersectionP);
      std::cout << intersectionP;
    } else
      std::cerr << "N1 is not a 2-manifold." << std::endl;

    std::cout << "Intersection has " << intersection.number_of_vertices() << " vertices\n";

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
      dune_static_assert((Dune::is_same<T,double>::value), "T must be 'double'");

      Dune::array<Dune::FieldVector<T,dim>, 4> duneVertexPos;

      typename Polyhedron_3::Point_iterator vIt = intersectionP.points_begin();

      for (int i=0; i<4; ++i, ++vIt)
        for (int j=0; j<3; j++)
          duneVertexPos[i][j] = CGAL::to_double((*vIt)[j]);

      assert(vIt==intersectionP.points_end());

      // ///////////////////////////////////////////////////
      // Output the tetrahedron
      // ///////////////////////////////////////////////////

      this->intersections_.push_back(RemoteSimplicialIntersection());

      // Compute local coordinates in the grid1 and grid2 elements
      for (int i=0; i<4; i++) {
        this->intersections_.back().grid1Local_[i] = grid1Geometry.local(duneVertexPos[i]);
        this->intersections_.back().grid2Local_[i] = grid2Geometry.local(duneVertexPos[i]);
      }

      // Set indices
      this->intersections_.back().grid1Entity_ = grid1Index;
      this->intersections_.back().grid2Entity_ = grid2Index;

      //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

    } else {

      //////////////////////////////////////////////////////////////////
      //   Compute the centroid
      //////////////////////////////////////////////////////////////////

      Dune::FieldVector<T,dim> centroid(0);

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

        std::vector<Dune::FieldVector<T,dim> > facetVertices(nVertices);

        typename Polyhedron_3::Facet::Halfedge_around_facet_circulator fVCirc = fIt->facet_begin();
        assert (fVCirc != 0);          // an empty circulator would be an error
        int i=0;           // array iterator

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
          dune_static_assert((Dune::is_same<T,double>::value), "T must be 'double'");

          Dune::FieldVector<T,dim> anchorFieldVector;
          Dune::FieldVector<T,dim> nextFieldVector;
          Dune::FieldVector<T,dim> nextNextFieldVector;

          for (int i=0; i<dim; i++) {
            anchorFieldVector[i]   = CGAL::to_double(anchor->vertex()->point()[i]);
            nextFieldVector[i]     = CGAL::to_double(next->vertex()->point()[i]);
            nextNextFieldVector[i] = CGAL::to_double(nextNext->vertex()->point()[i]);
          }

          // ////////////////////////////////////////////////////////////
          // Output the tetrahedron (anchor, next, nextNext, centroid)
          // ////////////////////////////////////////////////////////////

          this->intersections_.push_back(RemoteSimplicialIntersection());

          // Compute local coordinates in the grid1 element
          this->intersections_.back().grid1Local_[0] = grid1Geometry.local(anchorFieldVector);
          this->intersections_.back().grid1Local_[1] = grid1Geometry.local(nextFieldVector);
          this->intersections_.back().grid1Local_[2] = grid1Geometry.local(nextNextFieldVector);
          this->intersections_.back().grid1Local_[3] = grid1Geometry.local(centroid);

          // Compute local coordinates in the grid1 element
          this->intersections_.back().grid2Local_[0] = grid2Geometry.local(anchorFieldVector);
          this->intersections_.back().grid2Local_[1] = grid2Geometry.local(nextFieldVector);
          this->intersections_.back().grid2Local_[2] = grid2Geometry.local(nextNextFieldVector);
          this->intersections_.back().grid2Local_[3] = grid2Geometry.local(centroid);

          // Set indices
          this->intersections_.back().grid1Entity_ = grid1Index;
          this->intersections_.back().grid2Entity_ = grid2Index;

          //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

          // move to the next triangle
          ++next;
          ++nextNext;

        } while (nextNext != fIt->facet_begin());

      }

    }

  }
  break;
#endif

  default :
    DUNE_THROW(Dune::NotImplemented, "CGALMerge is not implemented for dim==" << dim << "!");

  }

}

template<int dim, typename T>
inline bool CGALMerge<dim, T>::grid1SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid1 entities
  return true;
}


template<int dim, typename T>
inline bool CGALMerge<dim, T>::grid2SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid2 entities
  return true;
}


// Explicit instantiation
#ifdef CGAL_EXTERN
#define DECL extern
#else
#define DECL
#endif
#define CGAL_INSTANTIATION(D,T)                                          \
  DECL template void CGALMerge<D, T>::computeIntersection(const Dune::GeometryType& grid1ElementType, \
                                                          const std::vector<Dune::FieldVector<T,D> >& grid1ElementCorners, \
                                                          unsigned int grid1Index, \
                                                          const Dune::GeometryType& grid2ElementType, \
                                                          const std::vector<Dune::FieldVector<T,D> >& grid2ElementCorners, \
                                                          unsigned int grid2Index); \
\
  DECL template bool CGALMerge<D, T>::grid1SimplexMatched(unsigned int idx) const; \
\
  DECL template bool CGALMerge<D, T>::grid2SimplexMatched(unsigned int idx) const

CGAL_INSTANTIATION(1, double);
CGAL_INSTANTIATION(2, double);
CGAL_INSTANTIATION(3, double);
// CGAL_INSTANTIATION(1, float);
// CGAL_INSTANTIATION(2, float);
// CGAL_INSTANTIATION(3, float);

#undef DECL
