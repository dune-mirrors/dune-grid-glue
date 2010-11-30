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

      std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

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

#if 0
    typedef CGAL::Homogeneous<CGAL::Gmpz>  Kernel;
    typedef CGAL::Polyhedron_3<Kernel>  Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron;
    typedef Kernel::Vector_3 Vector_3;
    typedef Kernel::Aff_transformation_3 Aff_transformation_3;

    Polyhedron P;
    std::cin >> P;
    if(P.is_closed()) {
      Nef_polyhedron N1(P);
      Nef_polyhedron N2(N1);
      Aff_transformation_3 aff(CGAL::TRANSLATION, Vector_3(2,2,0,1));
      N2.transform(aff);
      N1 += N2;

      if(N1.is_simple()) {
        N1.convert_to_polyhedron(P);
        std::cout << P;
      }
      else
        std::cerr << "N1 is not a 2-manifold." << std::endl;
    }

#endif
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
// CGAL_INSTANTIATION(3, double);
// CGAL_INSTANTIATION(1, float);
// CGAL_INSTANTIATION(2, float);
// CGAL_INSTANTIATION(3, float);

#undef DECL
