// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Implementation of the Merger concept using the CGAL library
 */

#ifndef CGAL_MERGE_HH
#define CGAL_MERGE_HH


#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/genericgeometry/geometry.hh>

#include <dune/grid-glue/merging/standardmerge.hh>

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


/** \brief Implementation of the Merger concept using the CGAL library

   \tparam dim Grid dimension of the coupling grids.  The world dimension is assumed to be the same.
   \tparam T Type used for coordinates
 */
template<int dim, typename T = double>
class CGALMerge
  : public StandardMerge<T,dim,dim,dim>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

protected:

  typedef typename StandardMerge<T,dim,dim,dim>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersection(const Dune::GeometryType& domainElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& domainElementCorners,
                           unsigned int domainIndex,
                           const Dune::GeometryType& targetElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& targetElementCorners,
                           unsigned int targetIndex);

public:

  /**
   * @brief check if given domain simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a domain simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "domainLocalToMerged" or "domainGlobalToMerged".
   * @param idx the index of the domain simplex
   * @return TRUE <=> refined in merged grid
   */
  bool domainSimplexMatched(unsigned int idx) const;

  /**
   * @brief check if given target simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a target simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "targetLocalToMerged" or "targetGlobalToMerged".
   * @param idx the index of the target simplex
   * @return TRUE <=> refined in merged grid
   */
  bool targetSimplexMatched(unsigned int idx) const;

};


/* IMPLEMENTATION */

template<int dim, typename T>
void CGALMerge<dim, T>::
computeIntersection(const Dune::GeometryType& domainElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& domainElementCorners,
                    unsigned int domainIndex,
                    const Dune::GeometryType& targetElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& targetElementCorners,
                    unsigned int targetIndex)
{

  // A few consistency checks
  assert((Dune::GenericReferenceElements<T,dim>::general(domainElementType).size(dim) == domainElementCorners.size()));
  assert((Dune::GenericReferenceElements<T,dim>::general(targetElementType).size(dim) == targetElementCorners.size()));

  // Make generic geometries representing the domain- and target element.
  // this eases computation of local coordinates.
  typedef Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<T,dim,dim> > Geometry;

  //     std::vector<Coords> corners(subEntities_[index].nCorners());
  //     for (int i = 0; i < subEntities_[index].nCorners(); ++i)
  //         corners[i] = this->coords_[this->subEntities_[index].corners[i].idx].coord;

#define DUNE_GRID_VERSION_NUMBER (DUNE_GRID_VERSION_MAJOR * 10 + DUNE_GRID_VERSION_MINOR)
#if DUNE_GRID_VERSION_NUMBER > 20
  Geometry domainGeometry(Dune::GenericGeometry::topologyId(domainElementType), domainElementCorners);
  Geometry targetGeometry(Dune::GenericGeometry::topologyId(targetElementType), targetElementCorners);
#else
  Geometry domainGeometry(domainElementType, domainElementCorners);
  Geometry targetGeometry(targetElementType, targetElementCorners);
#endif

  // /////////////////////////////////////////////////////////////////////////////////////
  //   Compute the intersection between the two elements.  The 1d case is implemented
  //   by hand; 2d and 3d use CGAL.
  // /////////////////////////////////////////////////////////////////////////////////////

  switch (dim) {
  case 1 : {

    // Check consistent orientation
    // \todo Reverse the orientation if this check fails
    assert(domainElementCorners[0][0] <= domainElementCorners[1][0]);
    assert(targetElementCorners[0][0] <= targetElementCorners[1][0]);

    T lowerBound = std::max(domainElementCorners[0][0], targetElementCorners[0][0]);
    T upperBound = std::min(domainElementCorners[1][0], targetElementCorners[1][0]);

    if (lowerBound <= upperBound) {      // Intersection is non-empty

      this->intersections_.push_back(RemoteSimplicialIntersection());

      // Compute local coordinates in the domain element
      this->intersections_.back().domainLocal_[0] = domainGeometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().domainLocal_[1] = domainGeometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Compute local coordinates in the target element
      this->intersections_.back().targetLocal_[0] = targetGeometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().targetLocal_[1] = targetGeometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Set indices
      this->intersections_.back().domainEntity_ = domainIndex;
      this->intersections_.back().targetEntity_ = targetIndex;

      std::cout << "Intersection between elements " << domainIndex << " and " << targetIndex << std::endl;

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
    if (domainElementType.isQuadrilateral()) {
      // Vertex renumbering Dune --> CGAL
      P.push_back( Point_2(domainElementCorners[0][0], domainElementCorners[0][1]));
      P.push_back( Point_2(domainElementCorners[1][0], domainElementCorners[1][1]));
      P.push_back( Point_2(domainElementCorners[3][0], domainElementCorners[3][1]));
      P.push_back( Point_2(domainElementCorners[2][0], domainElementCorners[2][1]));

    } else
      for (std::size_t i=0; i<domainElementCorners.size(); i++)
        P.push_back (Point_2 (domainElementCorners[i][0], domainElementCorners[i][1]));

    //std::cout << "P = "; print_polygon (P);

    Polygon_2 Q;
    if (targetElementType.isQuadrilateral()) {
      // Vertex renumbering Dune --> CGAL
      Q.push_back( Point_2(targetElementCorners[0][0], targetElementCorners[0][1]));
      Q.push_back( Point_2(targetElementCorners[1][0], targetElementCorners[1][1]));
      Q.push_back( Point_2(targetElementCorners[3][0], targetElementCorners[3][1]));
      Q.push_back( Point_2(targetElementCorners[2][0], targetElementCorners[2][1]));

    } else
      for (std::size_t i=0; i<targetElementCorners.size(); i++)
        Q.push_back (Point_2 (targetElementCorners[i][0], targetElementCorners[i][1]));

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

        // Compute local coordinates in the domain element
        this->intersections_.back().domainLocal_[0] = domainGeometry.local(anchorFieldVector);
        this->intersections_.back().domainLocal_[1] = domainGeometry.local(nextFieldVector);
        this->intersections_.back().domainLocal_[2] = domainGeometry.local(nextNextFieldVector);

        // Compute local coordinates in the domain element
        this->intersections_.back().targetLocal_[0] = targetGeometry.local(anchorFieldVector);
        this->intersections_.back().targetLocal_[1] = targetGeometry.local(nextFieldVector);
        this->intersections_.back().targetLocal_[2] = targetGeometry.local(nextNextFieldVector);

        // Set indices
        this->intersections_.back().domainEntity_ = domainIndex;
        this->intersections_.back().targetEntity_ = targetIndex;

        //std::cout << "Intersection between elements " << domainIndex << " and " << targetIndex << std::endl;

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
inline bool CGALMerge<dim, T>::domainSimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all domain entities
  return true;
}


template<int dim, typename T>
inline bool CGALMerge<dim, T>::targetSimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all target entities
  return true;
}

#endif // CGAL_MERGE_HH
