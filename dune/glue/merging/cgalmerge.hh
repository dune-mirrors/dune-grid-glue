// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Implementation of the Merger concept using the CGAL library
 */

#ifndef CGAL_MERGE_HH
#define CGAL_MERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <set>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/glue/merging/merger.hh>

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
  : public Merger<T,dim,dim,dim>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

private:

  struct RemoteSimplicialIntersection
  {
    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,dim>, dim+1> domainLocal_;

    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,dim>, dim+1> targetLocal_;

    //
    int domainEntity_;

    int targetEntity_;

  };

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersection(const Dune::GeometryType& domainElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& domainElementCorners,
                           unsigned int domainIndex,
                           const Dune::GeometryType& targetElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& targetElementCorners,
                           unsigned int targetIndex);

private:

  /*   M E M B E R   V A R I A B L E S   */

  /** \brief The computed intersections */
  std::vector<RemoteSimplicialIntersection> intersections_;

public:

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @brief builds the merged grid
   *
   * Note that the indices are used consequently throughout the whole class interface just like they are
   * introduced here.
   *
   * @param domain_coords the domain vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param domain_simplices array with all domain simplices represented as corner indices into @c domain_coords;
   * the simplices are just written to this array one after another
   * @param target_coords the target vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param target_simplices just like with the domain_simplices and domain_coords
   */
  void build(const std::vector<Dune::FieldVector<T,dim> >& domainCoords,
             const std::vector<unsigned int>& domain_elements,
             const std::vector<Dune::GeometryType>& domain_element_types,
             const std::vector<Dune::FieldVector<T,dim> >& targetCoords,
             const std::vector<unsigned int>& target_elements,
             const std::vector<Dune::GeometryType>& target_element_types
             );


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const;

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


  /*   M A P P I N G   O N   I N D E X   B A S I S   */

  /**
   * @brief get index of domain parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the domain parent simplex
   */
  unsigned int domainParent(unsigned int idx) const;

  /**
   * @brief get index of target parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the target parent simplex
   */
  unsigned int targetParent(unsigned int idx) const;

  /**
   * @brief get the merged grid simplices refining a given domain simplex
   * @param idx index of domain simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const;

  /**
   * @brief get the merged grid simplices refining a given target simplex
   * @param idx index of target simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const;


  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the domain parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "domainParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent domain simplex
   */
  LocalCoords domainParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent target simplex
   */
  LocalCoords targetParentLocal(unsigned int idx, unsigned int corner) const;

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
  //         corners[i] = this->_coords[this->subEntities_[index].corners[i].idx].coord;

  Geometry domainGeometry(domainElementType, domainElementCorners);
  Geometry targetGeometry(targetElementType, targetElementCorners);

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

      intersections_.push_back(RemoteSimplicialIntersection());

      // Compute local coordinates in the domain element
      intersections_.back().domainLocal_[0] = domainGeometry.local(Dune::FieldVector<T,dim>(lowerBound));
      intersections_.back().domainLocal_[1] = domainGeometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Compute local coordinates in the target element
      intersections_.back().targetLocal_[0] = targetGeometry.local(Dune::FieldVector<T,dim>(lowerBound));
      intersections_.back().targetLocal_[1] = targetGeometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Set indices
      intersections_.back().domainEntity_ = domainIndex;
      intersections_.back().targetEntity_ = targetIndex;

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
    for (std::size_t i=0; i<domainElementCorners.size(); i++)
      P.push_back (Point_2 (domainElementCorners[i][0], domainElementCorners[i][1]));

    //std::cout << "P = "; print_polygon (P);

    Polygon_2 Q;
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

        intersections_.push_back(RemoteSimplicialIntersection());

        // Compute local coordinates in the domain element
        intersections_.back().domainLocal_[0] = domainGeometry.local(anchorFieldVector);
        intersections_.back().domainLocal_[1] = domainGeometry.local(nextFieldVector);
        intersections_.back().domainLocal_[2] = domainGeometry.local(nextNextFieldVector);

        // Compute local coordinates in the domain element
        intersections_.back().targetLocal_[0] = targetGeometry.local(anchorFieldVector);
        intersections_.back().targetLocal_[1] = targetGeometry.local(nextFieldVector);
        intersections_.back().targetLocal_[2] = targetGeometry.local(nextNextFieldVector);

        // Set indices
        intersections_.back().domainEntity_ = domainIndex;
        intersections_.back().targetEntity_ = targetIndex;

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
void CGALMerge<dim, T>::build(const std::vector<Dune::FieldVector<T,dim> >& domainCoords,
                              const std::vector<unsigned int>& domain_elements,
                              const std::vector<Dune::GeometryType>& domain_element_types,
                              const std::vector<Dune::FieldVector<T,dim> >& targetCoords,
                              const std::vector<unsigned int>& target_elements,
                              const std::vector<Dune::GeometryType>& target_element_types
                              )
{

  std::cout << "CGALMerge building merged grid..." << std::endl;

  // /////////////////////////////////////////////////////////////////////
  //   Compute the intersection of all pairs of elements
  //   \todo This is only the naive quadratic algorithm
  // /////////////////////////////////////////////////////////////////////

  unsigned int domainCornerCounter = 0;

  for (std::size_t i=0; i<domain_element_types.size(); i++) {

    // Select vertices of the domain element
    int domainNumVertices = Dune::GenericReferenceElements<T,dim>::general(domain_element_types[i]).size(dim);
    std::vector<Dune::FieldVector<T,dim> > domainElementCorners(domainNumVertices);
    for (int j=0; j<domainNumVertices; j++)
      domainElementCorners[j] = domainCoords[domain_elements[domainCornerCounter++]];

    unsigned int targetCornerCounter = 0;

    for (std::size_t j=0; j<target_element_types.size(); j++) {

      // Select vertices of the domain element
      int targetNumVertices = Dune::GenericReferenceElements<T,dim>::general(target_element_types[j]).size(dim);
      std::vector<Dune::FieldVector<T,dim> > targetElementCorners(targetNumVertices);
      for (int k=0; k<targetNumVertices; k++)
        targetElementCorners[k] = targetCoords[target_elements[targetCornerCounter++]];

      // ///////////////////////////////////////////////////////
      //   Compute the intersection between the two elements
      // ///////////////////////////////////////////////////////

      computeIntersection(domain_element_types[i], domainElementCorners, i,
                          target_element_types[j], targetElementCorners, j);


    }

  }
}


template<int dim, typename T>
inline unsigned int CGALMerge<dim, T>::nSimplices() const
{
  return intersections_.size();
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


template<int dim, typename T>
inline unsigned int CGALMerge<dim, T>::domainParent(unsigned int idx) const
{
  return intersections_[idx].domainEntity_;
}


template<int dim, typename T>
inline unsigned int CGALMerge<dim, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return intersections_[idx].targetEntity_;
}


template<int dim, typename T>
bool CGALMerge<dim, T>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  // WARNING: stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if (intersections_[i].domainEntity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<int dim, typename T>
bool CGALMerge<dim, T>::targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  // WARNING: stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if (intersections_[i].targetEntity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<int dim, typename T>
typename CGALMerge<dim, T>::LocalCoords CGALMerge<dim, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].domainLocal_[corner];
}


template<int dim, typename T>
typename CGALMerge<dim, T>::LocalCoords CGALMerge<dim, T>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].targetLocal_[corner];
}

#endif // CONFORMING_MERGE_HH
