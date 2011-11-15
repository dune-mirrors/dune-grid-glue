// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef OVERLAPPING_MERGE_HH
#define OVERLAPPING_MERGE_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/geometry/genericgeometry/geometry.hh>

#include <dune/grid-glue/merging/standardmerge.hh>

/** \brief Computing overlapping grid intersections without using an external library

   Implementation of computeIntersection, described in 'An Algorithm for Non-Matching Grid Projections with Linear Complexity,
   M.J. Gander and C. Japhet, Domain Decomposition Methods in Science and Engineering XVIII, pp. 185--192, Springer-Verlag, 2009.'

   \tparam dim Grid dimension of the coupling grids.  The world dimension is assumed to be the same.
   \tparam T Type used for coordinates
 */
template<int dim, typename T = double>
class OverlappingMerge
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

private:

  typedef typename StandardMerge<T,dim,dim,dim>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersection(const Dune::GeometryType& grid1ElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& grid1ElementCorners,
                           unsigned int grid1Index,
                           std::bitset<(1<<dim)>& neighborIntersects1,
                           const Dune::GeometryType& grid2ElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& grid2ElementCorners,
                           unsigned int grid2Index,
                           std::bitset<(1<<dim)>& neighborIntersects2);

private:

  void EdgeIntersections( const std::vector<Dune::FieldVector<T,dim> >   X,
                          const std::vector<Dune::FieldVector<T,dim> >   Y,
                          std::vector<Dune::FieldVector<T,dim> > & P ) ;

  void PointsofXinY( const std::vector<Dune::FieldVector<T,dim> >   X,
                     const std::vector<Dune::FieldVector<T,dim> >   Y,
                     std::vector<Dune::FieldVector<T,dim> > & P ) ;

  void SortAndRemoveDoubles( std::vector<Dune::FieldVector<T,dim> > & P ) ;

};

#include "overlappingmerge.cc"

#endif // OVERLAPPING_MERGE_HH
