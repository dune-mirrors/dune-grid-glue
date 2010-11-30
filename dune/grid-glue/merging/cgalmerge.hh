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

private:

  typedef typename StandardMerge<T,dim,dim,dim>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersection(const Dune::GeometryType& grid1ElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& grid1ElementCorners,
                           unsigned int grid1Index,
                           const Dune::GeometryType& grid2ElementType,
                           const std::vector<Dune::FieldVector<T,dim> >& grid2ElementCorners,
                           unsigned int grid2Index);

private:

  /**
   * @brief check if given grid1 simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a grid1 simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "grid1LocalToMerged" or "grid1GlobalToMerged".
   * @param idx the index of the grid1 simplex
   * @return TRUE <=> refined in merged grid
   */
  bool grid1SimplexMatched(unsigned int idx) const;

  /**
   * @brief check if given grid2 simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a grid2 simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "grid2LocalToMerged" or "grid2GlobalToMerged".
   * @param idx the index of the grid2 simplex
   * @return TRUE <=> refined in merged grid
   */
  bool grid2SimplexMatched(unsigned int idx) const;

};

#ifdef CGAL_EXTRA_TYPES
#define CGAL_EXTERN
#include "cgalmerge.cc"
#undef CGAL_EXTERN
#endif

#endif // CGAL_MERGE_HH
