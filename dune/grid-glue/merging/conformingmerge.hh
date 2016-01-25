// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    conformingmerge.hh
 *  Version:     1.0
 *  Created on:  Sep 14, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: implementation of the Merger concept for conforming interfaces
 *
 */
/**
 * @file
 * @brief
 * Implementation of the Merger concept for conforming interfaces
 */

#ifndef DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH
#define DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH

#include <iomanip>
#include <vector>
#include <algorithm>
#include <bitset>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/geometry/referenceelements.hh>

#include <dune/grid-glue/merging/standardmerge.hh>

namespace Dune {

  namespace GridGlue {

/** \brief Implementation of the Merger concept for conforming interfaces

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class ConformingMerge
  : public StandardMerge<T,dim,dim,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  T tolerance_;

  typedef typename StandardMerge<T,dim,dim,dimworld>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  void computeIntersections(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects1,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   std::bitset<(1<<dim)>& neighborIntersects2,
                                   unsigned int grid2Index,
                                   std::vector<RemoteSimplicialIntersection>& intersections);

public:

  ConformingMerge(T tolerance = 1E-4) :
    tolerance_(tolerance)
  {}

private:

  /*   M A P P I N G   O N   I N D E X   B A S I S   */

  /**
   * @brief get index of grid1 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid1 parent simplex
   */
  unsigned int grid1Parent(unsigned int idx, unsigned int parId = 0) const;

  /**
   * @brief get index of grid2 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid2 parent simplex
   */
  unsigned int grid2Parent(unsigned int idx, unsigned int parId = 0) const;

  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1 simplex
   */
  LocalCoords grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2 simplex
   */
  LocalCoords grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

};


template<int dim, int dimworld, typename T>
void ConformingMerge<dim, dimworld, T>::computeIntersections(const Dune::GeometryType& grid1ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                                            std::bitset<(1<<dim)>& neighborIntersects1,
                                                            unsigned int grid1Index,
                                                            const Dune::GeometryType& grid2ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                                            std::bitset<(1<<dim)>& neighborIntersects2,
                                                            unsigned int grid2Index,
                                                            std::vector<RemoteSimplicialIntersection>& intersections)
{
  this->counter++;

  // A few consistency checks
  assert((unsigned int)(Dune::ReferenceElements<T,dim>::general(grid1ElementType).size(dim)) == grid1ElementCorners.size());
  assert((unsigned int)(Dune::ReferenceElements<T,dim>::general(grid2ElementType).size(dim)) == grid2ElementCorners.size());
  // any intersection we may find will be the entire elements.
  neighborIntersects1.reset();
  neighborIntersects2.reset();

  // the intersection is either conforming or empty, hence the GeometryTypes have to match
  if (grid1ElementType != grid2ElementType)
    return;

  // ////////////////////////////////////////////////////////////
  //   Find correspondences between the different corners
  // ////////////////////////////////////////////////////////////
  std::vector<int> other(grid1ElementCorners.size(), -1);

  for (unsigned int i=0; i<grid1ElementCorners.size(); i++) {

    for (unsigned int j=0; j<grid2ElementCorners.size(); j++) {

      if ( (grid1ElementCorners[i]-grid2ElementCorners[j]).two_norm() < tolerance_ ) {

        other[i] = j;
        break;

      }

    }

    // No corresponding grid2 vertex found for this grid1 vertex
    if (other[i] == -1)
      return;

  }

  // ////////////////////////////////////////////////////////////
  //   Set up the new remote intersection
  // ////////////////////////////////////////////////////////////

  const Dune::ReferenceElement<T,dim>& refElement = Dune::ReferenceElements<T,dim>::general(grid1ElementType);

  /** \todo Currently the RemoteIntersections have to be simplices */
  if (grid1ElementType.isSimplex()) {

    intersections.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

    for (int i=0; i<refElement.size(dim); i++) {
      intersections.back().grid1Local_[0][i] = refElement.position(i,dim);
      intersections.back().grid2Local_[0][i] = refElement.position(other[i],dim);
    }

  } else if (grid1ElementType.isQuadrilateral()) {

    // split the quadrilateral into two triangles
    const unsigned int subVertices[2][3] = {{0,1,3}, {0,3,2}};

    for (int i=0; i<2; i++) {

      RemoteSimplicialIntersection newSimplicialIntersection(grid1Index, grid2Index);

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.grid1Local_[0][j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.grid2Local_[0][j] = refElement.position(subVertices[i][other[j]],dim);
      }

      intersections.push_back(newSimplicialIntersection);

    }

  } else if (grid1ElementType.isHexahedron()) {

    // split the hexahedron into five tetrahedra
    // This can be removed if ever we allow RemoteIntersections that are not simplices
    const unsigned int subVertices[5][4] = {{0,1,3,5}, {0,3,2,6}, {4,5,0,6}, {6,7,6,3}, {6,0,5,3}};

    for (int i=0; i<5; i++) {

      RemoteSimplicialIntersection newSimplicialIntersection(grid1Index, grid2Index);

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.grid1Local_[0][j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.grid2Local_[0][j] = refElement.position(subVertices[i][other[j]],dim);
      }

      intersections.push_back(newSimplicialIntersection);

    }

  } else
    DUNE_THROW(Dune::GridError, "Unsupported element type");

}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::grid1Parent(unsigned int idx, unsigned int parId) const
{
    return this->intersections_[idx].grid1Entities_[parId];
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::grid2Parent(unsigned int idx, unsigned int parId) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after grid1 parent indices!!
    return this->intersections_[idx].grid2Entities_[parId];
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId) const
{
  return this->intersections_[idx].grid1Local_[parId][corner];
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId) const
{
  return this->intersections_[idx].grid2Local_[parId][corner];
}

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_MERGING_CONFORMINGMERGE_HH
