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
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief
 * Implementation of the Merger concept for conforming interfaces
 */

#ifndef CONFORMING_MERGE_HH
#define CONFORMING_MERGE_HH

#include <iomanip>
#include <vector>
#include <algorithm>
#include <bitset>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/grid-glue/merging/standardmerge.hh>


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

  typedef ::RemoteSimplicialIntersection<T,dim,dim,dimworld> RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  virtual void computeIntersection(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   unsigned int grid1Index,
                                   std::bitset<(1<<dim)>& neighborIntersects1,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   unsigned int grid2Index,
                                   std::bitset<(1<<dim)>& neighborIntersects2);

public:

  ConformingMerge(T tolerance = 1E-4) :
    tolerance_(tolerance)
  {}

private:

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

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


  /*   M A P P I N G   O N   I N D E X   B A S I S   */

  /**
   * @brief get index of grid1 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid1 parent simplex
   */
  unsigned int grid1Parent(unsigned int idx) const;

  /**
   * @brief get index of grid2 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid2 parent simplex
   */
  unsigned int grid2Parent(unsigned int idx) const;

  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1 simplex
   */
  LocalCoords grid1ParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2 simplex
   */
  LocalCoords grid2ParentLocal(unsigned int idx, unsigned int corner) const;

};


template<int dim, int dimworld, typename T>
void ConformingMerge<dim, dimworld, T>::computeIntersection(const Dune::GeometryType& grid1ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                                            unsigned int grid1Index,
                                                            std::bitset<(1<<dim)>& neighborIntersects1,
                                                            const Dune::GeometryType& grid2ElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                                            unsigned int grid2Index,
                                                            std::bitset<(1<<dim)>& neighborIntersects2)
{
  this->counter++;

  // A few consistency checks
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid1ElementType).size(dim)) == grid1ElementCorners.size());
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid2ElementType).size(dim)) == grid2ElementCorners.size());

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

  const Dune::GenericReferenceElement<T,dim>& refElement = Dune::GenericReferenceElements<T,dim>::general(grid1ElementType);

  /** \todo Currently the RemoteIntersections have to be simplices */
  if (grid1ElementType.isSimplex()) {

    this->intersections_.push_back(RemoteSimplicialIntersection());

    for (int i=0; i<refElement.size(dim); i++) {
      this->intersections_.back().grid1Local_[i] = refElement.position(i,dim);
      this->intersections_.back().grid2Local_[i] = refElement.position(other[i],dim);
    }

    this->intersections_.back().grid1Entity_ = grid1Index;
    this->intersections_.back().grid2Entity_ = grid2Index;

  } else if (grid1ElementType.isQuadrilateral()) {

    // split the quadrilateral into two triangles
    const unsigned int subVertices[2][3] = {{0,1,3}, {0,3,2}};

    for (int i=0; i<2; i++) {

      RemoteSimplicialIntersection newSimplicialIntersection;

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.grid1Local_[j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.grid2Local_[j] = refElement.position(subVertices[i][other[j]],dim);
      }

      newSimplicialIntersection.grid1Entity_ = grid1Index;
      newSimplicialIntersection.grid2Entity_ = grid2Index;

      this->intersections_.push_back(newSimplicialIntersection);

    }

  } else if (grid1ElementType.isHexahedron()) {

    // split the hexahedron into five tetrahedra
    // This can be removed if ever we allow RemoteIntersections that are not simplices
    const unsigned int subVertices[5][4] = {{0,1,3,5}, {0,3,2,6}, {4,5,0,6}, {6,7,6,3}, {6,0,5,3}};

    for (int i=0; i<5; i++) {

      RemoteSimplicialIntersection newSimplicialIntersection;

      for (int j=0; j<dim+1; j++) {
        newSimplicialIntersection.grid1Local_[j] = refElement.position(subVertices[i][j],dim);
        newSimplicialIntersection.grid2Local_[j] = refElement.position(subVertices[i][other[j]],dim);
      }

      newSimplicialIntersection.grid1Entity_ = grid1Index;
      newSimplicialIntersection.grid2Entity_ = grid2Index;

      this->intersections_.push_back(newSimplicialIntersection);

    }

  } else
    DUNE_THROW(Dune::GridError, "Unsupported element type");

}



template<int dim, int dimworld, typename T>
inline bool ConformingMerge<dim, dimworld, T>::grid1SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid1 entities
  return true;
}


template<int dim, int dimworld, typename T>
inline bool ConformingMerge<dim, dimworld, T>::grid2SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid2 entities
  return true;
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::grid1Parent(unsigned int idx) const
{
  return this->intersections_[idx].grid1Entity_;
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::grid2Parent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after grid1 parent indices!!
  return this->intersections_[idx].grid2Entity_;
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::grid1ParentLocal(unsigned int idx, unsigned int corner) const
{
  return this->intersections_[idx].grid1Local_[corner];
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::grid2ParentLocal(unsigned int idx, unsigned int corner) const
{
  return this->intersections_[idx].grid2Local_[corner];
}

#endif // CONFORMING_MERGE_HH
