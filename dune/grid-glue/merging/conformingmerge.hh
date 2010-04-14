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

  typedef typename StandardMerge<T,dim,dim,dimworld>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  virtual void computeIntersection(const Dune::GeometryType& domainElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& domainElementCorners,
                                   unsigned int domainIndex,
                                   const Dune::GeometryType& targetElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& targetElementCorners,
                                   unsigned int targetIndex);

public:

  ConformingMerge(T tolerance = 1E-4) :
    tolerance_(tolerance)
  {}


  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

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


template<int dim, int dimworld, typename T>
void ConformingMerge<dim, dimworld, T>::computeIntersection(const Dune::GeometryType& domainElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& domainElementCorners,
                                                            unsigned int domainIndex,
                                                            const Dune::GeometryType& targetElementType,
                                                            const std::vector<Dune::FieldVector<T,dimworld> >& targetElementCorners,
                                                            unsigned int targetIndex)
{
  // A few consistency checks
  assert((Dune::GenericReferenceElements<T,dim>::general(domainElementType).size(dim) == domainElementCorners.size()));
  assert((Dune::GenericReferenceElements<T,dim>::general(targetElementType).size(dim) == targetElementCorners.size()));

  // the intersection is either conforming or empty, hence the GeometryTypes have to match
  if (domainElementType != targetElementType)
    return;

  // ////////////////////////////////////////////////////////////
  //   Find correspondences between the different corners
  // ////////////////////////////////////////////////////////////
  std::vector<int> other(domainElementCorners.size(), -1);

  for (int i=0; i<domainElementCorners.size(); i++) {

    for (int j=0; j<targetElementCorners.size(); j++) {

      if ( (domainElementCorners[i]-targetElementCorners[j]).two_norm() < tolerance_ ) {

        other[i] = j;
        break;

      }

    }

    // No corresponding target vertex found for this domain vertex
    if (other[i] == -1)
      return;

  }

  // ////////////////////////////////////////////////////////////
  //   Set up the new remote intersection
  // ////////////////////////////////////////////////////////////

  /** \todo Currently the RemoteIntersections have to be simplices */
  if (domainElementType.isSimplex()) {

    this->intersections_.push_back(RemoteSimplicialIntersection());

    const Dune::GenericReferenceElement<T,dim>& refElement = Dune::GenericReferenceElements<T,dim>::general(domainElementType);

    for (int i=0; i<refElement.size(dim); i++) {
      this->intersections_.back().domainLocal_[i] = refElement.position(i,dim);
      this->intersections_.back().targetLocal_[i] = refElement.position(other[i],dim);
    }

    this->intersections_.back().domainEntity_ = domainIndex;
    this->intersections_.back().targetEntity_ = targetIndex;

  } else if (domainElementType.isQuadrilateral()) {

    // split the quadrilateral into two triangles
    const Dune::GenericReferenceElement<T,dim>& refElement = Dune::GenericReferenceElements<T,dim>::general(domainElementType);

    // triangle 1
    RemoteSimplicialIntersection newSimplicialIntersection1;

    for (int i=0; i<3; i++) {
      newSimplicialIntersection1.domainLocal_[i] = refElement.position(i,dim);
      newSimplicialIntersection1.targetLocal_[i] = refElement.position(other[i],dim);
    }

    newSimplicialIntersection1.domainEntity_ = domainIndex;
    newSimplicialIntersection1.targetEntity_ = targetIndex;

    this->intersections_.push_back(newSimplicialIntersection1);

    // triangle 2
    RemoteSimplicialIntersection newSimplicialIntersection2;

    newSimplicialIntersection2.domainLocal_[0] = refElement.position(2,dim);
    newSimplicialIntersection2.targetLocal_[0] = refElement.position(other[2],dim);
    newSimplicialIntersection2.domainLocal_[1] = refElement.position(1,dim);
    newSimplicialIntersection2.targetLocal_[1] = refElement.position(other[1],dim);
    newSimplicialIntersection2.domainLocal_[2] = refElement.position(3,dim);
    newSimplicialIntersection2.targetLocal_[2] = refElement.position(other[3],dim);

    newSimplicialIntersection2.domainEntity_ = domainIndex;
    newSimplicialIntersection2.targetEntity_ = targetIndex;

    this->intersections_.push_back(newSimplicialIntersection2);

  } else
    DUNE_THROW(Dune::GridError, "Unsupported element type");

}



template<int dim, int dimworld, typename T>
inline bool ConformingMerge<dim, dimworld, T>::domainSimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all domain entities
  return true;
}


template<int dim, int dimworld, typename T>
inline bool ConformingMerge<dim, dimworld, T>::targetSimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all target entities
  return true;
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::domainParent(unsigned int idx) const
{
  return this->intersections_[idx].domainEntity_;
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return this->intersections_[idx].targetEntity_;
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  return this->intersections_[idx].domainLocal_[corner];
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  return this->intersections_[idx].targetLocal_[corner];
}

#endif // CONFORMING_MERGE_HH
