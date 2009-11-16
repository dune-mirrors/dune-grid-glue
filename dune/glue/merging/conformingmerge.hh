// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    conformingmerge.hh
 *  Version:     1.0
 *  Created on:  Sep 14, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: implementation of the SurfaceMerge concept for conforming interfaces
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief
 * Standard implementation of the SurfaceMerge concept for the use in 2d and 3d.
 * Uses psurface routines to compute the merged grid  and provides access to it
 * via the interface specified by the concept.
 *
 * While this adapter is correctly implemented and thoroughly tested there are still some problems with the
 * implementation behind it. E.g. in a case where edges are mapped onto other edges the code is not stable
 * and the program might crash. Specifying a normal field in the 3d case does not guarantee the expected
 * improvements in the mapping either.
 */

#ifndef CONFORMING_MERGE_HH
#define CONFORMING_MERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <set>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/glue/misc/geometry.hh>
#include <dune/glue/merging/merger.hh>


/** \brief Implementation of the SurfaceMerge concept for conforming interfaces

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class ConformingMerge
  : public Merger<T,dim,dim,dimworld>
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

  struct ConformingRemoteIntersection
  {
    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,dim>, dim+1> domainLocal_;

    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,dim>, dim+1> targetLocal_;

    //
    int domainEntity_;

    int targetEntity_;

  };


private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  T tolerance_;

  /* geometric data for both domain and targt */

  /// @brief domain coordinates
  std::vector<Dune::FieldVector<double,dimworld> >   domainCoords_;

  /// @brief target coordinates
  std::vector<Dune::FieldVector<double,dimworld> >   targetCoords_;


  /* topologic information for domain and target */

  /// @ brief domain indices (internal copy)
  std::vector<Dune::array<int,dim+1> >         _domi;

  /// @brief target indices (internal copy)
  std::vector<Dune::array<int,dim+1> >         _tari;

  /** \brief The computed intersections */
  std::vector<ConformingRemoteIntersection> intersections_;

public:

  ConformingMerge(T tolerance = 1E-4) :
    tolerance_(tolerance)
  {}


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
  void build(const std::vector<T>& domain_coords,
             const std::vector<unsigned int>& domain_elements,
             const std::vector<Dune::GeometryType>& domain_element_types,
             const std::vector<T>& target_coords,
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


/* IMPLEMENTATION OF CLASS   C O N T A C T  M A P P I N G  S U R F A C E  M E R G E */


template<int dim, int dimworld, typename T>
void ConformingMerge<dim, dimworld, T>::build(
  const std::vector<T>& domain_coords,
  const std::vector<unsigned int>& domain_elements,
  const std::vector<Dune::GeometryType>& domain_element_types,
  const std::vector<T>& target_coords,
  const std::vector<unsigned int>& target_elements,
  const std::vector<Dune::GeometryType>& target_element_types
  )
{
  for (size_t i=0; i<domain_element_types.size(); i++)
    if (!domain_element_types[i].isSimplex())
      DUNE_THROW(Dune::GridError, "ConformingMerge is only implemented for simplicial elements!");

  for (size_t i=0; i<target_element_types.size(); i++)
    if (!target_element_types[i].isSimplex())
      DUNE_THROW(Dune::GridError, "ConformingMerge is only implemented for simplicial elements!");

  // copy domain and target elements to internal arrays
  // (cannot keep refs since outside modification would destroy information)
  this->_domi.resize(domain_elements.size()/(dim+1));
  this->_tari.resize(target_elements.size()/(dim+1));

  for (unsigned int i = 0; i < domain_elements.size()/(dim+1); ++i)
    for (int j=0; j<dim+1; j++)
      _domi[i][j] = domain_elements[i*(dim+1)+j];

  // dim==dimworld-1: just copy the two arrays
  for (unsigned int i = 0; i < target_elements.size()/(dim+1); ++i)
    for (int j=0; j<dim+1; j++)
      _tari[i][j] = target_elements[i*(dim+1)+j];

  // copy the coordinates to internal arrays of coordinates
  // (again cannot just keep refs and this representation has advantages)
  domainCoords_.resize(domain_coords.size() / dimworld);
  for (unsigned int i = 0; i < domainCoords_.size(); ++i)
    for (unsigned int j = 0; j < dimworld; ++j)
      domainCoords_[i][j] = domain_coords[i*dimworld + j];

  targetCoords_.resize(target_coords.size() / dimworld);
  for (unsigned int i = 0; i < targetCoords_.size(); ++i)
    for (unsigned int j = 0; j < dimworld; ++j)
      targetCoords_[i][j] = target_coords[i*dimworld + j];

  std::cout << "ConformingMerge building merged grid..." << std::endl;

  // /////////////////////////////////////////////////////////////////////
  //   Compute correspondences between vertices
  //   \todo This is only the naive quadratic algorithm
  // /////////////////////////////////////////////////////////////////////

  std::vector<int> pairs(domainCoords_.size(), -1);

  // Loop over first boundary
  for (int i=0; i<domainCoords_.size(); i++) {

    // Loop over second boundary
    for (int j=0; j<targetCoords_.size(); j++) {

      if ((domainCoords_[i]-targetCoords_[j]).two_norm() < tolerance_) {

        pairs[i] = j;
        break;

      }

    }

  }

  std::cout << "ConformingMerge: " << pairs.size() - std::count(pairs.begin(), pairs.end(), -1) << " node pairs found" << std::endl;

  // //////////////////////////////////////////////////////////////////////
  //   Make the actual intersections
  // //////////////////////////////////////////////////////////////////////

  intersections_.resize(0);

  for (size_t i=0; i<_domi.size(); i++) {

    // --------------------------------------------
    //  Find the corresponding target simplex
    // --------------------------------------------

    // Assemble a set of the vertices that the target simplex has to consist of
    std::set<unsigned int> domainSimplexSorted;
    for (int j=0; j<_domi[i].size(); j++)
      if (pairs[_domi[i][j]] != -1)
        domainSimplexSorted.insert(pairs[_domi[i][j]]);

    // Do nothing if not all of the vertices of this simplex could be matched on the other side
    if (domainSimplexSorted.size() != _domi[i].size())
      continue;

    // Find the index of this target simplex
    int targetSimplexIdx = -1;
    for (size_t j=0; j<_tari.size(); j++) {

      std::set<unsigned int> targetSimplexSorted(_tari[j].begin(), _tari[j].end());

      if (domainSimplexSorted == targetSimplexSorted) {
        targetSimplexIdx = j;
        break;
      }

    }

    // We have found a conforming intersection.  Let's store it!
    if (targetSimplexIdx != -1) {

      intersections_.push_back(ConformingRemoteIntersection());
      intersections_.back().domainEntity_ = i;
      intersections_.back().targetEntity_ = targetSimplexIdx;

      const Dune::GenericReferenceElement<T,dim>& refElement
        = Dune::GenericReferenceElements<T,dim>::general(Dune::GeometryType(Dune::GeometryType::simplex,dim));

      // Loop over the vertices of the intersection
      for (int j=0; j<refElement.size(dim); j++) {

        // local coordinates in the domain
        // \todo This is so trivial we may not even want to save it...
        intersections_.back().domainLocal_[j] = refElement.position(j,dim);

        // global vertex number on the target side
        int globalDomainNumber = _domi[i][j];
        int globalTargetNumber = pairs[globalDomainNumber];

        // local number in the target simplex
        assert(std::find(_tari[targetSimplexIdx].begin(),
                         _tari[targetSimplexIdx].end(), globalTargetNumber) != _tari[targetSimplexIdx].end());
        int localTargetNumber = std::find(_tari[targetSimplexIdx].begin(),
                                          _tari[targetSimplexIdx].end(), globalTargetNumber)
                                - _tari[targetSimplexIdx].begin();

        intersections_.back().targetLocal_[localTargetNumber] = refElement.position(j,dim);

      }

    }

  }

}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::nSimplices() const
{
  return intersections_.size();
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
  return intersections_[idx].domainEntity_;
}


template<int dim, int dimworld, typename T>
inline unsigned int ConformingMerge<dim, dimworld, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return intersections_[idx].targetEntity_;
}


template<int dim, int dimworld, typename T>
bool ConformingMerge<dim, dimworld, T>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  indices.resize(1);
  indices[0] = idx;

  // ...return success
  return true;
}


template<int dim, int dimworld, typename T>
bool ConformingMerge<dim, dimworld, T>::targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  indices.resize(1);
  indices[0] = intersections_[idx].targetEntity_;

  // ...return success
  return true;
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].domainLocal_[corner];
}


template<int dim, int dimworld, typename T>
typename ConformingMerge<dim, dimworld, T>::LocalCoords ConformingMerge<dim, dimworld, T>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].targetLocal_[corner];
}

#endif // CONFORMING_MERGE_HH
