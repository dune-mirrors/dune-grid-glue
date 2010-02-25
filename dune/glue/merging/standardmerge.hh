// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect
 */

#ifndef STANDARD_MERGE_HH
#define STANDARD_MERGE_HH


#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/glue/merging/merger.hh>



/** \brief Common base class for many merger implementations: produce pairs of entities that _may_ intersect

   Many merger algorithms consist of two parts: on the one hand there is a mechanism that produces pairs of
   elements that may intersect.  On the other hand there is an algorithm that computes the intersection of two
   given elements.  For the pairs-producing algorithm there appears to be a canonical choice, namely the algorithm
   by Gander and Japhet described in 'An Algorithm for Non-Matching Grid Projections with Linear Complexity,
   M.J. Gander and C. Japhet, Domain Decomposition Methods in Science and Engineering XVIII,
   pp. 185--192, Springer-Verlag, 2009.'  This class implements this algorithm, calling a pure virtual function
   computeIntersections() to compute the intersection between two elements.  Actual merger implementations
   can derive from this class and only implement computeIntersections().

   \tparam T The type used for coordinates (assumed to be the same for both grids)
   \tparam domainDim Dimension of the domain grid
   \tparam targetDim Dimension of the target grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template<class T, int domainDim, int targetDim, int dimworld>
class StandardMerge
  : public Merger<T,domainDim,targetDim,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief Type used for local coordinates on the domain side
  typedef typename Merger<T,domainDim,targetDim,dimworld>::DomainCoords DomainCoords;

  /// @brief Type used for local coordinates on the target side
  typedef typename Merger<T,domainDim,targetDim,dimworld>::TargetCoords TargetCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

protected:

  struct RemoteSimplicialIntersection
  {
    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,domainDim>, dimworld+1> domainLocal_;

    // Local coordinates in the domain entity
    Dune::array<Dune::FieldVector<T,targetDim>, dimworld+1> targetLocal_;

    //
    int domainEntity_;

    int targetEntity_;

  };

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  virtual void computeIntersection(const Dune::GeometryType& domainElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& domainElementCorners,
                                   unsigned int domainIndex,
                                   const Dune::GeometryType& targetElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& targetElementCorners,
                                   unsigned int targetIndex) = 0;


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
  void build(const std::vector<Dune::FieldVector<T,dimworld> >& domainCoords,
             const std::vector<unsigned int>& domain_elements,
             const std::vector<Dune::GeometryType>& domain_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& targetCoords,
             const std::vector<unsigned int>& target_elements,
             const std::vector<Dune::GeometryType>& target_element_types
             );


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const;

#if 0
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
#endif

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
  DomainCoords domainParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent target simplex
   */
  TargetCoords targetParentLocal(unsigned int idx, unsigned int corner) const;

};


/* IMPLEMENTATION */

template<typename T, int domainDim, int targetDim, int dimworld>
void StandardMerge<T,domainDim,targetDim,dimworld>::build(const std::vector<Dune::FieldVector<T,dimworld> >& domainCoords,
                                                          const std::vector<unsigned int>& domain_elements,
                                                          const std::vector<Dune::GeometryType>& domain_element_types,
                                                          const std::vector<Dune::FieldVector<T,dimworld> >& targetCoords,
                                                          const std::vector<unsigned int>& target_elements,
                                                          const std::vector<Dune::GeometryType>& target_element_types
                                                          )
{

  std::cout << "StandardMerge building merged grid..." << std::endl;

  // /////////////////////////////////////////////////////////////////////
  //   Compute the intersection of all pairs of elements
  //   \todo This is only the naive quadratic algorithm
  // /////////////////////////////////////////////////////////////////////

  unsigned int domainCornerCounter = 0;

  for (std::size_t i=0; i<domain_element_types.size(); i++) {

    // Select vertices of the domain element
    int domainNumVertices = Dune::GenericReferenceElements<T,domainDim>::general(domain_element_types[i]).size(domainDim);
    std::vector<Dune::FieldVector<T,dimworld> > domainElementCorners(domainNumVertices);
    for (int j=0; j<domainNumVertices; j++)
      domainElementCorners[j] = domainCoords[domain_elements[domainCornerCounter++]];

    unsigned int targetCornerCounter = 0;

    for (std::size_t j=0; j<target_element_types.size(); j++) {

      // Select vertices of the domain element
      int targetNumVertices = Dune::GenericReferenceElements<T,targetDim>::general(target_element_types[j]).size(targetDim);
      std::vector<Dune::FieldVector<T,dimworld> > targetElementCorners(targetNumVertices);
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


template<typename T, int domainDim, int targetDim, int dimworld>
inline unsigned int StandardMerge<T,domainDim,targetDim,dimworld>::nSimplices() const
{
  return intersections_.size();
}

template<typename T, int domainDim, int targetDim, int dimworld>
inline unsigned int StandardMerge<T,domainDim,targetDim,dimworld>::domainParent(unsigned int idx) const
{
  return intersections_[idx].domainEntity_;
}


template<typename T, int domainDim, int targetDim, int dimworld>
inline unsigned int StandardMerge<T,domainDim,targetDim,dimworld>::targetParent(unsigned int idx) const
{
  return intersections_[idx].targetEntity_;
}

template<typename T, int domainDim, int targetDim, int dimworld>
bool StandardMerge<T,domainDim,targetDim,dimworld>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  // WARNING: stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if (intersections_[i].domainEntity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<typename T, int domainDim, int targetDim, int dimworld>
bool StandardMerge<T,domainDim,targetDim,dimworld>::targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  // WARNING: stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if (intersections_[i].targetEntity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<typename T, int domainDim, int targetDim, int dimworld>
typename StandardMerge<T,domainDim,targetDim,dimworld>::DomainCoords StandardMerge<T,domainDim,targetDim,dimworld>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].domainLocal_[corner];
}


template<typename T, int domainDim, int targetDim, int dimworld>
typename StandardMerge<T,domainDim,targetDim,dimworld>::TargetCoords StandardMerge<T,domainDim,targetDim,dimworld>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  return intersections_[idx].targetLocal_[corner];
}


#endif // STANDARD_MERGE_HH
