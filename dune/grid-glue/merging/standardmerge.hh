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
#include <stack>
#include <set>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/merger.hh>



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

  /** \brief Test whether the two given elements intersect
      \note If they do intersect their intersection is automatically added to the list!
   */
  bool testIntersection(unsigned int candidate0, unsigned int candidate1,
                        const std::vector<Dune::FieldVector<T,dimworld> >& domainCoords,
                        const std::vector<Dune::GeometryType>& domain_element_types,
                        const std::vector<Dune::FieldVector<T,dimworld> >& targetCoords,
                        const std::vector<Dune::GeometryType>& target_element_types);


  /*   M E M B E R   V A R I A B L E S   */

  /** \brief The computed intersections */
  std::vector<RemoteSimplicialIntersection> intersections_;

  /** \brief Temporary internal data */
  std::vector<std::vector<unsigned int> > domainElementCorners_;
  std::vector<std::vector<unsigned int> > targetElementCorners_;
  std::vector<std::vector<int> > elementsPerVertex0_;
  std::vector<std::vector<int> > elementsPerVertex1_;


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
bool StandardMerge<T,domainDim,targetDim,dimworld>::testIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                     const std::vector<Dune::FieldVector<T,dimworld> >& domainCoords,
                                                                     const std::vector<Dune::GeometryType>& domain_element_types,
                                                                     const std::vector<Dune::FieldVector<T,dimworld> >& targetCoords,
                                                                     const std::vector<Dune::GeometryType>& target_element_types)
{
  unsigned int oldNumberOfIntersections = intersections_.size();

  // Select vertices of the domain element
  int domainNumVertices = domainElementCorners_[candidate0].size();
  std::vector<Dune::FieldVector<T,dimworld> > domainElementCorners(domainNumVertices);
  for (int i=0; i<domainNumVertices; i++)
    domainElementCorners[i] = domainCoords[domainElementCorners_[candidate0][i]];

  // Select vertices of the target element
  int targetNumVertices = targetElementCorners_[candidate1].size();
  std::vector<Dune::FieldVector<T,dimworld> > targetElementCorners(targetNumVertices);
  for (int i=0; i<targetNumVertices; i++)
    targetElementCorners[i] = targetCoords[targetElementCorners_[candidate1][i]];

  // ///////////////////////////////////////////////////////
  //   Compute the intersection between the two elements
  // ///////////////////////////////////////////////////////

  computeIntersection(domain_element_types[candidate0], domainElementCorners, candidate0,
                      target_element_types[candidate1], targetElementCorners, candidate1);

  // Have we found an intersection?
  return (intersections_.size() > oldNumberOfIntersections);

}


// /////////////////////////////////////////////////////////////////////
//   Compute the intersection of all pairs of elements
//   Linear algorithm by Gander and Japhet, Proc. of DD18
// /////////////////////////////////////////////////////////////////////

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
  //   Copy element corners into a data structure with block-structure.
  //   This is not as efficient but a lot easier to use.
  //   We may think about efficiency later.
  // /////////////////////////////////////////////////////////////////////

  // first the domain side
  domainElementCorners_.resize(domain_element_types.size());

  unsigned int domainCornerCounter = 0;

  for (std::size_t i=0; i<domain_element_types.size(); i++) {

    // Select vertices of the domain element
    int numVertices = Dune::GenericReferenceElements<T,domainDim>::general(domain_element_types[i]).size(domainDim);
    domainElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      domainElementCorners_[i][j] = domain_elements[domainCornerCounter++];

  }

  // then the target side
  targetElementCorners_.resize(target_element_types.size());

  unsigned int targetCornerCounter = 0;

  for (std::size_t i=0; i<target_element_types.size(); i++) {

    // Select vertices of the target element
    int numVertices = Dune::GenericReferenceElements<T,targetDim>::general(target_element_types[i]).size(targetDim);
    targetElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      targetElementCorners_[i][j] = target_elements[targetCornerCounter++];

  }

  // /////////////////////////////////////////////////////////////////////
  //   Set of lists of elements per vertex.  This is needed to find
  //   neighboring elements later.
  // /////////////////////////////////////////////////////////////////////

  // first the domain side
  elementsPerVertex0_.resize(domainCoords.size());

  for (std::size_t i=0; i<domainElementCorners_.size(); i++)
    for (int j=0; j<domainElementCorners_[i].size(); j++)
      elementsPerVertex0_[domainElementCorners_[i][j]].push_back(i);

  // then the target side
  elementsPerVertex1_.resize(targetCoords.size());

  for (std::size_t i=0; i<target_element_types.size(); i++)
    for (int j=0; j<targetElementCorners_[i].size(); j++)
      elementsPerVertex1_[targetElementCorners_[i][j]].push_back(i);

  std::stack<unsigned int> candidates0;
  std::stack<std::pair<unsigned int,unsigned int> > candidates1;

  // /////////////////////////////////////////////////////////////////////
  //   Do a brute-force search to find one pair of intersecting elements
  //   to start the advancing-front type algorithm with.
  // /////////////////////////////////////////////////////////////////////

  for (std::size_t i=0; i<domain_element_types.size(); i++) {

    for (std::size_t j=0; j<target_element_types.size(); j++) {

      if (testIntersection(i,j,domainCoords,domain_element_types,targetCoords,target_element_types)) {
        candidates1.push(std::make_pair(j,i));          // the candidate and a seed for the candidate
        break;
      }

    }

    // Have we found an intersection?
    if (intersections_.size() > 0)
      break;

  }

  // clear global intersection list again:  the main algorithm wants to start with
  // an empty list.
  intersections_.resize(0);

  // /////////////////////////////////////////////////////////////////////
  //   Main loop
  // /////////////////////////////////////////////////////////////////////

  Dune::BitSetVector<1> isHandled0(domain_element_types.size());
  Dune::BitSetVector<1> isHandled1(target_element_types.size());

  Dune::BitSetVector<1> isCandidate0(domain_element_types.size());
  Dune::BitSetVector<1> isCandidate1(target_element_types.size());

  while (!candidates1.empty()) {

    // Get the next element on the target side
    unsigned int currentCandidate1 = candidates1.top().first;
    unsigned int seed = candidates1.top().second;
    candidates1.pop();
    isHandled1[currentCandidate1] = true;

    // Start advancing front algorithm on the domain side from the 'seed' element that
    // we stored along with the current target element
    candidates0.push(seed);

    isHandled0.unsetAll();
    isCandidate0.unsetAll();


    while (!candidates0.empty()) {

      unsigned int currentCandidate0 = candidates0.top();
      candidates0.pop();
      isHandled0[currentCandidate0] = true;

      // Test whether there is an intersection between currentCandidate0 and currentCandidate1
      bool intersectionFound = testIntersection(currentCandidate0, currentCandidate1,
                                                domainCoords,domain_element_types,targetCoords,target_element_types);

      if (intersectionFound) {

        // get all neighbors of currentCandidate0, but not currentCandidate0 itself
        for (int i=0; i<domainElementCorners_[currentCandidate0].size(); i++) {
          unsigned int v = domainElementCorners_[currentCandidate0][i];

          // The neighbors of currentCandidate0 are all possible candidates
          for (typename std::vector<int>::iterator it = elementsPerVertex0_[v].begin();
               it != elementsPerVertex0_[v].end(); ++it) {

            if (!isHandled0[*it][0] && !isCandidate0[*it][0]) {
              candidates0.push(*it);
              isCandidate0[*it] = true;
            }

          }

        }

      }

    }

    // We have now found all intersections of elements in the domain side with currentCandidate1
    // Now we add all neighbors of currentCandidate1 that have not been treated yet as new
    // candidates.
    std::set<unsigned int> neighbors1;

    // get all neighbors of currentCandidate1, but not currentCandidate1 itself
    for (int i=0; i<targetElementCorners_[currentCandidate1].size(); i++) {
      unsigned int v = domainElementCorners_[currentCandidate1][i];
      neighbors1.insert(elementsPerVertex1_[v].begin(), elementsPerVertex1_[v].end());
    }

    // The unhandled neighbors of currentCandidate1 are all possible candidates
    for (typename std::set<unsigned int>::iterator it = neighbors1.begin();
         it != neighbors1.end(); ++it) {

      if (!isHandled1[*it][0] && !isCandidate1[*it][0]) {

        // Get a seed element for the new target element
        // Look for an element on the domain side that intersects the new target element.
        // Look only among the ones that have been tested during the last iteration.
        // Since currentCandidate1 is a neighbor of the previous element, there has to be one.
        int seed = -1;
        for (int i=0; i<isHandled0.size(); i++) {

          if (!isHandled0[i][0])
            continue;

          int oldSize = intersections_.size();
          bool intersectionFound = testIntersection(i, *it,
                                                    domainCoords,domain_element_types,
                                                    targetCoords,target_element_types);

          if (intersectionFound) {

            // i is our new seed candidate on the domain side
            seed = i;

            while (intersections_.size()> oldSize)
              intersections_.pop_back();

            break;

          }

        }

        if (seed < 0) {
          // The fast method didn't find a domain element that intersections with
          // the new target candidate.  We have to do a brute-force search.
          for (std::size_t i=0; i<domain_element_types.size(); i++) {

            int oldSize = intersections_.size();
            bool intersectionFound = testIntersection(i, *it,
                                                      domainCoords,domain_element_types,
                                                      targetCoords,target_element_types);

            if (intersectionFound) {

              // i is our new seed candidate on the domain side
              seed = i;

              while (intersections_.size()> oldSize)
                intersections_.pop_back();

              break;

            }

          }

        }

        // We have tried all we could: the candidate is 'handled' now
        isCandidate1[*it] = true;

        if (seed < 0)
          // still no seed?  Then the new target candidate isn't overlapped by anything
          continue;

        // we have a seed now
        candidates1.push(std::make_pair(*it,seed));

      }

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
