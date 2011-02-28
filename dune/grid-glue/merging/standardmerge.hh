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
   computeIntersection() to compute the intersection between two elements.  Actual merger implementations
   can derive from this class and only implement computeIntersection().

   \tparam T The type used for coordinates (assumed to be the same for both grids)
   \tparam grid1Dim Dimension of the grid1 grid
   \tparam grid2Dim Dimension of the grid2 grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template<class T, int grid1Dim, int grid2Dim, int dimworld>
class StandardMerge
  : public Merger<T,grid1Dim,grid2Dim,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief Type used for local coordinates on the grid1 side
  typedef typename Merger<T,grid1Dim,grid2Dim,dimworld>::Grid1Coords Grid1Coords;

  /// @brief Type used for local coordinates on the grid2 side
  typedef typename Merger<T,grid1Dim,grid2Dim,dimworld>::Grid2Coords Grid2Coords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

protected:

  bool valid;

  StandardMerge() : valid(false) {}

  struct RemoteSimplicialIntersection
  {
    // Local coordinates in the grid1 entity
    Dune::array<Dune::FieldVector<T,grid1Dim>, dimworld+1> grid1Local_;

    // Local coordinates in the grid1 entity
    Dune::array<Dune::FieldVector<T,grid2Dim>, dimworld+1> grid2Local_;

    //
    int grid1Entity_;

    int grid2Entity_;

  };

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.
   */
  virtual void computeIntersection(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   unsigned int grid2Index) = 0;

  /** \brief Test whether the two given elements intersect
      \note If they do intersect their intersection is automatically added to the list!
   */
  bool testIntersection(unsigned int candidate0, unsigned int candidate1,
                        const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                        const std::vector<Dune::GeometryType>& grid1_element_types,
                        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                        const std::vector<Dune::GeometryType>& grid2_element_types);


  /*   M E M B E R   V A R I A B L E S   */

  /** \brief The computed intersections */
  std::vector<RemoteSimplicialIntersection> intersections_;

  /** \brief Temporary internal data */
  std::vector<std::vector<unsigned int> > grid1ElementCorners_;
  std::vector<std::vector<unsigned int> > grid2ElementCorners_;
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
   * @param grid1_coords the grid1 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid1_simplices array with all grid1 simplices represented as corner indices into @c grid1_coords;
   * the simplices are just written to this array one after another
   * @param grid2_coords the grid2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid2_simplices just like with the grid1_simplices and grid1_coords
   */
  void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
             const std::vector<unsigned int>& grid1_elements,
             const std::vector<Dune::GeometryType>& grid1_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
             const std::vector<unsigned int>& grid2_elements,
             const std::vector<Dune::GeometryType>& grid2_element_types
             );


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const;

  void clear()
  {
    // Delete old internal data, from a possible previous run
    purge(intersections_);
    purge(grid1ElementCorners_);
    purge(grid2ElementCorners_);
    purge(elementsPerVertex0_);
    purge(elementsPerVertex1_);

    valid = false;
  }

private:

  /** clear arbitrary containers */
  template<typename V>
  static void purge(V & v)
  {
    v.clear();
    V v2(v);
    v.swap(v2);
  }

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

  /**
   * @brief get the merged grid simplices refining a given grid1 simplex
   * @param idx index of grid1 simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool grid1SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const;

  /**
   * @brief get the merged grid simplices refining a given grid2 simplex
   * @param idx index of grid2 simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool grid2SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const;


  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1 simplex
   */
  Grid1Coords grid1ParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2 simplex
   */
  Grid2Coords grid2ParentLocal(unsigned int idx, unsigned int corner) const;

};


/* IMPLEMENTATION */

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::testIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                   const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                   const std::vector<Dune::GeometryType>& grid2_element_types)
{
  unsigned int oldNumberOfIntersections = intersections_.size();

  // Select vertices of the grid1 element
  int grid1NumVertices = grid1ElementCorners_[candidate0].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid1ElementCorners(grid1NumVertices);
  for (int i=0; i<grid1NumVertices; i++)
    grid1ElementCorners[i] = grid1Coords[grid1ElementCorners_[candidate0][i]];

  // Select vertices of the grid2 element
  int grid2NumVertices = grid2ElementCorners_[candidate1].size();
  std::vector<Dune::FieldVector<T,dimworld> > grid2ElementCorners(grid2NumVertices);
  for (int i=0; i<grid2NumVertices; i++)
    grid2ElementCorners[i] = grid2Coords[grid2ElementCorners_[candidate1][i]];

  // ///////////////////////////////////////////////////////
  //   Compute the intersection between the two elements
  // ///////////////////////////////////////////////////////

  computeIntersection(grid1_element_types[candidate0], grid1ElementCorners, candidate0,
                      grid2_element_types[candidate1], grid2ElementCorners, candidate1);

  // Have we found an intersection?
  return (intersections_.size() > oldNumberOfIntersections);

}


// /////////////////////////////////////////////////////////////////////
//   Compute the intersection of all pairs of elements
//   Linear algorithm by Gander and Japhet, Proc. of DD18
// /////////////////////////////////////////////////////////////////////

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                        const std::vector<unsigned int>& grid1_elements,
                                                        const std::vector<Dune::GeometryType>& grid1_element_types,
                                                        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                        const std::vector<unsigned int>& grid2_elements,
                                                        const std::vector<Dune::GeometryType>& grid2_element_types
                                                        )
{

  std::cout << "StandardMerge building merged grid..." << std::endl;

  clear();

  // /////////////////////////////////////////////////////////////////////
  //   Copy element corners into a data structure with block-structure.
  //   This is not as efficient but a lot easier to use.
  //   We may think about efficiency later.
  // /////////////////////////////////////////////////////////////////////

  // first the grid1 side
  grid1ElementCorners_.resize(grid1_element_types.size());

  unsigned int grid1CornerCounter = 0;

  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    // Select vertices of the grid1 element
    int numVertices = Dune::GenericReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
    grid1ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid1ElementCorners_[i][j] = grid1_elements[grid1CornerCounter++];

  }

  // then the grid2 side
  grid2ElementCorners_.resize(grid2_element_types.size());

  unsigned int grid2CornerCounter = 0;

  for (std::size_t i=0; i<grid2_element_types.size(); i++) {

    // Select vertices of the grid2 element
    int numVertices = Dune::GenericReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
    grid2ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid2ElementCorners_[i][j] = grid2_elements[grid2CornerCounter++];

  }

  // /////////////////////////////////////////////////////////////////////
  //   Set of lists of elements per vertex.  This is needed to find
  //   neighboring elements later.
  // /////////////////////////////////////////////////////////////////////

  // first the grid1 side
  elementsPerVertex0_.resize(grid1Coords.size());

  for (std::size_t i=0; i<grid1ElementCorners_.size(); i++)
    for (std::size_t j=0; j<grid1ElementCorners_[i].size(); j++)
      elementsPerVertex0_[grid1ElementCorners_[i][j]].push_back(i);

  // then the grid2 side
  elementsPerVertex1_.resize(grid2Coords.size());

  for (std::size_t i=0; i<grid2_element_types.size(); i++)
    for (std::size_t j=0; j<grid2ElementCorners_[i].size(); j++)
      elementsPerVertex1_[grid2ElementCorners_[i][j]].push_back(i);

  std::stack<unsigned int> candidates0;
  std::stack<std::pair<unsigned int,unsigned int> > candidates1;

  // /////////////////////////////////////////////////////////////////////
  //   Do a brute-force search to find one pair of intersecting elements
  //   to start the advancing-front type algorithm with.
  // /////////////////////////////////////////////////////////////////////

  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    for (std::size_t j=0; j<grid2_element_types.size(); j++) {

      if (testIntersection(i,j,grid1Coords,grid1_element_types,grid2Coords,grid2_element_types)) {
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
  intersections_.clear();

  // /////////////////////////////////////////////////////////////////////
  //   Main loop
  // /////////////////////////////////////////////////////////////////////

  Dune::BitSetVector<1> isHandled0(grid1_element_types.size());
  Dune::BitSetVector<1> isHandled1(grid2_element_types.size());

  Dune::BitSetVector<1> isCandidate0(grid1_element_types.size());
  Dune::BitSetVector<1> isCandidate1(grid2_element_types.size());

  while (!candidates1.empty()) {

    // Get the next element on the grid2 side
    unsigned int currentCandidate1 = candidates1.top().first;
    unsigned int seed = candidates1.top().second;
    candidates1.pop();
    isHandled1[currentCandidate1] = true;

    // Start advancing front algorithm on the grid1 side from the 'seed' element that
    // we stored along with the current grid2 element
    candidates0.push(seed);

    isHandled0.unsetAll();
    isCandidate0.unsetAll();

    std::set<unsigned int> potentialSeeds;

    while (!candidates0.empty()) {

      unsigned int currentCandidate0 = candidates0.top();
      candidates0.pop();
      isHandled0[currentCandidate0] = true;
      potentialSeeds.insert(currentCandidate0);

      // Test whether there is an intersection between currentCandidate0 and currentCandidate1
      bool intersectionFound = testIntersection(currentCandidate0, currentCandidate1,
                                                grid1Coords,grid1_element_types,grid2Coords,grid2_element_types);

      if (intersectionFound) {

        // get all neighbors of currentCandidate0, but not currentCandidate0 itself
        for (std::size_t i=0; i<grid1ElementCorners_[currentCandidate0].size(); i++) {
          unsigned int v = grid1ElementCorners_[currentCandidate0][i];

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

    // We have now found all intersections of elements in the grid1 side with currentCandidate1
    // Now we add all neighbors of currentCandidate1 that have not been treated yet as new
    // candidates.

    // get all neighbors of currentCandidate1, but not currentCandidate1 itself
    for (std::size_t i=0; i<grid2ElementCorners_[currentCandidate1].size(); i++) {
      unsigned int v = grid2ElementCorners_[currentCandidate1][i];

      // The unhandled neighbors of currentCandidate1 are all possible candidates
      for (typename std::vector<int>::iterator it = elementsPerVertex1_[v].begin();
           it != elementsPerVertex1_[v].end(); ++it) {

        if (!isHandled1[*it][0] && !isCandidate1[*it][0]) {

          // Get a seed element for the new grid2 element
          // Look for an element on the grid1 side that intersects the new grid2 element.
          // Look only among the ones that have been tested during the last iteration.
          int seed = -1;

          for (typename std::set<unsigned int>::iterator seedIt = potentialSeeds.begin();
               seedIt != potentialSeeds.end(); ++seedIt) {

            std::size_t oldSize = intersections_.size();
            bool intersectionFound = testIntersection(*seedIt, *it,
                                                      grid1Coords,grid1_element_types,
                                                      grid2Coords,grid2_element_types);

            if (intersectionFound) {

              // i is our new seed candidate on the grid1 side
              seed = *seedIt;

              while (intersections_.size() > oldSize)
                intersections_.pop_back();

              break;

            }

          }

          if (seed < 0) {
            // The fast method didn't find a grid1 element that intersects with
            // the new grid2 candidate.  We have to do a brute-force search.
            for (std::size_t i=0; i<grid1_element_types.size(); i++) {

              std::size_t oldSize = intersections_.size();
              bool intersectionFound = testIntersection(i, *it,
                                                        grid1Coords,grid1_element_types,
                                                        grid2Coords,grid2_element_types);

              if (intersectionFound) {

                // i is our new seed candidate on the grid1 side
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
            // still no seed?  Then the new grid2 candidate isn't overlapped by anything
            continue;

          // we have a seed now
          candidates1.push(std::make_pair(*it,seed));

        }

      }

    }

  }

  valid = true;

}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::nSimplices() const
{
  assert(valid);
  return intersections_.size();
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1Parent(unsigned int idx) const
{
  assert(valid);
  return intersections_[idx].grid1Entity_;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2Parent(unsigned int idx) const
{
  assert(valid);
  return intersections_[idx].grid2Entity_;
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  assert(valid);

#warning stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if ((unsigned int) intersections_[i].grid1Entity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  assert(valid);

#warning stupid linear algorithm!
  indices.resize(0);

  for (size_t i=0; i<intersections_.size(); i++)
    if ((unsigned int) intersections_[i].grid2Entity_ == idx)
      indices.push_back(i);

  return indices.size() > 0;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid1Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1ParentLocal(unsigned int idx, unsigned int corner) const
{
  assert(valid);
  return intersections_[idx].grid1Local_[corner];
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid2Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2ParentLocal(unsigned int idx, unsigned int corner) const
{
  assert(valid);
  return intersections_[idx].grid2Local_[corner];
}

#define DECL extern
#define STANDARD_MERGE_INSTANTIATE(T,A,B,C) \
  DECL template \
  void StandardMerge<T,A,B,C>::build(const std::vector<Dune::FieldVector<T,C> >& grid1Coords, \
                                     const std::vector<unsigned int>& grid1_elements, \
                                     const std::vector<Dune::GeometryType>& grid1_element_types, \
                                     const std::vector<Dune::FieldVector<T,C> >& grid2Coords, \
                                     const std::vector<unsigned int>& grid2_elements, \
                                     const std::vector<Dune::GeometryType>& grid2_element_types \
                                     )

STANDARD_MERGE_INSTANTIATE(double,1,1,1);
STANDARD_MERGE_INSTANTIATE(double,2,2,2);
STANDARD_MERGE_INSTANTIATE(double,3,3,3);
#undef STANDARD_MERGE_INSTANTIATE
#undef DECL

#endif // STANDARD_MERGE_HH
