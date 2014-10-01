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
#include <map>
#include <algorithm>

#include <dune/common/fvector.hh>
#include <dune/common/bitsetvector.hh>
#include <dune/common/stdstreams.hh>
#include <dune/common/timer.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/merger.hh>
#include <dune/grid-glue/merging/computeintersection.hh>


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
    /** \brief Dimension of this intersection */
    enum {intersectionDim = (grid1Dim<grid2Dim) ? grid1Dim : grid2Dim};

    /** \brief Number of vertices of the intersection (it's a simplex) */
    enum {nVertices = intersectionDim + 1};

    /** \brief Default constructor */
    RemoteSimplicialIntersection()
    {
        grid1Entities_.resize(1);
        grid2Entities_.resize(1);
        grid1Local_.resize(1);
        grid2Local_.resize(1);
    }

    /** \brief Constructor for two given entity indices */
    RemoteSimplicialIntersection(int grid1Entity, int grid2Entity)
    {
        grid1Entities_.resize(1);
        grid2Entities_.resize(1);
        grid1Local_.resize(1);
        grid2Local_.resize(1);

        grid1Entities_[0] = grid1Entity;
        grid2Entities_[0] = grid2Entity;
    }

    // Local coordinates in the grid1 entity
    std::vector<Dune::array<Dune::FieldVector<T,grid1Dim>, nVertices> > grid1Local_;

    // Local coordinates in the grid1 entity
    std::vector<Dune::array<Dune::FieldVector<T,grid2Dim>, nVertices> > grid2Local_;

    //
    std::vector<int> grid1Entities_;

    std::vector<int> grid2Entities_;

  };

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices stored in the vector intersections.
   */
  virtual void computeIntersections(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                   unsigned int grid1Index,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                                   unsigned int grid2Index,
                                   std::vector<RemoteSimplicialIntersection>& intersections) = 0;

  /** \brief Compute the intersection between two overlapping elements
   * \return true if at least one intersection point was found
   */
  bool computeIntersection(unsigned int candidate0, unsigned int candidate1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                           const std::vector<Dune::GeometryType>& grid1_element_types,
                           std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                           const std::vector<Dune::GeometryType>& grid2_element_types,
                           std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                           bool insert = true);

  /*   M E M B E R   V A R I A B L E S   */

  /** \brief The computed intersections */
  std::vector<RemoteSimplicialIntersection> intersections_;

  /** \brief Temporary internal data */
  std::vector<std::vector<unsigned int> > grid1ElementCorners_;
  std::vector<std::vector<unsigned int> > grid2ElementCorners_;

  std::vector<std::vector<int> > elementNeighbors1_;
  std::vector<std::vector<int> > elementNeighbors2_;

public:

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @copydoc Merger<T,grid1Dim,grid2Dim,dimworld>::build
   */
  virtual void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1_Coords,
             const std::vector<unsigned int>& grid1_elements,
             const std::vector<Dune::GeometryType>& grid1_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& grid2_coords,
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
   * @brief get number of grid1 parents to the intersection idx
   * @param idx index of the merged grid simplex
   * @return amount of parent simplices
   */
  unsigned int grid1Parents(unsigned int idx) const;

  /**
   * @brief get number of grid2 parents to the intersection idx
   * @param idx index of the merged grid simplex
   * @return amount of parent simplices
   */
  unsigned int grid2Parents(unsigned int idx) const;

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
  Grid1Coords grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2 simplex
   */
  Grid2Coords grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

  /**
   * Do a brute-force search to find one pair of intersecting elements
   * to start or continue the advancing-front type algorithm.
   */
  void generateSeed(std::vector<int>& seeds,
                    Dune::BitSetVector<1>& isHandled2,
                    std::stack<unsigned>& candidates2,
                    const std::vector<Dune::FieldVector<T, dimworld> >& grid1Coords,
                    const std::vector<Dune::GeometryType>& grid1_element_types,
                    const std::vector<Dune::FieldVector<T, dimworld> >& grid2Coords,
                    const std::vector<Dune::GeometryType>& grid2_element_types);

  /**
    * Insert intersections into this->intersection_ and return index
    */
  int insertIntersections(unsigned int candidate1, unsigned int candidate2,std::vector<RemoteSimplicialIntersection>& intersections);

  /**
   * Find a grid2 element intersecting the candidate1 grid1 element by brute force search
   */
  int bruteForceSearch(int candidate1,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                       const std::vector<Dune::GeometryType>& grid1_element_types,
                       const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                       const std::vector<Dune::GeometryType>& grid2_element_types);

  /**
   * Get the index of the intersection in intersections_ (= size if it is a new intersection)
   */
  int intersectionIndex(unsigned int grid1Index, unsigned int grid2Index,
                                 RemoteSimplicialIntersection& intersection);

  /**
   * get the neighbor relations between the given elements
   */
  template <int gridDim>
  void computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                                  const std::vector<std::vector<unsigned int> >& gridElementCorners,
                                  std::vector<std::vector<int> >& elementNeighbors);
};


/* IMPLEMENTATION */

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
bool StandardMerge<T,grid1Dim,grid2Dim,dimworld>::computeIntersection(unsigned int candidate0, unsigned int candidate1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                      const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                      std::bitset<(1<<grid1Dim)>& neighborIntersects1,
                                                                      const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                      const std::vector<Dune::GeometryType>& grid2_element_types,
                                                                      std::bitset<(1<<grid2Dim)>& neighborIntersects2,
                                                                      bool insert)
{
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

  std::vector<RemoteSimplicialIntersection> intersections(0);

  // compute the intersections
  computeIntersections(grid1_element_types[candidate0], grid1ElementCorners,
                       neighborIntersects1, candidate0,
                       grid2_element_types[candidate1], grid2ElementCorners,
                       neighborIntersects2, candidate1,
                       intersections);

  // insert intersections if needed
  if(insert && intersections.size() > 0)
      insertIntersections(candidate0,candidate1,intersections);

  // Have we found an intersection?
  return (intersections.size() > 0 || neighborIntersects1.any() || neighborIntersects2.any());

}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::bruteForceSearch(int candidate1,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
                                                                  const std::vector<Dune::GeometryType>& grid1_element_types,
                                                                  const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
                                                                  const std::vector<Dune::GeometryType>& grid2_element_types)
{
  std::bitset<(1<<grid1Dim)> neighborIntersects1;
  std::bitset<(1<<grid2Dim)> neighborIntersects2;
  for (std::size_t i=0; i<grid1_element_types.size(); i++) {

    bool intersectionFound = computeIntersection(i, candidate1,
                                                 grid1Coords, grid1_element_types, neighborIntersects1,
                                                 grid2Coords, grid2_element_types, neighborIntersects2,
                                                 false);

    // if there is an intersection, i is our new seed candidate on the grid1 side
    if (intersectionFound)
      return i;

  }

  return -1;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
template<int gridDim>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::
computeNeighborsPerElement(const std::vector<Dune::GeometryType>& gridElementTypes,
                           const std::vector<std::vector<unsigned int> >& gridElementCorners,
                           std::vector<std::vector<int> >& elementNeighbors)
{
  typedef std::vector<unsigned int> FaceType;
  typedef std::map<FaceType, std::pair<unsigned int, unsigned int> > FaceSetType;

  ///////////////////////////////////////////////////////////////////////////////////////
  //  First: grid 1
  ///////////////////////////////////////////////////////////////////////////////////////
  FaceSetType faces;
  elementNeighbors.resize(gridElementTypes.size());

  for (size_t i=0; i<gridElementTypes.size(); i++)
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    elementNeighbors[i].resize(Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]).size(1), -1);
#else
    elementNeighbors[i].resize(Dune::GenericReferenceElements<T,gridDim>::general(gridElementTypes[i]).size(1), -1);
#endif

  for (size_t i=0; i<gridElementTypes.size(); i++) { //iterate over all elements

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    const Dune::ReferenceElement<T,gridDim>& refElement = Dune::ReferenceElements<T,gridDim>::general(gridElementTypes[i]);
#else
    const Dune::GenericReferenceElement<T,gridDim>& refElement = Dune::GenericReferenceElements<T,gridDim>::general(gridElementTypes[i]);
#endif

    for (size_t j=0; j<(size_t)refElement.size(1); j++) { // iterate over all faces of the element

      FaceType face;
      // extract element face
      for (size_t k=0; k<(size_t)refElement.size(j,1,gridDim); k++)
        face.push_back(gridElementCorners[i][refElement.subEntity(j,1,k,gridDim)]);

      // sort the face vertices to get rid of twists and other permutations
      std::sort(face.begin(), face.end());

      typename FaceSetType::iterator faceHandle = faces.find(face);

      if (faceHandle == faces.end()) {

        // face has not been visited before
        faces.insert(std::make_pair(face, std::make_pair(i,j)));

      } else {

        // face has been visited before: store the mutual neighbor information
        elementNeighbors[i][j] = faceHandle->second.first;
        elementNeighbors[faceHandle->second.first][faceHandle->second.second] = i;

        faces.erase(faceHandle);

      }

    }

  }
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
  Dune::Timer watch;

  clear();
  // clear global intersection list
  intersections_.clear();
  this->counter = 0;

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
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    int numVertices = Dune::ReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
#else
    int numVertices = Dune::GenericReferenceElements<T,grid1Dim>::general(grid1_element_types[i]).size(grid1Dim);
#endif
    grid1ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid1ElementCorners_[i][j] = grid1_elements[grid1CornerCounter++];

  }

  // then the grid2 side
  grid2ElementCorners_.resize(grid2_element_types.size());

  unsigned int grid2CornerCounter = 0;

  for (std::size_t i=0; i<grid2_element_types.size(); i++) {

    // Select vertices of the grid2 element
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    int numVertices = Dune::ReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
#else
    int numVertices = Dune::GenericReferenceElements<T,grid2Dim>::general(grid2_element_types[i]).size(grid2Dim);
#endif
    grid2ElementCorners_[i].resize(numVertices);
    for (int j=0; j<numVertices; j++)
      grid2ElementCorners_[i][j] = grid2_elements[grid2CornerCounter++];

  }

  ////////////////////////////////////////////////////////////////////////
  //  Compute the face neighbors for each element
  ////////////////////////////////////////////////////////////////////////

  computeNeighborsPerElement<grid1Dim>(grid1_element_types, grid1ElementCorners_, elementNeighbors1_);
  computeNeighborsPerElement<grid2Dim>(grid2_element_types, grid2ElementCorners_, elementNeighbors2_);

  std::cout << "setup took " << watch.elapsed() << " seconds." << std::endl;

  ////////////////////////////////////////////////////////////////////////
  //   Data structures for the advancing-front algorithm
  ////////////////////////////////////////////////////////////////////////

  std::stack<unsigned int> candidates1;
  std::stack<unsigned int> candidates2;

  std::vector<int> seeds(grid2_element_types.size(), -1);

  // /////////////////////////////////////////////////////////////////////
  //   Do a brute-force search to find one pair of intersecting elements
  //   to start the advancing-front type algorithm with.
  // /////////////////////////////////////////////////////////////////////

  // Set flag if element has been handled
  Dune::BitSetVector<1> isHandled2(grid2_element_types.size());

  // Set flag if the element has been entered in the queue
  Dune::BitSetVector<1> isCandidate2(grid2_element_types.size());

  generateSeed(seeds, isHandled2, candidates2, grid1Coords, grid1_element_types, grid2Coords, grid2_element_types);

  // /////////////////////////////////////////////////////////////////////
  //   Main loop
  // /////////////////////////////////////////////////////////////////////

  std::set<unsigned int> isHandled1;

  std::set<unsigned int> isCandidate1;

  while (!candidates2.empty()) {

    // Get the next element on the grid2 side
    unsigned int currentCandidate2 = candidates2.top();
    int seed = seeds[currentCandidate2];
    assert(seed >= 0);

    candidates2.pop();
    isHandled2[currentCandidate2] = true;

    // Start advancing front algorithm on the grid1 side from the 'seed' element that
    // we stored along with the current grid2 element
    candidates1.push(seed);

    isHandled1.clear();
    isCandidate1.clear();

    while (!candidates1.empty()) {

      unsigned int currentCandidate1 = candidates1.top();
      candidates1.pop();
      isHandled1.insert(currentCandidate1);

      // Test whether there is an intersection between currentCandidate0 and currentCandidate1
      std::bitset<(1<<grid1Dim)> neighborIntersects1;
      std::bitset<(1<<grid2Dim)> neighborIntersects2;
      bool intersectionFound = computeIntersection(currentCandidate1, currentCandidate2,
                                                   grid1Coords,grid1_element_types, neighborIntersects1,
                                                   grid2Coords,grid2_element_types, neighborIntersects2);

      for (size_t i=0; i<neighborIntersects2.size(); i++)
        if (neighborIntersects2[i] && elementNeighbors2_[currentCandidate2][i] != -1)
          seeds[elementNeighbors2_[currentCandidate2][i]] = currentCandidate1;

      // add neighbors of candidate0 to the list of elements to be checked
      if (intersectionFound) {

        for (size_t i=0; i<elementNeighbors1_[currentCandidate1].size(); i++) {

          int neighbor = elementNeighbors1_[currentCandidate1][i];

          if (neighbor == -1)            // do nothing at the grid boundary
            continue;

          if (isHandled1.find(neighbor) == isHandled1.end()
              && isCandidate1.find(neighbor) == isCandidate1.end()) {
            candidates1.push(neighbor);
            isCandidate1.insert(neighbor);
          }

        }

      }

    }

    // We have now found all intersections of elements in the grid1 side with currentCandidate2
    // Now we add all neighbors of currentCandidate2 that have not been treated yet as new
    // candidates.

    // Do we have an unhandled neighbor with a seed?
    bool seedFound = !candidates2.empty();
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      // Add all unhandled intersecting neighbors to the queue
      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0] && seeds[neighbor]>-1) {

        isCandidate2[neighbor][0] = true;
        candidates2.push(neighbor);
        seedFound = true;
      }
    }

    if (seedFound)
      continue;

    // There is no neighbor with a seed, so we need to be a bit more aggressive...
    // get all neighbors of currentCandidate2, but not currentCandidate2 itself
    for (size_t i=0; i<elementNeighbors2_[currentCandidate2].size(); i++) {

      int neighbor = elementNeighbors2_[currentCandidate2][i];

      if (neighbor == -1)         // do nothing at the grid boundary
        continue;

      if (!isHandled2[neighbor][0] && !isCandidate2[neighbor][0]) {

        // Get a seed element for the new grid2 element
        // Look for an element on the grid1 side that intersects the new grid2 element.
        int seed = -1;

        // Look among the ones that have been tested during the last iteration.
        for (typename std::set<unsigned int>::iterator seedIt = isHandled1.begin();
             seedIt != isHandled1.end(); ++seedIt) {

          std::bitset<(1<<grid1Dim)> neighborIntersects1;
          std::bitset<(1<<grid2Dim)> neighborIntersects2;
          bool intersectionFound = computeIntersection(*seedIt, neighbor,
                                                       grid1Coords, grid1_element_types, neighborIntersects1,
                                                       grid2Coords, grid2_element_types, neighborIntersects2,
                                                       false);

          // if the intersection is nonempty, *seedIt is our new seed candidate on the grid1 side
          if (intersectionFound) {
            seed = *seedIt;
            Dune::dwarn << "Algorithm entered first fallback method and found a new seed in the build algorithm." <<
                     "Probably, the neighborIntersects bitsets computed in computeIntersection specialization is wrong." << std::endl;
            break;
          }

        }

        if (seed < 0) {
          // The fast method didn't find a grid1 element that intersects with
          // the new grid2 candidate.  We have to do a brute-force search.
          seed = bruteForceSearch(neighbor,
                                  grid1Coords,grid1_element_types,
                                  grid2Coords,grid2_element_types);
          Dune::dwarn << "Algorithm entered second fallback method. This probably should not happen." << std::endl;

        }

        // We have tried all we could: the candidate is 'handled' now
        isCandidate2[neighbor] = true;

        // still no seed?  Then the new grid2 candidate isn't overlapped by anything
        if (seed < 0)
          continue;

        // we have a seed now
        candidates2.push(neighbor);
        seeds[neighbor] = seed;
        seedFound = true;

      }

    }

    /* Do a brute-force search if there is still no seed:
     * There might still be a disconnected region out there.
     */
    if (!seedFound && candidates2.empty()) {
      generateSeed(seeds, isHandled2, candidates2, grid1Coords, grid1_element_types, grid2Coords, grid2_element_types);
    }
  }

  valid = true;
  std::cout << "intersection construction took " << watch.elapsed() << " seconds." << std::endl;
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::nSimplices() const
{
  assert(valid);
  return intersections_.size();
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1Parents(unsigned int idx) const
{
  assert(valid);
  return (intersections_[idx].grid1Entities_).size();
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2Parents(unsigned int idx) const
{
  assert(valid);
  return (intersections_[idx].grid2Entities_).size();
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1Parent(unsigned int idx, unsigned int parId) const
{
  assert(valid);
  return intersections_[idx].grid1Entities_[parId];
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
inline unsigned int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2Parent(unsigned int idx, unsigned int parId) const
{
  assert(valid);
  return intersections_[idx].grid2Entities_[parId];
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid1Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId) const
{
  assert(valid);
  return intersections_[idx].grid1Local_[parId][corner];
}


template<typename T, int grid1Dim, int grid2Dim, int dimworld>
typename StandardMerge<T,grid1Dim,grid2Dim,dimworld>::Grid2Coords StandardMerge<T,grid1Dim,grid2Dim,dimworld>::grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId) const
{
  assert(valid);
  return intersections_[idx].grid2Local_[parId][corner];
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
void StandardMerge<T,grid1Dim,grid2Dim,dimworld>::generateSeed(std::vector<int>& seeds, Dune::BitSetVector<1>& isHandled2, std::stack<unsigned>& candidates2, const std::vector<Dune::FieldVector<T, dimworld> >& grid1Coords, const std::vector<Dune::GeometryType>& grid1_element_types, const std::vector<Dune::FieldVector<T, dimworld> >& grid2Coords, const std::vector<Dune::GeometryType>& grid2_element_types)
{
  for (std::size_t j=0; j<grid2_element_types.size(); j++) {

    if (seeds[j] > 0 || isHandled2[j][0])
      continue;

    int seed = bruteForceSearch(j,grid1Coords,grid1_element_types,grid2Coords,grid2_element_types);

    if (seed >= 0) {
      candidates2.push(j);        // the candidate and a seed for the candidate
      seeds[j] = seed;
      break;
    } else // If the brute force search did not find any intersection we can skip this element
      isHandled2[j] = true;
  }
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::insertIntersections(unsigned int candidate1, unsigned int candidate2,
                                                                        std::vector<RemoteSimplicialIntersection>& intersections)
{
    typedef typename std::vector<RemoteSimplicialIntersection>::size_type size_t;
    int count = 0;

    for (size_t i = 0; i < intersections.size(); ++i) {
        // get the intersection index of the current intersection from intersections in this->intersections
        int index = intersectionIndex(candidate1,candidate2,intersections[i]);

        if (index >= this->intersections_.size()) { //the intersection is not yet contained in this->intersections
            this->intersections_.push_back(intersections[i]);   // insert

            ++count;
        } else if (index > -1) {
            // insert each grid1 element and local representation of intersections[i] with parent candidate1
            for (size_t j = 0; j < intersections[i].grid1Entities_.size(); ++j) {
                this->intersections_[index].grid1Entities_.push_back(candidate1);
                this->intersections_[index].grid1Local_.push_back(intersections[i].grid1Local_[j]);
            }

            // insert each grid2 element and local representation of intersections[i] with parent candidate2
            for (size_t j = 0; j < intersections[i].grid2Entities_.size(); ++j) {
                this->intersections_[index].grid2Entities_.push_back(candidate2);
                this->intersections_[index].grid2Local_.push_back(intersections[i].grid2Local_[j]);
            }

            ++count;
        } else {
            Dune::dwarn << "Computed the same intersection twice!" << std::endl;
        }
    }
    return count;
}

template<typename T, int grid1Dim, int grid2Dim, int dimworld>
int StandardMerge<T,grid1Dim,grid2Dim,dimworld>::intersectionIndex(unsigned int grid1Index, unsigned int grid2Index,
                                                                            RemoteSimplicialIntersection& intersection) {


    // return index in intersections_ if at least one local representation of a Remote Simplicial Intersection (RSI)
    // of intersections_ is equal to the local representation of one element in intersections

    std::size_t n_intersections = this->intersections_.size();
    T eps = 1e-10;

    for (std::size_t i = 0; i < n_intersections; ++i) {

        // compare the local representation of the subelements of the RSI
        for (std::size_t ei = 0; ei < this->intersections_[i].grid1Entities_.size(); ++ei) // merger subelement
        {
            if (this->intersections_[i].grid1Entities_[ei] == grid1Index)
            {
                for (std::size_t er = 0; er < intersection.grid1Entities_.size(); ++er) // list subelement
                {
                    bool found_all = true;
                    // compare the local coordinate representations
                    for (std::size_t ci = 0; ci < this->intersections_[i].grid1Local_[ei].size(); ++ci)
                    {
                        Dune::FieldVector<T,grid1Dim> ni = this->intersections_[i].grid1Local_[ei][ci];
                        bool found_ni = false;
                        for (std::size_t cr = 0; cr < intersection.grid1Local_[er].size(); ++cr)
                        {
                            Dune::FieldVector<T,grid1Dim> nr = intersection.grid1Local_[er][cr];

                            found_ni = found_ni || ((ni-nr).infinity_norm() < eps);
                            if (found_ni)
                                break;
                        }
                        found_all = found_all && found_ni;

                        if (!found_ni)
                            break;
                    }

                    if (found_all && (this->intersections_[i].grid2Entities_[ei] != grid2Index))
                        return i;
                    else if (found_all)
                        return -1;
                }
            }
        }

        // compare the local representation of the subelements of the RSI
        for (std::size_t ei = 0; ei < this->intersections_[i].grid2Entities_.size(); ++ei) // merger subelement
        {
            if (this->intersections_[i].grid2Entities_[ei] == grid2Index)
            {
                for (std::size_t er = 0; er < intersection.grid2Entities_.size(); ++er) // list subelement
                {
                    bool found_all = true;
                    // compare the local coordinate representations
                    for (std::size_t ci = 0; ci < this->intersections_[i].grid2Local_[ei].size(); ++ci)
                    {
                        Dune::FieldVector<T,grid2Dim> ni = this->intersections_[i].grid2Local_[ei][ci];
                        bool found_ni = false;
                        for (std::size_t cr = 0; cr < intersection.grid2Local_[er].size(); ++cr)
                        {
                            Dune::FieldVector<T,grid2Dim> nr = intersection.grid2Local_[er][cr];
                            found_ni = found_ni || ((ni-nr).infinity_norm() < eps);

                            if (found_ni)
                                break;
                        }
                        found_all = found_all && found_ni;

                        if (!found_ni)
                            break;
                    }

                    if (found_all && (this->intersections_[i].grid1Entities_[ei] != grid1Index))
                        return i;
                    else if (found_all)
                        return -1;
                }
            }
        }
    }

    return n_intersections;
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
