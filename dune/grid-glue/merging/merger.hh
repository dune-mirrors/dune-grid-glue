// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_MERGING_MERGER_HH
#define DUNE_GRIDGLUE_MERGING_MERGER_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>

// forward declaration
template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
class Merger;

namespace {

  // forward calls to the (grid2/grid1)-calls depending on the grid number (0/1)
  template <class ctype, int grid1Dim, int grid2Dim, int dimworld, int n>
  struct MergerGridPolicy;

  template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
  struct MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 0>
  {
    typedef Merger<ctype, grid1Dim, grid2Dim, dimworld> Parent;

    /// @brief the local coordinate type for the grid1 coordinates
    typedef Dune::FieldVector<ctype, grid1Dim>  GridCoords;

    static
    unsigned int parents(const Parent & m, unsigned int idx)
    {
        return m.grid1Parents(idx);
    }

    static
    unsigned int parent(const Parent & m, unsigned int idx, unsigned int parId = 0)
    {
      return m.grid1Parent(idx, parId);
    }

    static
    GridCoords parentLocal(const Parent & m, unsigned int idx, unsigned int corner, unsigned int parId = 0)
    {
      return m.grid1ParentLocal(idx, corner, parId);
    }

  };

  template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
  struct MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 1>
  {
    typedef Merger<ctype, grid1Dim, grid2Dim, dimworld> Parent;

    /// @brief the local coordinate type for the grid2 coordinates
    typedef Dune::FieldVector<ctype, grid2Dim>  GridCoords;

    static
    unsigned int parents(const Parent & m, unsigned int idx)
    {
      return m.grid2Parents(idx);
    }

    static
    unsigned int parent(const Parent & m, unsigned int idx, unsigned int parId = 0)
    {
      return m.grid2Parent(idx, parId);
    }

    static
    GridCoords parentLocal(const Parent & m, unsigned int idx, unsigned int corner, unsigned int parId = 0)
    {
      return m.grid2ParentLocal(idx, corner, parId);
    }

  };

} // end empty namespace

/** \brief Abstract base for all classes that take extracted grids and build sets of intersections

   \tparam ctype The type used for coordinates (assumed to be the same for both grids)
   \tparam grid1Dim Dimension of the grid1 grid
   \tparam grid2Dim Dimension of the grid2 grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
class Merger
{

  // the policy class get's access to the Merger
  friend struct MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 0>;
  friend struct MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 1>;

public:

  /// @brief the local coordinate type for the grid1 coordinates
  typedef Dune::FieldVector<ctype, grid1Dim>  Grid1Coords;

  /// @brief the local coordinate type for the grid2 coordinates
  typedef Dune::FieldVector<ctype, grid2Dim>  Grid2Coords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<ctype, dimworld>  WorldCoords;

  template<int n>
  struct GridTraits
  {
    /// @brief the policy class for this grid number
    typedef MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, n> Policy;
    /// @brief the local coordinate type for the grid-n coordinates
    typedef typename MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, n>::GridCoords Coords;
  };

  /**
   * @brief builds the merged grid
   *
   * Note that the indices are used consequently throughout the whole class interface just like they are
   * introduced here.
   *
   * @param grid1_coords the grid1 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid1_elements array with all grid1 elements represented as corner indices into @c grid1_coords
   * @param grid1_element_types array with the GeometryType of the elements listed grid1_elements
   * @param grid2_coords the grid2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid2_elements just like with the grid1_elements and grid1_coords
   * @param grid2_element_types array with the GeometryType of the elements listed grid2_elements
   */
  virtual void build(const std::vector<Dune::FieldVector<ctype,dimworld> >& grid1_coords,
                     const std::vector<unsigned int>& grid1_elements,
                     const std::vector<Dune::GeometryType>& grid1_element_types,
                     const std::vector<Dune::FieldVector<ctype,dimworld> >& grid2_coords,
                     const std::vector<unsigned int>& grid2_elements,
                     const std::vector<Dune::GeometryType>& grid2_element_types) = 0;

  /** @brief get the number of simplices in the merged grid
      The indices are then in 0..nSimplices()-1
   */
  virtual unsigned int nSimplices() const = 0;

  virtual void clear() = 0;

   /**
    * doc me
    */
  template<int n>
  unsigned int parents(unsigned int idx) const {
    return GridTraits<n>::Policy::parents(*this, idx);
  }

  /**
   * @brief get index of grid-n's parent simplex for given merged grid simplex
   * @tparam n specify which grid
   * @param idx index of the merged grid simplex
   * @return index of the parent simplex
   */
  template<int n>
  unsigned int parent(unsigned int idx, unsigned int parId = 0) const
  {
    return GridTraits<n>::Policy::parent(*this, idx, parId);
  }

  /**
   * @brief get the merged grid simplices refining a given grid-n simplex
   * @tparam n specify which grid (grid1/grid2: 0/1)
   * @param idx index of grid-n simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  template<int n>
  bool simplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
  {
    return GridTraits<n>::Policy::simplexRefined(*this, idx, indices);
  }

  /**
   * @brief get the grid-n parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "parent<n>")
   * @tparam n specify which grid
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in grid-n grid1
   */
  template<int n>
  typename GridTraits<n>::Coords parentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const
  {
    return GridTraits<n>::Policy::parentLocal(*this, idx, corner, parId);
  }

  /** \brief Counts the number of times the computeIntersection method has been called
   *
   * Used temporarily to speed up the implementation
   */
  unsigned int counter;


private:

  virtual unsigned int grid1Parents(unsigned int idx) const = 0;

  virtual unsigned int grid2Parents(unsigned int idx) const = 0;

  /**
   * @brief get index of grid1 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid1 parent simplex
   */
  virtual unsigned int grid1Parent(unsigned int idx, unsigned int parId = 0) const = 0;

  /**
   * @brief get index of grid2 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid2 parent simplex
   */
  virtual unsigned int grid2Parent(unsigned int idx, unsigned int parId = 0) const = 0;

  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1
   */
  virtual Grid1Coords grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const = 0;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2
   */
  virtual Grid2Coords grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const = 0;

};

#endif
