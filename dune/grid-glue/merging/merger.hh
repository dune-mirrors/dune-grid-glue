// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MERGER_HH
#define DUNE_MERGER_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

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
    bool simplexMatched(const Parent & m, unsigned int idx)
    {
      return m.grid1SimplexMatched(idx);
    }

    static
    unsigned int parent(const Parent & m, unsigned int idx)
    {
      return m.grid1Parent(idx);
    }

    static
    GridCoords parentLocal(const Parent & m, unsigned int idx, unsigned int corner)
    {
      return m.grid1ParentLocal(idx, corner);
    }

    static
    bool simplexRefined(const Parent & m, unsigned int idx, std::vector<unsigned int>& indices)
    {
      return m.grid1SimplexRefined(idx, indices);
    }

  };

  template <class ctype, int grid1Dim, int grid2Dim, int dimworld>
  struct MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 1>
  {
    typedef Merger<ctype, grid1Dim, grid2Dim, dimworld> Parent;

    /// @brief the local coordinate type for the grid2 coordinates
    typedef Dune::FieldVector<ctype, grid2Dim>  GridCoords;

    static
    bool simplexMatched(const Parent & m, unsigned int idx)
    {
      return m.grid2SimplexMatched(idx);
    }

    static
    unsigned int parent(const Parent & m, unsigned int idx)
    {
      return m.grid2Parent(idx);
    }

    static
    GridCoords parentLocal(const Parent & m, unsigned int idx, unsigned int corner)
    {
      return m.grid2ParentLocal(idx, corner);
    }

    static
    bool simplexRefined(const Parent & m, unsigned int idx, std::vector<unsigned int>& indices)
    {
      return m.grid2SimplexRefined(idx, indices);
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
  friend class MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 0>;
  friend class MergerGridPolicy<ctype, grid1Dim, grid2Dim, dimworld, 1>;

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
   * @param grid1_simplices array with all grid1 simplices represented as corner indices into @c grid1_coords;
   * the simplices are just written to this array one after another
   * @param grid2_coords the grid2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param grid2_simplices just like with the grid1_simplices and grid1_coords
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
   * @brief check if given grid-n simplex could be matched in the merged grid
   *
   * @tparam n specify which grid
   *
   * The result of this member even is positive if a grid-n simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "localToMerged<n>" or "globalToMerged<n>".
   * @param idx the index of the grid simplex
   * @return TRUE <=> refined in merged grid
   */
  template<int n>
  bool simplexMatched(unsigned int idx) const
  {
    return GridTraits<n>::Policy::simplexMatched(*this, idx);
  }

  /**
   * @brief get index of grid-n's parent simplex for given merged grid simplex
   * @tparam n specify which grid
   * @param idx index of the merged grid simplex
   * @return index of the parent simplex
   */
  template<int n>
  unsigned int parent(unsigned int idx) const
  {
    return GridTraits<n>::Policy::parent(*this, idx);
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
  typename GridTraits<n>::Coords parentLocal(unsigned int idx, unsigned int corner) const
  {
    return GridTraits<n>::Policy::parentLocal(*this, idx, corner);
  }

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
  virtual bool grid1SimplexMatched(unsigned int idx) const = 0;


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
  virtual bool grid2SimplexMatched(unsigned int idx) const = 0;



  /**
   * @brief get index of grid1 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid1 parent simplex
   */
  virtual unsigned int grid1Parent(unsigned int idx) const = 0;

  /**
   * @brief get index of grid2 parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the grid2 parent simplex
   */
  virtual unsigned int grid2Parent(unsigned int idx) const = 0;


  /**
   * @brief get the merged grid simplices refining a given grid1 simplex
   * @param idx index of grid1 simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  virtual bool grid1SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const = 0;


  /**
   * @brief get the merged grid simplices refining a given grid2 simplex
   * @param idx index of grid2 simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  virtual bool grid2SimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const = 0;


  /**
   * @brief get the grid1 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid1Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid1
   */
  virtual Grid1Coords grid1ParentLocal(unsigned int idx, unsigned int corner) const = 0;

  /**
   * @brief get the grid2 parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "grid2Parent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent grid2
   */
  virtual Grid2Coords grid2ParentLocal(unsigned int idx, unsigned int corner) const = 0;

};

#endif
