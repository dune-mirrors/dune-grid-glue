// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_MERGER_HH
#define DUNE_MERGER_HH

#include <vector>

#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

/** \brief Abstract base for all classes that take extracted grids and build sets of intersections

   \tparam ctype The type used for coordinates (assumed to be the same for both grids)
   \tparam domainDim Dimension of the domain grid
   \tparam targetDim Dimension of the target grid
   \tparam dimworld Dimension of the world space where the coupling takes place
 */
template <class ctype, int domainDim, int targetDim, int dimworld>
class Merger
{
public:

  /// @brief the local coordinate type for the domain coordinates
  typedef Dune::FieldVector<ctype, domainDim>  DomainCoords;

  /// @brief the local coordinate type for the target coordinates
  typedef Dune::FieldVector<ctype, targetDim>  TargetCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<ctype, dimworld>  WorldCoords;


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
  virtual void build(const std::vector<Dune::FieldVector<ctype,dimworld> >& domain_coords,
                     const std::vector<unsigned int>& domain_elements,
                     const std::vector<Dune::GeometryType>& domain_element_types,
                     const std::vector<Dune::FieldVector<ctype,dimworld> >& target_coords,
                     const std::vector<unsigned int>& target_elements,
                     const std::vector<Dune::GeometryType>& target_element_types) = 0;

  /** @brief get the number of simplices in the merged grid
      The indices are then in 0..nSimplices()-1
   */
  virtual unsigned int nSimplices() const = 0;

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
  virtual bool domainSimplexMatched(unsigned int idx) const = 0;


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
  virtual bool targetSimplexMatched(unsigned int idx) const = 0;



  /**
   * @brief get index of domain parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the domain parent simplex
   */
  virtual unsigned int domainParent(unsigned int idx) const = 0;

  /**
   * @brief get index of target parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the target parent simplex
   */
  virtual unsigned int targetParent(unsigned int idx) const = 0;


  /**
   * @brief get the merged grid simplices refining a given domain simplex
   * @param idx index of domain simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  virtual bool domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const = 0;


  /**
   * @brief get the merged grid simplices refining a given target simplex
   * @param idx index of target simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  virtual bool targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const = 0;


  /**
   * @brief get the domain parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "domainParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent domain
   */
  virtual DomainCoords domainParentLocal(unsigned int idx, unsigned int corner) const = 0;

  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent target
   */
  virtual TargetCoords targetParentLocal(unsigned int idx, unsigned int corner) const = 0;

};

#endif
