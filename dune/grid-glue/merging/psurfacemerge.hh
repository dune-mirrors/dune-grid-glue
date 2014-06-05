// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file
   @author Gerrit Buse, Christian Engwer, Oliver Sander
   @brief implementation of the SurfaceMerge concept based on libpsurface

   This implementation of the SurfaceMerge concept can be used to
   compute 1d, 2d and 3d intersections. It uses psurface routines to
   compute the merged grid and provides access to it via the interface
   specified by the concept.
 */

#ifndef PSURFACEMERGE_HH
#define PSURFACEMERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/merger.hh>

#include <tr1/array>
#if HAVE_PSURFACE
#include <psurface/DirectionFunction.h>
#include <psurface/IntersectionPrimitive.h>
#else
// forward declaration of PSurface classes
template <int dim, typename ctype> class DirectionFunction;
// switch off the macro that contains (in certain versions) the psurface namespace prefix
#define PSURFACE_NAMESPACE
#endif


/** \brief Standard implementation of the SurfaceMerge concept using the psurface library.

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.  Must be equal to dim or to dim+1
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class PSurfaceMerge
  : public Merger<T,dim,dim,dimworld>
{
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,3,0)
  static_assert( dim==1 || dim==2,
                      "PSurface can only handle the cases dim==1 and dim==2!");

  static_assert( dim==dimworld || dim+1==dimworld,
                      "PSurface can only handle the cases dim==dimworld and dim+1==dimworld!");
#else
  dune_static_assert( dim==1 || dim==2,
                      "PSurface can only handle the cases dim==1 and dim==2!");

  dune_static_assert( dim==dimworld || dim+1==dimworld,
                      "PSurface can only handle the cases dim==dimworld and dim+1==dimworld!");
#endif

  // The psurface library itself always expects dimworld to be dim+1
  // To be able to handle the case dim==dimworld we keep an artificial world
  // around with the additional dimension
  enum {psurfaceDimworld = dimworld + (dim==dimworld)};

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;


private:

#if HAVE_PSURFACE

  /*   P R I V A T E   H E L P E R S   */

  /**
   * @brief transforms barycentric coordinates to Dune-style local coordinates
   */
  static Dune::FieldVector<ctype, dim> barycentricToReference(const Dune::FieldVector<ctype, dim+1>& bar)
  {
    Dune::FieldVector<ctype, dim> result;
    for (int i=0; i<dim; i++)
      result[i] = bar[i+1];

    return result;
  }


  /**
   * @brief transforms Dune-style local coordinates to barycentric coordinates.
   */
  static Dune::FieldVector<ctype, dim+1> referenceToBarycentric(const Dune::FieldVector<ctype, dim>& ref)
  {
    Dune::FieldVector<ctype, dim+1> result;
    result[0] = 1.0;
    for (int i=0; i<dim; i++) {
      result[i+1] = ref[i];
      result[0] -= ref[i];
    }

    return result;
  }


  /**
   * @brief compare function dep. on parent domain simplex
   * @param a first item
   * @param b second item
   * @return TRUE <=> parent domain simplex' index of @c a smaller
   */
  static bool domainParentSmaller(const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& a, const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& b)
  {
    return a.tris[0] < b.tris[0];
  }

  /**
   * @brief compare function dep. on parent target simplex
   * @param a first item
   * @param b second item
   * @return TRUE <=> parent target simplex' index of @c a smaller
   */
  static bool targetParentSmaller(const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>* a, const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>* b)
  {
    return a->tris[1] < b->tris[1];
  }


  /*   P R I V A T E   S U B C L A S S E S   */

  /**
   * @class OverlapManager
   * @brief provides convenient and efficient access to the merged grid's entities
   */
  struct OverlapManager
  {

    OverlapManager()
    {}


    /**
     * @brief sets up the internal structures and containers using the overlaps
     * from the given array
     * @param unordered array containing all overlap simplices
     */
    void setOverlaps(const std::vector<PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype> >& unordered);

    /**
     * @brief getter for the number of simplex overlaps in the merged grid
     */
    unsigned int nOverlaps() const
    {
      return this->domOrder.size();
    }

    /**
     * Getter for the index of the first simplex in the domain parent sorted
     * array that has the given target simplex as parent.
     * @param parent the index of the domain parent simplex
     * @return the index in the domain parent sorted array
     */
    unsigned int firstDomainParent(unsigned int parent) const;

    /**
     * Getter for the index of the first simplex in the target parent sorted
     * array that has the given target simplex as parent.
     * @param parent the index of the target parent simplex
     * @return the index in the target parent sorted array
     */
    unsigned int firstTargetParent(unsigned int parent) const;

    /**
     * @brief getter for an simplex in the array sorted after the domain parent indices
     * @param idx index into the array
     * @return the simplex
     */
    const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& domain(unsigned int idx) const
    {
      return this->domOrder[idx];
    }

    /**
     * @brief getter for an simplex in the array sorted after the target parent indices
     * @param idx index into the array
     * @return the simplex
     */
    const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& target(unsigned int idx) const
    {
      return *this->tarOrder[idx];
    }

    /**
     * @brief getter for the DOMINATING domain index of an simplex
     * @param idx the target index of an simplex
     * @return the actual index (used in the domain list)
     */
    unsigned int domainIndex(unsigned int idx) const
    {
      return (this->tarOrder[idx] - this->baseptr);
    }

    /** @brief clear the internal state
     */
    void clear()
    {
      purge(domOrder);
      purge(tarOrder);
    }

  private:

    template<typename V>
    static void purge(V & v)
    {
      v.clear();
      V v2(v);
      v.swap(v2);
    }

    /// @brief the computed merged grid simplices sorted after
    /// ascending domain parent simplex indices (see domi_)
    std::vector<PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype> > domOrder;

    /// @brief the computed merged grid simplices sorted after
    /// ascending target parent simplex indices (see tari_)
    std::vector<PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>*> tarOrder;

    /// @brief used to recompute original index of an overlap in "tarOrder"
    PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>*         baseptr;
  };

  /** \brief Return a constant direction field for PSurface

     This is used for overlapping couplings.  We need the direction fields (0,0,1) and (0,0,-1)
     \tparam dir Value of the last component: either -1 or 1.
   */
  template <int dir>
  struct ConstantDirection : public PSURFACE_NAMESPACE AnalyticDirectionFunction<psurfaceDimworld,ctype>
  {
    PSURFACE_NAMESPACE StaticVector<ctype,psurfaceDimworld> operator()(const PSURFACE_NAMESPACE StaticVector<ctype,psurfaceDimworld>& position) const
    {
      PSURFACE_NAMESPACE StaticVector<ctype,psurfaceDimworld> result;
      for (size_t i=0; i<psurfaceDimworld-1; i++)
        result[i] = 0;
      result[psurfaceDimworld-1] = dir;
      return result;
    }
  };

private:

  /*   M E M B E R   V A R I A B L E S   */

  /* geometric data for both domain and target */

  /// @brief domain coordinates
  std::vector<std::tr1::array<double,psurfaceDimworld> >   domc_;

  /// @brief target coordinates
  std::vector<std::tr1::array<double,psurfaceDimworld> >   tarc_;


  /* topologic information for domain and target */

  /// @ brief domain indices (internal copy)
  std::vector<std::tr1::array<int,dim+1> >         domi_;

  /// @brief target indices (internal copy)
  std::vector<std::tr1::array<int,dim+1> >         tari_;


  /* members associated with the merged grid */

  /// @brief provides an interface for convenient and efficient
  /// access to the set of computed simplex overlaps in the merged grid
  OverlapManager olm_;

  /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
   */
  const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* domainDirections_;

  /** \brief Vector field on the target surface which prescribes a 'forward'
      direction.

      PSurface uses the normals of the target side to increase projection
      robustness.  If these cannot be computed from the surface directly
      (e.g. because it is not properly oriented), they can be given
      explicitly through the targetDirections field.
   */
  const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* targetDirections_;

  bool valid;

#endif // if HAVE_PSURFACE

public:

  PSurfaceMerge(const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* domainDirections = NULL,
                const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* targetDirections = NULL);

  /*   M O D E L   S P E C I F I C   E X T E N D I N G   F U N C T I O N A L I T Y   */

  /**
   * @brief Set surface direction functions
   *
   * The matching of the geometries offers the possibility to specify a function for
   * the exact evaluation of domain surface normals. If no such function is specified
   * (default) normals are interpolated.
   * @param domainDirections the new function for the outer normal of grid0 (domain) (or NULL to unset the function)
   * @param targetDirections the new function for the outer normal of grid1 (domain) (or NULL to unset the function)
   */
  inline
  void setSurfaceDirections(const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* domainDirections,
                            const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* targetDirections);

  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @copydoc Merger<T,dim,dim,dimworld>::build
   */
  void build(const std::vector<Dune::FieldVector<ctype,dimworld> >& grid1_coords,
             const std::vector<unsigned int>& grid1_elements,
             const std::vector<Dune::GeometryType>& grid1_element_types,
             const std::vector<Dune::FieldVector<ctype,dimworld> >& grid2_coords,
             const std::vector<unsigned int>& grid2_elements,
             const std::vector<Dune::GeometryType>& grid2_element_types);


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const;

  /// @brief clear the internal state, so that we can run a new merging process
  void clear()
  {
#ifdef HAVE_PSURFACE
    // Delete old internal data, from a possible previous run
    purge(domi_);
    purge(tari_);
    purge(domc_);
    purge(domi_);

    olm_.clear();

    valid = false;
#endif
  }

private:

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
  unsigned int grid1Parent(unsigned int idx, unsigned int parId = 0) const;

  /**
   * @brief get index of target parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the target parent simplex
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
  LocalCoords grid1ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return local coordinates in parent target simplex
   */
  LocalCoords grid2ParentLocal(unsigned int idx, unsigned int corner, unsigned int parId = 0) const;

};

#if HAVE_PSURFACE

/* IMPLEMENTATION OF CLASS   C O N T A C T  M A P P I N G  S U R F A C E  M E R G E */

template<int dim, int dimworld, typename T>
PSurfaceMerge<dim, dimworld, T>::PSurfaceMerge(const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* domainDirections,
                                               const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* targetDirections)
  : domainDirections_(domainDirections), targetDirections_(targetDirections), valid(false)
{}


template<int dim, int dimworld, typename T>
inline void PSurfaceMerge<dim, dimworld, T>::setSurfaceDirections(const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* domainDirections,
                                                                  const PSURFACE_NAMESPACE DirectionFunction<psurfaceDimworld,ctype>* targetDirections)
{
  domainDirections_ = domainDirections;
  targetDirections_ = targetDirections;
  valid = false;
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::nSimplices() const
{
  assert(valid);
  return this->olm_.nOverlaps();
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::grid1Parent(unsigned int idx, unsigned int parId) const
{
  assert(valid);
  return this->olm_.domain(idx).tris[0];
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::grid2Parent(unsigned int idx, unsigned int parId) const
{
  assert(valid);
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return this->olm_.domain(idx).tris[1];
}

#endif // HAVE_PSURFACE

#ifdef PSURFACE_EXTRA_TYPES
#define PSURFACE_EXTERN
#include "psurfacemerge.cc"
#undef PSURFACE_EXTERN
#endif

#endif // PSURFACEMERGE_HH
