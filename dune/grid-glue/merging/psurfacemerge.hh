// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    PSurfaceMerge.hh
 *  Version:     1.0
 *  Created on:  Jan 25, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: standard implementation of the SurfaceMerge concept
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Standard implementation of the SurfaceMerge concept for the use in 2d and 3d.
 *
 * Uses psurface routines to compute the merged grid  and provides access to it
 * via the interface specified by the concept.
 *
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

#include <psurface/ContactMapping.h>
#include <psurface/DirectionFunction.h>


/** \brief Standard implementation of the SurfaceMerge concept using the psurface library.

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.  Must be equal to dim or to dim+1
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class PSurfaceMerge
  : public Merger<T,dim,dim,dimworld>
{
  dune_static_assert( dim==1 || dim==2,
                      "PSurface can only handle the cases dim==1 and dim==2!");

  dune_static_assert( dim==dimworld || dim+1==dimworld,
                      "PSurface can only handle the cases dim==dimworld and dim+1==dimworld!");

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
  static bool domainParentSmaller(const IntersectionPrimitive<dim,ctype>& a, const IntersectionPrimitive<dim,ctype>& b)
  {
    return a.tris[0] < b.tris[0];
  }

  /**
   * @brief compare function dep. on parent target simplex
   * @param a first item
   * @param b second item
   * @return TRUE <=> parent target simplex' index of @c a smaller
   */
  static bool targetParentSmaller(const IntersectionPrimitive<dim,ctype>* a, const IntersectionPrimitive<dim,ctype>* b)
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
    void setOverlaps(const std::vector<IntersectionPrimitive<dim,ctype> >& unordered);

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
    const IntersectionPrimitive<dim,ctype>& domain(unsigned int idx) const
    {
      return this->domOrder[idx];
    }

    /**
     * @brief getter for an simplex in the array sorted after the target parent indices
     * @param idx index into the array
     * @return the simplex
     */
    const IntersectionPrimitive<dim,ctype>& target(unsigned int idx) const
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

  private:

    /// @brief the computed merged grid simplices sorted after
    /// ascending domain parent simplex indices (see domi_)
    std::vector<IntersectionPrimitive<dim,ctype> > domOrder;

    /// @brief the computed merged grid simplices sorted after
    /// ascending target parent simplex indices (see tari_)
    std::vector<IntersectionPrimitive<dim,ctype>*> tarOrder;

    /// @brief used to recompute original index of an overlap in "tarOrder"
    IntersectionPrimitive<dim,ctype>*         baseptr;
  };

  /** \brief Return a constant direction field for PSurface

     This is used for overlapping couplings.  We need the direction fields (0,0,1) and (0,0,-1)
     \tparam dir Value of the last component: either -1 or 1.
   */
  template <int dir>
  struct ConstantDirection : public AnalyticDirectionFunction<psurfaceDimworld,ctype>
  {
    StaticVector<ctype,psurfaceDimworld> operator()(const StaticVector<ctype,psurfaceDimworld>& position) const
    {
      StaticVector<ctype,psurfaceDimworld> result;
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

  /// @brief make use of psurface contact mapping functionality
  ContactMapping<dim+1,ctype> cm_;

  /// @brief provides an interface for convenient and efficient
  /// access to the set of computed simplex overlaps in the merged grid
  OverlapManager olm_;

  /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
   */
  const DirectionFunction<psurfaceDimworld,ctype>* domainDirections_;

  /** \brief Vector field on the target surface which prescribes a 'forward'
      direction.

      PSurface uses the normals of the target side to increase projection
      robustness.  If these cannot be computed from the surface directly
      (e.g. because it is not properly oriented), they can be given
      explicitly through the targetDirections field.
   */
  const DirectionFunction<psurfaceDimworld,ctype>* targetDirections_;


public:

  PSurfaceMerge(const DirectionFunction<psurfaceDimworld,ctype>* domainDirections = NULL,
                const DirectionFunction<psurfaceDimworld,ctype>* targetDirections = NULL
                )
    : domainDirections_(domainDirections), targetDirections_(targetDirections)
  {}


  /*   M O D E L   S P E C I F I C   E X T E N D I N G   F U N C T I O N A L I T Y   */

  /**
   * @brief Set surface direction functions
   *
   * The matching of the geometries offers the possibility to specify a function for
   * the exact evaluation of domain surface normals. If no such function is specified
   * (default) normals are interpolated.
   * @param value the new function (or NULL to unset the function)
   */
  void setSurfaceDirections(const DirectionFunction<psurfaceDimworld,ctype>* domainDirections,
                            const DirectionFunction<psurfaceDimworld,ctype>* targetDirections)
  {
    domainDirections_ = domainDirections;
    targetDirections_ = targetDirections;
  }


  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @brief builds the merged grid
   *
   * Note that the indices are used consequently throughout the whole class interface just like they are
   * introduced here.
   *
   * @param domain_coords the domain vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param domain_elements array with all domain simplices represented as corner indices into @c domain_coords;
   * the simplices are just written to this array one after another
   * @param target_coords the target vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param target_elements just like with the domain_elements and domain_coords
   */
  void build(const std::vector<Dune::FieldVector<T,dimworld> >& domain_coords,
             const std::vector<unsigned int>& domain_elements,
             const std::vector<Dune::GeometryType>& domain_element_types,
             const std::vector<Dune::FieldVector<T,dimworld> >& target_coords,
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
void PSurfaceMerge<dim, dimworld, T>::build(const std::vector<Dune::FieldVector<T,dimworld> >& domain_coords,
                                            const std::vector<unsigned int>& domain_elements,
                                            const std::vector<Dune::GeometryType>& domain_element_types,
                                            const std::vector<Dune::FieldVector<T,dimworld> >& target_coords,
                                            const std::vector<unsigned int>& target_elements,
                                            const std::vector<Dune::GeometryType>& target_element_types
                                            )
{
  // //////////////////////////////////////////
  //   check input data for consistency
  // //////////////////////////////////////////
#ifndef NDEBUG
  // first the domain side
  unsigned int expectedCorners = 0;

  for (size_t i=0; i<domain_element_types.size(); i++) {

    // Check whether GeometryType has the correct dimension
    if (domain_element_types[i].dim() != dim)
      DUNE_THROW(Dune::GridError, "You cannot hand a " << domain_element_types[i]
                                                       << " to a " << dim << "-dimensional PSurfaceMerge!");

    expectedCorners += Dune::GenericReferenceElements<ctype,dim>::general(domain_element_types[i]).size(dim);

  }

  if (domain_elements.size() != expectedCorners)
    DUNE_THROW(Dune::GridError, domain_elements.size() << " element corners were handed over, "
               " but " << expectedCorners << " were expected!");

  // now the same for the target side
  expectedCorners = 0;

  for (size_t i=0; i<target_element_types.size(); i++) {

    // Check whether GeometryType has the correct dimension
    if (target_element_types[i].dim() != dim)
      DUNE_THROW(Dune::GridError, "You cannot hand a " << target_element_types[i]
                                                       << " to a " << dim << "-dimensional PSurfaceMerge!");

    expectedCorners += Dune::GenericReferenceElements<ctype,dim>::general(target_element_types[i]).size(dim);

  }

  if (target_elements.size() != expectedCorners)
    DUNE_THROW(Dune::GridError, target_elements.size() << " element corners were handed over, "
               " but " << expectedCorners << " were expected!");
#endif

  // //////////////////////////////////////////////////////////////////////////////////
  //   copy domain and target simplices to internal arrays
  //   split quadrilaterals to triangles because psurface can only handle triangles
  // //////////////////////////////////////////////////////////////////////////////////

  // First count the number of simplices we expect
  size_t numDomainSimplices = 0;
  size_t numTargetSimplices = 0;

  if (dim==1) {
    numDomainSimplices = domain_element_types.size();
    numTargetSimplices = target_element_types.size();
  } else {

    for (size_t i=0; i<domain_element_types.size(); i++)
      if (domain_element_types[i].isSimplex())
        numDomainSimplices++;
      else if (domain_element_types[i].isQuadrilateral())
        numDomainSimplices += 2;
      else
        DUNE_THROW(Dune::GridError, "Unknown element " << domain_element_types[i] << " found!");

    for (size_t i=0; i<target_element_types.size(); i++)
      if (target_element_types[i].isSimplex())
        numTargetSimplices++;
      else if (target_element_types[i].isQuadrilateral())
        numTargetSimplices += 2;
      else
        DUNE_THROW(Dune::GridError, "Unknown element " << target_element_types[i] << " found!");

  }

  this->domi_.resize(numDomainSimplices);
  this->tari_.resize(numTargetSimplices);

  // copy and split the domain surface
  Dune::BitSetVector<1> domainIsSecondTriangle(numDomainSimplices, false);
  std::vector<int> unsplitDomainElementNumbers(numDomainSimplices);

  size_t vertexCounter = 0;
  size_t simplexCounter = 0;

  for (size_t i = 0; i < domain_element_types.size(); ++i) {

    // dim is known at compile time, hence if dim==1 the following conditional is removed
    if (dim == 1 || domain_element_types[i].isSimplex()) {
      for (int j=0; j<dim+1; j++)
        domi_[simplexCounter][j] = domain_elements[vertexCounter++];

      unsplitDomainElementNumbers[simplexCounter] = i;

      simplexCounter++;

    } else {
      // quadrilateral: split it in two triangles
      domi_[simplexCounter][0] = domain_elements[vertexCounter];
      domi_[simplexCounter][1] = domain_elements[vertexCounter+1];
      domi_[simplexCounter][2] = domain_elements[vertexCounter+2];

      unsplitDomainElementNumbers[simplexCounter] = i;
      simplexCounter++;

      domi_[simplexCounter][0] = domain_elements[vertexCounter+3];
      domi_[simplexCounter][1] = domain_elements[vertexCounter+2];
      domi_[simplexCounter][2] = domain_elements[vertexCounter+1];

      domainIsSecondTriangle[simplexCounter] = true;
      unsplitDomainElementNumbers[simplexCounter] = i;

      simplexCounter++;
      vertexCounter += 4;
    }

  }

  // copy and split the target surface
  Dune::BitSetVector<1> targetIsSecondTriangle(numTargetSimplices, false);
  std::vector<int> unsplitTargetElementNumbers(numTargetSimplices);

  vertexCounter  = 0;
  simplexCounter = 0;

  for (size_t i = 0; i < target_element_types.size(); ++i) {

    // dim is known at compile time, hence if dim==1 the following conditional is removed
    if (dim==1 || target_element_types[i].isSimplex()) {
      for (size_t j=0; j<dim+1; j++)
        tari_[simplexCounter][j] = target_elements[vertexCounter++];

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;
    } else {
      // quadrilateral: split it in two triangles
      tari_[simplexCounter][0] = target_elements[vertexCounter];
      tari_[simplexCounter][1] = target_elements[vertexCounter+1];
      tari_[simplexCounter][2] = target_elements[vertexCounter+2];

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;

      tari_[simplexCounter][0] = target_elements[vertexCounter+3];
      tari_[simplexCounter][1] = target_elements[vertexCounter+2];
      tari_[simplexCounter][2] = target_elements[vertexCounter+1];

      targetIsSecondTriangle[simplexCounter] = true;

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;
      vertexCounter += 4;
    }

  }

  // ////////////////////////////////////////////////////////////
  //   copy the coordinates to internal arrays of coordinates
  // ////////////////////////////////////////////////////////////
  this->domc_.resize(domain_coords.size());
  for (size_t i = 0; i < this->domc_.size(); ++i)
    for (size_t j = 0; j < dimworld; ++j)
      this->domc_[i][j] = domain_coords[i][j];

  this->tarc_.resize(target_coords.size());
  for (size_t i = 0; i < this->tarc_.size(); ++i)
    for (size_t j = 0; j < dimworld; ++j)
      this->tarc_[i][j] = target_coords[i][j];

  // psurface doesn't actually support the case dim==dimworld.  Therefore we
  // use a trick: we just embed everything in a space of dimension dimworld+1
  if (dim==dimworld) {
    for (size_t i = 0; i < this->domc_.size(); ++i)
      domc_[i][dim] = 0;

    for (size_t i = 0; i < this->domc_.size(); ++i)
      tarc_[i][dim] = 1;
  }

  std::cout << "PSurfaceMerge building merged grid..." << std::endl;

  if (dim==dimworld) {

    ConstantDirection<+1> positiveDirection;
    ConstantDirection<-1> negativeDirection;

    // compute the merged grid using the psurface library
    this->cm_.build(domc_, domi_,tarc_, tari_,
                    &positiveDirection, &negativeDirection);
  } else {

    // compute the merged grid using the psurface library
    this->cm_.build(domc_, domi_,tarc_, tari_,
                    domainDirections_, targetDirections_);

  }

  std::cout << "Finished building merged grid!" << std::endl;

  // get the representation from the contact mapping object
  std::vector<IntersectionPrimitive<dim,ctype> > overlaps;
  this->cm_.getOverlaps(overlaps);

  // /////////////////////////////////////////////////////////////////////////////
  //  If overlaps refer to triangular elements that have been created
  //  by splitting a quadrilateral we have to reorder the element references
  //  and the local coordinates.
  // /////////////////////////////////////////////////////////////////////////////

  for (size_t i=0; i<overlaps.size(); i++) {

    //
    if (domainIsSecondTriangle[overlaps[i].tris[0]][0]) {

      assert(dim==2);

      // loop over the intersection corners
      for (int j=0; j<dim+1; j++) {

        // kinda umstaendlich: go from barycentric to reference coords, do the transformation
        // there, then go back to barycentric coordinates.  Streamlining this is left for later.
        Dune::FieldVector<double,dim+1> barycentric;
        barycentric[0] = overlaps[i].localCoords[0][j][0];
        barycentric[1] = overlaps[i].localCoords[0][j][1];
        barycentric[2] = 1 - barycentric[0] - barycentric[1];

        Dune::FieldVector<double,dim> ref = barycentricToReference(barycentric);

        // this is the actual transformation.  (0,0) becomes (1,1), (1,0) becomes (0,1),  (0,1) becomes (1,0)
        ref *= -1;
        ref +=  1;

        // back to barycentric
        barycentric = referenceToBarycentric(ref);
        overlaps[i].localCoords[0][j][0] = barycentric[0];
        overlaps[i].localCoords[0][j][1] = barycentric[1];

      }

    }

    if (targetIsSecondTriangle[overlaps[i].tris[1]][0]) {

      assert(dim==2);

      // loop over the intersection corners
      for (size_t j=0; j<dim+1; j++) {

        // kinda umstaendlich: go from barycentric to reference coords, do the transformation
        // there, then go back to barycentric coordinates.  Streamlining this is left for later.
        Dune::FieldVector<double,dim+1> barycentric;
        barycentric[0] = overlaps[i].localCoords[1][j][0];
        barycentric[1] = overlaps[i].localCoords[1][j][1];
        barycentric[2] = 1 - barycentric[0] - barycentric[1];

        Dune::FieldVector<double,dim> ref = barycentricToReference(barycentric);

        // this is the actual transformation.  (0,0) becomes (1,1), (1,0) becomes (0,1),  (0,1) becomes (1,0)
        ref *= -1;
        ref +=  1;

        // back to barycentric
        barycentric = referenceToBarycentric(ref);
        overlaps[i].localCoords[1][j][0] = barycentric[0];
        overlaps[i].localCoords[1][j][1] = barycentric[1];

      }


    }

    // The numbers in overlaps[].tri refer to the split elements.
    // Replace that with the actual numbers
    overlaps[i].tris[0] = unsplitDomainElementNumbers[overlaps[i].tris[0]];
    overlaps[i].tris[1] = unsplitTargetElementNumbers[overlaps[i].tris[1]];

  }

  // //////////////////////////////////////////////
  // initialize the merged grid overlap manager
  // //////////////////////////////////////////////

  this->olm_.setOverlaps(overlaps);

}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::nSimplices() const
{
  return this->olm_.nOverlaps();
}


template<int dim, int dimworld, typename T>
inline bool PSurfaceMerge<dim, dimworld, T>::domainSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstDomainParent" is in the valid range
  return (this->olm_.firstDomainParent(idx) < this->olm_.nOverlaps());
}


template<int dim, int dimworld, typename T>
inline bool PSurfaceMerge<dim, dimworld, T>::targetSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstTargetParent" is in the valid range
  return (this->olm_.firstTargetParent(idx) < this->olm_.nOverlaps());
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::domainParent(unsigned int idx) const
{
  return this->olm_.domain(idx).tris[0];
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return this->olm_.domain(idx).tris[1];
}


template<int dim, int dimworld, typename T>
bool PSurfaceMerge<dim, dimworld, T>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  unsigned int first = this->olm_.firstDomainParent(idx);
  unsigned int count = 0;
  // if the result is an index outside the range abort with error
  if (first >= this->olm_.nOverlaps())
  {
    indices.resize(0);
    return false;
  }
  // in case of a valid index:
  // ...count the number of simplex overlaps with given domain parent
  while (first + count < this->olm_.nOverlaps() && static_cast<unsigned int>(this->olm_.domain(first + count).tris[0]) == idx)
    count++;
  // ...fill the vector with the indices of the overlaps
  indices.resize(count);
  for (unsigned int i = 0; i < count; ++i)
    indices[i] = first + i;
  // ...return success
  return true;
}


template<int dim, int dimworld, typename T>
bool PSurfaceMerge<dim, dimworld, T>::targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  unsigned int first = this->olm_.firstTargetParent(idx);
  unsigned int count = 0;
  // if the result is an index outside the range abort with error
  if (first >= this->olm_.nOverlaps())
  {
    indices.resize(0);
    return false;
  }
  // in case of a valid index:
  // ...count the number of simplex overlaps with given domain parent
  while (first + count < this->olm_.nOverlaps() && static_cast<unsigned int>(this->olm_.target(first + count).tris[1]) == idx)
    count++;
  // ...fill the vector with the indices of the overlaps
  indices.resize(count);
  for (unsigned int i = 0; i < count; ++i)
    indices[i] = this->olm_.domainIndex(first + i);     // return the CORRECT INDEX here!!!
  // ...return success
  return true;
}


template<int dim, int dimworld, typename T>
typename PSurfaceMerge<dim, dimworld, T>::LocalCoords PSurfaceMerge<dim, dimworld, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const IntersectionPrimitive<dim,ctype>& ip = this->olm_.domain(idx);

  // read the local coordinates from the overlap's struct,
  LocalCoords result;

  if (dim==1) {
    result[0] = ip.localCoords[0][corner][0];
  } else {

    // The psurface local coordinates are strange: The coordinates x[0], x[1] are the first two
    // components of the three barycentric coordinates.  The third one is implied to be 1-x[0]-x[1].
    // Maybe I should change that some day.  OS
    assert(dim==2);
    Dune::FieldVector<double,dim+1> barycentric;
    barycentric[0] = ip.localCoords[0][corner][0];
    barycentric[1] = ip.localCoords[0][corner][1];
    barycentric[2] = 1 - barycentric[0] - barycentric[1];

    result = barycentricToReference(barycentric);
  }

  return result;
}


template<int dim, int dimworld, typename T>
typename PSurfaceMerge<dim, dimworld, T>::LocalCoords PSurfaceMerge<dim, dimworld, T>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const IntersectionPrimitive<dim,ctype>& ip = this->olm_.domain(idx);

  // read the local coordinates from the overlap's struct,
  LocalCoords result;
  if (dim==1) {
    result[0] = ip.localCoords[1][corner][0];
  } else {

    // The psurface local coordinates are strange: The coordinates x[0], x[1] are the first two
    // components of the three barycentric coordinates.  The third one is implied to be 1-x[0]-x[1].
    // Maybe I should change that some day.  OS
    assert(dim==2);
    Dune::FieldVector<double,dim+1> barycentric;
    barycentric[0] = ip.localCoords[1][corner][0];
    barycentric[1] = ip.localCoords[1][corner][1];
    barycentric[2] = 1 - barycentric[0] - barycentric[1];

    result = barycentricToReference(barycentric);
  }
  return result;
}


/* IMPLEMENTATION OF   O V E R L A P  M A N A G E R  SUBCLASS */

template<int dim, int dimworld, typename T>
void PSurfaceMerge<dim, dimworld, T>::OverlapManager::setOverlaps(const std::vector<IntersectionPrimitive<dim,ctype> >& unordered)
{
  this->domOrder.resize(unordered.size());
  this->tarOrder.resize(unordered.size(), NULL);

  // initialize the pointer arrays
  for (unsigned int i = 0; i < unordered.size(); ++i)
  {
    this->domOrder[i] = unordered[i];
  }

  // now sort to ascending domain parent index order
  std::sort(this->domOrder.begin(), this->domOrder.end(), &domainParentSmaller);
  // fill the 2nd array with the pointers to the elements from the 1st array
  for (unsigned int i = 0; i < this->domOrder.size(); ++i)
    this->tarOrder[i] = &this->domOrder[i];
  // now sort only the pointers in O(N*logN) comparisons
  std::sort(this->tarOrder.begin(), this->tarOrder.end(), &targetParentSmaller);

  // get the base pointer of the array
  this->baseptr = &this->domOrder[0];

}


template<int dim, int dimworld, typename T>
unsigned int PSurfaceMerge<dim, dimworld, T>::OverlapManager::firstDomainParent(unsigned int parent) const
{
  // perform a binary search for a simplex overlap with the given parent
  unsigned int first = 0, last = this->domOrder.size(), p = last/2;
  bool continuing = true;
  // search until either found or until not found with the search interval's size decreased to 1
  while (((p = this->domain((first + last)/2).tris[0]) != parent) && (continuing = (last > first + 1)))
  {
    if (p > parent)
      last = (first + last) / 2;
    else
      first = (first + last) / 2;
  }

  // if the parent is not found return a value too high
  if (!continuing)
    return this->domOrder.size();

  // else go back if the found element is not the first one
  // (note that p is now the INDEX of the overlap)
  p = (first + last) / 2;
  while (p > 0 && static_cast<unsigned int>(this->domain(p-1).tris[0]) == parent)
    p--;
  return p;
}


template<int dim, int dimworld, typename T>
unsigned int PSurfaceMerge<dim, dimworld, T>::OverlapManager::firstTargetParent(unsigned int parent) const
{
  // perform a binary search for a simplex overlap with the given parent
  unsigned int first = 0, last = this->domOrder.size(), p = last/2;
  bool continuing = true;
  // search until either found or until not found with the search interval's size decreased to 1
  while (((p = this->target((first + last)/2).tris[1]) != parent) && (continuing = (last > first + 1)))
  {
    if (p > parent)
      last = (first + last) / 2;
    else
      first = (first + last) / 2;
  }

  // if the parent is not found return a value too high
  if (!continuing)
    return this->domOrder.size();

  // else go back if the found element is not the first one
  // (note that p is now the INDEX of the overlap)
  p = (first + last) / 2;
  while (p > 0 && static_cast<unsigned int>(this->target(p-1).tris[1]) == parent)
    p--;
  return p;
}


#endif // PSURFACEMERGE_HH_
