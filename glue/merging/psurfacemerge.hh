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

#ifndef PSURFACEMERGE_HH
#define PSURFACEMERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/glue/misc/geometry.hh>
#include <dune/glue/merging/merger.hh>

#include <psurface/ContactMapping.h>


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
  // around which the additional dimension
  enum {psurfaceDimworld = dimworld + (dim==dimworld)};

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

  /// @brief a function (or functor) describing the domain normals
  typedef void (*SurfaceNormal)(const ctype* pos, ctype* dir);


private:

  /*   P R I V A T E   H E L P E R S   */

  /**
   * @brief compare function dep. on parent domain simplex
   * @param a first item
   * @param b second item
   * @return TRUE <=> parent domain simplex' index of @c a smaller
   */
  static bool domainParentSmaller(const IntersectionPrimitive<float>& a, const IntersectionPrimitive<float>& b)
  {
    return a.tris[0] < b.tris[0];
  }

  /**
   * @brief compare function dep. on parent target simplex
   * @param a first item
   * @param b second item
   * @return TRUE <=> parent target simplex' index of @c a smaller
   */
  static bool targetParentSmaller(const IntersectionPrimitive<float>* a, const IntersectionPrimitive<float>* b)
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
    void setOverlaps(const std::vector<IntersectionPrimitive<float> >& unordered);

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
    const IntersectionPrimitive<float>& domain(unsigned int idx) const
    {
      return this->domOrder[idx];
    }

    /**
     * @brief getter for an simplex in the array sorted after the target parent indices
     * @param idx index into the array
     * @return the simplex
     */
    const IntersectionPrimitive<float>& target(unsigned int idx) const
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
    /// ascending domain parent simplex indices (see _domi)
    std::vector<IntersectionPrimitive<float> > domOrder;

    /// @brief the computed merged grid simplices sorted after
    /// ascending target parent simplex indices (see _tari)
    std::vector<IntersectionPrimitive<float>*> tarOrder;

    /// @brief used to recompute original index of an overlap in "tarOrder"
    IntersectionPrimitive<float>*         baseptr;
  };


private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  T _maxdist;

  /* geometric data for both domain and targt */

  /// @brief domain coordinates
  std::vector<std::tr1::array<double,psurfaceDimworld> >   _domc;

  /// @brief target coordinates
  std::vector<std::tr1::array<double,psurfaceDimworld> >   _tarc;


  /* topologic information for domain and target */

  /// @ brief domain indices (internal copy)
  std::vector<std::tr1::array<int,dim+1> >         _domi;

  /// @brief target indices (internal copy)
  std::vector<std::tr1::array<int,dim+1> >         _tari;


  /* members associated with the merged grid */

  /// @brief make use of psurface contact mapping functionality
  ContactMapping<dim+1> _cm;

  /// @brief provides an interface for convenient and efficient
  /// access to the set of computed simplex overlaps in the merged grid
  OverlapManager _olm;

  /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
   */
  SurfaceNormal _domnormals;


public:

  PSurfaceMerge(T max_distance = 1E-4, const SurfaceNormal domain_normals = NULL) :
    _maxdist(max_distance), _domnormals(domain_normals)
  {}


  /*   M O D E L   S P E C I F I C   E X T E N D I N G   F U N C T I O N A L I T Y   */

  /**
   * @brief setter for the maximum distance between domain and target surface
   *
   * If the distance between a point on the target surface and a point on the domain
   * surface exceeds this value (measured along the search direction which is
   * usually the normal) the points will not be matched.
   * @param value the new value ( value > 0.0 )
   */
  void setMaxDistance(ctype value)
  {
    if (value > 0.0)
      this->_maxdist = value;
  }


  /**
   * @brief setter for the domain surface normal function
   *
   * The matching of the geometries offers the possibility to specify a function for
   * the exact evaluation of domain surface normals. If no such function is specified
   * (default) normals are interpolated.
   * @param value the new function (or NULL to unset the function)
   */
  void setDomainNormals(const SurfaceNormal value)
  {
    this->_domnormals = value;
  }


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
             const std::vector<unsigned int>& domain_simplices,
             const std::vector<T>& target_coords,
             const std::vector<unsigned int>& target_simplices
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
void PSurfaceMerge<dim, dimworld, T>::build(
  const std::vector<T>& domain_coords,
  const std::vector<unsigned int>& domain_simplices,
  const std::vector<T>& target_coords,
  const std::vector<unsigned int>& target_simplices
  )
{
  // copy domain and target simplices to internal arrays
  // (cannot keep refs since outside modification would destroy information)
  this->_domi.resize(domain_simplices.size()/(dim+1));
  this->_tari.resize(target_simplices.size()/(dim+1));

  for (unsigned int i = 0; i < domain_simplices.size()/(dim+1); ++i)
    for (int j=0; j<dim+1; j++)
      this->_domi[i][j] = domain_simplices[i*(dim+1)+j];

  if (dim!=dimworld) {

    // dim==dimworld-1: just copy the two arrays
    for (unsigned int i = 0; i < target_simplices.size()/(dim+1); ++i)
      for (int j=0; j<dim+1; j++)
        this->_tari[i][j] = target_simplices[i*(dim+1)+j];

  } else {

    // dim==dimworld: "hen grids are artificially embedded in a dim+1 space,
    // the second grid needs to have its orientation reversed.
    // That way, the 'surface normals' of the two grids point towards each other.
    for (unsigned int i = 0; i < target_simplices.size()/(dim+1); ++i) {
      _tari[i][0] = target_simplices[i*(dim+1)+1];
      _tari[i][1] = target_simplices[i*(dim+1)+0];
      if (dim==2)
        _tari[i][2] = target_simplices[i*(dim+1)+2];
    }

  }
  // copy the coordinates to internal arrays of coordinates
  // (again cannot just keep refs and this representation has advantages)
  this->_domc.resize(domain_coords.size() / dimworld);
  for (unsigned int i = 0; i < this->_domc.size(); ++i)
    for (unsigned int j = 0; j < dimworld; ++j)
      this->_domc[i][j] = domain_coords[i*dimworld + j];

  this->_tarc.resize(target_coords.size() / dimworld);
  for (unsigned int i = 0; i < this->_tarc.size(); ++i)
    for (unsigned int j = 0; j < dimworld; ++j)
      this->_tarc[i][j] = target_coords[i*dimworld + j];

  // psurface doesn't actually support the case dim==dimworld.  Therefore we
  // use a trick: we just embed everything in a space of dimension dimworld+1
  if (dim==dimworld) {
    for (unsigned int i = 0; i < this->_domc.size(); ++i)
      _domc[i][dim] = 0;

    for (unsigned int i = 0; i < this->_domc.size(); ++i)
      _tarc[i][dim] = 1;

    // Needs to be more than 1
    _maxdist = 2;
  }

  std::cout << "Building merged grid... (wait for finish!)" << std::endl;

  // compute the merged grid using the psurface library
  this->_cm.build(_domc, _domi,
                  _tarc, _tari,
                  this->_maxdist, this->_domnormals);

  std::cout << "Finished building merged grid!" << std::endl;

  // get the representation from the contact mapping object
  std::vector<IntersectionPrimitive<float> > overlaps;
  this->_cm.getOverlaps(overlaps);

  // initialize the merged grid overlap manager
  this->_olm.setOverlaps(overlaps);

}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::nSimplices() const
{
  return this->_olm.nOverlaps();
}


template<int dim, int dimworld, typename T>
inline bool PSurfaceMerge<dim, dimworld, T>::domainSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstDomainParent" is in the valid range
  return (this->_olm.firstDomainParent(idx) < this->_olm.nOverlaps());
}


template<int dim, int dimworld, typename T>
inline bool PSurfaceMerge<dim, dimworld, T>::targetSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstTargetParent" is in the valid range
  return (this->_olm.firstTargetParent(idx) < this->_olm.nOverlaps());
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::domainParent(unsigned int idx) const
{
  return this->_olm.domain(idx).tris[0];
}


template<int dim, int dimworld, typename T>
inline unsigned int PSurfaceMerge<dim, dimworld, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return this->_olm.domain(idx).tris[1];
}


template<int dim, int dimworld, typename T>
bool PSurfaceMerge<dim, dimworld, T>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
{
  unsigned int first = this->_olm.firstDomainParent(idx);
  unsigned int count = 0;
  // if the result is an index outside the range abort with error
  if (first >= this->_olm.nOverlaps())
  {
    indices.resize(0);
    return false;
  }
  // in case of a valid index:
  // ...count the number of simplex overlaps with given domain parent
  while (first + count < this->_olm.nOverlaps() && static_cast<unsigned int>(this->_olm.domain(first + count).tris[0]) == idx)
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
  unsigned int first = this->_olm.firstTargetParent(idx);
  unsigned int count = 0;
  // if the result is an index outside the range abort with error
  if (first >= this->_olm.nOverlaps())
  {
    indices.resize(0);
    return false;
  }
  // in case of a valid index:
  // ...count the number of simplex overlaps with given domain parent
  while (first + count < this->_olm.nOverlaps() && static_cast<unsigned int>(this->_olm.target(first + count).tris[1]) == idx)
    count++;
  // ...fill the vector with the indices of the overlaps
  indices.resize(count);
  for (unsigned int i = 0; i < count; ++i)
    indices[i] = this->_olm.domainIndex(first + i);             // return the CORRECT INDEX here!!!
  // ...return success
  return true;
}


template<int dim, int dimworld, typename T>
typename PSurfaceMerge<dim, dimworld, T>::LocalCoords PSurfaceMerge<dim, dimworld, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const IntersectionPrimitive<float>& ip = this->_olm.domain(idx);

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
  const IntersectionPrimitive<float>& ip = this->_olm.domain(idx);

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
void PSurfaceMerge<dim, dimworld, T>::OverlapManager::setOverlaps(const std::vector<IntersectionPrimitive<float> >& unordered)
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
