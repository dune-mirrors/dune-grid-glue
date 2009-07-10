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
 * @file PSurfaceMerge.hh
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

#ifndef PSURFACEMERGE_HH_
#define PSURFACEMERGE_HH_


// This flag can be set to avoid illegal barycentric coordinates.
// It may happen, that the sum of the first dim-1 coordinates
// exceeds 1.0 which results in the last coordinate being
// negative. It is set to 0.0 then instead.
#define CHECK_BARYCENTRIC_COORDINATES

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include "../misc/geometry.hh"
#include "../misc/multidimoctree.hh"
#include <psurface/ContactMapping.h>


/** \brief Standard implementation of the SurfaceMerge concept using the psurface library.

   \tparam dim Grid dimension
   \tparam T Type used for coordinates
 */
template<int dim, typename T = double>
class PSurfaceMerge
{
public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the dimension via compile time constant
  enum { dimw = dim };

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  Coords;

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



  /**
   * @class OverlapCorner
   * @brief identifies one corner point in a particular overlap simplex
   *
   * Since the number of corners is so small one "shared" int is enough.
   */
  struct OverlapCorner
  {
    unsigned int idx : 24;
    unsigned int corner : 8;
  };


  /**
   * @class OverlapCornerGeometry
   * @brief for an OverlapCorner struct provides the geometric information needed in the point locator
   */
  struct OverlapCornerGeometry
  {
    OverlapCornerGeometry(const OverlapManager& olm_) : olm(olm_)
    {}

    /// @brief mandatory for the point locator
    bool operator()(const Coords& lower, const Coords& upper, const OverlapCorner& item) const
    {
      const IntersectionPrimitive<float>& ip = this->olm.domain(item.idx);
      for (typename Coords::size_type i = 0; i < dim; ++i)
        if (lower[i] > ip.points[item.corner][i] || upper[i] <= ip.points[item.corner][i])
          return false;
      return true;
    }

    /// @brief the source of the geometric information
    const OverlapManager& olm;
  };

  /// @brief class used for domain and target vertex locating in the merged grid
  typedef MultiDimOctree<OverlapCorner, OverlapCornerGeometry, Coords, dim, true>  Locator;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  T _maxdist;

  /// @brief tolerance value for geometric comparisons with float and double precision
  T _eps;

  /* geometric data for both domain and targt */

  /// @brief domain coordinates
  std::vector<std::tr1::array<double,dim> >   _domc;

  /// @brief target coordinates
  std::vector<std::tr1::array<double,dim> >   _tarc;


  /* topologic information for domain and target */

  /// @ brief domain indices (internal copy)
  std::vector<std::tr1::array<int,dim> >         _domi;

  /// @brief target indices (internal copy)
  std::vector<std::tr1::array<int,dim> >         _tari;


  /* members associated with the merged grid */

  /// @brief make use of psurface contact mapping functionality
  ContactMapping<dim> _cm;

  /// @brief provides an interface for convenient and efficient
  /// access to the set of computed simplex overlaps in the merged grid
  OverlapManager _olm;


  /* members needed for the identification of domain and target vertices in the merged grid */

  /// @brief spatial map like an octree or quadtree for efficient search of point associated data
  mutable Locator _olclocator;

  /// @brief the data repository for the entities in the locator (no direct access to this)
  std::vector<OverlapCorner>  _olcorners;

  /// @brief the geometry accessing functor the entities in the locator
  OverlapCornerGeometry* _olcgeometry;

  SurfaceNormal _domnormals;


protected:

  /**
   * @brief get corners (orientation like domain) of simplex at given index
   * @param idx index of a merged grid simplex, user must ensure that index is valid
   * or else the result will be invalid
   * @return const ref to an array of ordered coordinates denoting the corners and the orientation
   * of the simplex IFF given index was valid
   */
  Dune::FieldVector<Coords, dim> simplexCorners(unsigned int idx) const;


public:

  PSurfaceMerge(T max_distance = 1E-4, const SurfaceNormal domain_normals = NULL) :
    _maxdist(max_distance), _eps(1E-10),
    _olclocator(typename Locator::BoxType(Coords(0.0), Coords(1.0))), _olcgeometry(NULL), _domnormals(domain_normals)
  {}


  ~PSurfaceMerge()
  {
    if (this->_olcgeometry != NULL)
      delete this->_olcgeometry;
  }


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
   * @param epsilon the estimate maximum deformation for the contact oracle
   * @param obsDirections If given this function is used to compute the normals for the vertices of domain and target.
   * Otherwise the default algorithm is used to compute them.
   * @return TRUE <=> build successful and merged grid not empty
   */
  bool build(
    const std::vector<T>& domain_coords,
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
   * @return barycentric coordinates in parent domain simplex
   */
  Coords domainParentLocal(unsigned int idx, unsigned int corner) const;

  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return barycentric coordinates in parent target simplex
   */
  Coords targetParentLocal(unsigned int idx, unsigned int corner) const;

};


/* IMPLEMENTATION OF CLASS   C O N T A C T  M A P P I N G  S U R F A C E  M E R G E */


template<int dim, typename T>
bool PSurfaceMerge<dim, T>::build(
  const std::vector<T>& domain_coords,
  const std::vector<unsigned int>& domain_simplices,
  const std::vector<T>& target_coords,
  const std::vector<unsigned int>& target_simplices
  )
{
  try
  {
    // copy domain and target simplices to internal arrays
    // (cannot keep refs since outside modification would destroy information)
    this->_domi.resize(domain_simplices.size()/dim);
    for (unsigned int i = 0; i < domain_simplices.size()/dim; ++i)
      for (int j=0; j<dim; j++)
        this->_domi[i][j] = domain_simplices[i*dim+j];

    this->_tari.resize(target_simplices.size()/dim);
    for (unsigned int i = 0; i < target_simplices.size()/dim; ++i)
      for (int j=0; j<dim; j++)
        this->_tari[i][j] = target_simplices[i*dim+j];

    // copy the coordinates to internal arrays of coordinates
    // (again cannot just keep refs and this representation has advantages)
    this->_domc.resize(domain_coords.size() / dim);
    for (unsigned int i = 0; i < this->_domc.size(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        this->_domc[i][j] = domain_coords[i*dim + j];

    this->_tarc.resize(target_coords.size() / dim);
    for (unsigned int i = 0; i < this->_tarc.size(); ++i)
      for (unsigned int j = 0; j < dim; ++j)
        this->_tarc[i][j] = target_coords[i*dim + j];

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

    // introduce a quadtree or octree for fast identification of domain and target
    // vertices in the merged grid

    // compute the total numbers of overlap simplex corners
    unsigned int num_olcorners = dim * this->_olm.nOverlaps();

    // re-initialize the point locator
    Coords lower(std::numeric_limits<ctype>::max()), upper(-std::numeric_limits<ctype>::max());
    // compute the bounding box of the merged grid
    typename Dune::FieldVector<Coords, dim>::size_type j;
    typename Coords::size_type k;
    for (unsigned int i = 0; i < this->_olm.nOverlaps(); ++i)
    {
      // can work with a reference here since no allocation is done in the next two loops
      // (should save some calls to constructors)
      const Dune::FieldVector<Coords, dim>& corners = this->simplexCorners(i);
      for (j = 0; j < dim; ++j)
        for (k = 0; k < dim; ++k)
        {
          lower[k] = std::min(lower[k], corners[j][k]);
          upper[k] = std::max(upper[k], corners[j][k]);
        }
    }

    // increase a little, since "upper" boundary not inside locator domain
    /** \todo This should be scaling invariant! */
    lower -= 0.1;
    upper += 0.1;

    // get an estimate on how many levels will be needed in the octree,
    // the default value for max_items/cell is 10, here it will be 8
    int max_items_per_cell = 8;
    int levels = 1;
    while ((num_olcorners / pow((float) (1 << dim), levels)) > max_items_per_cell)
      levels++;
    // ...and to be on the safe side
    levels += 4;
    this->_olclocator.init(typename Locator::BoxType(lower, upper), levels, max_items_per_cell);

    // now fill the locator...

    // first allocate the data "repository" for the locator's content is be stored outside
    this->_olcorners.resize(num_olcorners);

    // create geometry "oracle" for the overlap corner structs in the point locator
    this->_olcgeometry = new OverlapCornerGeometry(this->_olm);

    // insert the data
    unsigned int current_index = 0;
    for (unsigned int i = 0; i < this->_olm.nOverlaps(); ++i)
    {
      for (int corner = 0; corner < dim; ++corner)
      {
        this->_olcorners[current_index].idx = i;
        this->_olcorners[current_index].corner = corner;
        this->_olclocator.insert(&this->_olcorners[current_index], this->_olcgeometry);
        current_index++;
      }
    }

    return true;
  }
  catch (...)
  {
    std::cerr << "PSurfaceMerge: Unknown exception occurred!" << std::endl;
  }
  // only arriving here after exception
  return false;
}


template<int dim, typename T>
inline unsigned int PSurfaceMerge<dim, T>::nSimplices() const
{
  return this->_olm.nOverlaps();
}


template<int dim, typename T>
inline bool PSurfaceMerge<dim, T>::domainSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstDomainParent" is in the valid range
  return (this->_olm.firstDomainParent(idx) < this->_olm.nOverlaps());
}


template<int dim, typename T>
inline bool PSurfaceMerge<dim, T>::targetSimplexMatched(unsigned int idx) const
{
  // if the simplex was matched the result returned by "firstTargetParent" is in the valid range
  return (this->_olm.firstTargetParent(idx) < this->_olm.nOverlaps());
}


//template<int dim, typename T>
//bool PSurfaceMerge<dim, T>::domainVertexMatched(unsigned int idx) const
//{
//	// TODO think about implementation similar to targetVertexMatched below...
//	Coords lower(-this->_eps), upper(this->_eps);
//	lower += this->_domc[idx];
//	upper += this->_domc[idx];
//	typename Locator::ResultContainer res;
//	return (this->_olclocator.lookup(typename Locator::BoxType(lower, upper), res) > 0);
//}
//
//
//template<int dim, typename T>
//bool PSurfaceMerge<dim, T>::targetVertexMatched(unsigned int idx) const
//{
//	Coords lower(-this->_eps - this->_maxdist), upper(this->_eps + this->_maxdist);
//	lower += this->_tarc[idx];
//	upper += this->_tarc[idx];
//	typename Locator::ResultContainer res;
//	unsigned int count = this->_olclocator.lookup(typename Locator::BoxType(lower, upper), res);
//
//	// now check in the results if there is a target vertex there
//	for (unsigned int i = 0; i < count; ++i)
//	{
//		const IntersectionPrimitive<float>& ip = this->_olm.domain(res[i]->idx);
//		for (int j = 0; j < dim-1; ++j)
//		{
//			// for a vertex only 0 and 1 are possible bar. coord values
//			if (ip.localCoords[1][res[i]->corner][j] != 1.0 &&
//					ip.localCoords[1][res[i]->corner][j] != 0.0)
//				continue;
//		}
//		// if this point is reached, a valid vertex is found
//		return true;
//	}
//	// being here: no success!
//	return false;
//}


template<int dim, typename T>
inline unsigned int PSurfaceMerge<dim, T>::domainParent(unsigned int idx) const
{
  return this->_olm.domain(idx).tris[0];
}


template<int dim, typename T>
inline unsigned int PSurfaceMerge<dim, T>::targetParent(unsigned int idx) const
{
  // Warning: Be careful to use the ACTUAL indexing here defined in the array sorted after domain parent indices!!
  return this->_olm.domain(idx).tris[1];
}


template<int dim, typename T>
bool PSurfaceMerge<dim, T>::domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
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


template<int dim, typename T>
bool PSurfaceMerge<dim, T>::targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
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


template<int dim, typename T>
Dune::FieldVector<typename PSurfaceMerge<dim, T>::Coords, dim> PSurfaceMerge<dim, T>::simplexCorners(unsigned int idx) const
{
  // get the simplex overlap
  const IntersectionPrimitive<float>& ip = this->_olm.domain(idx);
  // copy the relevant fields from the points stored with the simplex overlap
  Dune::FieldVector<Coords, dim> result;
  for (int i = 0; i < dim; ++i)       // #simplex_corners equals dim here
    for (int j = 0; j < dim; ++j)             // #dimensions of coordinate space
      result[i][j] = ip.points[i][j];
  return result;
}


template<int dim, typename T>
typename PSurfaceMerge<dim, T>::Coords PSurfaceMerge<dim, T>::domainParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const IntersectionPrimitive<float>& ip = this->_olm.domain(idx);
  // read the local coordinates from the overlap's struct,
  // but note that the last coordinate is only stored implicitly
  // (sum must be 1.0, i.e. last is lin. dep.)

  Coords result(1.0);

  if (dim == 2)
  {
    // ContactMapping::getOverlaps fills the IntersectionPrimitive<float> data objects with
    // local coordinates in the domain parent triangle.
    // So in order to have a barycentric representation of local corners we just
    // reverse the order of the local coordinates' components.
    result[0] = 1.0 - ip.localCoords[0][corner][0];
    result[1] = ip.localCoords[0][corner][0];
  }
  else
  {
    // ContactMapping::getOverlaps fills the IntersectionPrimitive<float> data objects with
    // POSITIVELY oriented barycentric coordinates in the domain parent.
    // Since this is what we want we just stick to this order of corners.
    for (int j = 0; j < dim-1; ++j)             // #dimensions of coordinate space
    {
      result[j] = ip.localCoords[0][corner][j];
      // assemble the last coordinate
      result[dim-1] -= result[j];
    }
  }
#ifdef CHECK_BARYCENTRIC_COORDINATES
  // sometimes the numerics are a bit off...
  if (result[dim-1] < 0.0)
    result[dim-1] = 0.0;
#endif

  return result;
}


template<int dim, typename T>
typename PSurfaceMerge<dim, T>::Coords PSurfaceMerge<dim, T>::targetParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const IntersectionPrimitive<float>& ip = this->_olm.domain(idx);

  // read the local coordinates from the overlap's struct,
  // but note that the last coordinate is only stored implicitly
  // (sum must be 1.0, i.e. last is lin. dep.)
  Coords result(1.0);

  if (dim == 2)
  {
    // ContactMapping::getOverlaps fills the IntersectionPrimitive<float> data objects with
    // local coordinates in the target parent face.
    // So in order to have a barycentric representation of local corners we just
    // reverse the order of the local coordinates' components.
    result[0] = 1.0 - ip.localCoords[1][corner][0];
    result[1] = ip.localCoords[1][corner][0];
  }
  else       // dim == 3
  {
    // ContactMapping::getOverlaps fills the IntersectionPrimitive<float> data objects with
    // NEGATIVELY oriented barycentric coordinates in the target parent.
    for (int j = 0; j < dim-1; ++j)             // #dimensions of coordinate space
    {
      result[j] = ip.localCoords[1][corner][j];
      // assemble the last coordinate
      result[dim-1] -= result[j];
    }
  }
#ifdef CHECK_BARYCENTRIC_COORDINATES
  // sometimes the numerics are a bit off...
  if (result[dim-1] < 0.0)
    result[dim-1] = 0.0;
#endif
  return result;
}


/* IMPLEMENTATION OF   O V E R L A P  M A N A G E R  SUBCLASS */

template<int dim, typename T>
void PSurfaceMerge<dim, T>::OverlapManager::setOverlaps(const std::vector<IntersectionPrimitive<float> >& unordered)
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

  //	// This does not really work either!
  //	// If CHECK_BARYCENTRIC_COORDINATES is set this should ensure that the sum of the first n-1 barycentric
  //	// coordinates does not exceed 1.0 - but it does not work!!
  //	// Alternatively the result of the {domain,target}ParentLocal methods can be modified (which does not work either...)
  //	ofstream fosd, fost;
  //	fosd.open("/tmp/domain_overlaps.dat");
  //	fost.open("/tmp/target_overlaps.dat");
  //	for (unsigned int i = 0; i < this->domOrder.size(); ++i)
  //	{
  //		for (int k = 0; k < dim; ++k)
  //		{
  //			ctype domsum = 1.0, tarsum = 1.0;
  //			fosd << "l(";
  //			fost << "l(";
  //			for (int j = 0; j < dim-1; ++j)
  //			{
  //				fosd << this->domOrder[i].localCoords[0][k][j] << "  ";
  //				domsum -= this->domOrder[i].localCoords[0][k][j];
  //				fost << this->domOrder[i].localCoords[1][k][j] << "  ";
  //				tarsum -= this->domOrder[i].localCoords[1][k][j];
  //			}
  //			fosd << domsum << ")      g(";
  //			fost << tarsum << ")      g(";
  //			for (int j = 0; j < dim; ++j)
  //			{
  //				fosd << "  " << this->domOrder[i].points[k][j];
  //				fost << "  " << this->domOrder[i].points[k][j];
  //			}
  //			fosd << ")\n";
  //			fost << ")\n";
  //		}
  //		fosd << std::endl;
  //		fost << std::endl;
  //	}
  //	fosd.close();
  //	fost.close();

}


template<int dim, typename T>
unsigned int PSurfaceMerge<dim, T>::OverlapManager::firstDomainParent(unsigned int parent) const
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


template<int dim, typename T>
unsigned int PSurfaceMerge<dim, T>::OverlapManager::firstTargetParent(unsigned int parent) const
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
