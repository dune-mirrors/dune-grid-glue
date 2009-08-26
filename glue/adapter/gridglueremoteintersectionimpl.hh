// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueRemoteIntersectionImpl.hh
 *  Version:     1.0
 *  Created on:  Mar 2, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Model of the RemoteIntersection concept provided by GridGlue.
 *  subversion:  $Id$
 *
 */
/**
 * @file GridGlueRemoteIntersectionImpl.hh
 * @brief Model of the RemoteIntersection concept provided by GridGlue.
 */

#ifndef GRIDGLUEREMOTEINTERSECTIONIMPL_HH_
#define GRIDGLUEREMOTEINTERSECTIONIMPL_HH_


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   REMOTE INTERSECTION IMPL   */


template<typename GET1, typename GET2>
class GridGlue<GET1, GET2>::RemoteIntersectionImpl
{

private:

  typedef GridGlue<GET1, GET2> Parent;

  friend class Parent::RemoteIntersectionIteratorImpl;

  friend class Parent::DomainIntersectionIteratorImpl;

  friend class Parent::TargetIntersectionIteratorImpl;


public:


  enum { coorddim = Parent::dimw };

  enum { mydim = coorddim - 1 };

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename DomainGridView::Grid DomainGridType;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename TargetGridView::Grid TargetGridType;


  typedef typename DomainGridView::Traits::template Codim<0>::Entity DomainEntity;

  typedef typename DomainGridView::Traits::template Codim<0>::EntityPointer DomainEntityPointer;


  typedef typename TargetGridView::Traits::template Codim<0>::Entity TargetEntity;

  typedef typename TargetGridView::Traits::template Codim<0>::EntityPointer TargetEntityPointer;


  typedef LocalSimplexGeometry<mydim, DomainGridView::dimension, DomainGridType> DomainLocalGeometry;

  typedef SimplexGeometry<mydim, DomainGridType::dimensionworld, DomainGridType> DomainGeometry;

  typedef LocalSimplexGeometry<mydim, TargetGridView::dimension, TargetGridType> TargetLocalGeometry;

  typedef SimplexGeometry<mydim, TargetGridType::dimensionworld, TargetGridType> TargetGeometry;


  typedef typename DomainGridType::ctype ctype;


  typedef Dune::FieldVector<ctype, mydim>      LocalCoords;

  typedef Dune::FieldVector<ctype, coorddim>   Coords;

  typedef unsigned int IndexType;

private:

  /*   M E M B E R   V A R  I A B L E S   */

  /// @brief the grid glue entity this is built on
  const Parent*       _glue;

  /// @brief index of this intersection after GridGlue interface
  IndexType _mergeindex;

  IndexType _index;
  IndexType _domainindex;
  IndexType _targetindex;

  DomainLocalGeometry _domlgeom;

  DomainGeometry _domggeom;

  TargetLocalGeometry _tarlgeom;

  TargetGeometry _targgeom;

public:

  static const Dune::GeometryType geometrytype;

  /*   F U N C T I O N A L I T Y   */

  RemoteIntersectionImpl(const Parent* glue)
    : _glue(glue), _mergeindex((IndexType)-1), _index((IndexType)-1),
      _domainindex((IndexType)-1), _targetindex((IndexType)-1)
  {}

  RemoteIntersectionImpl(const Parent* glue, unsigned int mergeindex);

  /** \brief Copy constructor */
  RemoteIntersectionImpl(const RemoteIntersectionImpl & impl)
    : _glue(impl._glue), _mergeindex(impl._mergeindex),
      _index(impl._index), _domainindex(impl._domainindex),
      _targetindex(impl._targetindex),
      _domlgeom(impl._domlgeom),
      _domggeom(impl._domggeom),
      _tarlgeom(impl._tarlgeom),
      _targgeom(impl._targgeom)
  {}

  /** \brief Assignment operator */
  RemoteIntersectionImpl& operator=(const RemoteIntersectionImpl& impl)
  {
    this->_glue = impl._glue;
    this->_mergeindex = impl._mergeindex;
    this->_index = impl._index;
    this->_domainindex = impl._domainindex;
    this->_targetindex = impl._targetindex;
    this->_domlgeom = impl._domlgeom;
    this->_domggeom = impl._domggeom;
    this->_tarlgeom = impl._tarlgeom;
    this->_targgeom = impl._targgeom;
    return *this;
  }

  // return EntityPointer to the Entity on the inside of this intersection. That is the Entity where we started this .
  DomainEntityPointer entityDomain() const
  {
    return this->_glue->_domext.element(this->_glue->_merg->domainParent(this->_mergeindex));
  }


  // return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity.
  TargetEntityPointer entityTarget() const
  {
    return this->_glue->_tarext.element(this->_glue->_merg->targetParent(this->_mergeindex));
  }


  /** \brief Return true if intersection is conforming */
  bool conforming() const;


  // geometrical information about this intersection in local coordinates of the inside() entity.
  const DomainLocalGeometry& geometryInDomainEntity() const
  {
    return this->_domlgeom;
  }


  // geometrical information about this intersection in local coordinates of the outside() entity.
  const TargetLocalGeometry& geometryInTargetEntity() const
  {
    return this->_tarlgeom;
  }


  // geometrical information about this intersection in global coordinates in the domain grid.
  const DomainGeometry& geometryDomain() const
  {
    return this->_domggeom;
  }


  // geometrical information about this intersection in global coordinates in the target grid.
  const TargetGeometry& geometryTarget() const
  {
    return this->_targgeom;
  }


  bool hasTarget() const
  {
    unsigned int localindex;
    return this->_glue->_tarext.contains(this->_glue->_merg->targetParent(this->_mergeindex), localindex);
  }

  bool hasDomain() const
  {
    unsigned int localindex;
    return this->_glue->_domext.contains(this->_glue->_merg->domainParent(this->_mergeindex), localindex);
  }

  // obtain the type of reference element for this intersection
  Dune::GeometryType type() const
  {
    return geometrytype;
  }


  // Local number of codim 1 entity in the inside() Entity where intersection is contained in.
  int numberInDomainEntity() const
  {
    return this->_glue->_domext.indexInInside(this->_glue->_merg->domainParent(this->_mergeindex));
  }


  // Local number of codim 1 entity in outside() Entity where intersection is contained in.
  int numberInTargetEntity() const
  {
    return this->_glue->_tarext.indexInInside(this->_glue->_merg->targetParent(this->_mergeindex));
  }


  // Return an outer normal (length not necessarily 1).
  Coords outerNormalDomain(const Dune::FieldVector<ctype, mydim> &local) const
  {
    return this->_domggeom.outerNormal(local);
  }

  // Return an outer normal (length not necessarily 1).
  Coords outerNormalTarget(const Dune::FieldVector<ctype, mydim> &local) const
  {
    return this->_targgeom.outerNormal(local);
  }

#ifdef QUICKHACK_INDEX
  IndexType & index()
  {
    return this->_index;
  }

  IndexType & domainIndex()
  {
    return this->_domainindex;
  }

  IndexType & targetIndex()
  {
    return this->_targetindex;
  }

  IndexType index() const
  {
    return this->_index;
  }

  IndexType globalIndex() const
  {
    return this->_mergeindex;
  }

  IndexType domainIndex() const
  {
    assert(this->_domainindex != (IndexType)-1);
    return this->_domainindex;
  }

  IndexType targetIndex() const
  {
    assert(this->_targetindex != (IndexType)-1);
    return this->_targetindex;
  }

#endif
};


template<typename GET1, typename GET2>
const Dune::GeometryType GridGlue<GET1, GET2>::RemoteIntersectionImpl::geometrytype(Dune::GeometryType::simplex, Parent::dimw-1);


template<typename GET1, typename GET2>
bool GridGlue<GET1, GET2>::RemoteIntersectionImpl::conforming() const
{
  std::vector<unsigned int> results;
  // first check the domain side
  bool is_conforming =
    this->_glue->_merg->domainSimplexRefined(this->_glue->_merg->domainParent(this->_mergeindex), results) && results.size() == 1;
  results.resize(0);
  // now check the target side
  if (is_conforming)
    return this->_glue->_merg->targetSimplexRefined(this->_glue->_merg->targetParent(this->_mergeindex), results) && results.size() == 1;
  return false;
}


template<typename GET1, typename GET2>
GridGlue<GET1, GET2>::RemoteIntersectionImpl::RemoteIntersectionImpl(const Parent* glue, unsigned int mergeindex)
  : _glue(glue), _mergeindex(mergeindex), _index((IndexType)-1),
    _domainindex((IndexType)-1), _targetindex((IndexType)-1)
{
  // if an invalid index is given do not proceed!
  // (happens when the parent GridGlue initializes the "end"-Intersection)
  assert (0 <= mergeindex || mergeindex < glue->_index_sz);

  // initialize the local and the global geometry of the domain
  {
    // coordinates within the subentity that contains the remote intersection
    Dune::array<LocalCoords, coorddim> corners_subEntity_local;

    for (int i = 0; i < coorddim; ++i)
      corners_subEntity_local[i] = glue->_merg->domainParentLocal(mergeindex, i);

    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementcoorddim = DomainGridType::template Codim<0>::Geometry::mydimension;

    // a face number is only important if dealing with surfaces, not meshes

    // Coordinates of the remote intersection corners wrt the element coordinate system
    Dune::array<Dune::FieldVector<ctype, elementcoorddim>, coorddim> corners_element_local;

    // world coordinates of the remote intersection corners
    Dune::array<Dune::FieldVector<ctype, DomainGridType::dimensionworld>, coorddim> corners_global;

    unsigned int domainIndex = glue->_merg->domainParent(mergeindex);
    unsigned int unused;
    if (glue->_domext.contains(domainIndex, unused))
    {
      glue->_domext.localAndGlobalCoords(domainIndex, corners_subEntity_local, corners_element_local, corners_global, coorddim);

      // set the corners of the geometries
      this->_domlgeom.setup(geometrytype, corners_element_local);
      this->_domggeom.setup(geometrytype, corners_global);
    }
  }

  // do the same for the local and the global geometry of the target
  {
    // coordinates within the subentity that contains the remote intersection
    Dune::array<LocalCoords, coorddim> corners_subEntity_local;

    for (int i = 0; i < coorddim; ++i)
      corners_subEntity_local[i] = glue->_merg->targetParentLocal(mergeindex, i);

    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementcoorddim = TargetGridType::template Codim<0>::Geometry::mydimension;

    // a face number is only important if dealing with surfaces, not manifolds or meshes

    // Coordinates of the remote intersection corners wrt the element coordinate system
    Dune::array<Dune::FieldVector<ctype, elementcoorddim>, coorddim> corners_element_local;

    // world coordinates of the remote intersection corners
    Dune::array<Dune::FieldVector<ctype, TargetGridType::dimensionworld>, coorddim> corners_global;

    unsigned int targetIndex = glue->_merg->targetParent(mergeindex);
    unsigned int unused;
    if (glue->_tarext.contains(targetIndex, unused))
    {
      glue->_tarext.localAndGlobalCoords(targetIndex, corners_subEntity_local, corners_element_local, corners_global, coorddim);

      // set the corners of the geometries
      this->_tarlgeom.setup(geometrytype, corners_element_local);
      this->_targgeom.setup(geometrytype, corners_global);
    }
  }
}

#endif // GRIDGLUEREMOTEINTERSECTIONIMPL_HH_
