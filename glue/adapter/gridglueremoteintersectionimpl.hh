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


template<typename GET1, typename GET2, typename SM>
class GridGlue<GET1, GET2, SM>::RemoteIntersectionImpl
{

private:

  typedef GridGlue<GET1, GET2, SM> Parent;

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


private:

  /*   M E M B E R   V A R  I A B L E S   */

  /// @brief the grid glue entity this is built on
  const Parent*   _glue;

  /// @brief index of this intersection after GridGlue interface
  int _index;


  DomainLocalGeometry _domlgeom;

  DomainGeometry _domggeom;

  TargetLocalGeometry _tarlgeom;

  TargetGeometry _targgeom;


public:

  static const Dune::GeometryType geometrytype;


  /*   F U N C T I O N A L I T Y   */

  RemoteIntersectionImpl(const Parent* glue, int index);

  /** \brief Copy constructor */
  RemoteIntersectionImpl(const RemoteIntersectionImpl & impl)
    : _glue(impl._glue), _index(impl._index)
      ,
      _domlgeom(impl._domlgeom),
      _domggeom(impl._domggeom),
      _tarlgeom(impl._tarlgeom),
      _targgeom(impl._targgeom)
  {}

  /** \brief Assignment operator */
  RemoteIntersectionImpl& operator=(const RemoteIntersectionImpl& impl)
  {
    this->_glue = impl._glue;
    this->_index = impl._index;
    this->_domlgeom = impl._domlgeom;
    this->_domggeom = impl._domggeom;
    this->_tarlgeom = impl._tarlgeom;
    this->_targgeom = impl._targgeom;
    return *this;
  }


  // return EntityPointer to the Entity on the inside of this intersection. That is the Entity where we started this .
  DomainEntityPointer entityDomain() const
  {
    return this->_glue->_domext.element(this->_glue->_sm.domainParent(this->_index));
  }


  // return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity.
  TargetEntityPointer entityTarget() const
  {
    return this->_glue->_tarext.element(this->_glue->_sm.targetParent(this->_index));
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


  /** direct mapping of local coordinates, i.e. target local -> domain local */
  /** \todo Documentation does not match implementation */
  LocalCoords domainLocals(const LocalCoords &local) const
  {
    return barycentricToReference(this->_glue->_sm.domainLocals(this->_index, referenceToBarycentric(local)));
  }


  /** \brief direct mapping of local coordinates, i.e. domain local -> target local */
  /** \todo Documentation does not match implementation */
  LocalCoords targetLocals(const LocalCoords &local) const
  {
    return barycentricToReference(this->_glue->_sm.targetLocals(this->_index, referenceToBarycentric(local)));
  }


  // obtain the type of reference element for this intersection
  Dune::GeometryType type() const
  {
    return geometrytype;
  }


  // Local number of codim 1 entity in the inside() Entity where intersection is contained in.
  int numberInDomainEntity() const
  {
    return this->_glue->_domext.numberInSelf(this->_glue->_sm.domainParent(this->_index));
  }


  // Local number of codim 1 entity in outside() Entity where intersection is contained in.
  int numberInTargetEntity() const
  {
    return this->_glue->_tarext.numberInSelf(this->_glue->_sm.targetParent(this->_index));
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
  typedef int IndexType;

  IndexType index() const
  {
    return this->_index;
  }

#endif
};


template<typename GET1, typename GET2, typename SM>
const Dune::GeometryType GridGlue<GET1, GET2, SM>::RemoteIntersectionImpl::geometrytype(Dune::GeometryType::simplex, Parent::dimw-1);


template<typename GET1, typename GET2, typename SM>
bool GridGlue<GET1, GET2, SM>::RemoteIntersectionImpl::conforming() const
{
  std::vector<unsigned int> results;
  // first check the domain side
  bool is_conforming =
    this->_glue->_sm.domainSimplexRefined(this->_glue->_sm.domainParent(this->_index), results) && results.size() == 1;
  results.resize(0);
  // now check the target side
  if (is_conforming)
    return this->_glue->_sm.targetSimplexRefined(this->_glue->_sm.targetParent(this->_index), results) && results.size() == 1;
  return false;
}


template<typename GET1, typename GET2, typename SM>
GridGlue<GET1, GET2, SM>::RemoteIntersectionImpl::RemoteIntersectionImpl(const Parent* glue, int index)
  : _glue(glue), _index(index)
{
  // if an invalid index is given do not proceed!
  // (happens when the parent GridGlue initializes the "end"-Intersection)
  if (0 > index || index >= static_cast<int>(glue->_sm.nSimplices()))
    return;

  // the number of simplex corners as well as
  // the number of the coordinate space's dimensions
  Dune::array<Coords, coorddim> corners_barycentric;

  // initialize the local and the global geometry of the domain
  {
    for (int i = 0; i < coorddim; ++i)
      corners_barycentric[i] = glue->_sm.domainParentLocal(index, i);

    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementcoorddim = DomainGridType::template Codim<0>::Geometry::mydimension;

    // a face number is only important if dealing with surfaces, not meshes

    // compute the local coordinates with the correct dimension and the global coordinates
    Dune::array<Dune::FieldVector<ctype, elementcoorddim>, coorddim> corners_local;
    Dune::array<Dune::FieldVector<ctype, DomainGridType::dimensionworld>, coorddim> corners_global;
    glue->_domext.localAndGlobalCoords(glue->_sm.domainParent(index), corners_barycentric, corners_local, corners_global, coorddim);

    // set the corners of the geometries
    this->_domlgeom.setup(geometrytype, corners_local);
    this->_domggeom.setup(geometrytype, corners_global);
  }

  // do the same for the local and the global geometry of the target
  {
    for (int i = 0; i < coorddim; ++i)
      corners_barycentric[i] = glue->_sm.targetParentLocal(index, i);

    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementcoorddim = TargetGridType::template Codim<0>::Geometry::mydimension;

    // a face number is only important if dealing with surfaces, not manifolds or meshes

    // compute the local coordinates with the correct dimension and the global coordinates
    Dune::array<Dune::FieldVector<ctype, elementcoorddim>, coorddim> corners_local;
    Dune::array<Dune::FieldVector<ctype, TargetGridType::dimensionworld>, coorddim> corners_global;
    glue->_tarext.localAndGlobalCoords(glue->_sm.targetParent(index), corners_barycentric, corners_local, corners_global, coorddim);

    // set the corners of the geometries
    this->_tarlgeom.setup(geometrytype, corners_local);
    this->_targgeom.setup(geometrytype, corners_global);
  }
}

#endif // GRIDGLUEREMOTEINTERSECTIONIMPL_HH_
