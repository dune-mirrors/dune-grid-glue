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

#include <dune/grid-glue/adapter/simplexgeometry.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

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


  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename DomainGridView::Grid DomainGridType;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename TargetGridView::Grid TargetGridType;

  enum { coorddim = Parent::dimworld };

  dune_static_assert(DomainGridType::dimension - DomainExtractor::codim
                     == TargetGridType::dimension - TargetExtractor::codim,
                     "Currently both coupling extracts need to have the same dimension!");

  // Use domain for the definition of mydim because we assume that the value for target is the same
  enum { mydim = DomainGridType::dimension - DomainExtractor::codim };


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
  const Parent*       glue_;

  /// @brief index of this intersection after GridGlue interface
  IndexType mergeindex_;

  IndexType index_;
  IndexType domainindex_;
  IndexType targetindex_;

  DomainLocalGeometry domlgeom_;

  DomainGeometry domggeom_;

  TargetLocalGeometry tarlgeom_;

  TargetGeometry targgeom_;

public:

  /*   F U N C T I O N A L I T Y   */

  RemoteIntersectionImpl(const Parent* glue)
    : glue_(glue), mergeindex_((IndexType)-1), index_((IndexType)-1),
      domainindex_((IndexType)-1), targetindex_((IndexType)-1)
  {}

  RemoteIntersectionImpl(const Parent* glue, unsigned int mergeindex);

  /** \brief Copy constructor */
  RemoteIntersectionImpl(const RemoteIntersectionImpl & impl)
  {
    *this = impl;
  }

  /** \brief Assignment operator */
  RemoteIntersectionImpl& operator=(const RemoteIntersectionImpl& impl)
  {
    this->glue_ = impl.glue_;
    this->mergeindex_ = impl.mergeindex_;
    this->index_ = impl.index_;
    this->domainindex_ = impl.domainindex_;
    this->targetindex_ = impl.targetindex_;
    this->domlgeom_ = impl.domlgeom_;
    this->domggeom_ = impl.domggeom_;
    this->tarlgeom_ = impl.tarlgeom_;
    this->targgeom_ = impl.targgeom_;
    return *this;
  }

  // return EntityPointer to the Entity on the inside of this intersection. That is the Entity where we started this .
  DomainEntityPointer entityDomain() const
  {
    return this->glue_->domext_.element(this->glue_->merger_->template parent<0>(this->mergeindex_));
  }


  // return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity.
  TargetEntityPointer entityTarget() const
  {
    return this->glue_->tarext_.element(this->glue_->merger_->template parent<1>(this->mergeindex_));
  }


  /** \brief Return true if intersection is conforming */
  bool conforming() const;


  // geometrical information about this intersection in local coordinates of the inside() entity.
  const DomainLocalGeometry& geometryInDomainEntity() const
  {
    return this->domlgeom_;
  }


  // geometrical information about this intersection in local coordinates of the outside() entity.
  const TargetLocalGeometry& geometryInTargetEntity() const
  {
    return this->tarlgeom_;
  }


  // geometrical information about this intersection in global coordinates in the domain grid.
  const DomainGeometry& geometryDomain() const
  {
    return this->domggeom_;
  }


  // geometrical information about this intersection in global coordinates in the target grid.
  const TargetGeometry& geometryTarget() const
  {
    return this->targgeom_;
  }

  bool hasTarget() const
  {
    unsigned int localindex;
    return this->glue_->tarext_.contains(this->glue_->merger_->template parent<1>(this->mergeindex_), localindex);
  }

  bool hasDomain() const
  {
    unsigned int localindex;
    return this->glue_->domext_.contains(this->glue_->merger_->template parent<0>(this->mergeindex_), localindex);
  }

  // obtain the type of reference element for this intersection
  Dune::GeometryType type() const
  {
    return Dune::GeometryType(Dune::GeometryType::simplex, mydim);
  }


  // Local number of codim 1 entity in the inside() Entity where intersection is contained in.
  int numberInDomainEntity() const
  {
    return this->glue_->domext_.indexInInside(this->glue_->merger_->template parent<0>(this->mergeindex_));
  }


  // Local number of codim 1 entity in outside() Entity where intersection is contained in.
  int numberInTargetEntity() const
  {
    return this->glue_->tarext_.indexInInside(this->glue_->merger_->template parent<1>(this->mergeindex_));
  }


  // Return an outer normal (length not necessarily 1).
  Coords outerNormalDomain(const Dune::FieldVector<ctype, mydim> &local) const
  {
    return this->domggeom_.outerNormal(local);
  }

  // Return an outer normal (length not necessarily 1).
  Coords outerNormalTarget(const Dune::FieldVector<ctype, mydim> &local) const
  {
    return this->targgeom_.outerNormal(local);
  }

#ifdef QUICKHACK_INDEX
  IndexType & index()
  {
    return this->index_;
  }

  IndexType & domainIndex()
  {
    return this->domainindex_;
  }

  IndexType & targetIndex()
  {
    return this->targetindex_;
  }

  IndexType index() const
  {
    return this->index_;
  }

  IndexType globalIndex() const
  {
    return this->mergeindex_;
  }

  IndexType domainIndex() const
  {
    assert(this->domainindex_ != (IndexType)-1);
    return this->domainindex_;
  }

  IndexType targetIndex() const
  {
    assert(this->targetindex_ != (IndexType)-1);
    return this->targetindex_;
  }

#endif
};


template<typename GET1, typename GET2>
bool GridGlue<GET1, GET2>::RemoteIntersectionImpl::conforming() const
{
  std::vector<unsigned int> results;
  // first check the domain side
  bool is_conforming =
    this->glue_->merger_->template simplexRefined<0>(this->glue_->merger_->template parent<0>(this->mergeindex_), results) && results.size() == 1;
  results.resize(0);
  // now check the target side
  if (is_conforming)
    return this->glue_->merger_->template simplexRefined<1>(this->glue_->merger_->template parent<1>(this->mergeindex_), results) && results.size() == 1;
  return false;
}


template<typename GET1, typename GET2>
GridGlue<GET1, GET2>::RemoteIntersectionImpl::RemoteIntersectionImpl(const Parent* glue, unsigned int mergeindex)
  : glue_(glue), mergeindex_(mergeindex), index_((IndexType)-1),
    domainindex_((IndexType)-1), targetindex_((IndexType)-1)
{
  // if an invalid index is given do not proceed!
  // (happens when the parent GridGlue initializes the "end"-Intersection)
  assert (0 <= mergeindex || mergeindex < glue->index__sz);

  // initialize the local and the global geometry of the domain
  {
    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementdim = DomainGridType::template Codim<0>::Geometry::mydimension;

    const int nSimplexCorners = elementdim - Parent::DomainExtractor::codim + 1;

    // coordinates within the subentity that contains the remote intersection
    Dune::array<LocalCoords, nSimplexCorners> corners_subEntity_local;

    for (int i = 0; i < nSimplexCorners; ++i)
      corners_subEntity_local[i] = glue->merger_->template parentLocal<0>(mergeindex, i);

    // Coordinates of the remote intersection corners wrt the element coordinate system
    Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

    // world coordinates of the remote intersection corners
    Dune::array<Dune::FieldVector<ctype, DomainGridType::dimensionworld>, nSimplexCorners> corners_global;

    unsigned int domainIndex = glue->merger_->template parent<0>(mergeindex);
    unsigned int unused;
    if (glue->domext_.contains(domainIndex, unused))
    {
      typename DomainExtractor::Geometry domainWorldGeometry = glue->domext_.geometry(domainIndex);
      typename DomainExtractor::LocalGeometry domainLocalGeometry = glue->domext_.geometryLocal(domainIndex);

      for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
        corners_element_local[i] = domainLocalGeometry.global(corners_subEntity_local[i]);
        corners_global[i]        = domainWorldGeometry.global(corners_subEntity_local[i]);
      }

      // set the corners of the geometries
      this->domlgeom_.setup(type(), corners_element_local);
      this->domggeom_.setup(type(), corners_global);
    }
  }

  // do the same for the local and the global geometry of the target
  {
    // compute the coordinates of the subface's corners in codim 0 entity local coordinates
    const int elementdim = TargetGridType::template Codim<0>::Geometry::mydimension;

    const int nSimplexCorners = elementdim - Parent::TargetExtractor::codim + 1;

    // coordinates within the subentity that contains the remote intersection
    Dune::array<LocalCoords, nSimplexCorners> corners_subEntity_local;

    for (int i = 0; i < nSimplexCorners; ++i)
      corners_subEntity_local[i] = glue->merger_->template parentLocal<1>(mergeindex, i);

    // Coordinates of the remote intersection corners wrt the element coordinate system
    Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

    // world coordinates of the remote intersection corners
    Dune::array<Dune::FieldVector<ctype, TargetGridType::dimensionworld>, nSimplexCorners> corners_global;

    unsigned int targetIndex = glue->merger_->template parent<0>(mergeindex);
    unsigned int unused;
    if (glue->tarext_.contains(targetIndex, unused))
    {
      typename TargetExtractor::Geometry targetWorldGeometry = glue->tarext_.geometry(targetIndex);
      typename TargetExtractor::LocalGeometry targetLocalGeometry = glue->tarext_.geometryLocal(targetIndex);

      for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
        corners_element_local[i] = targetLocalGeometry.global(corners_subEntity_local[i]);
        corners_global[i]        = targetWorldGeometry.global(corners_subEntity_local[i]);
      }

      // set the corners of the geometries
      this->tarlgeom_.setup(type(), corners_element_local);
      this->targgeom_.setup(type(), corners_global);
    }
  }
}

#endif // GRIDGLUEREMOTEINTERSECTIONIMPL_HH_
