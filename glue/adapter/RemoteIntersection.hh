// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    RemoteIntersection.hh
 *  Version:     1.0
 *  Created on:  Feb 22, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: abstract interface of the RemoteInteresection concept
 *  subversion:  $Id$
 *
 */
/**
 * @file RemoteIntersection.hh
 * @brief abstract interface of the RemoteInteresection concept
 */

#ifndef REMOTEINTERSECTION_HH_
#define REMOTEINTERSECTION_HH_


#include <dune/common/fvector.hh>
#include <dune/common/geometrytype.hh>

#ifdef GRID_GLUE_USE_CONCEPTS
#include "../misc/conceptchecking.hh"
#include "RemoteIntersectionConcepts.hh"
#endif

namespace RemoteIntersectionInterface
{


  /**
   * @class RemoteIntersection
   * @brief Analogous to the Intersection interface this is an interface for
   * "intersections" between boundary faces of two different grids.
   *
   *
   */
  template<typename RemoteIntersectionImpl>
  class RemoteIntersection
  {
  private:
    /*   C H E C K   C O N C E P T S   */

    typedef RemoteIntersectionImpl RemoteIntersectionImplType;
#ifdef GRID_GLUE_USE_CONCEPTS
    CLASS_REQUIRE(RemoteIntersectionImplType, RemoteIntersectionConcept);
#endif

  private:

    /// @brief the actual implementation of the intersection
    RemoteIntersectionImpl realIntersection;

  protected:

    RemoteIntersectionImpl& getImpl()
    {
      return this->realIntersection;
    }

  public:

    enum { coorddim = RemoteIntersectionImpl::coorddim };

    enum { mydim = RemoteIntersectionImpl::mydim };


    typedef typename RemoteIntersectionImpl::DomainGridView DomainGridView;

    typedef typename DomainGridView::Grid DomainGridType;

    typedef typename RemoteIntersectionImpl::TargetGridView TargetGridView;

    typedef typename TargetGridView::Grid TargetGridType;


    typedef typename DomainGridView::Traits::template Codim<0>::Entity DomainEntity;

    typedef typename DomainGridView::Traits::template Codim<0>::EntityPointer DomainEntityPointer;

    typedef typename TargetGridView::Traits::template Codim<0>::Entity TargetEntity;

    typedef typename TargetGridView::Traits::template Codim<0>::EntityPointer TargetEntityPointer;


    typedef typename RemoteIntersectionImpl::DomainLocalGeometry DomainLocalGeometry;

    typedef typename RemoteIntersectionImpl::DomainGeometry DomainGeometry;

    typedef typename RemoteIntersectionImpl::TargetLocalGeometry TargetLocalGeometry;

    typedef typename RemoteIntersectionImpl::TargetGeometry TargetGeometry;


    typedef typename RemoteIntersectionImpl::ctype ctype;


    typedef typename RemoteIntersectionImpl::LocalCoords LocalCoords;

    typedef typename RemoteIntersectionImpl::Coords Coords;


    /*   C O N S T R U C T O R S   */
    /** \brief Copy construction from another RemoteIntersection */
    RemoteIntersection(const RemoteIntersection& impl)
      : realIntersection(impl.realIntersection)
    {}

    /** \brief Copy construction from an implementation class */
    RemoteIntersection(const RemoteIntersectionImpl& impl) : realIntersection(impl)
    {}

    /** \brief Assigment from another RemoteIntersection */
    RemoteIntersection& operator=(const RemoteIntersection& impl)
    {
      this->realIntersection = impl.realIntersection;
      return *this;
    }

    /** \brief Assignment from an implementation class */
    RemoteIntersection& operator=(const RemoteIntersectionImpl& impl)
    {
      this->realIntersection = impl;
      return *this;
    }


    /*   F U N C T I O N A L I T Y   */

    /** \brief return EntityPointer to the Entity on the inside of this intersection.
        That is the Entity where we started this. */
    DomainEntityPointer entityDomain() const
    {
      return this->realIntersection.entityDomain();
    }


    /** \brief return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity. */
    TargetEntityPointer entityTarget() const
    {
      return this->realIntersection.entityTarget();
    }


    /** \brief geometrical information about this intersection in local coordinates of the inside() entity.
        takes local domain intersection coords and maps them to domain parent element local coords */
    const DomainLocalGeometry& intersectionDomainLocal() const
    {
      return this->realIntersection.geometryInDomainEntity();
    }


    /** \brief geometrical information about this intersection in local coordinates of the outside() entity.
        takes local target intersection coords and maps them to target parent element local coords */
    const TargetLocalGeometry& intersectionTargetLocal() const
    {
      return this->realIntersection.geometryInTargetEntity();
    }


    /** \brief geometrical information about this intersection in global coordinates in the domain grid.
        takes local domain intersection coords and maps them to domain grid world coords */
    const DomainGeometry& intersectionDomainGlobal() const
    {
      return this->realIntersection.geometryDomain();
    }


    /** \brief geometrical information about this intersection in global coordinates in the target grid. */
    const TargetGeometry& intersectionTargetGlobal() const
    {
      return this->realIntersection.geometryTarget();
    }


    /*! @brief mapping of local coordinates, i.e. target -> domain local coords

       Maps a point in local intersection coords on target intersection to a
       point in domain intersection local coordinates.
     */
    FieldVector<ctype, mydim> domainLocals(const FieldVector<ctype, mydim> &local) const
    {
      return this->realIntersection.domainLocals(local);
    }


    /*! @brief mapping of local coordinates, i.e. domain -> target local coords

       Maps a point in local intersection coords on domain intersection to a
       point in target intersection local coordinates.
     */
    FieldVector<ctype, mydim> targetLocals(const FieldVector<ctype, mydim> &local) const
    {
      return this->realIntersection.targetLocals(local);
    }


    /** \brief obtain the type of reference element for this intersection */
    GeometryType type() const
    {
      return this->realIntersection.type();
    }


    /** \brief Local number of codim 1 entity in the inside() Entity where intersection is contained in. */
    int numberInDomainEntity() const
    {
      return this->realIntersection.numberInDomainEntity();
    }


    /** \brief Local number of codim 1 entity in outside() Entity where intersection is contained in. */
    int numberInTargetEntity() const
    {
      return this->realIntersection.numberInTargetEntity();
    }


    /** \brief Return an outer normal (length not necessarily 1). */
    FieldVector<ctype, coorddim>    outerNormalDomain(const FieldVector<ctype, mydim> &local) const
    {
      return this->realIntersection.outerNormalDomain(local);
    }


    /** \brief Return an outer normal */
    FieldVector<ctype, coorddim>    unitOuterNormalDomain(const FieldVector<ctype, mydim> &local) const
    {
      FieldVector<ctype, coorddim> normal = this->realIntersection.outerNormalDomain(local);
      normal /= normal.two_norm();
      return normal;
    }


    /** \brief Return an outer normal (length not necessarily 1) */
    FieldVector<ctype, coorddim>    integrationOuterNormalDomain(const FieldVector<ctype, mydim> &local) const
    {
      return this->unitOuterNormalDomain(local) *= this->realIntersection.geometryDomain().integrationElement(local);
    }


    /** \brief Return an outer normal (length not necessarily 1). */
    FieldVector<ctype, coorddim>    outerNormalTarget(const FieldVector<ctype, mydim> &local) const
    {
      return this->realIntersection.outerNormalTarget(local);
    }


    /** \brief Return a unit outer normal of the target intersection */
    FieldVector<ctype, coorddim>    unitOuterNormalTarget(const FieldVector<ctype, mydim> &local) const
    {
      FieldVector<ctype, coorddim> normal = this->realIntersection.outerNormalTarget(local);
      normal /= normal.two_norm();
      return normal;
    }


    /** \brief Return an outer normal (length not necessarily 1) */
    FieldVector<ctype, coorddim>    integrationOuterNormalTarget(const FieldVector<ctype, mydim> &local) const
    {
      return this->unitOuterNormalTarget(local) *= this->realIntersection.geometryTarget().integrationElement(local);
    }


  };

} // end namespace RemoteIntersectionInterface

#endif // REMOTEINTERSECTION_HH_
