// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    RemoteIntersectionConcepts.hh
 *  Version:     1.0
 *  Created on:  Feb 28, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Concept class for the RemoteIntersection concept
 *  subversion:  $Id$
 *
 */
/**
 * @file RemoteIntersectionConcepts.hh
 * @brief Class that can be used to syntactically check models of the
 * RemoteIntersection concept.
 */

#ifndef REMOTEINTERSECTIONCONCEPTS_HH_
#define REMOTEINTERSECTIONCONCEPTS_HH_

#include "../../dune/common/geometrytype.hh"
#include "../../dune/common/fvector.hh"


namespace RemoteIntersectionInterface
{


  /**
   * @class RemoteInterfaceConcept
   * @brief concept class for the implementation of the RemoteIntersection interface
   */
  template<typename T>
  class RemoteIntersectionConcept
  {
  public:
    /// @brief the method that calls to everything that is required of a model of this concept
    void constraints()
    {

      // return EntityPointer to the Entity on the inside of this intersection. That is the Entity where we started this .
      deptr = impl.entityDomain();

      // return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity.
      teptr = impl.entityTarget();

      // geometrical information about this intersection in local coordinates of the inside() entity.
      // takes local domain intersection coords and maps them to domain parent element local coords
      dlgeom = impl.geometryInDomainEntity();

      // geometrical information about this intersection in local coordinates of the outside() entity.
      // takes local target intersection coords and maps them to target parent element local coords
      tlgeom = impl.geometryInTargetEntity();

      // geometrical information about this intersection in global coordinates in the domain grid.
      // takes local domain intersection coords and maps them to domain grid world coords
      dggeom = impl.geometryDomain();

      // geometrical information about this intersection in global coordinates in the target grid.
      // takes local target intersection coords and maps them to target grid world coords
      tggeom = impl.geometryTarget();

      // a direct mapping from local target intersection coords to local domain intersection coords
      local = impl.domainLocals(local);

      // a direct mapping from local domain intersection coords to local target intersection coords
      local = impl.targetLocals(local);

      // obtain the type of reference element for this intersection
      gt = impl.type();

      // Local number of codim 1 entity in the inside() Entity where intersection is contained in.
      i = impl.numberInDomainEntity();

      // Local number of codim 1 entity in outside() Entity where intersection is contained in.
      i = impl.numberInTargetEntity();

      // Return an outer normal (length not necessarily 1).
      normal = impl.outerNormalDomain(local);

      // Return an outer normal (length not necessarily 1).
      normal = impl.outerNormalTarget(local);
    }

  private:


    typedef typename T::DomainGridView DomainGridView;

    typedef typename T::TargetGridView TargetGridView;


    enum { coorddim = T::coorddim };

    enum { mydim = T::mydim };


    typedef typename T::DomainLocalGeometry DomainLocalGeometry;

    typedef typename T::DomainGeometry DomainGeometry;


    typedef typename T::TargetLocalGeometry TargetLocalGeometry;

    typedef typename T::TargetGeometry TargetGeometry;


    typedef typename T::ctype ctype;


    typedef typename T::LocalCoords LocalCoords;

    typedef typename T::Coords Coords;


    T impl;

    typename DomainGridView::Traits::template Codim<0>::EntityPointer deptr;
    typename TargetGridView::Traits::template Codim<0>::EntityPointer teptr;

    DomainLocalGeometry dlgeom;
    DomainGeometry dggeom;
    TargetLocalGeometry tlgeom;
    TargetGeometry tggeom;

    Dune::FieldVector<ctype, coorddim> normal;
    Dune::FieldVector<ctype, mydim> local;

    Dune::GeometryType gt;

    bool b;
    int i;

  };


  /**
   * @class RemoteInterfaceIteratorConcept
   * @brief concept class for the implementations of the iterators in RemoteIntersectionIterators interface class
   */
  template<typename T>
  class RemoteIntersectionIteratorConcept
  {
  public:
    /// @brief the method that calls to everything that is required of a model of this concept
    void constraints()
    {
      // dereferencing the iterator to get a reference to the
      // RemoteIntersection itself
      ris = constiter.dereference();

      // comparison between iterators
      b = constiter.equals(iter);

      // moving the iterator
      iter.increment();
    }
  private:
    const T constiter;
    T iter;
    typename T::RemoteIntersection ris;
    bool b;
  };






}

#endif // REMOTEINTERSECTIONCONCEPTS_HH_
