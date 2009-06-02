// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    RemoteIntersectionIterators.hh
 *  Version:     1.0
 *  Created on:  Feb 24, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: abstract interface for the iterators on remote intersections
 *  subversion:  $Id$
 *
 */
/**
 * @file RemoteIntersectionIterators.hh
 * @brief Provide iterators over remote intersections
 */

#ifndef REMOTEINTERSECTIONITERATORS_HH_
#define REMOTEINTERSECTIONITERATORS_HH_

#ifdef GRID_GLUE_USE_CONCEPTS
#include "../misc/conceptchecking.hh"
#include "RemoteIntersectionConcepts.hh"
#endif
#include "RemoteIntersection.hh"


namespace RemoteIntersectionInterface
{

  /** \brief Iterator over remote intersections */
  template<typename RemoteIntersectionImpl,
      typename IndependentIteratorImpl,
      typename DomainIteratorImpl,
      typename TargetIteratorImpl
      >
  class RemoteIntersectionIterators
  {
  private:
    /*   C H E C K   C O N C E P T S   */

#ifdef GRID_GLUE_USE_CONCEPTS
    typedef RemoteIntersectionImpl RemoteIntersectionImplType;
    CLASS_REQUIRE(RemoteIntersectionImplType, RemoteIntersectionConcept);

    typedef IndependentIteratorImpl IndependentIteratorImplType;
    CLASS_REQUIRE(IndependentIteratorImplType, RemoteIntersectionIteratorConcept);

    typedef DomainIteratorImpl DomainIteratorImplType;
    CLASS_REQUIRE(DomainIteratorImplType, RemoteIntersectionIteratorConcept);

    typedef TargetIteratorImpl TargetIteratorImplType;
    CLASS_REQUIRE(TargetIteratorImplType, RemoteIntersectionIteratorConcept);
#else
    typedef RemoteIntersectionImpl RemoteIntersectionImplType;

    typedef IndependentIteratorImpl IndependentIteratorImplType;

    typedef DomainIteratorImpl DomainIteratorImplType;

    typedef TargetIteratorImpl TargetIteratorImplType;
#endif

  public:

    typedef RemoteIntersection<RemoteIntersectionImpl>  ThisRemoteIntersection;

    /** \todo Please doc me! */
    class IndependentIterator
    {

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


    private:

      IndependentIteratorImplType realIterator;


    protected:

      IndependentIteratorImplType& getImpl()
      {
        return this->realIterator;
      }


    public:

      /**
       * @brief copy constructor
       * @param iter the iterator to copy
       */
      IndependentIterator(const IndependentIterator& iter) : realIterator(iter.realIterator)
      {}


      /**
       * @brief only constructor to build a "new" iterator
       * @param iter the iterator to initialize from
       */
      IndependentIterator(const IndependentIteratorImplType& iter) : realIterator(iter)
      {}


      /** \brief Dereferencing operator. */
      const ThisRemoteIntersection& operator*() const
      {
        return this->realIterator.dereference();
      }


      /** \brief Pointer operator. */
      const ThisRemoteIntersection* operator->() const
      {
        return &this->realIterator.dereference();
      }


      /** \brief Increment operator. */
      IndependentIterator& operator++()
      {
        this->realIterator.increment();
        return *this;
      }


      bool operator==(const IndependentIterator& rhs) const
      {
        return rhs.realIterator.equals(this->realIterator);
      }


      bool operator!=(const IndependentIterator& rhs) const
      {
        return !rhs.realIterator.equals(this->realIterator);
      }

    };


    class DomainIterator
    {

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


    private:

      DomainIteratorImpl realIterator;


    protected:

      DomainIteratorImpl& getImpl()
      {
        return this->realIterator;
      }


    public:

      /**
       * @brief copy constructor
       * @param iter the iterator to copy
       */
      DomainIterator(const DomainIterator& iter) : realIterator(iter.realIterator)
      {}


      /**
       * @brief only constructor to build a "new" iterator
       * @param iter the iterator to initialize from
       */
      DomainIterator(const DomainIteratorImpl& iter) : realIterator(iter)
      {}


      /** \brief Dereferencing operator. */
      const ThisRemoteIntersection& operator*() const
      {
        return this->realIterator.dereference();
      }


      /** \brief Pointer operator. */
      const ThisRemoteIntersection* operator->() const
      {
        return &this->realIterator.dereference();
      }


      /** \brief Increment operator. */
      DomainIterator& operator++()
      {
        this->realIterator.increment();
        return *this;
      }


      bool operator==(const DomainIterator& rhs) const
      {
        return rhs.realIterator.equals(this->realIterator);
      }


      bool operator!=(const DomainIterator& rhs) const
      {
        return !rhs.realIterator.equals(this->realIterator);
      }


      /*! @brief return EntityPointer to the Entity on the inside of this
         intersection. That is the Entity where we started this .
       */
      DomainEntityPointer inside() const
      {
        return this->realIterator.dereference().entityDomain();
      }


      /*! @brief return EntityPointer to the Entity on the outside of this
         intersection. That is the neighboring Entity.

         @warning Don't call this method if there is no neighboring Entity
         (neighbor() returns false). In this case the result is undefined.
       */
      TargetEntityPointer outside() const
      {
        return this->realIterator.dereference().entityTarget();
      }


      /*! @brief mapping of local coordinates, i.e. selflocal -> neighborlocal

         Maps a point in local intersection coords on domain intersection to a
         point in target intersection local coordinates.
       */
      FieldVector<ctype, mydim> neighborLocals(const FieldVector<ctype, mydim> &selflocal) const
      {
        return this->realIterator.dereference().targetLocals(selflocal);
      }


      /*! @brief mapping of local coordinates, i.e. neighborlocal -> selflocal

         Maps a point in local intersection coords on target intersection to a
         point in domain intersection local coordinates.
       */
      FieldVector<ctype, mydim> selfLocals(const FieldVector<ctype, mydim> &neighborlocal) const
      {
        return this->realIterator.dereference().domainLocals(neighborlocal);
      }


      /** \brief obtain the type of reference element for this intersection */
      GeometryType type() const
      {
        return this->realIterator.dereference().type();
      }


      //! Local number of codim 1 entity in the inside() Entity where intersection is contained in
      int numberInSelf() const
      {
        return this->realIterator.dereference().numberInDomainEntity();
      }


      //! Local number of codim 1 entity in outside() Entity where intersection is contained in
      int numberInNeighbor() const
      {
        return this->realIterator.dereference().numberInTargetEntity();
      }

    };

    class TargetIterator
    {
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


    private:

      TargetIteratorImpl realIterator;


    protected:


      TargetIteratorImpl& getImpl()
      {
        return this->realIterator;
      }

    public:

      /**
       * @brief copy constructor
       * @param iter the iterator to copy
       */
      TargetIterator(const TargetIterator& iter) : realIterator(iter.realIterator)
      {}


      /**
       * @brief only constructor to build a "new" iterator
       * @param iter the iterator to initialize from
       */
      TargetIterator(const TargetIteratorImpl& iter) : realIterator(iter)
      {}


      /** \brief Dereferencing operator. */
      const ThisRemoteIntersection& operator*() const
      {
        return this->realIterator.dereference();
      }


      /** \brief Pointer operator. */
      const ThisRemoteIntersection* operator->() const
      {
        return &this->realIterator.dereference();
      }


      /** \brief Increment operator. */
      TargetIterator& operator++()
      {
        this->realIterator.increment();
        return *this;
      }


      bool operator==(const TargetIterator& rhs) const
      {
        return rhs.realIterator.equals(this->realIterator);
      }


      bool operator!=(const TargetIterator& rhs) const
      {
        return !rhs.realIterator.equals(this->realIterator);
      }


      /*! @brief return EntityPointer to the Entity on the inside of this
         intersection. That is the Entity where we started this .
       */
      TargetEntityPointer inside() const
      {
        return this->realIterator.dereference().entityTarget();
      }


      /*! @brief return EntityPointer to the Entity on the outside of this
         intersection. That is the neighboring Entity.

         @warning Don't call this method if there is no neighboring Entity
         (neighbor() returns false). In this case the result is undefined.
       */
      DomainEntityPointer outside() const
      {
        return this->realIterator.dereference().entityDomain();
      }


      /*! @brief mapping of local coordinates, i.e. selflocal -> neighborlocal

         Maps a point in local intersection coords on target intersection to a
         point in domain intersection local coordinates.
       */
      FieldVector<ctype, mydim> neighborLocals(const FieldVector<ctype, mydim> &selflocal) const
      {
        return this->realIterator.dereference().domainLocals(selflocal);
      }


      /*! @brief mapping of local coordinates, i.e. neighborlocal -> selflocal

         Maps a point in local intersection coords on domain intersection to a
         point in target intersection local coordinates.
       */
      FieldVector<ctype, mydim> selfLocals(const FieldVector<ctype, mydim> &neighborlocal) const
      {
        return this->realIterator.dereference().targetLocals(neighborlocal);
      }


      /** \brief obtain the type of reference element for this intersection */
      GeometryType type() const
      {
        return this->realIterator.dereference().type();
      }

    };


  public:

    typedef DomainIterator RemoteIntersectionDomainIterator;

    typedef TargetIterator RemoteIntersectionTargetIterator;

    typedef IndependentIterator RemoteIntersectionIterator;

  };

} // end namespace RemoteIntersectionInterface

#endif // REMOTEINTERSECTIONITERATORS_HH_
