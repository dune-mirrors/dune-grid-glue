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
#include <dune/glue/common/conceptchecking.hh>
#include "remoteintersectionconcepts.hh"
#endif
#include "remoteintersection.hh"

namespace RemoteIntersectionInterface
{

  /** \brief Iterator over remote intersections */
  template<typename RemoteIntersectionImpl,
      typename RemoteIntersectionIteratorImpl>
  class RemoteIntersectionIterator
  {
  private:
    /*   C H E C K   C O N C E P T S   */

#ifdef GRID_GLUE_USE_CONCEPTS
    typedef RemoteIntersectionImpl RemoteIntersectionImplType;
    CLASS_REQUIRE(RemoteIntersectionImplType, RemoteIntersectionConcept);

    typedef RemoteIntersectionIteratorImpl RemoteIntersectionIteratorImplType;
    CLASS_REQUIRE(RemoteIntersectionIteratorImplType, RemoteIntersectionIteratorConcept);
#else
    typedef RemoteIntersectionImpl RemoteIntersectionImplType;

    typedef RemoteIntersectionIteratorImpl RemoteIntersectionIteratorImplType;
#endif

  public:

    typedef RemoteIntersection<RemoteIntersectionImpl>  ThisRemoteIntersection;

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


    typedef typename RemoteIntersectionImpl::ctype ctype;


  private:

    RemoteIntersectionIteratorImplType realIterator;


  protected:

    RemoteIntersectionIteratorImplType& getImpl()
    {
      return this->realIterator;
    }


  public:

    /**
     * @brief copy constructor
     * @param iter the iterator to copy
     */
    RemoteIntersectionIterator(const RemoteIntersectionIterator& iter) : realIterator(iter.realIterator)
    {}


    /**
     * @brief only constructor to build a "new" iterator
     * @param iter the iterator to initialize from
     */
    RemoteIntersectionIterator(const RemoteIntersectionIteratorImplType& iter) : realIterator(iter)
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
    RemoteIntersectionIterator& operator++()
    {
      this->realIterator.increment();
      return *this;
    }


    bool operator==(const RemoteIntersectionIterator& rhs) const
    {
      return rhs.realIterator.equals(this->realIterator);
    }


    bool operator!=(const RemoteIntersectionIterator& rhs) const
    {
      return !rhs.realIterator.equals(this->realIterator);
    }

  };

} // end namespace RemoteIntersectionInterface

#endif // REMOTEINTERSECTIONITERATORS_HH_
