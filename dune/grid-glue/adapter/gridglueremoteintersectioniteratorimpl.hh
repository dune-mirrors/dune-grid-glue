// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueRemoteIntersectionIteratorImpl.hh
 *  Version:     1.0
 *  Created on:  Mar 2, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Implementations of iterators over remote intersections provided by GridGlue.
 *  subversion:  $Id$
 *
 */
/**
 * @file GridGlueRemoteIntersectionIteratorImpl.hh
 * @brief Implementations of iterators over remote intersections provided by GridGlue.
 */

#ifndef GRIDGLUEREMOTEINTERSECTIONITERATORIMPL_HH_
#define GRIDGLUEREMOTEINTERSECTIONITERATORIMPL_HH_

#include <dune/grid-glue/adapter/gridglue.hh>

/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   REMOTE INTERSECTION ITERATOR IMPL   */

template<typename GET1, typename GET2>
class GridGlue<GET1, GET2>::RemoteIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2> Parent;

public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              glue_;

  RemoteIntersection intersection_;

  unsigned int index_;


public:

  RemoteIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : glue_(intersectionImpl_.glue_),
      intersection_(intersectionImpl_),
      index_(intersectionImpl_.index())
  {}


  const RemoteIntersection& dereference() const
  {
    if (this->index_ == this->glue_->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->intersection_;
  }


  void increment()
  {
    if (++this->index_ < this->glue_->index__sz)
      this->intersection_ = this->glue_->intersections_[this->index_];
    else
      this->index_ = this->glue_->NULL_INTERSECTION.index();
  }


  bool equals(const RemoteIntersectionIteratorImpl& iter) const
  {
    return iter.index_ == this->index_;
  }

};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   DOMAIN INTERSECTION ITERATOR IMPL   */

template<typename GET1, typename GET2>
class GridGlue<GET1, GET2>::DomainIntersectionIteratorImpl
{
private:
  typedef GridGlue<GET1, GET2> Parent;

public:
  typedef typename Parent::RemoteIntersection RemoteIntersection;

private:
  const Parent*              glue_;
  RemoteIntersection intersection_;
  unsigned int index_;
  unsigned int parent_id_;
  unsigned int current_;
  std::vector<unsigned int>  parts_;

public:

  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : glue_(intersectionImpl_.glue_),
      intersection_(intersectionImpl_),
      index_(intersectionImpl_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<0>(this->index_);
    this->parts_.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->parts_.begin());
  }


  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : glue_(intersectionImpl_.glue_),
      intersection_(intersectionImpl_),
      index_(intersectionImpl_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<0>(this->index_);
  }


  const RemoteIntersection& dereference() const
  {
    if (this->index_ == this->glue_->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->intersection_;
  }


  void increment()
  {
    if (++this->current_ < this->parts_.size())
    {
      const RemoteIntersectionImpl& next = this->glue_->intersections_[this->parts_[this->current_]];
      this->index_ = next.index();
      this->intersection_ = next;
    }
    else
      this->index_ = this->glue_->NULL_INTERSECTION.index();
  }


  bool equals(const DomainIntersectionIteratorImpl& iter) const
  {
    return iter.index_ == this->index_;
  }

};



/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   TARGET INTERSECTION ITERATOR IMPL   */



template<typename GET1, typename GET2>
class GridGlue<GET1, GET2>::TargetIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2> Parent;


public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              glue_;

  RemoteIntersection intersection_;

  unsigned int index_;

  unsigned int parent_id_;

  unsigned int current_;

  std::vector<unsigned int>  parts_;


public:

  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : glue_(intersectionImpl_.glue_),
      intersection_(intersectionImpl_),
      index_(intersectionImpl_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<1>(this->index_);
    this->parts_.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->parts_.begin());
  }


  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : glue_(intersectionImpl_.glue_),
      intersection_(intersectionImpl_),
      index_(intersectionImpl_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<1>(this->index_);
  }

public:

  const RemoteIntersection& dereference() const
  {
    if (this->index_ == this->glue_->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->intersection_;
  }


  void increment()
  {
    if (++this->current_ < this->parts_.size())
    {
      const RemoteIntersectionImpl& next = this->glue_->intersections_[this->parts_[this->current_]];
      this->index_ = next.index();
      this->intersection_ = next;
    }
    else
      this->index_ = this->glue_->NULL_INTERSECTION.index();
  }


  bool equals(const TargetIntersectionIteratorImpl& iter) const
  {
    return iter.index_ == this->index_;
  }

};


#endif // GRIDGLUEREMOTEINTERSECTIONITERATORIMPL_HH_
