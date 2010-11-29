// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueRemoteIntersectionIterator.hh
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
 * @file GridGlueRemoteIntersectionIterator.hh
 * @brief Implementations of iterators over remote intersections provided by GridGlue.
 */

#ifndef DUNE_GRIDGLUE_REMOTEINTERSECTIONITERATOR_HH
#define DUNE_GRIDGLUE_REMOTEINTERSECTIONITERATOR_HH

#include <dune/grid-glue/adapter/gridglue.hh>

/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   REMOTE INTERSECTION ITERATOR IMPL   */

template<typename P0, typename P1>
class GridGlue<P0, P1>::RemoteIntersectionIterator :
  public Dune::ForwardIteratorFacade< RemoteIntersectionIterator,
      const Dune::GridGlue::AbsoluteIntersection<P0,P1> >
{
private:

  typedef GridGlue<P0, P1> Parent;

public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              glue_;

  RemoteIntersection intersection_;

  unsigned int index_;


public:

  RemoteIntersectionIterator(RemoteIntersection& intersection_)
    : glue_(intersection_.glue_),
      intersection_(intersection_),
      index_(intersection_.index())
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


  bool equals(const RemoteIntersectionIterator& iter) const
  {
    return iter.index_ == this->index_;
  }

};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   DOMAIN INTERSECTION ITERATOR IMPL   */

template<typename P0, typename P1>
class GridGlue<P0, P1>::DomainIntersectionIterator :
  public Dune::ForwardIteratorFacade< DomainIntersectionIterator, const RemoteIntersection>
{
private:
  typedef GridGlue<P0, P1> Parent;

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

  DomainIntersectionIterator(RemoteIntersection& intersection_, const std::vector<unsigned int>& parts_)
    : glue_(intersection_.glue_),
      intersection_(intersection_),
      index_(intersection_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<0>(this->index_);
    this->parts_.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->parts_.begin());
  }


  DomainIntersectionIterator(RemoteIntersection& intersection_)
    : glue_(intersection_.glue_),
      intersection_(intersection_),
      index_(intersection_.index()),
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
      const RemoteIntersection& next = this->glue_->intersections_[this->parts_[this->current_]];
      this->index_ = next.index();
      this->intersection_ = next;
    }
    else
      this->index_ = this->glue_->NULL_INTERSECTION.index();
  }


  bool equals(const DomainIntersectionIterator& iter) const
  {
    return iter.index_ == this->index_;
  }

};



/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   TARGET INTERSECTION ITERATOR IMPL   */



template<typename P0, typename P1>
class GridGlue<P0, P1>::TargetIntersectionIterator :
  public Dune::ForwardIteratorFacade< TargetIntersectionIterator, const RemoteIntersection>
{
private:

  typedef GridGlue<P0, P1> Parent;


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

  TargetIntersectionIterator(RemoteIntersection& intersection_, const std::vector<unsigned int>& parts_)
    : glue_(intersection_.glue_),
      intersection_(intersection_),
      index_(intersection_.index()),
      parent_id_(0),
      current_(0)
  {
    if (this->index_ < 0 || this->glue_->index__sz <= this->index_)
      return;

    this->parent_id_ = this->glue_->merger_->template parent<1>(this->index_);
    this->parts_.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->parts_.begin());
  }


  TargetIntersectionIterator(RemoteIntersection& intersection_)
    : glue_(intersection_.glue_),
      intersection_(intersection_),
      index_(intersection_.index()),
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
      const RemoteIntersection& next = this->glue_->intersections_[this->parts_[this->current_]];
      this->index_ = next.index();
      this->intersection_ = next;
    }
    else
      this->index_ = this->glue_->NULL_INTERSECTION.index();
  }


  bool equals(const TargetIntersectionIterator& iter) const
  {
    return iter.index_ == this->index_;
  }

};


#endif // DUNE_GRIDGLUE_REMOTEINTERSECTIONITERATOR_HH
