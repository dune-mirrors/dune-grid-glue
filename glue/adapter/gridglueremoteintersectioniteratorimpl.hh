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


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   REMOTE INTERSECTION ITERATOR IMPL   */

template<typename GET1, typename GET2>
class GridGlue<GET1, GET2>::RemoteIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2> Parent;

public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              _glue;

  RemoteIntersection _intersection;

  unsigned int _index;


public:

  RemoteIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_.index())
  {}


  const RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_index < this->_glue->_index_sz)
      this->_intersection = this->_glue->_intersections[this->_index];
    else
      this->_index = this->_glue->NULL_INTERSECTION.index();
  }


  bool equals(const RemoteIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
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

  const Parent*              _glue;

  RemoteIntersection _intersection;

  unsigned int _index;

  unsigned int _domain_parent;

  unsigned int _current;

  std::vector<unsigned int>  _parts;


public:

  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_.index()),
      _domain_parent(0),
      _current(0)
  {
    if (this->_index < 0 || this->_glue->_index_sz <= this->_index)
      return;

    this->_domain_parent = this->_glue->_merg->domainParent(this->_index);
    this->_parts.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->_parts.begin());
  }


  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_.index()),
      _domain_parent(0),
      _current(0)
  {
    if (this->_index < 0 || this->_glue->_index_sz <= this->_index)
      return;

    this->_domain_parent = this->_glue->_merg->domainParent(this->_index);
  }


  const RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_current < this->_parts.size())
    {
      const RemoteIntersectionImpl& next = this->_glue->_intersections[this->_parts[this->_current]];
      this->_index = next.index();
      this->_intersection = next;
    }
    else
      this->_index = this->_glue->NULL_INTERSECTION.index();
  }


  bool equals(const DomainIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
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

  const Parent*              _glue;

  RemoteIntersection _intersection;

  unsigned int _index;

  unsigned int _target_parent;

  unsigned int _current;

  std::vector<unsigned int>  _parts;


public:

  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_.index()),
      _target_parent(0),
      _current(0)
  {
    if (this->_index < 0 || this->_glue->_index_sz <= this->_index)
      return;

    this->_target_parent = this->_glue->_merg->targetParent(this->_index);
    this->_parts.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->_parts.begin());
  }


  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_.index()),
      _target_parent(0),
      _current(0)
  {
    if (this->_index < 0 || this->_glue->_index_sz <= this->_index)
      return;

    this->_target_parent = this->_glue->_merg->targetParent(this->_index);
  }

public:

  const RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION.index())
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_current < this->_parts.size())
    {
      const RemoteIntersectionImpl& next = this->_glue->_intersections[this->_parts[this->_current]];
      this->_index = next.index();
      this->_intersection = next;
    }
    else
      this->_index = this->_glue->NULL_INTERSECTION.index();
  }


  bool equals(const TargetIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
  }

};


#endif // GRIDGLUEREMOTEINTERSECTIONITERATORIMPL_HH_
