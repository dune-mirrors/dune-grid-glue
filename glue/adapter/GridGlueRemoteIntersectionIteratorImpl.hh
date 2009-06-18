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

template<typename GET1, typename GET2, typename SM>
class GridGlue<GET1, GET2, SM>::RemoteIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2, SM> Parent;

  //	typedef RemoteIntersectionTempl<Parent>  RemoteIntersectionImpl;


public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              _glue;

  mutable RemoteIntersection _intersection;

  int _index;


public:

  RemoteIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_._index)
  {}


  //	GridGlueDomainIntersectionIterator& operator=(RemoteIntersectionImpl& intersectionImpl_)
  //	{
  //		this->_index = intersectionImpl_._index;
  //		this->_glue = intersectionImpl_._glue;
  //		this->_intersection = RemoteIntersection(intersectionImpl_);
  //	}

  RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION._index)
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_index < static_cast<int>(this->_glue->_sm.nSimplices()))
      this->_intersection = this->_glue->_intersections[this->_index];
    else
      this->_index = this->_glue->NULL_INTERSECTION._index;
  }


  bool equals(const RemoteIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
  }

};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   DOMAIN INTERSECTION ITERATOR IMPL   */

template<typename GET1, typename GET2, typename SM>
class GridGlue<GET1, GET2, SM>::DomainIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2, SM> Parent;

  //	typedef RemoteIntersectionTempl<Parent>  RemoteIntersectionImpl;


public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              _glue;

  mutable RemoteIntersection _intersection;

  int _index;

  unsigned int _domain_parent;

  unsigned int _current;

  std::vector<unsigned int>  _parts;


public:

  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_._index),
      _domain_parent(0),
      _current(0)
  {
    if (this->_index < 0 || static_cast<int>(this->_glue->_sm.nSimplices()) <= this->_index)
      return;

    this->_domain_parent = this->_glue->_sm.domainParent(this->_index);
    this->_parts.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->_parts.begin());
  }


  DomainIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_._index),
      _domain_parent(0),
      _current(0)
  {
    if (this->_index < 0 || static_cast<int>(this->_glue->_sm.nSimplices()) <= this->_index)
      return;

    this->_domain_parent = this->_glue->_sm.domainParent(this->_index);
  }


  //	GridGlueDomainIntersectionIterator& operator=(RemoteIntersectionImpl& intersectionImpl_)
  //	{
  //		this->_index = intersectionImpl_._index;
  //		this->_glue = intersectionImpl_._glue;
  //		this->_intersection = RemoteIntersection(intersectionImpl_);
  //	}

  RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION._index)
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_current < this->_parts.size())
    {
      const RemoteIntersectionImpl& next = this->_glue->_intersections[this->_parts[this->_current]];
      this->_index = next._index;
      this->_intersection = next;
    }
    else
      this->_index = this->_glue->NULL_INTERSECTION._index;
  }


  bool equals(const DomainIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
  }

};



/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   TARGET INTERSECTION ITERATOR IMPL   */



template<typename GET1, typename GET2, typename SM>
class GridGlue<GET1, GET2, SM>::TargetIntersectionIteratorImpl
{
private:

  typedef GridGlue<GET1, GET2, SM> Parent;


public:

  typedef typename Parent::RemoteIntersection RemoteIntersection;


private:

  const Parent*              _glue;

  mutable RemoteIntersection _intersection;

  int _index;

  unsigned int _target_parent;

  unsigned int _current;

  std::vector<unsigned int>  _parts;


public:

  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_, const std::vector<unsigned int>& parts_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_._index),
      _target_parent(0),
      _current(0)
  {
    if (this->_index < 0 || static_cast<int>(this->_glue->_sm.nSimplices()) <= this->_index)
      return;

    this->_target_parent = this->_glue->_sm.targetParent(this->_index);
    this->_parts.resize(parts_.size());
    copy(parts_.begin(), parts_.end(), this->_parts.begin());
  }


  TargetIntersectionIteratorImpl(RemoteIntersectionImpl& intersectionImpl_)
    : _glue(intersectionImpl_._glue),
      _intersection(intersectionImpl_),
      _index(intersectionImpl_._index),
      _target_parent(0),
      _current(0)
  {
    if (this->_index < 0 || static_cast<int>(this->_glue->_sm.nSimplices()) <= this->_index)
      return;

    this->_target_parent = this->_glue->_sm.targetParent(this->_index);
  }

  //	GridGlueTargetIntersectionIterator& operator=(RemoteIntersectionImpl& intersectionImpl_)
  //	{
  //		this->_index = intersectionImpl_._index;
  //		this->_glue = intersectionImpl_._glue;
  //		this->_intersection = RemoteIntersection(intersectionImpl_);
  //	}

public:

  RemoteIntersection& dereference() const
  {
    if (this->_index == this->_glue->NULL_INTERSECTION._index)
      DUNE_THROW(Dune::GridError, "dereferencing end iterator");
    return this->_intersection;
  }


  void increment()
  {
    if (++this->_current < this->_parts.size())
    {
      const RemoteIntersectionImpl& next = this->_glue->_intersections[this->_parts[this->_current]];
      this->_index = next._index;
      this->_intersection = next;
    }
    else
      this->_index = this->_glue->NULL_INTERSECTION._index;
  }


  bool equals(const TargetIntersectionIteratorImpl& iter) const
  {
    return iter._index == this->_index;
  }

};


#endif // GRIDGLUEREMOTEINTERSECTIONITERATORIMPL_HH_
