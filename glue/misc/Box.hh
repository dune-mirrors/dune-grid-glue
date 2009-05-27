// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    Box.hh
 *  Version:     1.0
 *  Created on:  Jan 13, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     liboctree
 *  Description: dimension independent boxes (cubes) used in MultiDimOctree
 *  subversion:  $Id$
 *
 */
/**
 * @file Box.hh
 * @brief dimension independent boxes (cubes) used in MultiDimOctree
 */

#ifndef BOX_HH_
#define BOX_HH_


template<typename C, int dim>
class Box
{
public:

  ~Box()
  {}

  Box(const C& lower, const C& upper) : _lower(lower), _upper(upper)
  {
    for (int i = 0; i < dim; ++i)
      _center[i] = 0.5*(upper[i]+lower[i]);
  }

  Box(const Box& b) : _lower(b._lower), _upper(b._upper)
  {
    for (int i = 0; i < dim; ++i)
      _center[i] = 0.5*(_upper[i]+_lower[i]);
  }

  bool contains(const C& c) const
  {
    for (int i = 0; i < dim; ++i)
      if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
        return false;
    return true;
  }

  bool intersects(const Box& b)
  {
    for (int i = 0; i < dim; ++i)
      if (this->_lower[i] >= b._upper[i] || b._lower[i] >= this->_upper[i])
        return false;
    return true;
  }

  const C& center() const
  {
    return this->_center;
  }


  double size(int i)
  {
    return _upper[i]-_lower[i];
  }

  const C& lower() const
  {
    return _lower;
  }

  const C& upper() const
  {
    return _upper;
  }

private:
  C _lower;
  C _upper;
  C _center;
};

#ifdef Vector3r_h

template<>
class Box<hxa7241_graphics::Vector3r, 3>
{
private:

  typedef hxa7241_graphics::Vector3r C;

public:
  Box(const C& lower, const C& upper) :
    _lower(lower), _upper(upper), _center(0.5*(upper[0]+lower[0]), 0.5*(upper[1]+lower[1]), 0.5*(upper[2]+lower[2]))
  {}

  Box(const Box& b) :
    _lower(b._lower), _upper(b._upper), _center(0.5*(upper[0]+lower[0]), 0.5*(upper[1]+lower[1]), 0.5*(upper[2]+lower[2]))
  {}

  bool contains(const C& c) const
  {
    for (int i = 0; i < 3; ++i)
      if (c[i] < this->_lower[i] || c[i] >= this->_upper[i])
        return false;
    return true;
  }

  bool intersects(const Box& b)
  {
    for (int i = 0; i < 3; ++i)
      if (this->_lower[i] >= b._upper[i] || b._lower[i] >= this->_upper[i])
        return false;
    return true;
  }

  const C& center() const
  {
    return this->_center;
  }

  double size(int i) const
  {
    return _upper[i]-_lower[i];
  }

  const C& _lower() const
  {
    return this->_lower;
  }

  const C& _upper() const
  {
    return this->_upper;
  }

private:
  C _lower;
  C _upper;
  C _center;
};
#endif


#endif // BOX_HH_
