// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    box.hh
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
 * @file box.hh
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

  Box(const C& lower, const C& upper) : lower_(lower), upper_(upper)
  {
    for (int i = 0; i < dim; ++i)
      center_[i] = 0.5*(upper[i]+lower[i]);
  }

  Box(const Box& b) : lower_(b.lower_), upper_(b.upper_)
  {
    for (int i = 0; i < dim; ++i)
      center_[i] = 0.5*(upper_[i]+lower_[i]);
  }

  bool contains(const C& c) const
  {
    for (int i = 0; i < dim; ++i)
      if (c[i] < this->lower_[i] || c[i] >= this->upper_[i])
        return false;
    return true;
  }

  bool intersects(const Box& b)
  {
    for (int i = 0; i < dim; ++i)
      if (this->lower_[i] >= b.upper_[i] || b.lower_[i] >= this->upper_[i])
        return false;
    return true;
  }

  const C& center() const
  {
    return this->center_;
  }


  double size(int i)
  {
    return upper_[i]-lower_[i];
  }

  const C& lower() const
  {
    return lower_;
  }

  const C& upper() const
  {
    return upper_;
  }

private:
  C lower_;
  C upper_;
  C center_;
};

#ifdef Vector3r_h

template<>
class Box<hxa7241_graphics::Vector3r, 3>
{
private:

  typedef hxa7241_graphics::Vector3r C;

public:
  Box(const C& lower, const C& upper) :
    lower_(lower), upper_(upper), center_(0.5*(upper[0]+lower[0]), 0.5*(upper[1]+lower[1]), 0.5*(upper[2]+lower[2]))
  {}

  Box(const Box& b) :
    lower_(b.lower_), upper_(b.upper_), center_(0.5*(upper[0]+lower[0]), 0.5*(upper[1]+lower[1]), 0.5*(upper[2]+lower[2]))
  {}

  bool contains(const C& c) const
  {
    for (int i = 0; i < 3; ++i)
      if (c[i] < this->lower_[i] || c[i] >= this->upper_[i])
        return false;
    return true;
  }

  bool intersects(const Box& b)
  {
    for (int i = 0; i < 3; ++i)
      if (this->lower_[i] >= b.upper_[i] || b.lower_[i] >= this->upper_[i])
        return false;
    return true;
  }

  const C& center() const
  {
    return this->center_;
  }

  double size(int i) const
  {
    return upper_[i]-lower_[i];
  }

  const C& lower_() const
  {
    return this->lower_;
  }

  const C& upper_() const
  {
    return this->upper_;
  }

private:
  C lower_;
  C upper_;
  C center_;
};
#endif


#endif // BOX_HH_
