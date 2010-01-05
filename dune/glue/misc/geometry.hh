// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    geometry.h
 *  Version:     1.0
 *  Created on:  Jan 28, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: collection of useful geometry-related functions
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief header with helpful geometric computation related functions
 */

#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>


template<typename K>
Dune::FieldVector<K, 3> computeNormal(const Dune::array<Dune::FieldVector<K, 3>, 2>& v)
{
  Dune::FieldVector<K, 3> result(0.0);
  result[0] = v[0][1]*v[1][2] - v[0][2]*v[1][1];
  result[1] = v[0][2]*v[1][0] - v[0][0]*v[1][2];
  result[2] = v[0][0]*v[1][1] - v[0][1]*v[1][0];
  return result;
}


template<typename GEO>
void printGeometry(const GEO &geo, const char *name)
{
  std::cout << name << " :";
  for (int i = 0; i < geo.corners(); ++i)
    STDOUT(" (" << geo.corner(i) << ")");
  std::cout << std::endl;
}


#endif // GEOMETRY_H_
