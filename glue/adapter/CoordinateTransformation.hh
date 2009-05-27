// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    CoordinateTransformation.hh
 *  Version:     1.0
 *  Created on:  Mar 12, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Simple class for uniform transformations of vertices.
 *  subversion:  $Id$
 *
 */
/**
 * @file CoordinateTransformation.hh
 * @brief Simple class for uniform transformations of vertices.
 */

#ifndef COORDINATETRANSFORMATION_HH_
#define COORDINATETRANSFORMATION_HH_


template<int dim, typename ctype = double>
class CoordinateTransformation
{
public:
  typedef Dune::FieldVector<ctype, dim> Coords;

  virtual const Coords operator()(const Coords& c) const = 0;

};

#endif // COORDINATETRANSFORMATION_HH_
