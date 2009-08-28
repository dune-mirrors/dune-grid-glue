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
 * @brief Simple class for transformations of vertices.
 */

#ifndef COORDINATETRANSFORMATION_HH
#define COORDINATETRANSFORMATION_HH

/** \brief Base class for vertex transformations
    \tparam gridDimworld The dimension of the grid world
    \tparam couplingDimworld The dimension of the coupling world.  This may be larger than gridDimworld
    \tparam ctype The type used for coordinates
 */
template<int gridDimworld, int couplingDimworld, typename ctype = double>
class CoordinateTransformation
{
public:

  /** \brief Map a point to a new position */
  virtual Dune::FieldVector<ctype, couplingDimworld> operator()(const Dune::FieldVector<ctype, gridDimworld>& c) const = 0;

};

#endif // COORDINATETRANSFORMATION_HH_
