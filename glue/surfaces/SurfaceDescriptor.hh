// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SurfaceDescriptor.hh
 *  Version:     1.0
 *  Created on:  Mar 10, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: simple uniform descriptor for surface or mesh parts
 *  subversion:  $Id$
 *
 */
/**
 * @file SurfaceDescriptor.hh
 * @brief simple uniform descriptor for surface or mesh parts
 */

#ifndef SURFACEDESCRIPTOR_HH_
#define SURFACEDESCRIPTOR_HH_


template<typename GV>
class ElementDescriptor
{
public:

  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr) const = 0;

  virtual ~ElementDescriptor() {}
};


template<typename GV>
class FaceDescriptor
{
public:

  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr, unsigned int face) const = 0;

  virtual ~FaceDescriptor() {}
};


#endif // SURFACEDESCRIPTOR_HH_
