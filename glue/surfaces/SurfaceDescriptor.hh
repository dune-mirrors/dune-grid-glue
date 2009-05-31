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


/** \brief Base class for element-selecting predicates
    \tparam GV GridView that the elements are selected from
 */
template<typename GV>
class ElementDescriptor
{
public:
  /** \brief Return true if the element pointed to be eptr should be extracted */
  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr) const = 0;

  /** \brief Dummy virtual destructor */
  virtual ~ElementDescriptor() {}
};


/** \brief Base class for face-selecting predicates
    \tparam GV GridView that the faces are selected from
 */
template<typename GV>
class FaceDescriptor
{
public:

  /** \brief Return true if a face should be extracted.  Faces are specified through an element and a face number
     \param eptr An element
     \param face Face number
   */
  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr, unsigned int face) const = 0;

  /** \brief Dummy virtual destructor */
  virtual ~FaceDescriptor() {}
};


#endif // SURFACEDESCRIPTOR_HH_
