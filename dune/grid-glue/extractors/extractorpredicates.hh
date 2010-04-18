// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    extractorpredicates.hh
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
 * @file
 * @brief Base class for predicates selecting the part of a grid to be extracted
 */

#ifndef EXTRACTOR_PREDICATES_HH
#define EXTRACTOR_PREDICATES_HH

/** \brief Base class for face-selecting predicates
    \tparam GV GridView that the faces are selected from
 */
template<typename GV, int codim>
class ExtractorPredicate
{
public:

  /** \brief Return true if a subentity should be extracted.
      \param eptr An element
      \param subentity Subentity number
   */
  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr, unsigned int subentity) const = 0;

  /** \brief Dummy virtual destructor */
  virtual ~ExtractorPredicate() {}
};

/** \brief Base class for element-selecting predicates
    \tparam GV GridView that the elements are selected from
 */
template<typename GV>
class ExtractorPredicate<GV,0>
{
public:

  /** \brief Return true if a subentity should be extracted.
      \param eptr An element
      \param subentity Subentity number
   */
  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& eptr) const = 0;

  /** \brief Dummy virtual destructor */
  virtual ~ExtractorPredicate() {}
};

#endif // EXTRACTOR_PREDICATES_HH
