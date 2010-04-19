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

/** \brief Base class for subentity-selecting predicates
    \tparam GV GridView that the subentities are extracted from
 */
template<typename GV, int codim>
class ExtractorPredicate
{
public:

  /** \brief Return true if a subentity should be extracted.
      \param element An element
      \param subentity Subentity number
   */
  virtual bool contains(const typename GV::Traits::template Codim<0>::EntityPointer& element, unsigned int subentity) const = 0;

  /** \brief Dummy virtual destructor */
  virtual ~ExtractorPredicate() {}
};

#endif // EXTRACTOR_PREDICATES_HH
