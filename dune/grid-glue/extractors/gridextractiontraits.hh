// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridExtractionTraits.hh
 *  Version:     1.0
 *  Created on:  Feb 19, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: traits classes used to describe compile-time properties of a coupling scenario
 *  subversion:  $Id$
 *
 */
/**
 * @file GridExtractionTraits.hh
 * @brief traits classes used to describe compile-time properties of a coupling scenario
 */


#ifndef GRIDEXTRACTIONTRAITS_HH_
#define GRIDEXTRACTIONTRAITS_HH_

/**
 * @class GridExtractionTraitsConcept
 * @brief concept check class for extraction traits classes
 */
template<typename T>
struct GridExtractionTraitsConcept
{
public:
  void constraints()
  {
    typedef typename T::GridView TempGridView;
    enum { codimension = T::codimension };
    this->par  = T::parallel;
  }

private:

  bool par;
};


/**
 * @brief default configuration of extraction traits,
 * deduced from type of grid and extractor codimension.
 *
 * This can be overloaded in order to change single components.
 * But actually such a small class is rewritten pretty quickly...
 * (efficiency!)
 */
template<typename GV, int extractorCodim, bool par=false>
struct DefaultExtractionTraits
{
  /// @brief the grid's type
  typedef GV GridView;

  /// @brief Extractor codimension
  enum { codimension = extractorCodim };

  /// @brief is this a parallel mesh?
  static const bool parallel = par;
};


#endif // GRIDEXTRACTIONTRAITS_HH_
