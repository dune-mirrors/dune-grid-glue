// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    ExtractorSelector.hh
 *  Version:     1.0
 *  Created on:  Feb 19, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: classes implementing automatic selection of suitable grid extractors for different grids through specialization
 *  subversion:  $Id$
 *
 */
/**
 * @file ExtractorSelector.hh
 * @brief classes implementing automatic selection of suitable grid extractors for different grids through specialization
 */

#ifndef EXTRACTORSELECTOR_HH
#define EXTRACTORSELECTOR_HH

#include "gridextractiontraits.hh"
#include "codim0extractor.hh"
#include "codim1extractor.hh"
#include "parallelextractor.hh"


/**
 * @brief automatic deduction of a suitable extractor for the configuration
 * given in the extraction traits class
 *
 * \tparam GSET Traits class which is a model of
 * SurfaceExtractionTraitsConcept. Given the configuration in this traits class
 * a suitable extractor is chosen and exported in the typedef ExtractorType.
 */
template<typename GSET>
struct ExtractorSelector
{
private:

  typedef ExtractorSelector<GSET> This;

  /**
   * @class Helper
   * @brief helper class to ExtractorSelector
   *
   * Actually the whole selection mechanism is implemented in this class which is only
   * accessible by the selector class.
   * Selection is executed by partial specialization. All valid combinations of the
   * world dimension (coordinate dimension), the grid's dimension and the dimension of
   * the extracted surface are treated. Only thing left is to implement classes that
   * can do the extraction for all of these combinations.
   */
  template<typename LSET, int codimension>
  struct Helper
  {};

  template<typename EXTR, bool parallel>
  struct ParallelHelper
  {};

public:

  /// @brief export the type of the deduced extractor
  typedef
  typename ParallelHelper<
      typename Helper<GSET, GSET::codimension>::ExtractorType,
      GSET::parallel
      >::ExtractorType
  ExtractorType;
};

/*   S P E C I A L I Z A T I O N   F O R   P A R A L L E L    E X T R A C T O R   */
/*   --------------------------------------------------------------------------   */

template<typename GSET>
template<typename EXTR>
struct ExtractorSelector<GSET>::ParallelHelper<EXTR, false>
{
  typedef EXTR ExtractorType;
};

template<typename GSET>
template<typename EXTR>
struct ExtractorSelector<GSET>::ParallelHelper<EXTR, true>
{
  typedef ParallelExtractor<EXTR> ExtractorType;
};

/*   S P E C I A L I Z A T I O N   F O R   C O D I M   1   E X T R A C T I O N   */
/*   --------------------------------------------------------------------------   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 1>
{
  typedef Codim1Extractor<typename LSET::GridView>  ExtractorType;
};



/*   S P E C I A L I Z A T I O N   F O R   C O D I M   0   E X T R A C T I O N   */
/*   --------------------------------------------------------------------   */


template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 0>
{
  typedef Codim0Extractor<typename LSET::GridView>  ExtractorType;
};

#endif // EXTRACTORSELECTOR_HH
