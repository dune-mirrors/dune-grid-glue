// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    codim1extractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: base class for grid extractors extracting surface grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Surface grid extractor base class
 */

#ifndef DUNE_CODIM_1_EXTRACTOR_HH
#define DUNE_CODIM_1_EXTRACTOR_HH

#include "extractor.hh"

template<typename GV>
class Codim1Extractor : public Extractor<GV,1>
{
public:
  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  Codim1Extractor(const GV& gv)
    :  Extractor<GV,1>(gv)
  {}
};

#endif // DUNE_CODIM_1_EXTRACTOR_HH_
