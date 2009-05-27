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

#include "GridMeshTraits.hh"



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
    enum { dimensionsurface = T::dimS };
    this->mesh = T::mesh;
  }

private:

  MeshClassification::MeshType mesh;
};


/**
 * @class DefaultExtractionTraits
 * @brief default configuration of surface extraction traits,
 * deduced from type of grid and descriptor
 *
 * This can be overloaded in order to change single components.
 * But actually such a small class is rewritten pretty quickly...
 * (efficieny!)
 * !!IMPORTANT!!
 * Use the third parameter with care! E.g. specifying "simplex" for a given grid for which a simplicial
 * extractor would have been chosen anyway results in ambiguous template specializations and thus
 * deducing the "right" specialization of the GridMeshTraits struct is not possible for the compiler!
 * Remove the explicit MeshType argument in that case and you're fine.
 */
template<typename GV, bool extract_mesh = false, MeshClassification::MeshType mtype = MeshClassification::hybrid>
struct DefaultExtractionTraits
{
  /// @brief the grid's type
  typedef GV GridView;

  /// @brief default is the extraction of a manifold
  enum { dimS = GridView::dimension-static_cast<int>(!extract_mesh)       };

  /// @brief determines the type of the extractor used, defaults
  /// to automatic determination through traits class (recommended).
  static const MeshClassification::MeshType mesh = GridMeshTraits<typename GridView::Grid, mtype>::mesh;
};


#endif // GRIDEXTRACTIONTRAITS_HH_
