// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridMeshTraits.hh
 *  Version:     1.0
 *  Created on:  Feb 3, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: helper classes for the grid extractor selection mechanism
 *  subversion:  $Id$
 *
 */
/**
 * @file GridSurfaceExtractionTraits.hh
 * @brief helper classes for the grid extractor selection mechanism
 */

#ifndef GRIDMESHTRAITS_HH_
#define GRIDMESHTRAITS_HH_

// always include those that are shipped with dune
#include <dune/grid/sgrid.hh>
#include <dune/grid/onedgrid.hh>
#include <dune/grid/yaspgrid.hh>

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#if HAVE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#if HAVE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif





/// @brief enum helper to identify a grid's mesh type

namespace MeshClassification
{
  enum MeshType
  {
    hybrid  = 0,
    simplex = 1,
    cube    = 2,
    rectangular = 3
  };
}


/**
 * @class GeneralSurfaceExtractionTraits
 * @brief traits class helping to find an appropriate surface extractor for a particular grid
 */
template<typename G, MeshClassification::MeshType mtype = MeshClassification::hybrid>
struct GridMeshTraits
{
  static const MeshClassification::MeshType mesh = mtype;
};


/*  S P E C I A L I Z A T I O N S   F O R   F O R C E D   S I M P L I C I A L   M E S H E S  */

/**
 * @class GridMeshTraits
 * Forces to treat the given grid as a simplicial mesh. This can be used
 * e.g. for potential multi-element meshes like UG grids.
 */
template<typename G>
struct GridMeshTraits<G, MeshClassification::simplex>
{
  static const MeshClassification::MeshType mesh = MeshClassification::simplex;
};


/*  S P E C I A L I Z A T I O N S   F O R   F O R C E D   C U B E   M E S H E S  */

/**
 * @class GridMeshTraits
 * Forces to treat the given grid as a simplicial mesh. This can be used
 * e.g. for potential multi-element meshes like UG grids.
 */
template<typename G>
struct GridMeshTraits<G, MeshClassification::cube>
{
  static const MeshClassification::MeshType mesh = MeshClassification::cube;
};


/*  S P E C I A L I Z A T I O N S   F O R   F O R C E D   R E C T A N G U L A R   C U B E   M E S H E S  */

/**
 * @class GridMeshTraits
 * Forces to treat the given grid as a simplicial mesh. This can be used
 * e.g. for potential multi-element meshes like UG grids.
 */
template<typename G>
struct GridMeshTraits<G, MeshClassification::rectangular>
{
  static const MeshClassification::MeshType mesh = MeshClassification::rectangular;
};


/*  S P E C I A L I Z A T I O N S   F O R   A L U  G R I D  */

#if HAVE_ALUGRID

/**
 * @class GridMeshTraits
 * @brief for ALU grid triangle meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::ALUSimplexGrid<2, 2>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::simplex;
};


/**
 * @class GridMeshTraits
 * @brief for ALU grid tetrahedral meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::ALUSimplexGrid<3, 3>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::simplex;
};


/**
 * @class GridMeshTraits
 * @brief for ALU grid quadrilateral meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::ALUCubeGrid<2, 2>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::cube;
};


/**
 * @class GridMeshTraits
 * @brief for ALU grid hexahedral meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::ALUCubeGrid<3, 3>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::cube;
};

#endif


/*  S P E C I A L I Z A T I O N S   F O R   A L B E R T A  G R I D  */

#if HAVE_ALBERTA

/**
 * @class GridMeshTraits
 * @brief for Alberta grid triangle meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::AbertaGrid<2, 2>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::simplex;
};


/**
 * @class GridMeshTraits
 * @brief for Alberta grid tetrahedral meshes
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::AbertaGrid<3, 3>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::simplex;
};

#endif


/*  S P E C I A L I Z A T I O N   F O R   S  G R I D  */

/**
 * @class GridMeshTraits
 * @brief for dune's SGrid which uses quadrilaterals
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::SGrid<2, 2>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::rectangular;
};


/**
 * @class GridMeshTraits
 * @brief for dune's SGrid which uses hexahedra
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::SGrid<3, 3>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::rectangular;
};


/*  S P E C I A L I Z A T I O N   F O R   Y A S P  G R I D  */

/**
 * @class GridMeshTraits
 * @brief for YaspGrid which uses quadrilaterals
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::YaspGrid<2>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::cube;
};


/**
 * @class GridMeshTraits
 * @brief for YaspGrid which uses hexahedra
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::YaspGrid<3>, mtype>
{
  static const MeshClassification::MeshType mesh = MeshClassification::cube;
};


/*  S P E C I A L I Z A T I O N   F O R   O N E  D  G R I D  */

/**
 * @class GridMeshTraits
 * @brief for the one dimension grid surface extraction is trivial
 */
template<MeshClassification::MeshType mtype>
struct GridMeshTraits<Dune::OneDGrid, mtype>
{};


#endif // GRIDMESHTRAITS_HH_
