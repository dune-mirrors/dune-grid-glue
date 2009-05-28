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

#ifndef EXTRACTORSELECTOR_HH_
#define EXTRACTORSELECTOR_HH_

#ifdef GRID_GLUE_USE_CONCEPTS
#include "../misc/conceptchecking.hh"
#include "GridExtractor.hh"
#endif
#include "GridExtractionTraits.hh"
#include "SimplicialMeshExtractor.hh"
#include "SimplicialManifoldExtractor.hh"
#include "SimplicialSurfaceExtractor.hh"
#include "CubeSurfaceExtractor.hh"
#include "CubeManifoldExtractor.hh"
#include "CubeMeshExtractor.hh"


/**
 * @class ExtractorSelector
 * @brief automatic deduction of a suitable extractor for the configuration
 * given in the extraction traits class
 *
 * The template parameter SET stands for a traits class which is a model of
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
   * Note:
   * So far only the simplicial extractors for surfaces and meshes are fully implemented.
   * For all other configurations the helper class simply does not export an extractor
   * type, so the compiler will notice and fail.
   * Also note that only using the specialized versions of this helper struct will lead to success.
   */
  template<typename LSET, int dimW, int dimG, int dimS, MeshClassification::MeshType mtype>
  struct Helper
  {};


public:

  /// @brief export the type of the deduced extractor
#ifdef GRID_GLUE_USE_CONCEPTS
  typedef GridExtractor<typename Helper<GSET, GSET::GridView::dimensionworld, GSET::GridView::dimension, GSET::dimS, GSET::mesh>::ExtractorType>  ExtractorType;
#else
  typedef typename Helper<GSET, GSET::GridView::dimensionworld, GSET::GridView::dimension, GSET::dimS, GSET::mesh>::ExtractorType ExtractorType;
#endif
};


/*   S P E C I A L I Z A T I O N   F O R   S U R F A C E   E X T R A C T I O  N   */
/*   --------------------------------------------------------------------------   */

/*   Surface in 3D is a 2D manifold   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 3, 2, MeshClassification::simplex>
{
  typedef SimplicialSurfaceExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 3, 2, MeshClassification::cube>
{
  typedef CubeSurfaceExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 3, 2, MeshClassification::rectangular>
{
  typedef CubeSurfaceExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 3, 2, MeshClassification::hybrid>
{
  //	typedef GeneralSurfaceExtractor<typename LSET::GridView>  ExtractorType;
};


/*   Surface in 2D is a 1D manifold   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 1, MeshClassification::simplex>
{
  typedef SimplicialSurfaceExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 1, MeshClassification::cube>
{
  typedef CubeSurfaceExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 1, MeshClassification::rectangular>
{
  typedef CubeSurfaceExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 1, MeshClassification::hybrid>
{
  //	typedef GeneralSurfaceExtractor<typename LSET::GridView>  ExtractorType;
};



/*   S P E C I A L I Z A T I O N   F O R   M E S H   E X T R A C T I O  N   */
/*   --------------------------------------------------------------------   */


/*   2D mesh is extracted and interpreted as 2D manifold in 3D   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 2, MeshClassification::simplex>
{
  typedef SimplicialMeshExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 2, MeshClassification::cube>
{
  typedef CubeMeshExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 2, MeshClassification::rectangular>
{
  typedef CubeMeshExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 2, 2, MeshClassification::hybrid>
{
  //		typedef GeneralMeshExtractor<typename LSET::GridView>  ExtractorType;
};


/*   2D mesh is extracted and interpreted as 2D manifold in 3D   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 1, 1, 1, MeshClassification::simplex>
{
  typedef SimplicialMeshExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 1, 1, 1, MeshClassification::cube>
{
  typedef CubeMeshExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 1, 1, 1, MeshClassification::rectangular>
{
  typedef CubeMeshExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 1, 1, 1, MeshClassification::hybrid>
{
  //		typedef GeneralMeshExtractor<typename LSET::GridView>  ExtractorType;
};



/*   S P E C I A L I Z A T I O N   F O R   M A N I F O L D   E X T R A C T I O  N   */
/*   ----------------------------------------------------------------------------   */


/*   2D manifold in 3D is directly extracted   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 2, 2, MeshClassification::simplex>
{
  typedef SimplicialManifoldExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 2, 2, MeshClassification::cube>
{
  typedef CubeManifoldExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 2, 2, MeshClassification::rectangular>
{
  typedef CubeManifoldExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 3, 2, 2, MeshClassification::hybrid>
{
  //		typedef GeneralManifoldExtractor<typename LSET::GridView>  ExtractorType;
};


/*   1D manifold in 2D is directly extracted   */

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 1, 1, MeshClassification::simplex>
{
  typedef SimplicialManifoldExtractor<typename LSET::GridView>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 1, 1, MeshClassification::cube>
{
  typedef CubeManifoldExtractor<typename LSET::GridView, false>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 1, 1, MeshClassification::rectangular>
{
  typedef CubeManifoldExtractor<typename LSET::GridView, true>  ExtractorType;
};

template<typename GSET>
template<typename LSET>
struct ExtractorSelector<GSET>::Helper<LSET, 2, 1, 1, MeshClassification::hybrid>
{
  //		typedef GeneralManifoldExtractor<typename LSET::GridView>  ExtractorType;
};


#endif // EXTRACTORSELECTOR_HH_
