// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridExtractor.hh
 *  Version:     1.0
 *  Created on:  Apr 27, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: abstract interface for grid extractor implementations of any kind
 *  subversion:  $Id$
 *
 */
/**
 * @file GridExtractor.hh
 * @brief abstract interface for grid extractor implementations of any kind
 */

#ifndef GRIDEXTRACTOR_HH_
#define GRIDEXTRACTOR_HH_


#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include "../../dune/common/exceptions.hh"
#include "../../dune/common/fvector.hh"
#include "../../dune/common/fixedarray.hh"
#include "../../dune/grid/common/geometry.hh"
#include "SurfaceDescriptor.hh"
#include "GridMeshTraits.hh"

#include "GridExtractorConcept.hh"
#include "../misc/conceptchecking.h"


using namespace Dune;
using namespace std;


namespace ExtractorClassification
{
  enum ExtractorType
  {
    surface = 0,
    mesh = 1,
    manifold = 2
  };
}


template<int dimG, int dimGW, int dimEW>
struct ExtractorClassifier {};

template<int dim>
struct ExtractorClassifier<dim, dim, dim>
{
  static const ExtractorClassification::ExtractorType type = ExtractorClassification::surface;
};

template<int hyperdim, int dimEW>
struct ExtractorClassifier<hyperdim, hyperdim, dimEW>
{
  static const ExtractorClassification::ExtractorType type = ExtractorClassification::mesh;
};

template<int dimGW, int hyperdim>
struct ExtractorClassifier<hyperdim, dimGW, hyperdim>
{
  static const ExtractorClassification::ExtractorType type = ExtractorClassification::manifold;
};



/**
 * @class GridExtractor
 * @brief provides static methods for grid surface extraction
 *
 * Provides methods that build topology information for given grids.
 * The template parameters
 * @li Implementation the actual extractor implementation
 */
template<typename Implementation>
class GridExtractor
{
private:

  typedef Implementation Impl;

  /*   C H E C K   C O N C E P T S   */

  CLASS_REQUIRE(Impl, GridExtractorConcept);


public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum { dimw = Impl::dimw };

  enum { dim = Impl::dim };

  /// @brief compile time number of corners of surface simplices
  enum { simplex_corners = Impl::simplex_corners };


  typedef typename Impl::GridView GridView;

  typedef typename Impl::ctype ctype;
  typedef typename Impl::Coords Coords;
  typedef typename Impl::SimplexTopology SimplexTopology;

  typedef typename Impl::VertexPtr VertexPtr;
  typedef typename Impl::Vertex Vertex;
  typedef typename Impl::VertexIter VertexIter;

  typedef typename Impl::ElementPtr ElementPtr;
  typedef typename Impl::Element Element;
  typedef typename Impl::ElementIter ElementIter;

  typedef typename Impl::IsIter IsIter;

  // index sets and index types
  typedef typename Impl::IndexSet IndexSet;
  typedef typename Impl::IndexType IndexType;


  /// @brief for barycentric coordinates
  typedef Dune::FieldVector<ctype, dimw>                      BCoords;

  /// @brief for local element coordinates
  typedef Dune::FieldVector<ctype, dim>                       ECoords;

  /// @brief for grid's world coordinates
  typedef Dune::FieldVector<ctype, GridView::dimensionworld>  WCoords;


  /// @brief a constant classifying this extractor
  static const ExtractorClassification::ExtractorType type;


private:
  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /// @brief the grid object to extract the surface from
  const GridView&                         gv;

  Impl imp;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  GridExtractor(const GridView& gv_) :
    gv(gv), imp(gv_)
  {
    STDOUTLN("This is GridExtractor on a <" << GridView::dimension << "," << GridView::dimensionworld << "> grid working in " << dimw << " space!");
  }


  /**
   * @brief default destructor, frees memory
   */
  ~GridExtractor() {};


  /*  F U N C T I O N A L I T Y  */

  /**
   * ASSUMPTION:
   * dim == dimw
   *
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a decider function
   * that returns a boolean value for each codim 0 entity. If an entity is "in" that
   * means that all of its boundary faces are considered in the surface.
   * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are simplices.
   * @param descr a descriptor that "selects" the elements whose faces are to add to the surface
   */
  template<typename Descriptor>
  void update(const Descriptor& descr)
  {
    imp.update(descr);
  }


  //	/**
  //	 * ASSUMPTION:
  //	 * dim == dimw
  //	 *
  //	 * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
  //	 * with the topology of the surface written to them. The description of the
  //	 * surface part that is to be extracted is given in form of a decider function
  //	 * that returns a boolean value for each codim 0 entity's every face.
  //	 * Note that the function is only called for faces that are part of the domain boundary.
  //	 * Inner faces cannot be part of the surface.
  //	 * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
  //	 *
  //	 * Assumed that we are in 2D the coords array will have the structure
  //	 * x0 y0 x1 y1 ... x(n-1) y(n-1)
  //	 * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
  //	 * index 1 is associated with the position x1. If we the surface consists of triangles
  //	 * we have always groups of 3 indices describing one triangle.
  //	 *
  //	 * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
  //	 * segments have are simplices.
  //	 * @param descr a descriptor that "selects" the faces to add to the surface
  //	 */
  //	void update(const FaceDescriptor<GV>& descr)
  //	{
  //		imp.update(descr);
  //	}


  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear()
  {
    imp.clear();
  }


  /*  S E T T E R S  */


  /*  G E T T E R S  */

  /**
   * @brief getter for the coordinates array
   * @param coords a vector that will be resized (!) and filled with the coordinates,
   * note that the single components are written consecutively
   */
  void getCoords(vector<FieldVector<ctype, dimw> >& coords) const
  {
    imp.getCoords(coords);
  }


  /**
   * @brief getter for the count of coordinates
   * @return the count
   */
  unsigned int nCoords() const
  {
    return imp.nCoords();
  }


  /**
   * @brief getter for the indices array
   * It is strongly recommended not to modify its contents.
   * Deallocation is done in this class.
   * @return the _indices array
   */
  void getFaces(vector<SimplexTopology>& faces) const
  {
    imp.getFaces(faces);
  }


  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return this->gv.indexSet();
  }


  /**
   * @brief getter for the index of an entity of codim cc
   * @return the index specified by the grid's index set
   */
  template<int cc>
  IndexType index(const typename GridView::Traits::template Codim<cc>::Entity& e) const
  {
    return this->indexSet().template index<cc>(e);
  }


  /**
   * @brief gets index of coordinate in _coords associated with given vertex
   * @return the index if possible, -1 else
   */
  int coordinateIndex(const Vertex& v) const
  {
    return imp.coordinateIndex(v);
  }


  /**
   * @brief gets index of first face as well as the total number of faces that
   * were extracted from this element
   * @param e the element
   * @param first will contain the first index if found, else -1
   * @param count will contain the number of faces if found, else 0
   * @return success
   */
  bool faceIndices(const Element& e, int& first, int& count) const
  {
    return imp.faceIndices(e, first, count);
  }


  /**
   * @brief gets the number face in the parent element
   * @param index the index of the face
   * @return if failed -1, else the index
   */
  int numberInSelf(unsigned int index) const
  {
    return imp.numberInSelf(index);
  }


  /**
   * @brief for given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  void globalCoords(unsigned int index, const BCoords &bcoords, WCoords &wcoords) const
  {
    imp.globalCoords(index, bcoords, wcoords);
  }


  /**
   * @brief for given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  void localCoords(unsigned int index, const BCoords &bcoords, ECoords &ecoords) const
  {
    imp.localCoords(index, bcoords, ecoords);
  }


  /**
   * @brief for given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const BCoords &bcoords, ECoords &ecoords, WCoords &wcoords) const
  {
    imp.localAndGlobalCoords(index, bcoords, ecoords, wcoords);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  template<typename CoordContainerB, typename CoordContainerW>
  void globalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerW &wcoords, int size) const
  {
    imp.globalCoords(index, bcoords, wcoords, size);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  template<typename CoordContainerB, typename CoordContainerE>
  void localCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, int size) const
  {
    imp.localCoords(index, bcoords, ecoords, size);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   * @return
   */
  template<typename CoordContainerB, typename CoordContainerE, typename CoordContainerW>
  void localAndGlobalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, CoordContainerW &wcoords, int size) const
  {
    imp.localAndGlobalCoords(index, bcoords, ecoords, wcoords, size);
  }


  /**
   * @brief gets for each vertex corner of given face (by index) the number of
   * the vertex in parent element's ordering
   * @param index the face's index
   * @param corner the index of the corner
   * @return -1 <=> index invalid or array not filled, else index
   */
  int numCornerInParent(unsigned int index, unsigned int corner) const
  {
    return imp.numCornerInParent(index, corner);
  }


  /**
   * @brief gets the parent element for a given face index,
   * throws an exception if index not valid
   * @param index the index of the face
   * @return a reference to the element's stored pointer
   */
  const ElementPtr& element(unsigned int index) const
  {
    return imp.element(index);
  }


  /**
   * @brief gets the vertex for a given coordinate index
   * throws an exception if index not valid
   * @param index the index of the coordinate
   * @return a reference to the vertex' stored pointer
   */
  const VertexPtr& vertex(unsigned int index) const
  {
    return imp.vertex(index);
  }

}; // end of class GridExtractor

template<typename Implementation>
const ExtractorClassification::ExtractorType GridExtractor<Implementation>::type = ExtractorClassifier<GridView::dimension, GridView::dimensionworld, dimw>::type;


#endif // GRIDEXTRACTOR_HH_
