// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridExtractorConcept.hh
 *  Version:     1.0
 *  Created on:  Apr 27, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Concept class for syntactical concept checkig of models of GridExtractor concept.
 *  subversion:  $Id$
 *
 */
/**
 * @file GridExtractorConcept.hh
 * @brief
 */

#ifndef GRIDEXTRACTORCONCEPT_HH_
#define GRIDEXTRACTORCONCEPT_HH_


#include <vector>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/genericgeometry/misc.hh>
#include "surfacedescriptor.hh"



/*
 * @class GridExtractorConcept
 * @brief checks a class' conformity with the extractor concept
 */
template<typename T>
class GridExtractorConcept
{
public:

  /*  R E Q U E S T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum { dimw = T::dimw };

  enum { dim = T::dim };

  enum { simplex_corners = T::simplex_corners };


  typedef typename T::GridView GridView;

  typedef typename T::ctype ctype;
  typedef typename T::Coords Coords;
  typedef typename T::SimplexTopology SimplexTopology;

  typedef typename T::VertexPtr VertexPtr;
  typedef typename T::Vertex Vertex;
  typedef typename T::VertexIter VertexIter;

  typedef typename T::ElementPtr ElementPtr;
  typedef typename T::Element Element;
  typedef typename T::ElementIter ElementIter;

  typedef typename T::IsIter IsIter;

  // index sets and index types
  typedef typename T::IndexSet IndexSet;
  typedef typename T::IndexType IndexType;


  void constraints()
  {

    /*
     * constructor based on given grid view
     * @param gv the grid view object to extract from
     */
    extractor = new T(gv);

    /*
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
    extractor->update(*ed);


    /*
     * ASSUMPTION:
     * dim == dimw
     *
     * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
     * with the topology of the surface written to them. The description of the
     * surface part that is to be extracted is given in form of a decider function
     * that returns a boolean value for each codim 0 entity's every face.
     * Note that the function is only called for faces that are part of the domain boundary.
     * Inner faces cannot be part of the surface.
     * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
     *
     * Assumed that we are in 2D the coords array will have the structure
     * x0 y0 x1 y1 ... x(n-1) y(n-1)
     * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
     * index 1 is associated with the position x1. If we the surface consists of triangles
     * we have always groups of 3 indices describing one triangle.
     *
     * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
     * segments have are simplices.
     * @param descr a descriptor that "selects" the faces to add to the surface
     */
    fdt.constraints(*extractor, *fd);


    /*
     * @brief delete everything build up so far and free the memory
     */
    extractor->clear();


    /*  G E T T E R S  */

    /*
     * @brief getter for the coordinates array
     * @param coords a vector that will be resized (!) and filled with the coordinates,
     * note that the single components are written consecutively
     */
    extractor->getCoords(coords);

    /*
     * @brief getter for the count of coordinates
     * @return the count
     */
    u = extractor->nCoords();

    /*
     * @brief getter for the indices array
     * It is strongly recommended not to modify its contents.
     * Deallocation is done in this class.
     * @return the _indices array
     */
    extractor->getFaces(faces);


    /*
     * @brief gets index of coordinate in _coords associated with given vertex
     * @return the index if possible, -1 else
     */
    i = extractor->coordinateIndex(v);



    /*
     * @brief gets index of first face as well as the total number of faces that
     * were extracted from this element
     * @param e the element
     * @param first will contain the first index if found, else -1
     * @param count will contain the number of faces if found, else 0
     * @return success
     */
    b = extractor->faceIndices(e, i, i);


    /*
     * @brief gets the number face in the parent element
     * @param index the index of the face
     * @return if failed -1, else the index
     */
    i = extractor->indexInInside(u);


    /*
     * @brief for given barycentric coords in a simplex compute world coordinates
     *
     * If both are to be computed, element and world coordinates, then use the
     * combined method for efficiency!
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param wcoords to be filled with world coordinates
     */
    extractor->globalCoords(u, bc, wc);


    /*
     * @brief for given barycentric coords in a simplex compute element coordinates
     *
     * If both are to be computed, element and world coordinates, then use the
     * combined method for efficiency!
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param ecoords to be filled with element coordinates
     */
    extractor->localCoords(u, bc, lc);


    /*
     * @brief for given barycentric coords in a simplex compute element and world coordinates
     *
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param ecoords to be filled with element coordinates
     * @param wcoords to be filled with world coordinates
     */
    extractor->localAndGlobalCoords(u, bc, lc, wc);


    /*
     * @brief for several given barycentric coords in a simplex compute world coordinates
     *
     * If both are to be computed, element and world coordinates, then use the
     * combined method for efficiency!
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param wcoords to be filled with world coordinates
     */
    extractor->globalCoords(u, bcs, wcs, u);


    /*
     * @brief for several given barycentric coords in a simplex compute element coordinates
     *
     * If both are to be computed, element and world coordinates, then use the
     * combined method for efficiency!
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param ecoords to be filled with element coordinates
     */
    extractor->localCoords(u, bcs, lcs, u);


    /*
     * @brief for several given barycentric coords in a simplex compute element and world coordinates
     *
     * @param index the index of the simplex
     * @param bcoords the barycentric coordinates
     * @param ecoords to be filled with element coordinates
     * @param wcoords to be filled with world coordinates
     * @return
     */
    extractor->localAndGlobalCoords(u, bcs, lcs, wcs, u);


    /*
     * @brief gets for each vertex corner of given face (by index) the number of
     * the vertex in parent element's ordering
     * @param index the face's index
     * @param corner the index of the corner
     * @return -1 <=> index invalid or array not filled, else index
     */
    i = extractor->numCornerInParent(u, u);

    /*
     * @brief gets the parent element for a given face index,
     * throws an exception if index not valid
     * @param index the index of the face
     * @return a reference to the element's stored pointer
     */
    eptr = extractor->element(u);


    /*
     * @brief gets the vertex for a given coordinate index
     * throws an exception if index not valid
     * @param index the index of the coordinate
     * @return a reference to the vertex' stored pointer
     */
    vptr = extractor->vertex(u);


    /*
     * @brief default destructor
     */
    delete extractor;
  }


private:


  /* ------------------------------------------------------------------------------------------------- */
  /* There is a little variation between extractors: FaceDescriptor only works with SurfaceExtractors! */
  /* So to pack it all in one test case...                                                             */
  template<bool>
  struct Dummy
  {
    void constraints(T &extractor, FaceDescriptor<GridView> &fd) {}
  };

  template<bool>
  struct TesterTemplate
  {
    void constraints(T &extractor, FaceDescriptor<GridView> &fd)
    {
      extractor.update(fd);
    }
  };

  typedef Dune::GenericGeometry::ProtectedIf<static_cast<int>(dim) == static_cast<int>(dimw), TesterTemplate, Dummy>  FaceDescriptorTester;
  /* ------------------------------------------------------------------------------------------------- */



  typedef Dune::FieldVector<ctype, dim>                      LocalCoords;
  typedef Dune::FieldVector<ctype, GridView::dimensionworld> WorldCoords;
  typedef Dune::FieldVector<ctype, dimw>                     BarycentricCoords;

  T                           *extractor;
  GridView gv;
  unsigned int u;
  int i;
  bool b;
  LocalCoords lc;
  WorldCoords wc;
  BarycentricCoords bc;
  ElementDescriptor<GridView> *ed;
  FaceDescriptor<GridView>    *fd;
  FaceDescriptorTester fdt;
  std::vector<BarycentricCoords>   coords;
  std::vector<SimplexTopology>     faces;
  Vertex v;
  VertexPtr vptr;
  Element e;
  ElementPtr eptr;
  Dune::array<LocalCoords, simplex_corners>       lcs;
  Dune::array<WorldCoords, simplex_corners>       wcs;
  Dune::array<BarycentricCoords, simplex_corners> bcs;


public:


}; // end of class GridExtractorConcept




#endif // GRIDEXTRACTORCOCPET_HH_
