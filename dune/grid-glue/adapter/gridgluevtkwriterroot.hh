// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueVtkWriterRoot.hh
 *  Version:     1.0
 *  Created on:  June, 2014
 *  Author:      Natalie Schroeder
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Class thought to make graphical debugging of couplings easier.
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Write all remote intersections to a vtk file for debugging
 */

#ifndef GRIDGLUEVTKWRITERROOT_HH
#define GRIDGLUEVTKWRITERROOT_HH


#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

#include <dune/common/classname.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/scsgmapper.hh>


/** \brief Write remote intersections to a vtk file for debugging purposes
 */
class GridGlueVtkWriterRoot
{

  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeExtractedPart(const Glue& glue, const std::string& filename)
  {
    dune_static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fgrid;

    fgrid.open(filename.c_str());

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename GridView::Traits::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPointer;
    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;

    typedef typename GridView::ctype ctype;

    const int dim = GridView::dimension;
    const int domdimw = GridView::dimensionworld;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
      coordinatePadding += " 0";

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << filename << "\nASCII" << std::endl;

    // WRITE POINTS
    // ----------------
    typedef typename Glue::Grid0Patch Extractor;
    std::vector<typename Extractor::Coords> coords;
    glue.template patch<side>().getCoords(coords);

    fgrid << ((dim==3) ? "DATASET UNSTRUCTURED_GRID" : "DATASET POLYDATA") << std::endl;
    fgrid << "POINTS " << coords.size() << " " << Dune::className<ctype>() << std::endl;

    for (size_t i=0; i<coords.size(); i++)
      fgrid << coords[i] << coordinatePadding << std::endl;

    fgrid << std::endl;

    // WRITE POLYGONS
    // ----------------

    std::vector<typename Extractor::VertexVector> faces;
    std::vector<Dune::GeometryType> geometryTypes;
    glue.template patch<side>().getFaces(faces);
    glue.template patch<side>().getGeometryTypes(geometryTypes);

    unsigned int faceCornerCount = 0;
    for (size_t i=0; i<faces.size(); i++)
      faceCornerCount += faces[i].size();

    fgrid << ((dim==3) ? "CELLS " : "LINES ")
          << geometryTypes.size() << " " << geometryTypes.size() + faceCornerCount << std::endl;

    for (size_t i=0; i<faces.size(); i++) {

      fgrid << faces[i].size();
      // vtk expects the vertices to by cyclically ordered
      // therefore unfortunately we have to deal with several element types on a case-by-case basis
      if (geometryTypes[i].isSimplex()) {
        for (int j=0; j<=dim; j++)
          fgrid << " " << faces[i][j];

      } else if (geometryTypes[i].isQuadrilateral()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2];

      } else if (geometryTypes[i].isPyramid()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2] << " " << faces[i][4];

      } else if (geometryTypes[i].isPrism()) {
        fgrid << " " << faces[i][0] << " " << faces[i][2] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][5] << " " << faces[i][4];

      } else if (geometryTypes[i].isHexahedron()) {
        fgrid << " " << faces[i][0] << " " << faces[i][1]
              << " " << faces[i][3] << " " << faces[i][2]
              << " " << faces[i][4] << " " << faces[i][5]
              << " " << faces[i][7] << " " << faces[i][6];

      } else {
        DUNE_THROW(Dune::NotImplemented, "Geometry type " << geometryTypes[i] << " not supported yet");
      }

      fgrid << std::endl;
    }

    fgrid << std::endl;

    // 3d VTK files need an extra section specifying the CELL_TYPES aka GeometryTypes
    if (dim==3) {

      fgrid << "CELL_TYPES " << geometryTypes.size() << std::endl;

      for (size_t i=0; i<geometryTypes.size(); i++) {
        if (geometryTypes[i].isSimplex())
          fgrid << "10" << std::endl;
        else if (geometryTypes[i].isHexahedron())
          fgrid << "12" << std::endl;
        else if (geometryTypes[i].isPrism())
          fgrid << "13" << std::endl;
        else if (geometryTypes[i].isPyramid())
          fgrid << "14" << std::endl;
        else
          DUNE_THROW(Dune::NotImplemented, "Geometry type " << geometryTypes[i] << " not supported yet");

      }

    }

#if 0
    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (gridSubEntityData.size()-1);

    fgrid << "CELL_DATA " << gridSubEntityData.size() << std::endl;
    fgrid << "SCALARS property_coding " << Dune::className<ctype>() << " 1" << std::endl;
    fgrid << "LOOKUP_TABLE default" << std::endl;

    for (typename GridSubEntityData::const_iterator sEIt = gridSubEntityData.begin();
         sEIt != gridSubEntityData.end();
         ++sEIt, accum += delta)
    {
      // "encode" the parent with one color...
      fgrid << accum << std::endl;
    }
#endif
    fgrid.close();
  }


  /** \brief Write either the grid0 or the grid1-side into streams
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeIntersections(const Glue& glue, const std::string& filename)
  {
    dune_static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    std::ofstream fmerged;

    fmerged.open(filename.c_str());

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    static const Dune::PartitionIteratorType PType = Dune::All_Partition;
    typedef typename GridView::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;

    typedef typename GridView::ctype ctype;

    const int dim = GridView::dimension;
    const int domdimw = GridView::dimensionworld;

    // coordinates have to be in R^3 in the VTK format
    std::string coordinatePadding;
    for (int i=domdimw; i<3; i++)
      coordinatePadding += " 0";

    int overlaps = glue.size();

    // WRITE POINTS
    // ----------------
    typedef typename Glue::Grid0Patch Extractor;
    std::vector<typename Extractor::Coords> coords;
    glue.template patch<side>().getCoords(coords);

    RemoteIntersectionIterator isIt = glue.template ibegin<side>();
    unsigned int elementPoints = isIt->inside()->geometry().corners();

    // the merged grid (i.e. the set of remote intersections
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << filename << "\nASCII" << std::endl;
    fmerged << ((dim==3) ? "DATASET UNSTRUCTURED_GRID" : "DATASET POLYDATA") << std::endl;
    fmerged << "POINTS " << elementPoints*overlaps << " " << Dune::className<ctype>() << std::endl;

    for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
         isIt != glue.template iend<side>();
         ++isIt)

    {
      for (int i = 0; i < isIt->inside()->geometry().corners(); ++i){
        if (side == 0){
          //std::cout <<"isIt->inside()->geometry().corner(i)" <<isIt->inside()->geometry().corner(i) << std::endl;
          fmerged << isIt->inside()->geometry().corner(i) << coordinatePadding << std::endl;
        }
        if (side == 1){
          fmerged << isIt->geometry().corner(i) << coordinatePadding << std::endl;
         }
      }
    }

    // WRITE POLYGONS
    // ----------------

    std::vector<typename Extractor::VertexVector> faces;
    std::vector<Dune::GeometryType> geometryTypes;
    glue.template patch<side>().getFaces(faces);
    glue.template patch<side>().getGeometryTypes(geometryTypes);

    unsigned int faceCornerCount = 0;
    for (size_t i=0; i<faces.size(); i++){
      faceCornerCount += faces[i].size();
    }

    //int domainSimplexCorners = dim-Glue::Grid0Patch::codim+1;
    fmerged << ((dim==3) ? "CELLS " : "LINES ")
            << overlaps << " " << (elementPoints+1)*overlaps << std::endl;

    unsigned int i = 0;
    for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
         isIt != glue.template iend<side>();
         ++isIt)
      {
        fmerged << elementPoints;
        for (int j=0; j<elementPoints; j++)
          fmerged << " " <<elementPoints*i+j;
        fmerged << std::endl;
        i++;
      }

    // 3d VTK files need an extra section specifying the CELL_TYPES aka GeometryTypes
    if (dim==3) {

      fmerged << "CELL_TYPES " << overlaps << std::endl;
      for (size_t i=0; i<overlaps; i++)
        fmerged << "11" << std::endl;
    }

#if 0
    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (gridSubEntityData.size()-1);

    fmerged << "CELL_DATA " << overlaps << std::endl;
    fmerged << "SCALARS property_coding " << Dune::className<ctype>() << " 1" << std::endl;
    fmerged << "LOOKUP_TABLE default" << std::endl;

    for (typename GridSubEntityData::const_iterator sEIt = gridSubEntityData.begin();
         sEIt != gridSubEntityData.end();
         ++sEIt, accum += delta)
    {
      // ...and mark all of its merged grid parts with the same color
      for (int j = 0; j < sEIt->first.second; ++j)
        fmerged << accum << std::endl;
    }
#endif
    fmerged.close();
  }

 /** \brief Write multi_map for mapping elements to intersections
   * \tparam side Write the grid0-side if this is 0, and grid1 if it is 1.
   */
  template <class Glue, int side>
  static void writeElementIntersectionMap(const Glue& glue , std::multimap<int,int> mymm2)
  {
    dune_static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    static const Dune::PartitionIteratorType PType = Dune::All_Partition;
    typedef typename GridView::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;

    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ::SingleCodimSingleGeomTypeMapper Mapper;

     Mapper mapper(glue.template patch<side>().gridView());
     std::multimap<int,RemoteIntersectionIterator> mymm;

    // iteration over all remote intersections
    for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
         isIt != glue.template iend<side>();
         ++isIt)
      {
        // get element ID from inside-element of the remote intersection
        int eId = mapper.map(*isIt->inside());

        // get Intersection ID from remote intersection
        int interId = isIt->index();

        // insert element ID and RemoteIntersection iterator in multi-map
        mymm.insert (std::pair<int,RemoteIntersectionIterator>(eId,isIt));

        // insert element ID and intersection ID in multi-map
        mymm2.insert (std::pair<int,int>(eId,interId));
      }

    // test multi-maps by iteration over all elements
    ElementIter eIt =  glue.template patch<side>().gridView().template begin<0>();
    ElementIter eEndIt = glue.template patch<side>().gridView().template end<0>();
    for (; eIt != eEndIt; ++eIt)
      {
        int eId = mapper.map(*eIt);
        std::cout << "There are " << mymm2.count(eId) << " elements with key " << eId << ":";

        typename std::multimap<int,int>::iterator it;
        for (it=mymm2.equal_range(eId).first; it!=mymm2.equal_range(eId).second; ++it)
          std::cout << ' ' << (*it).second;
        std::cout << '\n';
      }
  }

public:
  template<typename Glue>
  static void write(const Glue& glue, const std::string& filenameTrunk)
  {

    // Write extracted grid and remote intersection on the grid0-side
    writeExtractedPart<Glue,0>(glue,
                               filenameTrunk + "-domain.vtk");

    writeIntersections<Glue,0>(glue,
                               filenameTrunk + "-intersections-domain.vtk");

    // Write extracted grid and remote intersection on the grid1-side
    writeExtractedPart<Glue,1>(glue,
                               filenameTrunk + "-target.vtk");

    writeIntersections<Glue,1>(glue,
                               filenameTrunk + "-intersections-target.vtk");

    std::multimap<int,int> mymm0;
    writeElementIntersectionMap<Glue,0>(glue, mymm0);

    std::multimap<int,int> mymm1;
    writeElementIntersectionMap<Glue,1>(glue, mymm1);
  }

};

#endif // GRIDGLUEVTKWRITER_ROOT_HH
