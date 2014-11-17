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

struct data2test{

  /*typedef typename Dune::conditional<(side==0),
                                     typename Glue::Grid0View,
                                     typename Glue::Grid1View>::type GridView;
  typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename Dune::conditional<(side==0),
                                     typename Glue::Grid0IntersectionIterator,
                                     typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;
  */
  int elm;
  int inter;
  //ElementPtr eIt;
  //RemoteIntersectionIterator rIt;
  double frac;
  double data;

  data2test()
    {

    }


};
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
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

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
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

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
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename  Glue::Grid0View GridView0;
    typedef typename  Glue::Grid1View GridView1;

    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    static const Dune::PartitionIteratorType PType = Dune::All_Partition;
    typedef typename GridView::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;

    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ::SingleCodimSingleGeomTypeMapper Mapper;

     std::multimap<int,RemoteIntersectionIterator> mymm0;
     std::multimap<int,int> mymm1;

     Mapper mapper(glue.template patch<side>().gridView());

    // iteration over all remote intersections
     for (RemoteIntersectionIterator isIt = glue.template ibegin<side>();
            isIt != glue.template iend<side>();
            ++isIt)
       {
         int eId = 0;
         int eId2 = 0;

         mymm0.insert (std::pair<int,int>(eId,isIt));
         mymm0.insert (std::pair<int,int>(eId,eId2));
       }

     // test multi-maps by iteration over all elements
     ElementIter eIt0 =  glue.template patch<side>().gridView().template begin<0>();
     ElementIter eEndIt0 = glue.template patch<side>().gridView().template end<0>();
     for (; eIt0 != eEndIt0; ++eIt0)
       {
         int eId0 = mapper.map(*eIt0);
         std::cout << "There are " << mymm0.count(eId0) << " elements with key " << eId0 << ":";

         typename std::multimap<int,int>::iterator it0;
         for (it0=mymm0.equal_range(eId0).first; it0!=mymm0.equal_range(eId0).second; ++it0)
           std::cout << ' ' << (*it0).second;
         std::cout << '\n';
       }
  }

  template <class Glue, int side>
  static void writeElementIntersectionMap0(const Glue& glue, std::multimap<int,data2test> &mymm3)
  {
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename  Glue::Grid0View GridView0;
    typedef typename  Glue::Grid1View GridView1;

    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    static const Dune::PartitionIteratorType PType = Dune::All_Partition;
    typedef typename GridView::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;
    typedef typename GridView0::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter0;
    typedef typename GridView1::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter1;

    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;
    typedef typename Glue::Grid0IntersectionIterator RemoteIntersectionIterator0;
    typedef typename Glue::Grid1IntersectionIterator RemoteIntersectionIterator1;

    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ::SingleCodimSingleGeomTypeMapper Mapper;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView0, 0> ::SingleCodimSingleGeomTypeMapper Mapper0;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView1, 0> ::SingleCodimSingleGeomTypeMapper Mapper1;

    std::multimap<int,RemoteIntersectionIterator0> mymm0;
    //std::multimap<int,data2test<Glue,side> > mymm3;

    Mapper mapper(glue.template patch<side>().gridView());
    Mapper0 mapper0(glue.template patch<0>().gridView());
    Mapper1 mapper1(glue.template patch<1>().gridView());

    int gridsize1d = glue.template patch<1>().gridView().size(0);
    std::vector<double > frac1d(gridsize1d);

    for (RemoteIntersectionIterator0 isIt = glue.template ibegin<0>();
          isIt != glue.template iend<0>();
            ++isIt)
       {
         int eId = mapper0.map(*isIt->inside());
         int eId2 = mapper1.map(*isIt->outside());
         int interId = isIt->index();
         data2test mapstruct;

         mapstruct.elm = eId2;
         mapstruct.inter = interId;
         frac1d[eId2] += isIt->geometry().volume()/ isIt->outside()->geometry().volume();

         mymm0.insert (std::pair<int,RemoteIntersectionIterator0>(eId,isIt));
         mymm3.insert (std::pair<int, data2test>(eId,mapstruct));
       }

     // test multi-maps by iteration over all elements
     ElementIter0 eIt1 =  glue.template patch<0>().gridView().template begin<0>();
     ElementIter0 eEndIt1 = glue.template patch<0>().gridView().template end<0>();
     for (; eIt1 != eEndIt1; ++eIt1)
       {
         int eId1 = mapper0.map(*eIt1);
         double interVolsum = 0.;
         double fracTot = 0;
         std::vector<double > frac(mymm3.count(eId1));

         if (mymm3.count(eId1) != 0)
           {
             double interVol = 0.;
             int k = 0;
             double eVol = 0.;

             typename std::multimap<int,RemoteIntersectionIterator0>::iterator it1;
             typename std::multimap<int,data2test>::iterator it3;
             it3=mymm3.equal_range(eId1).first;
             for (it1=mymm0.equal_range(eId1).first; it1!=mymm0.equal_range(eId1).second; ++it1){
               interVol = (*it1).second->geometry().volume();
               eVol = (*it1).second->outside()->geometry().volume();

               (*it3).second.frac = interVol/eVol/frac1d[(*it3).second.elm];
               k++;
               ++it3;
             }
           }
       }

     eIt1 =  glue.template patch<0>().gridView().template begin<0>();
     eEndIt1 = glue.template patch<0>().gridView().template end<0>();

     for (; eIt1 != eEndIt1; ++eIt1)
       {
         int eId1 = mapper0.map(*eIt1);
         std::cout << "There are " << mymm3.count(eId1) << " elements with key " << eId1 << ":"<< '\n';

         typename std::multimap<int,data2test >::iterator it3;
         for (it3=mymm3.equal_range(eId1).first; it3!=mymm3.equal_range(eId1).second; ++it3){
           std::cout << "elem: " << (*it3).second.elm << " inter: " << (*it3).second.inter << " frac: " << (*it3).second.frac<< '\n';
         }
       }
  }

  template <class Glue, int side>
  static void writeElementIntersectionMap1(const Glue& glue, std::multimap<int,data2test > &mymm3)
  {
    static_assert((side==0 || side==1), "'side' can only be 0 or 1");

    typedef typename Dune::conditional<(side==0), typename Glue::Grid0View, typename Glue::Grid1View>::type GridView;
    typedef typename  Glue::Grid0View GridView0;
    typedef typename  Glue::Grid1View GridView1;

    typedef typename GridView::Traits::template Codim<0>::EntityPointer ElementPtr;
    typedef typename GridView::Traits::template Codim<0>::Entity Element;

    static const Dune::PartitionIteratorType PType = Dune::All_Partition;
    typedef typename GridView::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;
    typedef typename GridView1::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter1;

    typedef typename Dune::conditional<(side==0),
        typename Glue::Grid0IntersectionIterator,
        typename Glue::Grid1IntersectionIterator>::type RemoteIntersectionIterator;
    typedef typename Glue::Grid0IntersectionIterator RemoteIntersectionIterator0;
    typedef typename Glue::Grid1IntersectionIterator RemoteIntersectionIterator1;

    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView, 0> ::SingleCodimSingleGeomTypeMapper Mapper;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView0, 0> ::SingleCodimSingleGeomTypeMapper Mapper0;
    typedef typename Dune::SingleCodimSingleGeomTypeMapper<GridView1, 0> ::SingleCodimSingleGeomTypeMapper Mapper1;

     std::multimap<int,RemoteIntersectionIterator1> mymm0;
     //std::multimap<int,data2test<Glue,side> > mymm3;

     Mapper mapper(glue.template patch<side>().gridView());
     Mapper0 mapper0(glue.template patch<0>().gridView());
     Mapper1 mapper1(glue.template patch<1>().gridView());

     for (RemoteIntersectionIterator1 isIt = glue.template ibegin<1>();
          isIt != glue.template iend<1>();
            ++isIt)
       {
         int eId = mapper1.map(*isIt->inside());
         int eId2 = mapper0.map(*isIt->outside());
         int interId = isIt->index();
         data2test mapstruct;

         mapstruct.elm = eId2;
         mapstruct.inter = interId;

         mymm0.insert (std::pair<int,RemoteIntersectionIterator1>(eId,isIt));
         mymm3.insert (std::pair<int, data2test>(eId,mapstruct));
       }

     // test multi-maps by iteration over all elements
     ElementIter1 eIt1 =  glue.template patch<1>().gridView().template begin<0>();
     ElementIter1 eEndIt1 = glue.template patch<1>().gridView().template end<0>();

     for (; eIt1 != eEndIt1; ++eIt1)
       {
         int eId1 = mapper1.map(*eIt1);
         double interVolsum = 0.;
         double fracTot = 0;
         std::vector<double > frac(mymm3.count(eId1));

         if (mymm3.count(eId1) != 0)
           {
             double eVol = eIt1->geometry().volume();
             double interVol = 0.;
             double interVolsum = 0.;
             int k = 0;

             typename std::multimap<int,RemoteIntersectionIterator1>::iterator it1;
             typename std::multimap<int,data2test>::iterator it3;
             for (it1=mymm0.equal_range(eId1).first; it1!=mymm0.equal_range(eId1).second; ++it1){
               interVol = (*it1).second->geometry().volume();
               interVolsum += (*it1).second->geometry().volume();
               frac[k] = interVol/eVol;
               k++;
             }

             it3=mymm3.equal_range(eId1).first;
             for (int m = 0; m < frac.size(); m++){
               frac[m] = frac[m]/(interVolsum/eVol);
               fracTot +=  frac[m];
               (*it3).second.frac =  frac[m];
               ++it3;
             }
             std::cout << "frac total: " << fracTot <<  std::endl;
           }
       }

     eIt1 =  glue.template patch<1>().gridView().template begin<0>();
     eEndIt1 = glue.template patch<1>().gridView().template end<0>();

     for (; eIt1 != eEndIt1; ++eIt1)
       {
         int eId1 = mapper1.map(*eIt1);
         std::cout << "There are " << mymm3.count(eId1) << " elements with key " << eId1 << ":"<< '\n';

         typename std::multimap<int,data2test>::iterator it3;
         for (it3=mymm3.equal_range(eId1).first; it3!=mymm3.equal_range(eId1).second; ++it3){
           std::cout << "elem: " << (*it3).second.elm << " inter: " << (*it3).second.inter << " frac: " << (*it3).second.frac<< '\n';
         }
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
  }

template<typename Glue>
static void writeMaps(const Glue& glue, std::multimap<int,data2test> &mymm0, std::multimap<int,data2test> &mymm1)
  {
    //writeElementIntersectionMap<Glue,0>(glue, mymm0);

    //writeElementIntersectionMap<Glue,1>(glue, mymm1);

    writeElementIntersectionMap0<Glue,0>(glue, mymm0);

    writeElementIntersectionMap1<Glue,1>(glue, mymm1);

  }

};

#endif // GRIDGLUEVTKWRITER_ROOT_HH
