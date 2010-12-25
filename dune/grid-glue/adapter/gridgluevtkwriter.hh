// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueVtkWriter.hh
 *  Version:     1.0
 *  Created on:  Mar 5, 2009
 *  Author:      Gerrit Buse
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

#ifndef GRIDGLUEVTKWRITER_HH_
#define GRIDGLUEVTKWRITER_HH_


#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/genericreferenceelements.hh>


const char TypeNames[][8] = { "double", "float  ", "int    " };

// double is default
template<typename T>
struct Nametraits
{
  enum { nameidx = 0 };
};


template<>
struct Nametraits<double>
{
  enum { nameidx = 0 };
};

template<>
struct Nametraits<float>
{
  enum { nameidx = 1 };
};

template<>
struct Nametraits<int>
{
  enum { nameidx = 2 };
};


/** \brief Write remote intersections to a vtk file for debugging purposes
 */
class GridGlueVtkWriter
{

public:
  template<typename Glue>
  static void write(const Glue& glue, const char* filename_trunk)
  {
    const int dimw = Glue::dimworld;

    typedef typename Glue::ctype ctype;

    typedef std::list<int> FaceList;
    FaceList faces, parts, face_corners;
    typename FaceList::iterator faceit, facecornerit;

    std::ofstream fgrid, fmerged;

    std::string fngrid(filename_trunk);
    fngrid += "-domain.vtk";
    std::string fnmerged(filename_trunk);
    fnmerged += "-overlaps-domain.vtk";

    fgrid.open(fngrid.c_str());
    fmerged.open(fnmerged.c_str());

    typedef typename Glue::Grid0View Grid0View;
    typedef typename Grid0View::Traits::template Codim<0>::Iterator DomainIter;
    typedef typename Grid0View::Traits::template Codim<0>::EntityPointer DomainEPtr;
    const int domainDim = Grid0View::dimension;
    const int domdimw = Grid0View::dimensionworld;

    std::string domainCoordinatePadding;
    for (int i=domdimw; i<dimw; i++)
      domainCoordinatePadding += " 0";

    const Grid0View& domgv = glue.template gridView<0>();

    // remember the entities that have been mapped
    std::list<DomainEPtr> domeptrs(0, (DomainEPtr) domgv.template begin<0>());

    int overlaps = 0, face_corner_count = 0;

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << fngrid << "\nASCII" << std::endl;
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << fnmerged << "\nASCII" << std::endl;

    // count the points and the polygons
    for (DomainIter pit = domgv.template begin<0>(); pit != domgv.template end<0>(); ++pit)
    {
      int face = 0;
      while ((face = glue.domainEntityNextFace(*pit, face)) != -1)
      {
        faces.push_back(face);
        domeptrs.push_back(pit);

        // count and remember the corners of this subEntity
        const Dune::GenericReferenceElement<ctype, domainDim>& refElement =
          Dune::GenericReferenceElements<ctype, domainDim>::general(pit->type());
        int size = refElement.size(face, Glue::Grid0Patch::codim, domainDim);
        face_corners.push_back(size);
        face_corner_count += size;

        int num_parts = 0;
        for (typename Glue::DomainIntersectionIterator domisit = glue.idomainbegin(*pit, face); domisit != glue.idomainend(); ++domisit)
        {
          num_parts++;
          overlaps++;
        }
        parts.push_back(num_parts);

        face++;                 // move to next face
      }
    }


    // WRITE POINTS
    // ----------------

    fgrid << "DATASET POLYDATA\nPOINTS " << face_corner_count << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;
    fmerged << "DATASET POLYDATA\nPOINTS " << overlaps*3 << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<DomainEPtr>::const_iterator pit = domeptrs.begin(); pit != domeptrs.end(); ++pit, ++faceit)
    {

      const Dune::GenericReferenceElement<ctype, domainDim>& refElement =
        Dune::GenericReferenceElements<ctype, domainDim>::general((*pit)->type());

      // Write the current subentity into the fgrid file
      for (int i=0; i<refElement.size(*faceit, Glue::Grid0Patch::codim, domainDim); i++)
        fgrid << (*pit)->geometry().corner(refElement.subEntity(*faceit, Glue::Grid0Patch::codim, i, domainDim))
              << domainCoordinatePadding
              << std::endl;

      // write the remote intersections of the current subentity into the fmerged file
      for (typename Glue::DomainIntersectionIterator domisit = glue.idomainbegin(**pit, *faceit); domisit != glue.idomainend(); ++domisit)
      {
        // Later down we assume that all remote intersections are simplices
        assert(domisit->intersectionDomainGlobal().type().isSimplex());

        for (int i = 0; i < domisit->intersectionDomainGlobal().corners(); ++i)
          fmerged << domisit->intersectionDomainGlobal().corner(i) << domainCoordinatePadding << std::endl;
      }
    }

    // WRITE POLYGONS
    // ----------------

    fgrid << "POLYGONS " << faces.size() << " " << faces.size() + face_corner_count << std::endl;

    facecornerit = face_corners.begin();
    face_corner_count = 0;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<DomainEPtr>::const_iterator pit = domeptrs.begin();
         pit != domeptrs.end();
         ++pit, ++faceit) {

      const Dune::GenericReferenceElement<ctype, domainDim>& refElement =
        Dune::GenericReferenceElements<ctype, domainDim>::general((*pit)->type());

      int size = refElement.size(*faceit, Glue::Grid0Patch::codim, domainDim);

      fgrid << size;

      // vtk expects the vertices to by cyclically ordered
      // therefore unfortunately we have to deal with several element types on a case-by-case basis
      if (refElement.type(*faceit,Glue::Grid0Patch::codim).isQuadrilateral()) {
        fgrid << " " << face_corner_count << " " << face_corner_count+1
              << " " << face_corner_count+3 << " " << face_corner_count+2;

        face_corner_count += 4;
      } else {
        for (int j = 0; j <size; ++j)
          fgrid << " " << face_corner_count++;
      }

      fgrid << std::endl;
    }
    fgrid << std::endl;

    int domainSimplexCorners = domainDim-Glue::Grid0Patch::codim+1;
    fmerged << "POLYGONS " << overlaps << " " << (domainSimplexCorners+1)*overlaps << std::endl;

    for (int i = 0; i < overlaps; ++i) {
      fmerged << domainSimplexCorners;
      for (int j=0; j<domainSimplexCorners; j++)
        fmerged << " " << domainSimplexCorners*i+j;
      fmerged << std::endl;
    }


    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (faces.size()-1);

    fgrid << "CELL_DATA " << faces.size() << std::endl;
    fgrid << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fgrid << "LOOKUP_TABLE default" << std::endl;

    fmerged << "CELL_DATA " << overlaps << std::endl;
    fmerged << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fmerged << "LOOKUP_TABLE default" << std::endl;

    for (faceit = parts.begin(); faceit != parts.end(); ++faceit, accum += delta)
    {
      // "encode" the parent with one color...
      fgrid << accum << std::endl;
      // ...and mark all of its merged grid parts with the same color
      for (int j = 0; j < *faceit; ++j)
        fmerged << accum << std::endl;
    }

    fgrid.close();
    fmerged.close();
    domeptrs.clear();


    // NOW FOR THE TARGET SIDE
    // ***********************

    fngrid = filename_trunk;
    fngrid += "-target.vtk";
    fnmerged = filename_trunk;
    fnmerged += "-overlaps-target.vtk";

    fgrid.open(fngrid.c_str());
    fmerged.open(fnmerged.c_str());

    typedef typename Glue::Grid1View Grid1View;
    typedef typename Grid1View::Traits::template Codim<0>::Iterator TargetIter;
    typedef typename Grid1View::Traits::template Codim<0>::EntityPointer TargetEPtr;
    const int targetDim = Grid1View::dimension;
    const int tardimw = Grid1View::dimensionworld;

    std::string targetCoordinatePadding;
    for (int i=tardimw; i<dimw; i++)
      targetCoordinatePadding += " 0";

    const Grid1View& targv = glue.template gridView<1>();

    // remember the entities that have been mapped
    std::list<TargetEPtr> tareptrs(0, targv.template begin<0>());

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << fngrid << "\nASCII" << std::endl;
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << fnmerged << "\nASCII" << std::endl;

    // reset some of the variables
    overlaps = 0;
    face_corner_count = 0;
    faces.clear();
    face_corners.clear();
    parts.clear();

    // count the points and the polygons
    for (TargetIter pit = targv.template begin<0>(); pit != targv.template end<0>(); ++pit)
    {
      int face = 0;
      while ((face = glue.targetEntityNextFace(*pit, face)) != -1)
      {
        faces.push_back(face);
        tareptrs.push_back(pit);

        // count and remember the corners of this face
        const Dune::GenericReferenceElement<ctype, targetDim>& refElement =
          Dune::GenericReferenceElements<ctype, targetDim>::general(pit->type());
        int size = refElement.size(face, Glue::Grid1Patch::codim, targetDim);
        face_corners.push_back(size);
        face_corner_count += size;

        int num_parts = 0;
        for (typename Glue::TargetIntersectionIterator tarisit = glue.itargetbegin(*pit, face); tarisit != glue.itargetend(); ++tarisit)
        {
          num_parts++;
          overlaps++;
        }
        parts.push_back(num_parts);

        face++;                 // move to next face
      }
    }


    // WRITE POINTS
    // ----------------

    fgrid << "DATASET POLYDATA\nPOINTS " << face_corner_count << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;
    fmerged << "DATASET POLYDATA\nPOINTS " << overlaps*3 << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<TargetEPtr>::const_iterator pit = tareptrs.begin(); pit != tareptrs.end(); ++pit, ++faceit)
    {

      const Dune::GenericReferenceElement<ctype, targetDim>& refElement =
        Dune::GenericReferenceElements<ctype, targetDim>::general((*pit)->type());

      // Write the current subentity into the fgrid file
      for (int i=0; i<refElement.size(*faceit, Glue::Grid1Patch::codim, targetDim); i++)
        fgrid << (*pit)->geometry().corner(refElement.subEntity(*faceit, Glue::Grid1Patch::codim, i, targetDim))
              << targetCoordinatePadding
              << std::endl;

      // write the merged grid refinement of this surface part
      for (typename Glue::TargetIntersectionIterator tarisit = glue.itargetbegin(**pit, *faceit); tarisit != glue.itargetend(); ++tarisit)
      {
        // Later down we assume that all remote intersections are simplices
        assert(tarisit->intersectionTargetGlobal().type().isSimplex());

        for (int i = 0; i < tarisit->intersectionTargetGlobal().corners(); ++i)
          fmerged << tarisit->intersectionTargetGlobal().corner(i) << targetCoordinatePadding << std::endl;
      }
    }

    // WRITE POLYGONS
    // ----------------

    fgrid << "POLYGONS " << faces.size() << " " << faces.size() + face_corner_count << std::endl;

    facecornerit = face_corners.begin();
    face_corner_count = 0;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<TargetEPtr>::const_iterator pit = tareptrs.begin();
         pit != tareptrs.end();
         ++pit, ++faceit) {

      const Dune::GenericReferenceElement<ctype, targetDim>& refElement =
        Dune::GenericReferenceElements<ctype, targetDim>::general((*pit)->type());

      int size = refElement.size(*faceit, Glue::Grid1Patch::codim, targetDim);

      fgrid << size;

      // vtk expects the vertices to by cyclically ordered
      // therefore unfortunately we have to deal with several element types on a case-by-case basis
      if (refElement.type(*faceit,Glue::Grid1Patch::codim).isQuadrilateral()) {
        fgrid << " " << face_corner_count << " " << face_corner_count+1
              << " " << face_corner_count+3 << " " << face_corner_count+2;

        face_corner_count += 4;
      } else {
        for (int j = 0; j <size; ++j)
          fgrid << " " << face_corner_count++;
      }

      fgrid << std::endl;
    }
    fgrid << std::endl;

    int targetSimplexCorners = targetDim-Glue::Grid1Patch::codim+1;
    fmerged << "POLYGONS " << overlaps << " " << (targetSimplexCorners+1)*overlaps << std::endl;

    for (int i = 0; i < overlaps; ++i) {
      fmerged << targetSimplexCorners;
      for (int j=0; j<targetSimplexCorners; j++)
        fmerged << " " << targetSimplexCorners*i+j;
      fmerged << std::endl;
    }

    // WRITE CELL DATA
    // ---------------
    accum = 0.0;
    delta = 1.0 / (ctype) (faces.size()-1);

    fgrid << "CELL_DATA " << faces.size() << std::endl;
    fgrid << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fgrid << "LOOKUP_TABLE default" << std::endl;

    fmerged << "CELL_DATA " << overlaps << std::endl;
    fmerged << "SCALARS property_coding " << TypeNames[Nametraits<ctype>::nameidx] << " 1" << std::endl;
    fmerged << "LOOKUP_TABLE default" << std::endl;

    for (faceit = parts.begin(); faceit != parts.end(); ++faceit, accum += delta)
    {
      // "encode" the parent with one color...
      fgrid << accum << std::endl;
      // ...and mark all of its merged grid parts with the same color
      for (int j = 0; j < *faceit; ++j)
        fmerged << accum << std::endl;
    }

    fgrid.close();
    fmerged.close();


  }

};



#endif // GRIDGLUEVTKWRITER_HH_
