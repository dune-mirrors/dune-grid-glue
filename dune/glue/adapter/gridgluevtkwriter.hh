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
 * @file GridGlueVtkWriter.hh
 * @brief Class thought to make graphical debugging of couplings easier.
 */

#ifndef GRIDGLUEVTKWRITER_HH_
#define GRIDGLUEVTKWRITER_HH_


#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/genericreferenceelements.hh>

#include <dune/glue/misc/orientedsubface.hh>


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

    typedef typename Glue::DomainGridView DomainGridView;
    typedef typename Glue::DomainGridType DomainGridType;
    typedef typename DomainGridView::Traits::template Codim<0>::Iterator DomainIter;
    typedef typename DomainGridView::Traits::template Codim<0>::EntityPointer DomainEPtr;
    const int domdimw = DomainGridType::dimensionworld;

    bool pad = dimw > domdimw;
    std::string coordinatePadding;
    for (int i=1; i<dimw; i++)
      coordinatePadding += " 0";

    bool hyper = dimw > DomainGridType::dimension;

    const DomainGridView& domgv = glue.domainGridView();

    // remember the entities that have been mapped
    std::list<DomainEPtr> domeptrs(0, (DomainEPtr) domgv.template begin<0>());

    int overlaps = 0, parents = 0, face_corner_count = 0;

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

        // count and remember the corners of this face
        const Dune::GenericReferenceElement<ctype, dimw>& refElement =
          Dune::GenericReferenceElements<ctype, dimw>::general(pit->type());
        int size = refElement.size(face, 1, dimw);
        face_corners.push_back(size);
        face_corner_count += size;

        parents++;

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

    fgrid << "DATASET POLYDATA\nPOINTS " << (dimw == 2 ? parents*4 : face_corner_count) << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;
    fmerged << "DATASET POLYDATA\nPOINTS " << overlaps*3 << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<DomainEPtr>::const_iterator pit = domeptrs.begin(); pit != domeptrs.end(); ++pit, ++faceit)
    {
      // write the domain element
      if (dimw == 2)
      {
        if (hyper)
        {
          fgrid << (*pit)->geometry().corner(0) << (pad ? " 0" : "  ") << " 0\n"
                << (*pit)->geometry().corner(0) << (pad ? " 0" : "  ") << " 0.01\n"
                << (*pit)->geometry().corner(1) << (pad ? " 0" : "  ") << " 0\n"
                << (*pit)->geometry().corner(1) << (pad ? " 0" : "  ") << " 0.01" << std::endl;
        }
        else                         // full-dimensional grid
        {
          for (int i=0; i<(*pit)->geometry().corners(); i++)
            fgrid << (*pit)->geometry().corner(i) << coordinatePadding << std::endl;
        }
      }
      else                   // dimw == 3
      {
        if (hyper)
        {
          // write 3 or 4 points, depending on face geometry
          for (int i = 0; i < 3; ++i)
          {
            if (i == 2 && (*facecornerit) == 4)
              fgrid << (*pit)->geometry().corner(3) << (pad ? " 0" : "  ") << std::endl;
            fgrid << (*pit)->geometry().corner(i) << (pad ? " 0" : "  ") << std::endl;
          }
        }
        else                         // full-dimensional grid
        {
          Dune::GeometryType temp_gt = (*pit)->type();
          // write 3 or 4 points, depending on face geometry
          for (int i = 0; i < 3; ++i)
          {
            if (i == 2 && (*facecornerit) == 4)
            {
              int corner = orientedSubface<DomainGridType::dimension>(temp_gt, *faceit, 3);
              fgrid << (*pit)->template subEntity<DomainGridType::dimension>(corner)->geometry().corner(0) << std::endl;
            }
            int corner = orientedSubface<DomainGridType::dimension>(temp_gt, *faceit, i);
            fgrid << (*pit)->template subEntity<DomainGridType::dimension>(corner)->geometry().corner(0) << std::endl;
          }
        }
        facecornerit++;
      }

      // write the merged grid refinement of this surface part
      for (typename Glue::DomainIntersectionIterator domisit = glue.idomainbegin(**pit, *faceit); domisit != glue.idomainend(); ++domisit)
      {
        for (int i = 0; i < domisit->intersectionDomainGlobal().corners(); ++i)
          fmerged << domisit->intersectionDomainGlobal().corner(i) << coordinatePadding << std::endl;
      }
    }

    // WRITE POLYGONS
    // ----------------

    if (dimw == 2)
    {
      fgrid << "POLYGONS " << parents << " " << 5*parents << std::endl;
      fmerged << "POLYGONS " << overlaps << " " << 4*overlaps << std::endl;

      for (int i = 0; i < 4*parents; i += 4)
        fgrid << "4 " << i << " " << i+1 << " " << i+3 << " " << i+2 << std::endl;
      for (int i = 0; i < 3*overlaps; i += 3)
        fmerged << "3 " << i << " " << i+1 << " " << i+2 << std::endl;
    }
    else             // dimw == 3
    {
      fgrid << "POLYGONS " << parents << " " << parents + face_corner_count;
      fmerged << "POLYGONS " << overlaps << " " << 4*overlaps << std::endl;

      facecornerit = face_corners.begin();
      face_corner_count = 0;
      for (int i = 0; i < parents; ++i, ++facecornerit)
      {
        fgrid << "\n" << (*facecornerit);
        for (int i = 0; i < (*facecornerit); ++i, ++face_corner_count)
          fgrid << " " << face_corner_count;
      }
      fgrid << std::endl;
      for (int i = 0; i < overlaps; ++i)
        fmerged << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << std::endl;
    }

    // WRITE CELL DATA
    // ---------------
    ctype accum = 0.0, delta = 1.0 / (ctype) (parents-1);

    fgrid << "CELL_DATA " << parents << std::endl;
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

    typedef typename Glue::TargetGridView TargetGridView;
    typedef typename Glue::TargetGridType TargetGridType;
    typedef typename TargetGridView::Traits::template Codim<0>::Iterator TargetIter;
    typedef typename TargetGridView::Traits::template Codim<0>::EntityPointer TargetEPtr;
    const int tardimw = TargetGridType::dimensionworld;

    const TargetGridView& targv = glue.targetGridView();

    pad = dimw > tardimw;
    hyper = dimw > TargetGridType::dimension;

    // remember the entities that have been mapped
    std::list<TargetEPtr> tareptrs(0, targv.template begin<0>());

    fgrid << "# vtk DataFile Version 2.0\nFilename: " << fngrid << "\nASCII" << std::endl;
    fmerged << "# vtk DataFile Version 2.0\nFilename: " << fnmerged << "\nASCII" << std::endl;

    // reset some of the variables
    parents = 0;
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
        const Dune::GenericReferenceElement<ctype, dimw>& refElement =
          Dune::GenericReferenceElements<ctype, dimw>::general(pit->type());
        int size = refElement.size(face, 1, dimw);
        face_corners.push_back(size);
        face_corner_count += size;

        parents++;

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

    fgrid << "DATASET POLYDATA\nPOINTS " << (dimw == 2 ? parents*4 : face_corner_count) << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;
    fmerged << "DATASET POLYDATA\nPOINTS " << overlaps*(dimw == 2 ? 4 : 3) << " " << TypeNames[Nametraits<ctype>::nameidx] << std::endl;

    faceit = faces.begin();
    facecornerit = face_corners.begin();
    for (typename std::list<TargetEPtr>::const_iterator pit = tareptrs.begin(); pit != tareptrs.end(); ++pit, ++faceit)
    {
      // write the target element
      if (dimw == 2)
      {
        if (hyper)
        {
          fgrid << (*pit)->geometry().corner(0) << (pad ? " 0" : "  ") << " 0\n"
                << (*pit)->geometry().corner(0) << (pad ? " 0" : "  ") << " 0.01\n"
                << (*pit)->geometry().corner(1) << (pad ? " 0" : "  ") << " 0\n"
                << (*pit)->geometry().corner(1) << (pad ? " 0" : "  ") << " 0.01" << std::endl;
        }
        else                         // full-dimensional grid
        {
          Dune::GeometryType temp_gt = (*pit)->type();
          for (int i = 0; i < 2; ++i)
          {
            int corner = orientedSubface<TargetGridType::dimension>(temp_gt, *faceit, i);
            fgrid << (*pit)->template subEntity<TargetGridType::dimension>(corner)->geometry().corner(0) << " 0\n"
                  << (*pit)->template subEntity<TargetGridType::dimension>(corner)->geometry().corner(0) << " 0.01" << std::endl;
          }
        }
      }
      else                   // dimw == 3
      {
        if (hyper)
        {
          // write 3 or 4 points, depending on face geometry
          for (int i = 0; i < 3; ++i)
          {
            if (i == 2 && (*facecornerit) == 4)
              fgrid << (*pit)->geometry().corner(3) << (pad ? " 0" : "  ") << std::endl;
            fgrid << (*pit)->geometry().corner(i) << (pad ? " 0" : "  ") << std::endl;
          }
        }
        else                         // full-dimensional grid
        {
          Dune::GeometryType temp_gt = (*pit)->type();
          // write 3 or 4 points, depending on face geometry
          for (int i = 0; i < 3; ++i)
          {
            if (i == 2 && (*facecornerit) == 4)
            {
              int corner = orientedSubface<TargetGridType::dimension>(temp_gt, *faceit, 3);
              fgrid << (*pit)->template subEntity<TargetGridType::dimension>(corner)->geometry().corner(0) << std::endl;
            }
            int corner = orientedSubface<TargetGridType::dimension>(temp_gt, *faceit, i);
            fgrid << (*pit)->template subEntity<TargetGridType::dimension>(corner)->geometry().corner(0) << std::endl;
          }

        }
        facecornerit++;
      }

      // write the merged grid refinement of this surface part
      for (typename Glue::TargetIntersectionIterator tarisit = glue.itargetbegin(**pit, *faceit); tarisit != glue.itargetend(); ++tarisit)
      {
        if (dimw == 2)
        {
          fmerged << tarisit->intersectionTargetGlobal().corner(0) << (pad ? " 0" : "  ") << " 0\n"
                  << tarisit->intersectionTargetGlobal().corner(0) << (pad ? " 0" : "  ") << " 0.01\n"
                  << tarisit->intersectionTargetGlobal().corner(1) << (pad ? " 0" : "  ") << " 0\n"
                  << tarisit->intersectionTargetGlobal().corner(1) << (pad ? " 0" : "  ") << " 0.01" << std::endl;
        }
        else                         // dimw == 3
          for (int i = 0; i < 3; ++i)
            fmerged << tarisit->intersectionTargetGlobal().corner(i) << (pad ? " 0" : "  ") << std::endl;
      }
    }

    // WRITE POLYGONS
    // ----------------

    if (dimw == 2)
    {
      fgrid << "POLYGONS " << parents << " " << 5*parents << std::endl;
      fmerged << "POLYGONS " << overlaps << " " << 5*overlaps << std::endl;

      for (int i = 0; i < 2*parents; i += 2)
        fgrid << "4 " << 2*i << " " << 2*(i+1) << " " << 2*(i+1)+1 << " " << 2*i+1 << std::endl;
      for (int i = 0; i < 2*overlaps; i += 2)
        fmerged << "4 " << 2*i << " " << 2*(i+1) << " " << 2*(i+1)+1 << " " << 2*i+1 << std::endl;
    }
    else             // dimw == 3
    {
      fgrid << "POLYGONS " << parents << " " << parents + face_corner_count;
      fmerged << "POLYGONS " << overlaps << " " << 4*overlaps << std::endl;

      facecornerit = face_corners.begin();
      face_corner_count = 0;
      for (int i = 0; i < parents; ++i, ++facecornerit)
      {
        fgrid << "\n" << (*facecornerit);
        for (int i = 0; i < (*facecornerit); ++i, ++face_corner_count)
          fgrid << " " << face_corner_count;
      }
      fgrid << std::endl;
      for (int i = 0; i < overlaps; ++i)
        fmerged << "3 " << 3*i << " " << 3*i+1 << " " << 3*i+2 << std::endl;
    }

    // WRITE CELL DATA
    // ---------------
    accum = 0.0;
    delta = 1.0 / (ctype) (parents-1);

    fgrid << "CELL_DATA " << parents << std::endl;
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
