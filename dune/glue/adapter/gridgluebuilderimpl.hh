// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueBuilderImpl.hh
 *  Version:     1.0
 *  Created on:  Mar 10, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: configurable unit to set up a coupling scenario for GridGlue
 *  subversion:  $Id$
 *
 */
/**
 * @file GridGlueBuilderImpl.hh
 * @brief Configurable unit to set up a coupling scenario for GridGlue.
 */

#ifndef GRIDGLUEBUILDERIMPL_HH_
#define GRIDGLUEBUILDERIMPL_HH_

#include <dune/glue/extractors/surfacedescriptor.hh>
#include <dune/glue/adapter/gridglue.hh>
#include <dune/glue/extractors/vtksurfacewriter.hh>

template <typename GV, int codim>
struct ExtractorTypeTraits {};

template <typename GV>
struct ExtractorTypeTraits<GV, 0>
{
  typedef ElementDescriptor<GV> Type;
};

template <typename GV>
struct ExtractorTypeTraits<GV, 1>
{
  typedef FaceDescriptor<GV> Type;
};

/** \brief   I M P L E M E N T A T I O N   O F   S U B C L A S S   BUILDER IMPL
   (Specialization for two "non-surface" grids, i.e. meshes, manifolds) */

template<typename GET1, typename GET2>
template<typename BGET1, typename BGET2, int codim1, int codim2>
class GridGlue<GET1, GET2>::BuilderImpl
{
private:

  typedef GridGlue<GET1, GET2>  Parent;

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename Parent::ctype ctype;


private:

  Parent& glue_;

  typedef typename ExtractorTypeTraits<DomainGridView, codim1>::Type DomainDescriptor;

  typedef typename ExtractorTypeTraits<TargetGridView, codim2>::Type TargetDescriptor;

  const DomainDescriptor*                   domelmntdescr_;

  const TargetDescriptor*                   tarelmntdescr_;

  const typename Parent::DomainTransformation*    domtrafo_;

  const typename Parent::TargetTransformation*    tartrafo_;

  template<typename Extractor>
  void extractGrid (const Extractor & extractor,
                    std::vector<Dune::FieldVector<typename Parent::ctype, Parent::dimworld> > & coords,
                    std::vector<unsigned int> & faces,
                    std::vector<Dune::GeometryType>& geometryTypes,
                    const CoordinateTransformation<Extractor::dimworld, Parent::dimworld, typename Parent::ctype>* trafo) const
  {
    std::vector<typename Extractor::Coords> tempcoords;
    std::vector<typename Extractor::VertexVector> tempfaces;

    extractor.getCoords(tempcoords);
    coords.clear();
    coords.reserve(tempcoords.size());

    if (trafo != NULL)
    {
      std::cout << "GridGlue::Builder : apply trafo\n";
      for (size_t i = 0; i < tempcoords.size(); ++i)
      {
        typename Parent::Coords temp = (*trafo)(tempcoords[i]);
        coords.push_back(temp);
        // coords.push_back(Dune::FieldVector<typename Parent::ctype, Parent::dimworld>());
        // for (size_t j = 0; j < Parent::dimworld; ++j)
        //     coords.back()[j] = temp[j];
      }
    }
    else
    {
      for (unsigned int i = 0; i < tempcoords.size(); ++i)
      {
        assert(Parent::dimworld == Extractor::dimworld);
        coords.push_back(Dune::FieldVector<typename Parent::ctype, Parent::dimworld>());
        for (size_t j = 0; j < Parent::dimworld; ++j)
          coords.back()[j] = tempcoords[i][j];
      }
    }

    extractor.getFaces(tempfaces);
    faces.clear();

    for (unsigned int i = 0; i < tempfaces.size(); ++i) {
      for (int j = 0; j < tempfaces[i].size(); ++j)
        faces.push_back(tempfaces[i][j]);
    }

    // get the list of geometry types from the extractor
    extractor.getGeometryTypes(geometryTypes);

  }

public:

  BuilderImpl(Parent& glue_)
    : glue_(glue_),
      domelmntdescr_(NULL),
      tarelmntdescr_(NULL),
      domtrafo_(NULL),tartrafo_(NULL)
  {}


  void setDomainDescriptor(const DomainDescriptor& descr)
  {
    this->domelmntdescr_ = &descr;
  }


  void setTargetDescriptor(const TargetDescriptor& descr)
  {
    this->tarelmntdescr_ = &descr;
  }


  void setDomainTransformation(const typename Parent::DomainTransformation* trafo)
  {
    this->domtrafo_ = trafo;
  }


  void setTargetTransformation(const typename Parent::TargetTransformation* trafo)
  {
    this->tartrafo_ = trafo;
  }


  void build()
  {
    // setup the domain surface extractor
    if (this->domelmntdescr_ != NULL)
      this->glue_.domext_.update(*this->domelmntdescr_);
    else
      DUNE_THROW(Dune::Exception, "GridGlue::Builder : no domain surface descriptor set");

    // setup the target surface extractor
    if (this->tarelmntdescr_ != NULL)
      this->glue_.tarext_.update(*this->tarelmntdescr_);
    else
      DUNE_THROW(Dune::Exception, "GridGlue::Builder : no target surface descriptor set");

    // clear the contents from the current intersections array
    {
      std::vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->glue_.NULL_INTERSECTION);
      this->glue_.intersections_.swap(dummy);
    }

    std::vector<Dune::FieldVector<typename Parent::ctype, Parent::dimworld> > domcoords;
    std::vector<unsigned int> domfaces;
    std::vector<Dune::GeometryType> domainElementTypes;
    std::vector<Dune::FieldVector<typename Parent::ctype, Parent::dimworld> > tarcoords;
    std::vector<unsigned int> tarfaces;
    std::vector<Dune::GeometryType> targetElementTypes;

    /*
     * extract global surface grids
     */

    // retrieve the coordinate and topology information from the extractors
    // and apply transformations if necessary
    extractGrid(this->glue_.domext_, domcoords, domfaces, domainElementTypes, domtrafo_);
    extractGrid(this->glue_.tarext_, tarcoords, tarfaces, targetElementTypes, tartrafo_);

#ifdef WRITE_TO_VTK
    const int dimw = Parent::dimworld;
    const char prefix[] = "GridGlue::Builder::build() : ";
    char domainsurf[256];
    sprintf(domainsurf, "/tmp/vtk-domain-test-%i", mpi_rank);
    char targetsurf[256];
    sprintf(targetsurf, "/tmp/vtk-target-test-%i", mpi_rank);

    std::cout << prefix << "Writing domain surface to '" << domainsurf << ".vtk'...\n";
    VtkSurfaceWriter vtksw(domainsurf);
    vtksw.writeSurface(domcoords, domfaces, dimw, dimw);
    std::cout << prefix << "Done writing domain surface!\n";

    std::cout << prefix << "Writing target surface to '" << targetsurf << ".vtk'...\n";
    vtksw.setFilename(targetsurf);
    vtksw.writeSurface(tarcoords, tarfaces, dimw, dimw);
    std::cout << prefix << "Done writing target surface!\n";
#endif // WRITE_TO_VTK


    // start the actual build process
    this->glue_.merg_->build(domcoords, domfaces, domainElementTypes,
                             tarcoords, tarfaces, targetElementTypes);

    // the intersections need to be recomputed
    this->glue_.updateIntersections();
  }
};

#endif // GRIDGLUEBUILDERIMPL_HH_
