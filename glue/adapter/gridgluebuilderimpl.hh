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

  Parent& _glue;

  typedef typename ExtractorTypeTraits<DomainGridView, codim1>::Type DomainDescriptor;

  typedef typename ExtractorTypeTraits<TargetGridView, codim2>::Type TargetDescriptor;

  const DomainDescriptor*                   _domelmntdescr;

  const TargetDescriptor*                   _tarelmntdescr;

  const typename Parent::DomainTransformation*    _domtrafo;

  const typename Parent::TargetTransformation*    _tartrafo;

  template<typename Extractor>
  void extractGrid (const Extractor & ext,
                    std::vector<typename Parent::ctype> & coords,
                    std::vector<unsigned int> & faces,
                    //const typename Parent::Transformation* trafo) const
                    const CoordinateTransformation<Extractor::dimworld, Parent::dimworld, typename Parent::ctype>* trafo) const
  {
    std::vector<typename Extractor::Coords> tempcoords;
    std::vector<typename Parent::DomainExtractor::SimplexTopology> tempfaces;

    ext.getCoords(tempcoords);
    coords.clear();
    coords.reserve(Parent::dimworld*tempcoords.size());

    if (trafo != NULL)
    {
      std::cout << "GridGlue::Builder : apply trafo\n";
      for (size_t i = 0; i < tempcoords.size(); ++i)
      {
        typename Parent::Coords temp = (*trafo)(tempcoords[i]);
        for (size_t j = 0; j < Parent::dimworld; ++j)
          coords.push_back(temp[j]);
      }
    }
    else
    {
      for (unsigned int i = 0; i < tempcoords.size(); ++i)
      {
        for (size_t j = 0; j < Parent::dimworld; ++j)
          coords.push_back(tempcoords[i][j]);
      }
    }

    ext.getFaces(tempfaces);
    faces.clear();
    faces.reserve(Parent::DomainExtractor::simplex_corners*tempfaces.size());
    for (unsigned int i = 0; i < tempfaces.size(); ++i)
      for (int j = 0; j < Parent::DomainExtractor::simplex_corners; ++j)
        faces.push_back(tempfaces[i][j]);
  }

public:

  BuilderImpl(Parent& glue_)
    : _glue(glue_),
      _domelmntdescr(NULL),
      _tarelmntdescr(NULL),
      _domtrafo(NULL),_tartrafo(NULL)
  {}


  void setDomainDescriptor(const DomainDescriptor& descr)
  {
    this->_domelmntdescr = &descr;
  }


  void setTargetDescriptor(const TargetDescriptor& descr)
  {
    this->_tarelmntdescr = &descr;
  }


  void setDomainTransformation(const typename Parent::DomainTransformation* trafo)
  {
    this->_domtrafo = trafo;
  }


  void setTargetTransformation(const typename Parent::TargetTransformation* trafo)
  {
    this->_tartrafo = trafo;
  }


  void build()
  {
    // setup the domain surface extractor
    if (this->_domelmntdescr != NULL)
      this->_glue._domext.update(*this->_domelmntdescr);
    else
      DUNE_THROW(Dune::Exception, "GridGlue::Builder : no domain surface descriptor set");

    // setup the target surface extractor
    if (this->_tarelmntdescr != NULL)
      this->_glue._tarext.update(*this->_tarelmntdescr);
    else
      DUNE_THROW(Dune::Exception, "GridGlue::Builder : no target surface descriptor set");

    // clear the contents from the current intersections array
    {
      std::vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->_glue.NULL_INTERSECTION);
      this->_glue._intersections.swap(dummy);
    }

    std::vector<typename Parent::ctype> domcoords;
    std::vector<unsigned int> domfaces;
    std::vector<typename Parent::ctype> tarcoords;
    std::vector<unsigned int> tarfaces;

    /*
     * extract global surface grids
     */

    // retrieve the coordinate and topology information from the extractors
    // and apply transformations if necessary
    extractGrid(this->_glue._domext, domcoords, domfaces, _domtrafo);
    extractGrid(this->_glue._tarext, tarcoords, tarfaces, _tartrafo);

#ifdef WRITE_TO_VTK
    const int dimw = Parent::dimw;
    const char prefix[] = "GridGlue::Builder::build() : ";
    const char domainsurf[] = "/tmp/vtk-domain-test";
    const char targetsurf[] = "/tmp/vtk-target-test";

    STDOUTLN(prefix << "Writing domain surface to '" << domainsurf << ".vtk'...");
    VtkSurfaceWriter vtksw(domainsurf);
    vtksw.writeSurface(domcoords, domfaces, dimw, dimw);
    STDOUTLN(prefix << "Done writing domain surface!");

    STDOUTLN(prefix << "Writing target surface to '" << targetsurf << ".vtk'...");
    vtksw.setFilename(targetsurf);
    vtksw.writeSurface(tarcoords, tarfaces, dimw, dimw);
    STDOUTLN(prefix << "Done writing target surface!");
#endif // WRITE_TO_VTK

    // start the actual build process
    this->_glue._merg->build(domcoords, domfaces, tarcoords, tarfaces);

    // the intersections need to be recomputed
    this->_glue.updateIntersections();
  }
};

#endif // GRIDGLUEBUILDERIMPL_HH_
