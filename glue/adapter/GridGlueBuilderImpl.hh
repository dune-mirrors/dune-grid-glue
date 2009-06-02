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


using namespace std;
using namespace ExtractorClassification;


/** \brief   I M P L E M E N T A T I O N   O F   S U B C L A S S   BUILDER IMPL
   (Specialization for two "non-surface" grids, i.e. meshes, manifolds) */

template<typename GET1, typename GET2, typename SM>
template<typename BGET1, typename BGET2, typename BSM, ExtractorType type1, ExtractorType type2>
class GridGlue<GET1, GET2, SM>::BuilderImpl
{
private:

  typedef GridGlue<GET1, GET2, SM>  Parent;

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename Parent::ctype ctype;


private:

  Parent& _glue;

  const ElementDescriptor<DomainGridView>*  _domelmntdescr;

  const ElementDescriptor<TargetGridView>*  _tarelmntdescr;

  const typename Parent::Transformation*    _domtrafo;

  const typename Parent::Transformation*    _tartrafo;


public:

  BuilderImpl(Parent& glue_)
    : _glue(glue_),
      _domelmntdescr(NULL),
      _tarelmntdescr(NULL),
      _domtrafo(NULL),_tartrafo(NULL)
  {}


  void setDomainElementDescriptor(const ElementDescriptor<DomainGridView>& descr)
  {
    this->_domelmntdescr = &descr;
  }


  void setTargetElementDescriptor(const ElementDescriptor<TargetGridView>& descr)
  {
    this->_tarelmntdescr = &descr;
  }


  void setDomainTransformation(const typename Parent::Transformation* trafo)
  {
    this->_domtrafo = trafo;
  }


  void setTargetTransformation(const typename Parent::Transformation* trafo)
  {
    this->_tartrafo = trafo;
  }


  bool build()
  {
    try
    {
      // extract the domain surface
      if (this->_domelmntdescr != NULL)
        this->_glue._domext.update(*this->_domelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no domain surface descriptor set" << endl;
        return false;
      }

      // extract the target surface
      if (this->_tarelmntdescr != NULL)
        this->_glue._tarext.update(*this->_tarelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no target surface descriptor set" << endl;
        return false;
      }

      // clear the contents from the current intersections array
      {
        vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->_glue.NULL_INTERSECTION);
        this->_glue._intersections.swap(dummy);
      }

      vector<typename Parent::ctype> domcoords;
      vector<unsigned int> domfaces;
      vector<typename Parent::ctype> tarcoords;
      vector<unsigned int> tarfaces;

      // retrieve the coordinate and topology information from the extractors
      // and apply transformations if necessary
      {
        vector<typename Parent::Coords> tempcoords;
        vector<typename Parent::DomainExtractor::SimplexTopology> tempfaces;

        this->_glue._domext.getCoords(tempcoords);
        domcoords.resize(Parent::dimw*tempcoords.size());
        typename Parent::ctype* currentc = &domcoords[0];
        if (this->_domtrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_domtrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._domext.getFaces(tempfaces);
        domfaces.resize(Parent::DomainExtractor::simplex_corners*tempfaces.size());
        unsigned int* currentf = &domfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::DomainExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];

        this->_glue._tarext.getCoords(tempcoords);
        tarcoords.resize(Parent::dimw*tempcoords.size());
        currentc = &tarcoords[0];
        if (this->_tartrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_tartrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._tarext.getFaces(tempfaces);
        tarfaces.resize(Parent::TargetExtractor::simplex_corners*tempfaces.size());
        currentf = &tarfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::TargetExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];
      }


#ifdef WRITE_TO_VTK
      const int dimw = Parent::dimw;
      const char prefix[] = "Builder<X, X>::build() : ";
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
      this->_glue._sm.build(domcoords, domfaces, tarcoords, tarfaces);

      // the intersections need to be recomputed
      this->_glue.updateIntersections();

      // success depends on whether the merged grid is empty or not
      return this->_glue._sm.nSimplices() != 0;
    }
    catch (Dune::Exception &e)
    {
      std::cerr << "GridGlue::Builder::build : Dune reported error: " << e << std::endl;
    }
    catch (...)
    {
      std::cerr << "GridGlue::Builder::build : Unknown exception occurred!" << endl;

    }
    // reaching this point is only possible after an exception has been thrown
    return false;
  }
};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   BUILDER IMPL   */
/*   (Specialization for two "surface" grids)                             */

template<typename GET1, typename GET2, typename SM>
template<typename BGET1, typename BGET2, typename BSM>
class GridGlue<GET1, GET2, SM>::BuilderImpl<
    BGET1,
    BGET2,
    BSM,
    ExtractorClassification::surface,
    ExtractorClassification::surface
    >
{
private:

  typedef GridGlue<GET1, GET2, SM>  Parent;

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename Parent::ctype ctype;


  Parent& _glue;

  const ElementDescriptor<DomainGridView>*  _domelmntdescr;

  const FaceDescriptor<DomainGridView>*     _domfacedescr;

  const ElementDescriptor<TargetGridView>*  _tarelmntdescr;

  const FaceDescriptor<TargetGridView>*     _tarfacedescr;


  const typename Parent::Transformation*    _domtrafo;

  const typename Parent::Transformation*    _tartrafo;


public:

  BuilderImpl(Parent& glue_)
    : _glue(glue_),
      _domelmntdescr(NULL), _domfacedescr(NULL),
      _tarelmntdescr(NULL), _tarfacedescr(NULL),
      _domtrafo(NULL),_tartrafo(NULL)
  {}


  void setDomainFaceDescriptor(const FaceDescriptor<DomainGridView>& descr)
  {
    this->_domelmntdescr = NULL;
    this->_domfacedescr = &descr;
  }


  void setDomainElementDescriptor(const ElementDescriptor<DomainGridView>& descr)
  {
    this->_domfacedescr = NULL;
    this->_domelmntdescr = &descr;
  }


  void setTargetFaceDescriptor(const FaceDescriptor<TargetGridView>& descr)
  {
    this->_tarelmntdescr = NULL;
    this->_tarfacedescr = &descr;
  }


  void setTargetElementDescriptor(const ElementDescriptor<TargetGridView>& descr)
  {
    this->_tarfacedescr = NULL;
    this->_tarelmntdescr = &descr;
  }


  void setDomainTransformation(const typename Parent::Transformation* trafo)
  {
    this->_domtrafo = trafo;
  }


  void setTargetTransformation(const typename Parent::Transformation* trafo)
  {
    this->_tartrafo = trafo;
  }


  bool build()
  {
    try
    {
      // extract the domain surface
      if (this->_domfacedescr != NULL)
        this->_glue._domext.update(*this->_domfacedescr);
      else if (this->_domelmntdescr != NULL)
        this->_glue._domext.update(*this->_domelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no domain surface descriptor set" << endl;
        return false;
      }

      // extract the target surface
      if (this->_tarfacedescr != NULL)
        this->_glue._tarext.update(*this->_tarfacedescr);
      else if (this->_tarelmntdescr != NULL)
        this->_glue._tarext.update(*this->_tarelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no target surface descriptor set" << endl;
        return false;
      }

      // clear the contents from the current intersections array
      {
        vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->_glue.NULL_INTERSECTION);
        this->_glue._intersections.swap(dummy);
      }

      vector<typename Parent::ctype> domcoords;
      vector<unsigned int> domfaces;
      vector<typename Parent::ctype> tarcoords;
      vector<unsigned int> tarfaces;


      // retrieve the coordinate and topology information from the extractors
      // and apply transformations if necessary
      {
        vector<typename Parent::Coords> tempcoords;
        vector<typename Parent::DomainExtractor::SimplexTopology> tempfaces;

        this->_glue._domext.getCoords(tempcoords);
        domcoords.resize(Parent::dimw*tempcoords.size());
        typename Parent::ctype* currentc = &domcoords[0];
        if (this->_domtrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_domtrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._domext.getFaces(tempfaces);
        domfaces.resize(Parent::DomainExtractor::simplex_corners*tempfaces.size());
        unsigned int* currentf = &domfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::DomainExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];

        this->_glue._tarext.getCoords(tempcoords);
        tarcoords.resize(Parent::dimw*tempcoords.size());
        currentc = &tarcoords[0];
        if (this->_tartrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_tartrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._tarext.getFaces(tempfaces);
        tarfaces.resize(Parent::TargetExtractor::simplex_corners*tempfaces.size());
        currentf = &tarfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::TargetExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];
      }

#ifdef WRITE_TO_VTK
      const int dimw = Parent::dimw;
      const char prefix[] = "Builder<surface, surface::build() : ";
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
      this->_glue._sm.build(domcoords, domfaces, tarcoords, tarfaces);

      // the intersections need to be recomputed
      this->_glue.updateIntersections();

      // success depends on whether the merged grid is empty or not
      return this->_glue._sm.nSimplices() != 0;
    }
    catch (Dune::Exception &e)
    {
      std::cerr << "GridGlue::Builder::build : Dune reported error: " << e << std::endl;
    }
    catch (...)
    {
      std::cerr << "GridGlue::Builder::build : Unknown exception occurred!" << endl;

    }
    // reaching this point is only possible after an exception has been thrown
    return false;
  }
};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   BUILDER IMPL    */
/*   (Specialization for a "surface" and a "non-surface" grid (order!)     */

template<typename GET1, typename GET2, typename SM>
template<typename BGET1, typename BGET2, typename BSM, ExtractorType type2>
class GridGlue<GET1, GET2, SM>::BuilderImpl<BGET1, BGET2, BSM, ExtractorClassification::surface, type2>
{
private:

  typedef GridGlue<GET1, GET2, SM>  Parent;

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename Parent::ctype ctype;


  Parent& _glue;

  const ElementDescriptor<DomainGridView>*  _domelmntdescr;

  const FaceDescriptor<DomainGridView>*     _domfacedescr;

  const ElementDescriptor<TargetGridView>*  _tarelmntdescr;

  const typename Parent::Transformation*    _domtrafo;

  const typename Parent::Transformation*    _tartrafo;


public:

  BuilderImpl(Parent& glue_)
    : _glue(glue_),
      _domelmntdescr(NULL), _domfacedescr(NULL),
      _tarelmntdescr(NULL),
      _domtrafo(NULL),_tartrafo(NULL)
  {}


  void setDomainFaceDescriptor(const FaceDescriptor<DomainGridView>& descr)
  {
    this->_domelmntdescr = NULL;
    this->_domfacedescr = &descr;
  }


  void setDomainElementDescriptor(const ElementDescriptor<DomainGridView>& descr)
  {
    this->_domfacedescr = NULL;
    this->_domelmntdescr = &descr;
  }


  void setTargetElementDescriptor(const ElementDescriptor<TargetGridView>& descr)
  {
    this->_tarelmntdescr = &descr;
  }


  void setDomainTransformation(const typename Parent::Transformation* trafo)
  {
    this->_domtrafo = trafo;
  }


  void setTargetTransformation(const typename Parent::Transformation* trafo)
  {
    this->_tartrafo = trafo;
  }


  bool build()
  {
    try
    {
      // extract the domain surface
      if (this->_domfacedescr != NULL)
        this->_glue._domext.update(*this->_domfacedescr);
      else if (this->_domelmntdescr != NULL)
        this->_glue._domext.update(*this->_domelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no domain surface descriptor set" << endl;
        return false;
      }

      // extract the target surface
      if (this->_tarelmntdescr != NULL)
        this->_glue._tarext.update(*this->_tarelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no target surface descriptor set" << endl;
        return false;
      }

      // clear the contents from the current intersections array
      {
        vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->_glue.NULL_INTERSECTION);
        this->_glue._intersections.swap(dummy);
      }

      vector<typename Parent::ctype> domcoords;
      vector<unsigned int> domfaces;
      vector<typename Parent::ctype> tarcoords;
      vector<unsigned int> tarfaces;

      // retrieve the coordinate and topology information from the extractors
      // and apply transformations if necessary
      {
        vector<typename Parent::Coords> tempcoords;
        vector<typename Parent::DomainExtractor::SimplexTopology> tempfaces;

        this->_glue._domext.getCoords(tempcoords);
        domcoords.resize(Parent::dimw*tempcoords.size());
        typename Parent::ctype* currentc = &domcoords[0];
        if (this->_domtrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_domtrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._domext.getFaces(tempfaces);
        domfaces.resize(Parent::DomainExtractor::simplex_corners*tempfaces.size());
        unsigned int* currentf = &domfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::DomainExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];

        this->_glue._tarext.getCoords(tempcoords);
        tarcoords.resize(Parent::dimw*tempcoords.size());
        currentc = &tarcoords[0];
        if (this->_tartrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_tartrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._tarext.getFaces(tempfaces);
        tarfaces.resize(Parent::TargetExtractor::simplex_corners*tempfaces.size());
        currentf = &tarfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::TargetExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];
      }


#ifdef WRITE_TO_VTK
      const int dimw = Parent::dimw;
      const char prefix[] = "Builder<surface, X>::build() : ";
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
      this->_glue._sm.build(domcoords, domfaces, tarcoords, tarfaces);

      // the intersections need to be recomputed
      this->_glue.updateIntersections();

      // success depends on whether the merged grid is empty or not
      return this->_glue._sm.nSimplices() != 0;
    }
    catch (Dune::Exception &e)
    {
      std::cerr << "GridGlue::Builder::build : Dune reported error: " << e << std::endl;
    }
    catch (...)
    {
      std::cerr << "GridGlue::Builder::build : Unknown exception occurred!" << endl;

    }
    // reaching this point is only possible after an exception has been thrown
    return false;
  }
};


/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   BUILDER IMPL   */
/*  (Specialization for a "non-surface" and a "surface" grid (order!)     */

template<typename GET1, typename GET2, typename SM>
template<typename BGET1, typename BGET2, typename BSM, ExtractorType type1>
class GridGlue<GET1, GET2, SM>::BuilderImpl<BGET1, BGET2, BSM, type1, ExtractorClassification::surface>
{
private:

  typedef GridGlue<GET1, GET2, SM>  Parent;

  typedef typename Parent::DomainGridView DomainGridView;

  typedef typename Parent::TargetGridView TargetGridView;

  typedef typename Parent::ctype ctype;


  Parent& _glue;

  const ElementDescriptor<DomainGridView>*  _domelmntdescr;

  const ElementDescriptor<TargetGridView>*  _tarelmntdescr;

  const FaceDescriptor<TargetGridView>*     _tarfacedescr;

  const typename Parent::Transformation*    _domtrafo;

  const typename Parent::Transformation*    _tartrafo;


public:

  BuilderImpl(Parent& glue_)
    : _glue(glue_),
      _domelmntdescr(NULL),
      _tarelmntdescr(NULL), _tarfacedescr(NULL),
      _domtrafo(NULL),_tartrafo(NULL)
  {}


  void setDomainElementDescriptor(const ElementDescriptor<DomainGridView>& descr)
  {
    this->_domelmntdescr = &descr;
  }


  void setTargetFaceDescriptor(const FaceDescriptor<TargetGridView>& descr)
  {
    this->_tarelmntdescr = NULL;
    this->_tarfacedescr = &descr;
  }


  void setTargetElementDescriptor(const ElementDescriptor<TargetGridView>& descr)
  {
    this->_tarfacedescr = NULL;
    this->_tarelmntdescr = &descr;
  }


  void setDomainTransformation(const typename Parent::Transformation* trafo)
  {
    this->_domtrafo = trafo;
  }


  void setTargetTransformation(const typename Parent::Transformation* trafo)
  {
    this->_tartrafo = trafo;
  }


  bool build()
  {
    try
    {
      // extract the domain surface
      if (this->_domelmntdescr != NULL)
        this->_glue._domext.update(*this->_domelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no domain surface descriptor set" << endl;
        return false;
      }

      // extract the target surface
      if (this->_tarfacedescr != NULL)
        this->_glue._tarext.update(*this->_tarfacedescr);
      else if (this->_tarelmntdescr != NULL)
        this->_glue._tarext.update(*this->_tarelmntdescr);
      else
      {
        cerr << "GridGlue::Builder : no target surface descriptor set" << endl;
        return false;
      }

      // clear the contents from the current intersections array
      {
        vector<typename Parent::RemoteIntersectionImpl> dummy(0, this->_glue.NULL_INTERSECTION);
        this->_glue._intersections.swap(dummy);
      }

      vector<typename Parent::ctype> domcoords;
      vector<unsigned int> domfaces;
      vector<typename Parent::ctype> tarcoords;
      vector<unsigned int> tarfaces;

      // retrieve the coordinate and topology information from the extractors
      // and apply transformations if necessary
      {
        vector<typename Parent::Coords> tempcoords;
        vector<typename Parent::DomainExtractor::SimplexTopology> tempfaces;

        this->_glue._domext.getCoords(tempcoords);
        domcoords.resize(Parent::dimw*tempcoords.size());
        typename Parent::ctype* currentc = &domcoords[0];
        if (this->_domtrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_domtrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._domext.getFaces(tempfaces);
        domfaces.resize(Parent::DomainExtractor::simplex_corners*tempfaces.size());
        unsigned int* currentf = &domfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::DomainExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];

        this->_glue._tarext.getCoords(tempcoords);
        tarcoords.resize(Parent::dimw*tempcoords.size());
        currentc = &tarcoords[0];
        if (this->_tartrafo != NULL)
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
          {
            typename Parent::Coords temp = (*this->_tartrafo)(tempcoords[i]);
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = temp[j];
          }
        }
        else
        {
          for (unsigned int i = 0; i < tempcoords.size(); ++i)
            for (int j = 0; j < Parent::dimw; ++j, ++currentc)
              *currentc = tempcoords[i][j];
        }

        this->_glue._tarext.getFaces(tempfaces);
        tarfaces.resize(Parent::TargetExtractor::simplex_corners*tempfaces.size());
        currentf = &tarfaces[0];
        for (unsigned int i = 0; i < tempfaces.size(); ++i)
          for (int j = 0; j < Parent::TargetExtractor::simplex_corners; ++j, ++currentf)
            *currentf = tempfaces[i][j];
      }


#ifdef WRITE_TO_VTK
      const int dimw = Parent::dimw;
      const char prefix[] = "Builder<X, surface>::build() : ";
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
      this->_glue._sm.build(domcoords, domfaces, tarcoords, tarfaces);

      // the intersections need to be recomputed
      this->_glue.updateIntersections();

      // success depends on whether the merged grid is empty or not
      return this->_glue._sm.nSimplices() != 0;
    }
    catch (Dune::Exception &e)
    {
      std::cerr << "GridGlue::Builder::build : Dune reported error: " << e << std::endl;
    }
    catch (...)
    {
      std::cerr << "GridGlue::Builder::build : Unknown exception occurred!" << endl;

    }
    // reaching this point is only possible after an exception has been thrown
    return false;
  }
};


#endif // GRIDGLUEBUILDERIMPL_HH_
