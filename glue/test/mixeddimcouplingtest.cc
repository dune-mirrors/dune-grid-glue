// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/extractors/surfacedescriptor.hh>
#include <dune/glue/extractors/gridextractiontraits.hh>
#include <dune/glue/merging/psurfacemerge.hh>
#include <dune/glue/adapter/gridglue.hh>

#include <dune/glue/test/couplingtest.hh>

using namespace Dune;

template <class GridView>
class HorizontalFaceDescriptor
  : public FaceDescriptor<GridView>
{
public:
  HorizontalFaceDescriptor(double sliceCoord)
    : sliceCoord_(sliceCoord)
  {}

  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& eptr,
                        unsigned int face) const
  {
    const int dim = GridView::dimension;
    const Dune::GenericReferenceElement<double,dim>& refElement = Dune::GenericReferenceElements<double, dim>::general(eptr->type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(eptr->geometry().corner(refElement.subEntity(face,1,i,dim))[1] - sliceCoord_) > 1e-6 )
        return false;

    return true;
  }

private:
  double sliceCoord_;
};

/** \brief Returns always true */
template <class GridView>
class AllElementsDescriptor
  : public ElementDescriptor<GridView>
{
public:
  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& eptr) const
  {
    return true;
  }
};



template <int dim, MeshClassification::MeshType ExtractorClassification>
void test1d2dCouplingMatchingDimworld()
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid0(elements, lower, upper);

  typedef SGrid<dim-1,dim> GridType1d;

  FieldVector<int, dim-1> elements1d(1);
  FieldVector<double,dim-1> lower1d(0);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d cubeGrid1(elements1d, lower1d, upper1d);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  typedef typename GridType1d::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,1,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,ExtractorClassification> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  HorizontalFaceDescriptor<DomGridView> domdesc(0);
  AllElementsDescriptor<TarGridView>  tardesc;

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification>
void test2d1dCouplingMatchingDimworld()
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  typedef SGrid<dim-1,dim> GridType1d;

  FieldVector<int, dim-1> elements1d(1);
  FieldVector<double,dim-1> lower1d(0);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d cubeGrid0(elements1d, lower1d, upper1d);

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid1(elements, lower, upper);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType1d::LevelGridView DomGridView;
  typedef typename GridType2d::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(0);

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification>
void test1d2dCoupling()
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid0(elements, lower, upper);

  typedef SGrid<dim-1,dim-1> GridType1d;

  FieldVector<int, dim-1> elements1d(1);
  FieldVector<double,dim-1> lower1d(0);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d cubeGrid1(elements1d, lower1d, upper1d);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  typedef typename GridType1d::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,1,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,ExtractorClassification> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  HorizontalFaceDescriptor<DomGridView> domdesc(0);
  AllElementsDescriptor<TarGridView>  tardesc;

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification>
void test2d1dCoupling()
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  typedef SGrid<dim-1,dim-1> GridType1d;

  FieldVector<int, dim-1> elements1d(1);
  FieldVector<double,dim-1> lower1d(0);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d cubeGrid0(elements1d, lower1d, upper1d);

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid1(elements, lower, upper);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType1d::LevelGridView DomGridView;
  typedef typename GridType2d::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), 0 /*&merger*/);

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(0);

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}




int main(int argc, char *argv[]) try
{
  // /////////////////////////////////////////////////////////////
  //   First set of tests: the grid have different dimensions,
  //   but the world dimension is the same for both of them.
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  test1d2dCouplingMatchingDimworld<2,MeshClassification::cube>();
  test1d2dCouplingMatchingDimworld<2,MeshClassification::simplex>();

  // Test a unit square versus a grid one dimension lower
  test2d1dCouplingMatchingDimworld<2,MeshClassification::cube>();
  test2d1dCouplingMatchingDimworld<2,MeshClassification::simplex>();

  // Test a unit cube versus a grid one dimension lower
  test2d1dCouplingMatchingDimworld<3,MeshClassification::cube>();

  // /////////////////////////////////////////////////////////////
  //   Second set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  // /////////////////////////////////////////////////////////////

#if 0
  // Test a unit square versus a grid one dimension lower
  test1d2dCoupling<2,MeshClassification::cube>();
  test1d2dCoupling<2,MeshClassification::simplex>();

  // Test a unit square versus a grid one dimension lower
  test2d1dCoupling<2,MeshClassification::cube>();
  test2d1dCoupling<2,MeshClassification::simplex>();

  // Test a unit cube versus a grid one dimension lower
  test2d1dCoupling<3,MeshClassification::cube>();
#endif

}
catch (Exception e) {
  std::cout << e << std::endl;
}
