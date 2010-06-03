// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>
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



/** \brief trafo from dim to dim+1 */
template<int dim, int dimw, class ctype>
class MixedDimTrafo : public CoordinateTransformation<dim,dimw,ctype>
{
  dune_static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  double yOffset_;
public:
  MixedDimTrafo(double yOffset) : yOffset_(yOffset) {}
  virtual Dune::FieldVector<ctype, dimw> operator()(const Dune::FieldVector<ctype, dim>& c) const
  {
    Dune::FieldVector<ctype, dimw> x(yOffset_);
    x[0] = c[0];
    for (int i=2; i<dimw ; i++)
      x[i] = c[i-1];
    return x;
  }
};



template <int dim, MeshClassification::MeshType ExtractorClassification>
void test1d2dCouplingMatchingDimworld()
{
  double slice = 0.0;

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
  glue.targetExtractor().positiveNormalDirection() = (slice == 0.0);

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
  double slice = 0.0;

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
  glue.domainExtractor().positiveNormalDirection() = (slice == 0.0);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification, bool par=false>
void test1d2dCoupling(double slice=0.0)
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

  typedef DefaultExtractionTraits<DomGridView,1,MeshClassification::cube, par> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,ExtractorClassification, par> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  HorizontalFaceDescriptor<DomGridView> domdesc(slice);
  AllElementsDescriptor<TarGridView>  tardesc;

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim
  glue.builder().setTargetTransformation(&trafo);
  glue.targetExtractor().positiveNormalDirection() = (slice == 0.0);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue, (CoordinateTransformation<dim,dim,double>*)nullptr, &trafo);

}


template <int dim, MeshClassification::MeshType ExtractorClassification, bool par=false>
void test2d1dCoupling(double slice=0.0)
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

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::cube, par> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification, par> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim-1,dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(slice);

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim
  glue.builder().setDomainTransformation(&trafo);
  glue.domainExtractor().positiveNormalDirection() = (slice == 0.0);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue, &trafo, (CoordinateTransformation<dim,dim,double>*)nullptr);

}

int main(int argc, char *argv[]) try
{
  // /////////////////////////////////////////////////////////////
  //   First set of tests: the grid have different dimensions,
  //   but the world dimension is the same for both of them.
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d cube ===== matching =============================\n";
  test1d2dCouplingMatchingDimworld<2,MeshClassification::cube>();
  std::cout << "==== 1d 2d simplex == matching =============================\n";
  test1d2dCouplingMatchingDimworld<2,MeshClassification::simplex>();
  std::cout << "============================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d cube ===== matching =============================\n";
  test2d1dCouplingMatchingDimworld<2,MeshClassification::cube>();
  std::cout << "==== 2d 1d simplex == matching =============================\n";
  test2d1dCouplingMatchingDimworld<2,MeshClassification::simplex>();
  std::cout << "============================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d simplex == matching =============================\n";
  test2d1dCouplingMatchingDimworld<3,MeshClassification::cube>();
  std::cout << "============================================================\n";

  // /////////////////////////////////////////////////////////////
  //   Second set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==0.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d cube ===== nonmatching ==========================\n";
  test1d2dCoupling<2,MeshClassification::cube>();
  test1d2dCoupling<2,MeshClassification::cube, true>();
  std::cout << "==== 1d 2d simplex=== nonmatching ==========================\n";
  test1d2dCoupling<2,MeshClassification::simplex>();
  test1d2dCoupling<2,MeshClassification::simplex, true>();
  std::cout << "============================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d cube ===== nonmatching ==========================\n";
  test2d1dCoupling<2,MeshClassification::cube>();
  test2d1dCoupling<2,MeshClassification::cube, true>();
  std::cout << "==== 2d 1d simplex == nonmatching ==========================\n";
  test2d1dCoupling<2,MeshClassification::simplex>();
  test2d1dCoupling<2,MeshClassification::simplex, true>();
  std::cout << "============================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d cube ===== nonmatching ==========================\n";
  test2d1dCoupling<3,MeshClassification::cube>();
  test2d1dCoupling<3,MeshClassification::cube, true>();
  std::cout << "============================================================\n";

  // /////////////////////////////////////////////////////////////
  //   Third set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==1.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d cube ===== nonmatching top ======================\n";
  test1d2dCoupling<2,MeshClassification::cube>(1.0);
  test1d2dCoupling<2,MeshClassification::cube, true>(1.0);
  std::cout << "==== 1d 2d simplex == nonmatching top ======================\n";
  test1d2dCoupling<2,MeshClassification::simplex>(1.0);
  test1d2dCoupling<2,MeshClassification::simplex, true>(1.0);
  std::cout << "============================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d cube ===== nonmatching top ======================\n";
  test2d1dCoupling<2,MeshClassification::cube>(1.0);
  test2d1dCoupling<2,MeshClassification::cube, true>(1.0);
  std::cout << "============================================================\n";
  std::cout << "==== 2d 1d simplex == nonmatching top ======================\n";
  test2d1dCoupling<2,MeshClassification::simplex>(1.0);
  test2d1dCoupling<2,MeshClassification::simplex, true>(1.0);
  std::cout << "============================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d cube ===== nonmatching top ======================\n";
  test2d1dCoupling<3,MeshClassification::cube>(1.0);
  test2d1dCoupling<3,MeshClassification::cube, true>(1.0);
  std::cout << "============================================================\n";

}
catch (Exception e) {
  std::cout << e << std::endl;
}
