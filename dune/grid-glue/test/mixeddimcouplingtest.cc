// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/extractors/parallelextractor.hh>
#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;

template <class GridView>
class HorizontalFaceDescriptor
  : public ExtractorPredicate<GridView,1>
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
  : public ExtractorPredicate<GridView,0>
{
public:
  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& element, unsigned int subentity) const
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



template <int dim>
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

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  HorizontalFaceDescriptor<DomGridView> domdesc(0);
  AllElementsDescriptor<TarGridView>  tardesc;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
  tarEx.positiveNormalDirection() = (slice == 0.0);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}


template <int dim>
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

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(0);

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  domEx.positiveNormalDirection() = (slice == 0.0);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}


template <int dim, bool par=false>
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

  // typedef DefaultExtractionTraits<DomGridView,1, par> DomTraits;
  // typedef DefaultExtractionTraits<TarGridView,0, par> TarTraits;
  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  HorizontalFaceDescriptor<DomGridView> domdesc(slice);
  AllElementsDescriptor<TarGridView>  tardesc;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
  tarEx.positiveNormalDirection() = (slice == 0.0);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim
  glue.setTargetTransformation(&trafo);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue, (CoordinateTransformation<dim,dim,double>*)nullptr, &trafo);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}


template <int dim, bool par=false>
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

  // typedef DefaultExtractionTraits<DomGridView,0, par> DomTraits;
  // typedef DefaultExtractionTraits<TarGridView,1, par> TarTraits;
  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(slice);

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  domEx.positiveNormalDirection() = (slice == 0.0);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim
  glue.setDomainTransformation(&trafo);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue, &trafo, (CoordinateTransformation<dim,dim,double>*)nullptr);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}

int main(int argc, char *argv[]) try
{
  // /////////////////////////////////////////////////////////////
  //   First set of tests: the grid have different dimensions,
  //   but the world dimension is the same for both of them.
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d == matching =============================\n";
  test1d2dCouplingMatchingDimworld<2>();
  std::cout << "====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == matching =============================\n";
  test2d1dCouplingMatchingDimworld<2>();
  std::cout << "====================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == matching =============================\n";
  test2d1dCouplingMatchingDimworld<3>();
  std::cout << "====================================================\n";

  // /////////////////////////////////////////////////////////////
  //   Second set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==0.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d === nonmatching ==========================\n";
  test1d2dCoupling<2>();
  test1d2dCoupling<2, true>();
  std::cout << "=====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching ==========================\n";
  test2d1dCoupling<2>();
  test2d1dCoupling<2, true>();
  std::cout << "====================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching ==========================\n";
  test2d1dCoupling<3>();
  test2d1dCoupling<3, true>();
  std::cout << "====================================================\n";

  // /////////////////////////////////////////////////////////////
  //   Third set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==1.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d == nonmatching top ======================\n";
  test1d2dCoupling<2>(1.0);
  test1d2dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching top ======================\n";
  test2d1dCoupling<2>(1.0);
  test2d1dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching top ======================\n";
  test2d1dCoupling<3>(1.0);
  test2d1dCoupling<3, true>(1.0);
  std::cout << "====================================================\n";

}
catch (Exception e) {
  std::cout << e << std::endl;
}
