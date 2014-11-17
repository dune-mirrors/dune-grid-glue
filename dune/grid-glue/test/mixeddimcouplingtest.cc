// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

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
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
    const Dune::ReferenceElement<double,dim>& refElement = Dune::ReferenceElements<double, dim>::general(eptr->type());
#else
    const Dune::GenericReferenceElement<double,dim>& refElement = Dune::GenericReferenceElements<double, dim>::general(eptr->type());
#endif

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
class MixedDimTrafo
  : public AnalyticalCoordFunction< ctype, dim, dimw, MixedDimTrafo<dim,dimw,ctype> >
{
  static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  double yOffset_;
public:
  MixedDimTrafo(double yOffset) : yOffset_(yOffset) {}

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dimw> &y ) const
  {
    y = yOffset_;
    y[0] = x[0];
    for (int i=2; i<dimw; i++)
      y[i] = x[i-1];
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

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
  DomExtractor domEx(cubeGrid0.levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
#endif
  tarEx.positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
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

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
  DomExtractor domEx(cubeGrid0.levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
#endif
  domEx.positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
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

  typedef GeometryGrid<GridType1d, MixedDimTrafo<dim-1,dim,double> > LiftedGridType;

  GridType1d cubeGrid1_in(elements1d, lower1d, upper1d);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim

  LiftedGridType cubeGrid1(cubeGrid1_in, trafo);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  // typedef typename GridType1d::LevelGridView TarGridView;
  typedef typename LiftedGridType::LevelGridView TarGridView;

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  HorizontalFaceDescriptor<DomGridView> domdesc(slice);
  AllElementsDescriptor<TarGridView>  tardesc;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
  tarEx.positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
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

  typedef GeometryGrid<GridType1d, MixedDimTrafo<dim-1,dim,double> > LiftedGridType;

  GridType1d cubeGrid0_in(elements1d, lower1d, upper1d);

  MixedDimTrafo<dim-1,dim,double> trafo(slice);   // transform dim-1 to dim

  LiftedGridType cubeGrid0(cubeGrid0_in, trafo);

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid1(elements, lower, upper);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename LiftedGridType::LevelGridView DomGridView;
  typedef typename GridType2d::LevelGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView>  domdesc;
  HorizontalFaceDescriptor<TarGridView> tardesc(slice);

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
  DomExtractor domEx(cubeGrid0.levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
#endif
  domEx.positiveNormalDirection() = (slice == 0.0);

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  PSurfaceMerge<dim-1,dim,double> merger;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#endif
}

int main(int argc, char *argv[]) try
{
#if !HAVE_PSURFACE
  exit(77); // Test is skipped, if PSurface is not present
#endif

  Dune::MPIHelper::instance(argc, argv);

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

#define _3DTEST 0
#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == matching =============================\n";
  test2d1dCouplingMatchingDimworld<3>();
  std::cout << "====================================================\n";
#endif

  // /////////////////////////////////////////////////////////////
  //   Second set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==0.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d === nonmatching ==========================\n";
  // test1d2dCoupling<2>();
  // test1d2dCoupling<2, true>();
  std::cout << "=====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching ==========================\n";
  test2d1dCoupling<2>();
  test2d1dCoupling<2, true>();
  std::cout << "====================================================\n";

#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching ==========================\n";
  test2d1dCoupling<3>();
  test2d1dCoupling<3, true>();
  std::cout << "====================================================\n";
#endif

  // /////////////////////////////////////////////////////////////
  //   Third set of tests: the grid have different dimensions,
  //   and the world dimension is different as well
  //   -- grids match at y==1.0
  // /////////////////////////////////////////////////////////////

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 1d 2d == nonmatching top ======================\n";
  // test1d2dCoupling<2>(1.0);
  // test1d2dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

  // Test a unit square versus a grid one dimension lower
  std::cout << "==== 2d 1d == nonmatching top ======================\n";
  test2d1dCoupling<2>(1.0);
  test2d1dCoupling<2, true>(1.0);
  std::cout << "====================================================\n";

#if _3DTEST
  // Test a unit cube versus a grid one dimension lower
  std::cout << "==== 3d 2d == nonmatching top ======================\n";
  test2d1dCoupling<3>(1.0);
  test2d1dCoupling<3, true>(1.0);
  std::cout << "====================================================\n";
#endif

}
catch (Exception e) {
  std::cout << e << std::endl;
  return 1;
}
