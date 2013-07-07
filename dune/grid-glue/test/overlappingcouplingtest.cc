// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/sgrid.hh>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/common/mpihelper.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/merging/conformingmerge.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

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


template <int dim>
void testCubeGrids(Merger<double,dim,dim,dim>& merger, const FieldVector<double,dim>& gridOffset)
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(elements, lower, upper);

  lower += gridOffset;
  upper += gridOffset;

  GridType grid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  DomExtractor domEx(grid0.leafView(), domdesc);
  TarExtractor tarEx(grid1.leafView(), tardesc);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim>
void testSimplexGrids(Merger<double,dim,dim,dim>& merger, const FieldVector<double,dim>& gridOffset)
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(elements, lower, upper);

  lower += gridOffset;
  upper += gridOffset;

  GridType grid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  DomExtractor domEx(grid0.leafView(), domdesc);
  TarExtractor tarEx(grid1.leafView(), tardesc);

  ::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


#if HAVE_UG
template <int dim>
void testSimplexGridsUG(Merger<double,dim,dim,dim>& merger, const FieldVector<double,dim>& gridOffset)
{
  // /////////////////////////////////////////////////////////////////
  //   Make two triangle grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  FieldVector<double,dim> lowerLeft(0);
  FieldVector<double,dim> upperRight(1);
  array<unsigned int, dim> elements;
  std::fill(elements.begin(), elements.end(), 10);

  StructuredGridFactory<GridType> factory;
  shared_ptr<GridType> grid0 = factory.createSimplexGrid(lowerLeft, upperRight, elements);

  lowerLeft  += gridOffset;
  upperRight += gridOffset;
  shared_ptr<GridType> grid1 = factory.createSimplexGrid(lowerLeft, upperRight, elements);

  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  DomExtractor domEx(grid0->leafView(), domdesc);
  TarExtractor tarEx(grid1->leafView(), tardesc);

  ::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim>
void testHybridGridsUG(Merger<double,dim,dim,dim>& merger, const FieldVector<double,dim>& gridOffset)
{
  // /////////////////////////////////////////////////////////////////////////
  //   Create the hybrid test grid from dune-grid twice and shift it once
  //   wrt the other grid.
  // /////////////////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  std::auto_ptr<Dune::UGGrid<dim> > grid0(make2DHybridTestGrid<Dune::UGGrid<dim> >());
  std::auto_ptr<Dune::UGGrid<dim> > grid1(make2DHybridTestGrid<Dune::UGGrid<dim> >());


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  DomExtractor domEx(grid0->leafView(), domdesc);
  TarExtractor tarEx(grid1->leafView(), tardesc);

  ::GridGlue<DomExtractor,TarExtractor> glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}
#endif


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  // //////////////////////////////////////////////////////////
  //   Test with the OverlappingMerge implementation
  // //////////////////////////////////////////////////////////
  OverlappingMerge<1,double> overlappingMerge1d;
  OverlappingMerge<2,double> overlappingMerge2d;

  testCubeGrids<1>(overlappingMerge1d, FieldVector<double,1>(0.05));
  testCubeGrids<2>(overlappingMerge2d, FieldVector<double,2>(0.05));

  testSimplexGrids<1>(overlappingMerge1d, FieldVector<double,1>(0.05));
#if HAVE_UG
  testSimplexGridsUG(overlappingMerge2d, FieldVector<double,2>(0.05));
  testHybridGridsUG<2>(overlappingMerge2d, FieldVector<double,2>(0.05));
#endif

  // //////////////////////////////////////////////////////////
  //   Test with the PSurfaceMerge implementation
  // //////////////////////////////////////////////////////////
#if HAVE_PSURFACE
  PSurfaceMerge<1,1,double> psurfaceMerge1d;
  PSurfaceMerge<2,2,double> psurfaceMerge2d;

  testCubeGrids<1>(psurfaceMerge1d, FieldVector<double,1>(0.05));
  testCubeGrids<2>(psurfaceMerge2d, FieldVector<double,2>(0.05));

  testSimplexGrids<1>(psurfaceMerge1d, FieldVector<double,1>(0.05));
#if HAVE_UG
  testSimplexGridsUG(psurfaceMerge2d, FieldVector<double,2>(0.05));
  testHybridGridsUG<2>(psurfaceMerge2d, FieldVector<double,2>(0.05));
#endif
#endif


  // //////////////////////////////////////////////////////////
  //   Test with the ConformingMerge implementation
  // //////////////////////////////////////////////////////////

  ConformingMerge<1,1,double> conformingMerge1d;
  ConformingMerge<2,2,double> conformingMerge2d;
  ConformingMerge<3,3,double> conformingMerge3d;

  testCubeGrids<1>(conformingMerge1d, FieldVector<double,1>(0));
  testCubeGrids<2>(conformingMerge2d, FieldVector<double,2>(0));

  testSimplexGrids<1>(conformingMerge1d, FieldVector<double,1>(0));
#if HAVE_UG
  testSimplexGridsUG(conformingMerge2d, FieldVector<double,2>(0));

  testHybridGridsUG<2>(conformingMerge2d, FieldVector<double,2>(0));
#endif
}
