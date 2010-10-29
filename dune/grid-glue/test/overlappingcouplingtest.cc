// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/sgrid.hh>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/common/quadraturerules.hh>
#include <doc/grids/gridfactory/hybridtestgrids.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/gridextractiontraits.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

#if HAVE_PSURFACE
#include <dune/grid-glue/merging/psurfacemerge.hh>
#endif
#if HAVE_CGAL
#include <dune/grid-glue/merging/cgalmerge.hh>
#endif
#include <dune/grid-glue/merging/conformingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;

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

  typedef DefaultExtractionTraits<DomGridView,0> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  GlueType glue(grid0.leafView(), grid1.leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  glue.setDomainDescriptor(domdesc);
  glue.setTargetDescriptor(tardesc);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

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

  typedef DefaultExtractionTraits<DomGridView,0> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0> TarTraits;

  GridGlue<DomTraits,TarTraits> glue(grid0.leafView(), grid1.leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  glue.setDomainDescriptor(domdesc);
  glue.setTargetDescriptor(tardesc);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


#if HAVE_UG
void testTriangleGridsUG(Merger<double,2,2,2>& merger, const FieldVector<double,2>& gridOffset)
{
  const int dim = 2;

  // /////////////////////////////////////////////////////////////////
  //   Make two triangle grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef UGGrid<dim> GridType;

  GridFactory<GridType> factory0, factory1;

  // insert vertices
  for (int i=0; i<11; i++) {
    for (int j=0; j<11; j++) {

      FieldVector<double,2> pos;
      pos[0] = i * 0.1;
      pos[1] = j * 0.1;

      factory0.insertVertex(pos);

      pos += gridOffset;

      factory1.insertVertex(pos);

    }

  }


  // insert elements
  std::vector<unsigned int> vertices(3);

  for (int i=0; i<10; i++) {

    for (int j=0; j<10; j++) {

      unsigned int base = j*11 + i;

      vertices[0] = base;   vertices[1] = base+1;   vertices[2] = base + 11;
      factory0.insertElement(GeometryType(GeometryType::simplex,dim), vertices);
      factory1.insertElement(GeometryType(GeometryType::simplex,dim), vertices);

      vertices[0] = base+1;   vertices[1] = base+12;   vertices[2] = base + 11;
      factory0.insertElement(GeometryType(GeometryType::simplex,dim), vertices);
      factory1.insertElement(GeometryType(GeometryType::simplex,dim), vertices);

    }

  }

  GridType* grid0 = factory0.createGrid();
  GridType* grid1 = factory1.createGrid();


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef GridType::LeafGridView DomGridView;
  typedef GridType::LeafGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,0> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0> TarTraits;

  GridGlue<DomTraits,TarTraits> glue(grid0->leafView(), grid1->leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  glue.setDomainDescriptor(domdesc);
  glue.setTargetDescriptor(tardesc);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

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

  typedef DefaultExtractionTraits<DomGridView,0> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0> TarTraits;

  GridGlue<DomTraits,TarTraits> glue(grid0->leafView(), grid1->leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  glue.setDomainDescriptor(domdesc);
  glue.setTargetDescriptor(tardesc);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}
#endif


int main()
{
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
  testTriangleGridsUG(psurfaceMerge2d, FieldVector<double,2>(0.05));
  testHybridGridsUG<2>(psurfaceMerge2d, FieldVector<double,2>(0.05));
#endif
#endif

  // //////////////////////////////////////////////////////////
  //   Test with the CGALMerge implementation
  // //////////////////////////////////////////////////////////

#if HAVE_CGAL
  CGALMerge<1,double> cgalMerge1d;
  CGALMerge<2,double> cgalMerge2d;

  testCubeGrids<1>(cgalMerge1d, FieldVector<double,1>(0.05));
  testCubeGrids<2>(cgalMerge2d, FieldVector<double,2>(0.05));

  testSimplexGrids<1>(cgalMerge1d, FieldVector<double,1>(0.05));
#if HAVE_UG
  testTriangleGridsUG(cgalMerge2d, FieldVector<double,2>(0.05));

  testHybridGridsUG<2>(cgalMerge2d, FieldVector<double,2>(0.05));
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
  testTriangleGridsUG(conformingMerge2d, FieldVector<double,2>(0));

  testHybridGridsUG<2>(conformingMerge2d, FieldVector<double,2>(0));
#endif
}
