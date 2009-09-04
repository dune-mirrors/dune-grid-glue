// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/grid/sgrid.hh>
#ifdef HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/extractors/surfacedescriptor.hh>
#include <dune/glue/extractors/gridextractiontraits.hh>
#include <dune/glue/merging/psurfacemerge.hh>
#include <dune/glue/adapter/gridglue.hh>

#include <dune/glue/test/couplingtest.hh>

using namespace Dune;

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


template <int dim>
void testCubeGrids()
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(elements, lower, upper);

  lower += 0.05;
  upper += 0.05;

  GridType grid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,MeshClassification::cube> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim,dim,double> merger;

  GlueType glue(grid0.leafView(), grid1.leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

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


template <int dim>
void testSimplexGrids()
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(10);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType grid0(elements, lower, upper);

  lower += 0.05;
  upper += 0.05;

  GridType grid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType::LeafGridView DomGridView;
  typedef typename GridType::LeafGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::simplex> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,MeshClassification::simplex> TarTraits;

  PSurfaceMerge<dim,dim,double> merger;

  GridGlue<DomTraits,TarTraits> glue(grid0.leafView(), grid1.leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

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


void testTriangleGridsUG()
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

      pos += 0.05;

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

  typedef DefaultExtractionTraits<DomGridView,0,MeshClassification::simplex> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,MeshClassification::simplex> TarTraits;

  PSurfaceMerge<dim,dim,double> merger;

  GridGlue<DomTraits,TarTraits> glue(grid0->leafView(), grid1->leafView(), &merger);

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

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


int main()
{

  testCubeGrids<1>();
  testCubeGrids<2>();

  testSimplexGrids<1>();
#if HAVE_UG
  testTriangleGridsUG();
#endif
}
