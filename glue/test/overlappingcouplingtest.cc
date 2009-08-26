// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <dune/grid/sgrid.hh>
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



int main()
{

  // /////////////////////////////////////////////////////////////////
  //   Make two cube grids that are slightly shifted wrt each other
  // /////////////////////////////////////////////////////////////////

  const int dim = 2;

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

  typedef GridType::LeafGridView DomGridView;
  typedef GridType::LeafGridView TarGridView;

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
