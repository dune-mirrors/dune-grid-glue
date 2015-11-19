// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/geometrygrid.hh>

#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

/** \brief Returns always true */
template<typename GridView>
typename Codim0Extractor<GridView>::Predicate
makeTruePredicate()
{
  using Element = typename GridView::Traits::template Codim<0>::Entity;
  auto predicate = [](const Element&, unsigned int) -> bool {
    return true;
  };
  return predicate;
}

/** \brief trafo from dim to dim+1 */
template<int dim, int dimw, class ctype>
class MixedDimTrafo
  : public AnalyticalCoordFunction< ctype, dim, dimw, MixedDimTrafo<dim,dimw,ctype> >
{
  static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  static_assert(dim==1, "MixedDimTrafo currently assumes dim==1");
public:

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dimw> &y ) const
  {
    y[0] = x[0]+0.2;
    y[1] = x[0]+0.1;
  }
};


int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  OverlappingMerge<2,1,2> overlappingMerge;

  // /////////////////////////////////////////////////////////
  //   Make a 2d unit cube grid and a 1d grid embedded in 2d
  // /////////////////////////////////////////////////////////

  typedef SGrid<2,2> GridType2d;

  FieldVector<int, 2> elements(1);
  FieldVector<double,2> lower(0);
  FieldVector<double,2> upper(1);

  GridType2d grid0(elements, lower, upper);

  typedef SGrid<1,1> GridType1d;

  FieldVector<int, 1> elements1d(1);
  FieldVector<double,1> lower1d(0);
  FieldVector<double,1> upper1d(1);

  typedef GeometryGrid<GridType1d, MixedDimTrafo<1,2,double> > LiftedGridType;

  GridType1d cubeGrid1_in(elements1d, lower1d, upper1d);

  MixedDimTrafo<1,2,double> trafo;   // transform dim-1 to dim

  LiftedGridType grid1(cubeGrid1_in, trafo);

  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType2d::LeafGridView DomGridView;
  typedef typename LiftedGridType::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  const DomExtractor::Predicate domdesc = makeTruePredicate<DomGridView>();
  const TarExtractor::Predicate tardesc = makeTruePredicate<TarGridView>();

#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
  DomExtractor domEx(grid0.leafGridView(), domdesc);
  TarExtractor tarEx(grid1.leafGridView(), tardesc);
#else
  DomExtractor domEx(grid0.leafView(), domdesc);
  TarExtractor tarEx(grid1.leafView(), tardesc);
#endif
  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  // The following code is out-commented, because the test functionality
  // doesn't actually work yet.
  OverlappingMerge<2,1,2> merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}
