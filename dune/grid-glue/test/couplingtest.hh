// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef GRIDGLUE_COUPLINGTEST_HH
#define GRIDGLUE_COUPLINGTEST_HH

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

template <class IntersectionIt>
void testIntersection(const IntersectionIt & rIIt)
{
  typedef typename IntersectionIt::value_type Intersection;
  // Dimension of the intersection
  const int dim = Intersection::mydim;

  // Dimension of world coordinates
  // not needed atm
  // const int coorddim = Intersection::coorddim;

  // Create a set of test points
  const Dune::QuadratureRule<double, dim>& quad = Dune::QuadratureRules<double, dim>::rule(rIIt->type(), 3);

  for (unsigned int l=0; l<quad.size(); l++) {

    Dune::FieldVector<double, dim> quadPos = quad[l].position();

    // Test whether local domain position is consistent with global domain position
    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> localDomainPos =
      rIIt->inside()->geometry().global(rIIt->geometryInInside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> globalDomainPos =
      rIIt->geometry().global(quadPos);

    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> localTargetPos =
      rIIt->outside()->geometry().global(rIIt->geometryInOutside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> globalTargetPos =
      rIIt->geometryOutside().global(quadPos);

    // Test whether local domain position is consistent with global domain position
    assert( (localDomainPos-globalDomainPos).two_norm() < 1e-6 );

    // Test whether local target position is consistent with global target position
    assert( (localTargetPos-globalTargetPos).two_norm() < 1e-6 );

    // Here we assume that the two interface match geometrically:
    if ( (globalDomainPos-globalTargetPos).two_norm() >= 1e-6 )
    {
      std::cout << __FILE__ << ":" << __LINE__ << ": error: assert( (globalDomainPos-globalTargetPos).two_norm() < 1e-6 ) failed\n";
      std::cerr << "localDomainPos  = " << localDomainPos << "\n";
      std::cerr << "globalDomainPos = " << globalDomainPos << "\n";
      std::cerr << "localTargetPos  = " << localTargetPos << "\n";
      std::cerr << "globalTargetPos = " << globalTargetPos << "\n";
    }
    //assert( (globalDomainPos-globalTargetPos).two_norm() < 1e-6 );

  }
}

template<int dim>
struct CellLayout
{
  bool contains (Dune::GeometryType gt)
  {
    return gt.dim()==dim;
  }
};


template <class GlueType>
void testCoupling(const GlueType& glue)
{
  typedef typename GlueType::ctype ctype;

  int dim = GlueType::domdim;
  dim = GlueType::tardim;

  typedef Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::Grid0View, CellLayout > View0Mapper;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::Grid1View, CellLayout > View1Mapper;
  View0Mapper view0mapper(glue.template gridView<0>());
  View1Mapper view1mapper(glue.template gridView<1>());

  std::vector<unsigned int> countInside0(view0mapper.size());
  std::vector<unsigned int> countOutside1(view1mapper.size());
  std::vector<unsigned int> countInside1(view1mapper.size(), 0);
  std::vector<unsigned int> countOutside0(view0mapper.size(), 0);

  // ///////////////////////////////////////
  //   MergedGrid centric Grid0->Grid1
  // ///////////////////////////////////////

  {
    typename GlueType::Grid0IntersectionIterator rIIt    = glue.template ibegin<0>();
    typename GlueType::Grid0IntersectionIterator rIEndIt = glue.template iend<0>();
    for (; rIIt!=rIEndIt; ++rIIt)
    {
      countInside0[view0mapper.map(*rIIt->inside())]++;
      countOutside1[view1mapper.map(*rIIt->outside())]++;
      testIntersection(rIIt);
    }
  }

  // ///////////////////////////////////////
  //   MergedGrid centric Grid1->Grid0
  // ///////////////////////////////////////

  {
    typename GlueType::Grid1IntersectionIterator rIIt    = glue.template ibegin<1>();
    typename GlueType::Grid1IntersectionIterator rIEndIt = glue.template iend<1>();
    for (; rIIt!=rIEndIt; ++rIIt)
    {
      countInside1[view1mapper.map(*rIIt->inside())]++;
      countOutside0[view0mapper.map(*rIIt->outside())]++;
      testIntersection(rIIt);
    }
  }

  // ///////////////////////////////////////
  //   Grid0 Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::Grid0View::template Codim<0>::Iterator DomainIterator;
  DomainIterator dit = glue.template gridView<0>().template begin<0>();
  DomainIterator dend = glue.template gridView<0>().template end<0>();
  for (; dit != dend; ++dit)
  {
    unsigned int icount = 0;

    // intersection iterators
    typename GlueType::Grid0CellIntersectionIterator rIIt    = glue.template ibegin<0>(*dit);
    typename GlueType::Grid0CellIntersectionIterator rIEndIt = glue.template iend<0>(*dit);
    for (; rIIt!=rIEndIt; ++rIIt) {
      // as we only have a single grid, even when testing the
      // parallel extractor, this assertion should be true
      assert (rIIt->inside() == dit);
      testIntersection(rIIt);
      icount++;
    }

    assert(icount == countInside0[view0mapper.map(*dit)]);
    assert(icount == countOutside0[view0mapper.map(*dit)]);
  }

  // ///////////////////////////////////////
  //   Grid1 Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::Grid1View::template Codim<0>::Iterator TargetIterator;
  TargetIterator tit = glue.template gridView<1>().template begin<0>();
  TargetIterator tend = glue.template gridView<1>().template end<0>();
  for (; tit != tend; ++tit)
  {
    unsigned int icount = 0;

    // intersection iterators
    typename GlueType::Grid1CellIntersectionIterator rIIt    = glue.template ibegin<1>(*tit);
    typename GlueType::Grid1CellIntersectionIterator rIEndIt = glue.template iend<1>(*tit);
    for (; rIIt!=rIEndIt; ++rIIt) {
      assert (rIIt->inside() == tit);
      testIntersection(rIIt);
      icount++;
    }

    assert(icount == countInside1[view1mapper.map(*tit)]);
    assert(icount == countOutside1[view1mapper.map(*tit)]);
  }

}

#endif // GRIDGLUE_COUPLINGTEST_HH
