// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef GRIDGLUE_COUPLINGTEST_HH
#define GRIDGLUE_COUPLINGTEST_HH

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/common/mcmgmapper.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/gridglue.hh>

template <class IntersectionIt>
void testIntersection(const IntersectionIt & rIIt)
{
  typedef typename IntersectionIt::value_type Intersection;
  // Dimension of the intersection
  const int dim = Intersection::mydim;

  // Dimension of world coordinates
  const int coorddim = Intersection::coorddim;

  // Create a set of test points
  const Dune::QuadratureRule<double, dim>& quad = Dune::QuadratureRules<double, dim>::rule(rIIt->type(), 3);

  for (unsigned int l=0; l<quad.size(); l++) {

    Dune::FieldVector<double, dim> quadPos = quad[l].position();

    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> localGrid0Pos =
      rIIt->inside()->geometry().global(rIIt->geometryInInside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::InsideGridView::dimensionworld> globalGrid0Pos =
      rIIt->geometry().global(quadPos);

    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> localGrid1Pos =
      rIIt->outside()->geometry().global(rIIt->geometryInOutside().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::OutsideGridView::dimensionworld> globalGrid1Pos =
      rIIt->geometryOutside().global(quadPos);

    // Test whether local grid0 position is consistent with global grid0 position
    assert( (localGrid0Pos-globalGrid0Pos).two_norm() < 1e-6 );

    // Test whether local grid1 position is consistent with global grid1 position
    assert( (localGrid1Pos-globalGrid1Pos).two_norm() < 1e-6 );

    // Here we assume that the two interfaces match geometrically:
    if ( (globalGrid0Pos-globalGrid1Pos).two_norm() >= 1e-4 )
    {
      std::cout << __FILE__ << ":" << __LINE__ << ": error: assert( (globalGrid0Pos-globalGrid1Pos).two_norm() < 1e-4 ) failed\n";
      std::cerr << "localGrid0Pos  = " << localGrid0Pos << "\n";
      std::cerr << "globalGrid0Pos = " << globalGrid0Pos << "\n";
      std::cerr << "localGrid1Pos  = " << localGrid1Pos << "\n";
      std::cerr << "globalGrid1Pos = " << globalGrid1Pos << "\n";
    }
    //assert( (globalGrid0Pos-globalGrid1Pos).two_norm() < 1e-6 );

    // Test the normal vector methods.  At least test whether they don't crash
    if (coorddim - dim == 1) // only test for codim 1
    {
      rIIt->outerNormal(quadPos);
      rIIt->unitOuterNormal(quadPos);
      rIIt->integrationOuterNormal(quadPos);
      rIIt->centerUnitOuterNormal();
    }
  }
}


template <class GlueType>
void testCoupling(const GlueType& glue)
{
  typedef Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::Grid0View, Dune::MCMGElementLayout > View0Mapper;
  typedef Dune::MultipleCodimMultipleGeomTypeMapper< typename GlueType::Grid1View, Dune::MCMGElementLayout > View1Mapper;
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
      if (rIIt->self() && rIIt->neighbor())
      {
        countInside0[view0mapper.map(*rIIt->inside())]++;
        countOutside1[view1mapper.map(*rIIt->outside())]++;
        testIntersection(rIIt);
      }
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
      if (rIIt->self() && rIIt->neighbor())
      {
        countInside1[view1mapper.map(*rIIt->inside())]++;
        countOutside0[view0mapper.map(*rIIt->outside())]++;
        testIntersection(rIIt);
      }
    }
  }
}

#endif // GRIDGLUE_COUPLINGTEST_HH
