// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef GRIDGLUE_COUPLINGTEST_HH
#define GRIDGLUE_COUPLINGTEST_HH

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/quadraturerules.hh>

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
      std::cerr << "localDomainPos  = " << localDomainPos << "\n";
      std::cerr << "globalDomainPos = " << globalDomainPos << "\n";
      std::cerr << "localTargetPos  = " << localTargetPos << "\n";
      std::cerr << "globalTargetPos = " << globalTargetPos << "\n";
    }
    assert( (globalDomainPos-globalTargetPos).two_norm() < 1e-6 );

  }
}


template <class GlueType>
void testCoupling(const GlueType& glue)
{
  typedef typename GlueType::ctype ctype;

  int dim = GlueType::domdim;
  dim = GlueType::tardim;

  // ///////////////////////////////////////
  //   MergedGrid centric
  // ///////////////////////////////////////

  typename GlueType::IntersectionIterator rIIt    = glue.ibegin();
  typename GlueType::IntersectionIterator rIEndIt = glue.iend();
  for (; rIIt!=rIEndIt; ++rIIt)
    testIntersection(rIIt);

  // ///////////////////////////////////////
  //   Domain Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::Grid0View::template Codim<0>::Iterator DomainIterator;
  DomainIterator dit = glue.domainGridView().template begin<0>();
  DomainIterator dend = glue.domainGridView().template end<0>();
  for (; dit != dend; ++dit)
  {
    int icount = 0;
    int icount2 = 0;

    // intersection iterators
    typename GlueType::Grid0IntersectionIterator rIIt    = glue.idomainbegin(*dit);
    typename GlueType::Grid0IntersectionIterator rIEndIt = glue.idomainend();
    for (; rIIt!=rIEndIt; ++rIIt) {
      // as we only have a single grid, even when testing the
      // parallel extractor, this assertion should be true
      assert (rIIt->inside() == dit);
      testIntersection(rIIt);
      icount++;
    }

    // face intersection iterators
    const int num_faces = dit->template count<1>();
    for (int f=0; f<num_faces; ++f)
    {
      typename GlueType::Grid0IntersectionIterator rIIt    = glue.idomainbegin(*dit, f);
      typename GlueType::Grid0IntersectionIterator rIEndIt = glue.idomainend();
      for (; rIIt!=rIEndIt; ++rIIt) {
        assert (rIIt->inside() == dit);
        testIntersection(rIIt);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

  // ///////////////////////////////////////
  //   Target Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::Grid1View::template Codim<0>::Iterator TargetIterator;
  TargetIterator tit = glue.targetGridView().template begin<0>();
  TargetIterator tend = glue.targetGridView().template end<0>();
  for (; tit != tend; ++tit)
  {
    int icount = 0;
    int icount2 = 0;

    // intersection iterators
    typename GlueType::Grid1IntersectionIterator rIIt    = glue.itargetbegin(*tit);
    typename GlueType::Grid1IntersectionIterator rIEndIt = glue.itargetend();
    for (; rIIt!=rIEndIt; ++rIIt) {
      assert (rIIt->inside() == tit);
      testIntersection(rIIt);
      icount++;
    }

    // face intersection iterators
    const int num_faces = tit->template count<1>();
    for (int f=0; f<num_faces; ++f)
    {
      typename GlueType::Grid1IntersectionIterator rIIt    = glue.itargetbegin(*tit, f);
      typename GlueType::Grid1IntersectionIterator rIEndIt = glue.itargetend();
      for (; rIIt!=rIEndIt; ++rIIt) {
        assert (rIIt->inside() == tit);
        testIntersection(rIIt);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

}

#endif // GRIDGLUE_COUPLINGTEST_HH
