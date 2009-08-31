// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/extractors/surfacedescriptor.hh>
#include <dune/glue/extractors/gridextractiontraits.hh>
#include <dune/glue/merging/psurfacemerge.hh>
#include <dune/glue/adapter/gridglue.hh>

template <class IntersectionIt>
void testIntersection(const IntersectionIt & rIIt)
{
  // Dimension of the intersection
  const int dim = IntersectionIt::mydim;

  // Dimension of world coordinates
  const int coorddim = IntersectionIt::coorddim;

  const Dune::QuadratureRule<double, dim>& quad = Dune::QuadratureRules<double, dim>::rule(rIIt->type(), 3);

  for (unsigned int l=0; l<quad.size(); l++) {

    Dune::FieldVector<double, dim> quadPos = quad[l].position();

    // Test whether local domain position is consistent with global domain position
    Dune::FieldVector<double,coorddim> localDomainPos =
      rIIt->entityDomain()->geometry().global(rIIt->intersectionDomainLocal().global(quadPos));
    Dune::FieldVector<double,coorddim> globalDomainPos = rIIt->intersectionDomainGlobal().global(quadPos);
    Dune::FieldVector<double,coorddim> localTargetPos =
      rIIt->entityTarget()->geometry().global(rIIt->intersectionTargetLocal().global(quadPos));
    Dune::FieldVector<double,coorddim> globalTargetPos = rIIt->intersectionTargetGlobal().global(quadPos);

    // Test whether local domain position is consistent with global domain position
    assert( (localDomainPos-globalDomainPos).two_norm() < 1e-6 );

    // Test whether local target position is consistent with global target position
    assert( (localTargetPos-globalTargetPos).two_norm() < 1e-6 );

    // Here we assume that the two interface match geometrically:
    assert( (globalDomainPos-globalTargetPos).two_norm() < 1e-6 );
  }
}

template <class GlueType>
void testCoupling(const GlueType& glue)
{
  //dune_static_assert(GlueType::domdim == GlueType::tardim, "For this test domain and target must have the same dimension");

  const int dim = GlueType::domdim;

  // ///////////////////////////////////////
  //   MergedGrid centric
  // ///////////////////////////////////////

  typename GlueType::RemoteIntersectionIterator rIIt    = glue.iremotebegin();
  typename GlueType::RemoteIntersectionIterator rIEndIt = glue.iremoteend();
  for (; rIIt!=rIEndIt; ++rIIt)
    testIntersection(rIIt);

  // ///////////////////////////////////////
  //   Domain Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::DomainGridView::template Codim<0>::Iterator DomainIterator;
  DomainIterator dit = glue.domainGridView().template begin<0>();
  DomainIterator dend = glue.domainGridView().template end<0>();
  for (; dit != dend; ++dit)
  {
    int icount = 0;
    int icount2 = 0;

    // intersection iterators
    typename GlueType::DomainIntersectionIterator rIIt    = glue.idomainbegin(*dit);
    typename GlueType::DomainIntersectionIterator rIEndIt = glue.idomainend();
    for (; rIIt!=rIEndIt; ++rIIt) {
      assert (rIIt->entityDomain() == dit);
      testIntersection(rIIt);
      icount++;
    }

    // face intersection iterators
    const int num_faces = dit->template count<1>();
    for (int f=0; f<num_faces; ++f)
    {
      typename GlueType::DomainIntersectionIterator rIIt    = glue.idomainbegin(*dit, f);
      typename GlueType::DomainIntersectionIterator rIEndIt = glue.idomainend();
      for (; rIIt!=rIEndIt; ++rIIt) {
        assert (rIIt->entityDomain() == dit);
        testIntersection(rIIt);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

  // ///////////////////////////////////////
  //   Target Entity centric
  // ///////////////////////////////////////

  typedef typename GlueType::TargetGridView::template Codim<0>::Iterator TargetIterator;
  TargetIterator tit = glue.targetGridView().template begin<0>();
  TargetIterator tend = glue.targetGridView().template end<0>();
  for (; tit != tend; ++tit)
  {
    int icount = 0;
    int icount2 = 0;

    // intersection iterators
    typename GlueType::TargetIntersectionIterator rIIt    = glue.itargetbegin(*tit);
    typename GlueType::TargetIntersectionIterator rIEndIt = glue.itargetend();
    for (; rIIt!=rIEndIt; ++rIIt) {
      assert (rIIt->entityTarget() == tit);
      testIntersection(rIIt);
      icount++;
    }

    // face intersection iterators
    const int num_faces = tit->template count<1>();
    for (int f=0; f<num_faces; ++f)
    {
      typename GlueType::TargetIntersectionIterator rIIt    = glue.itargetbegin(*tit, f);
      typename GlueType::TargetIntersectionIterator rIEndIt = glue.itargetend();
      for (; rIIt!=rIEndIt; ++rIIt) {
        assert (rIIt->entityTarget() == tit);
        testIntersection(rIIt);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

}
