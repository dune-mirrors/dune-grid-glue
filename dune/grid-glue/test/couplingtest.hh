// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef GRIDGLUE_COUPLINGTEST_HH
#define GRIDGLUE_COUPLINGTEST_HH

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

template <class IntersectionIt, int domdimw, int tardimw, class ctype>
void testIntersection(const IntersectionIt & rIIt,
                      CoordinateTransformation<domdimw, IntersectionIt::value_type::coorddim, ctype> & domTrafo,
                      CoordinateTransformation<tardimw, IntersectionIt::value_type::coorddim, ctype> & tarTrafo)
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
    Dune::FieldVector<double, Intersection::DomainGridView::dimensionworld> localDomainPos =
      rIIt->entityDomain()->geometry().global(rIIt->intersectionDomainLocal().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::DomainGridView::dimensionworld> globalDomainPos =
      rIIt->intersectionDomainGlobal().global(quadPos);

    Dune::FieldVector<double, Intersection::TargetGridView::dimensionworld> localTargetPos =
      rIIt->entityTarget()->geometry().global(rIIt->intersectionTargetLocal().global(quadPos));

    // currently the intersection maps to the GV::dimworld, this will hopefully change soon
    Dune::FieldVector<double, Intersection::TargetGridView::dimensionworld> globalTargetPos =
      rIIt->intersectionTargetGlobal().global(quadPos);

    // Test whether local domain position is consistent with global domain position
    assert( (localDomainPos-globalDomainPos).two_norm() < 1e-6 );

    // Test whether local target position is consistent with global target position
    assert( (localTargetPos-globalTargetPos).two_norm() < 1e-6 );

    // Here we assume that the two interface match geometrically:
    if ( (domTrafo(globalDomainPos)-tarTrafo(globalTargetPos)).two_norm() >= 1e-6 )
    {
      std::cerr << "globalDomainPos = " << globalDomainPos << "\n";
      std::cerr << "globalTargetPos = " << globalTargetPos << "\n";
      std::cerr << "domTrafo(globalDomainPos) = " << domTrafo(globalDomainPos) << "\n";
      std::cerr << "tarTrafo(globalTargetPos) = " << tarTrafo(globalTargetPos) << "\n";
    }
    assert( (domTrafo(globalDomainPos)-tarTrafo(globalTargetPos)).two_norm() < 1e-6 );

  }
}


template<int dim, int dimw, typename ctype>
class NoTransformation : public CoordinateTransformation<dim, dimw, ctype>
{
public:
  static
  void assignIfNull(CoordinateTransformation<dim, dimw, ctype>* & ptr)
  {
    assert(ptr != 0);
  }

  /** \brief Map a point to a new position */
  // this works only for dim == dimw
  virtual Dune::FieldVector<ctype, dim> operator()(const Dune::FieldVector<ctype, dimw>& c) const
  {
    assert(false);
  }
};

template<int dim, typename ctype>
class NoTransformation<dim, dim, ctype> : public CoordinateTransformation<dim, dim, ctype>
{
public:
  static
  void assignIfNull(CoordinateTransformation<dim, dim, ctype>* & ptr)
  {
    if (ptr == 0)
      ptr = new NoTransformation<dim, dim, ctype>;
  }

  /** \brief Map a point to a new position */
  // this works only for dim == dimw
  virtual Dune::FieldVector<ctype, dim> operator()(const Dune::FieldVector<ctype, dim>& c) const
  {
    return c;
  }
};

template <class GlueType>
void testCoupling(const GlueType& glue,
                  CoordinateTransformation<GlueType::Grid0Patch::dimworld, GlueType::dimworld, typename GlueType::ctype> * domTrafo = 0,
                  CoordinateTransformation<GlueType::Grid1Patch::dimworld, GlueType::dimworld, typename GlueType::ctype> * tarTrafo = 0 )
{
  typedef typename GlueType::ctype ctype;

  int dim = GlueType::domdim;
  dim = GlueType::tardim;

  // ///////////////////////////////////////
  //   set Identity trafo if necessary
  // ///////////////////////////////////////
  NoTransformation<GlueType::Grid0Patch::dimworld, GlueType::dimworld, ctype>::assignIfNull(domTrafo);
  NoTransformation<GlueType::Grid1Patch::dimworld, GlueType::dimworld, ctype>::assignIfNull(tarTrafo);

  // ///////////////////////////////////////
  //   MergedGrid centric
  // ///////////////////////////////////////

  typename GlueType::RemoteIntersectionIterator rIIt    = glue.iremotebegin();
  typename GlueType::RemoteIntersectionIterator rIEndIt = glue.iremoteend();
  for (; rIIt!=rIEndIt; ++rIIt)
    testIntersection(rIIt, *domTrafo, *tarTrafo);

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
      // as we only have a single grid, even when testing the
      // parallel extractor, this assertion should be true
      assert (rIIt->entityDomain() == dit);
      testIntersection(rIIt, *domTrafo, *tarTrafo);
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
        testIntersection(rIIt, *domTrafo, *tarTrafo);
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
      testIntersection(rIIt, *domTrafo, *tarTrafo);
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
        testIntersection(rIIt, *domTrafo, *tarTrafo);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

}

#endif // GRIDGLUE_COUPLINGTEST_HH
