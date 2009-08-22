// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/surfaces/surfacedescriptor.hh>
#include <dune/glue/surfaces/gridextractiontraits.hh>
#include <dune/glue/merging/psurfacemerge.hh>
#include <dune/glue/adapter/gridglue.hh>

using namespace Dune;

template <int dim, class IntersectionIt>
void testIntersection(const IntersectionIt & rIIt)
{
  const QuadratureRule<double, dim-1>& quad = QuadratureRules<double, dim-1>::rule(rIIt->type(), 3);

  for (unsigned int l=0; l<quad.size(); l++) {

    Dune::FieldVector<double, dim-1> quadPos = quad[l].position();

    // Test whether local domain position is consistent with global domain position
    FieldVector<double,dim> localDomainPos =
      rIIt->entityDomain()->geometry().global(rIIt->intersectionDomainLocal().global(quadPos));
    FieldVector<double,dim> globalDomainPos = rIIt->intersectionDomainGlobal().global(quadPos);
    FieldVector<double,dim> localTargetPos =
      rIIt->entityTarget()->geometry().global(rIIt->intersectionTargetLocal().global(quadPos));
    FieldVector<double,dim> globalTargetPos = rIIt->intersectionTargetGlobal().global(quadPos);

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
  for (; rIIt!=rIEndIt; ++rIIt) {
    // Create a set of test points
    testIntersection<dim, typename GlueType::RemoteIntersectionIterator>(rIIt);
  }

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
      testIntersection<dim>(rIIt);
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
        testIntersection<dim>(rIIt);
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
      testIntersection<dim>(rIIt);
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
        testIntersection<dim>(rIIt);
        icount2++;
      }
    }

    assert(icount == icount2);
  }

}

template <class GridView>
class VerticalFaceDescriptor
  : public FaceDescriptor<GridView>
{
public:
  VerticalFaceDescriptor(double sliceCoord)
    : sliceCoord_(sliceCoord)
  {}

  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& eptr,
                        unsigned int face) const
  {
    const int dim = GridView::dimension;
    const Dune::GenericReferenceElement<double,dim>& refElement = Dune::GenericReferenceElements<double, dim>::general(eptr->type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(eptr->geometry().corner(refElement.subEntity(face,1,i,dim))[0] - sliceCoord_) > 1e-6 )
        return false;

    return true;
  }

private:
  double sliceCoord_;
};

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

template <int dim, MeshClassification::MeshType ExtractorClassification>
void testMatchingCubeGrids()
{

  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType cubeGrid0(elements, lower, upper);

  lower[0] += 1;
  upper[0] += 1;

  GridType cubeGrid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType::LevelGridView DomGridView;
  typedef typename GridType::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,1,ExtractorClassification> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification> TarTraits;

  typedef PSurfaceMerge<dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  merger.setMaxDistance(0.01);

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification>
void testNonMatchingCubeGrids()
{

  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType;

  FieldVector<int, dim> elements(2);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType cubeGrid0(elements, lower, upper);

  elements = 4;
  lower[0] += 1;
  upper[0] += 1;

  GridType cubeGrid1(elements, lower, upper);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType::LevelGridView DomGridView;
  typedef typename GridType::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,1,ExtractorClassification> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification> TarTraits;

  typedef PSurfaceMerge<dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  merger.setMaxDistance(0.01);

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}


template <int dim, MeshClassification::MeshType ExtractorClassification>
void test1d2dCoupling()
{

  // ///////////////////////////////////////
  //   Make a cube grid and a 1d grid
  // ///////////////////////////////////////

  typedef SGrid<dim,dim> GridType2d;

  FieldVector<int, dim> elements(1);
  FieldVector<double,dim> lower(0);
  FieldVector<double,dim> upper(1);

  GridType2d cubeGrid0(elements, lower, upper);

  typedef SGrid<dim-1,dim> GridType1d;

  FieldVector<int, dim-1> elements1d(1);
  FieldVector<double,dim-1> lower1d(0);
  FieldVector<double,dim-1> upper1d(1);

  GridType1d cubeGrid1(elements1d, lower1d, upper1d);


  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef typename GridType2d::LevelGridView DomGridView;
  typedef typename GridType1d::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,1,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,0,ExtractorClassification> TarTraits;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  PSurfaceMerge<dim,double> merger;
  merger.setMaxDistance(0.01);

  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  AllElementsDescriptor<TarGridView>  tardesc;

  glue.builder().setDomainDescriptor(domdesc);
  glue.builder().setTargetDescriptor(tardesc);

  glue.builder().build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}




int main(int argc, char *argv[]) try
{

  // Test two unit squares, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::cube>();
  testNonMatchingCubeGrids<2,MeshClassification::cube>();

  // Test two unit squares, extract boundaries using the SimplexSurfaceExtractor
  // Should work, because the boundary consists of 1d simplices
  testMatchingCubeGrids<2,MeshClassification::simplex>();
  testNonMatchingCubeGrids<2,MeshClassification::simplex>();

  // Test two unit squares, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<2,MeshClassification::hybrid>();

  // Test two unit cubes, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::cube>();
  testNonMatchingCubeGrids<3,MeshClassification::cube>();

  // Test two unit cubes, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<3,MeshClassification::hybrid>();

  // Test a unit cube versus a grid one dimension lower
  test1d2dCoupling<2,MeshClassification::cube>();
  test1d2dCoupling<2,MeshClassification::simplex>();

}
catch (Exception e) {
  std::cout << e << std::endl;
}
