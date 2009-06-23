// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/surfaces/surfacedescriptor.hh>
#include <dune/glue/surfaces/gridextractiontraits.hh>
#include <dune/glue/merging/ContactMappingSurfaceMerge.hh>
#include <dune/glue/adapter/GridGlue.hh>

using namespace Dune;

template <class GlueType>
void testCoupling(const GlueType& glue)
{
  dune_static_assert(GlueType::domdim == GlueType::tardim, "For this test domain and target must have the same dimension");

  const int dim = GlueType::domdim;

  typename GlueType::RemoteIntersectionIterator rIIt    = glue.iremotebegin();
  typename GlueType::RemoteIntersectionIterator rIEndIt = glue.iremoteend();

  for (; rIIt!=rIEndIt; ++rIIt) {

    // Create a set of test points
    const QuadratureRule<double, GlueType::domdim-1>& quad = QuadratureRules<double, dim-1>::rule(rIIt->type(), 3);

    for (unsigned int l=0; l<quad.size(); l++) {

      Dune::FieldVector<double, dim-1> quadPos = quad[l].position();

      // Test whether local domain position is consistent with global domain position
      FieldVector<double,dim> localDomainPos
        =  rIIt->entityDomain()->geometry().global(rIIt->intersectionDomainLocal().global(quadPos));
      FieldVector<double,dim> globalDomainPos = rIIt->intersectionDomainGlobal().global(quadPos);
      FieldVector<double,dim> localTargetPos
        = rIIt->entityTarget()->geometry().global(rIIt->intersectionTargetLocal().global(quadPos));
      FieldVector<double,dim> globalTargetPos = rIIt->intersectionTargetGlobal().global(quadPos);

      // Test whether local domain position is consistent with global domain position
      assert( (localDomainPos-globalDomainPos).two_norm() < 1e-6 );

      // Test whether local target position is consistent with global target position
      assert( (localTargetPos-globalTargetPos).two_norm() < 1e-6 );

      // Here we assume that the two interface match geometrically:
      assert( (globalDomainPos-globalTargetPos).two_norm() < 1e-6 );

    }

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
    const Dune::ReferenceElement<double,dim>& refElement = Dune::ReferenceElements<double, dim>::general(eptr->type());

    int numVertices = refElement.size(face, 1, dim);

    for (int i=0; i<numVertices; i++)
      if ( std::abs(eptr->geometry().corner(refElement.subEntity(face,1,i,dim))[0] - sliceCoord_) > 1e-6 )
        return false;

    return true;
  }

private:
  double sliceCoord_;
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

  typedef DefaultExtractionTraits<DomGridView,false,ExtractorClassification> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,false,ExtractorClassification> TarTraits;

  typedef ContactMappingSurfaceMerge<dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits, SurfaceMergeImpl> GlueType;

  SurfaceMergeImpl matcher;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), matcher);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  glue.matcher().setMaxDistance(0.01);

  glue.builder().setDomainFaceDescriptor(domdesc);
  glue.builder().setTargetFaceDescriptor(tardesc);

  if (!glue.builder().build())
    DUNE_THROW(Dune::GridError, "Gluing failed!");

  std::cout << "Gluing successful!" << std::endl;

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);

}

int main(int argc, char *argv[]) try
{

  // Test two unit squares, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::cube>();

  // Test two unit squares, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::hybrid>();

  // Test two unit cubes, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::cube>();

  // Test two unit cubes, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::hybrid>();

}
catch (Exception e) {
  std::cout << e << std::endl;
}
