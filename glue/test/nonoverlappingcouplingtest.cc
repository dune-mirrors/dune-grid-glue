// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/extractors/surfacedescriptor.hh>
#include <dune/glue/extractors/gridextractiontraits.hh>
#include <dune/glue/merging/psurfacemerge.hh>
#include <dune/glue/adapter/gridglue.hh>

#include <dune/glue/test/couplingtest.hh>

using namespace Dune;


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

  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  merger.setMaxDistance(0.01);

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


template <int dim, MeshClassification::MeshType ExtractorClassification, bool domParallel = false, bool tarParallel = false>
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

  typedef DefaultExtractionTraits<DomGridView,1,ExtractorClassification,domParallel> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification,tarParallel> TarTraits;

  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  merger.setMaxDistance(0.01);

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


int main(int argc, char *argv[]) try
{

  // Test two unit squares, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::cube>();
  testNonMatchingCubeGrids<2,MeshClassification::cube>();
  testNonMatchingCubeGrids<2,MeshClassification::cube,true,true>();

  // Test two unit squares, extract boundaries using the SimplexSurfaceExtractor
  // Should work, because the boundary consists of 1d simplices
  testMatchingCubeGrids<2,MeshClassification::simplex>();
  testNonMatchingCubeGrids<2,MeshClassification::simplex>();
  testNonMatchingCubeGrids<2,MeshClassification::simplex,true,true>();

  // Test two unit squares, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<2,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<2,MeshClassification::hybrid,true,true>();

  // Test two unit cubes, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::cube>();
  testNonMatchingCubeGrids<3,MeshClassification::cube>();
  testNonMatchingCubeGrids<3,MeshClassification::cube,true,true>();

  // Test two unit cubes, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<3,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<3,MeshClassification::hybrid,true,true>();

}
catch (Exception e) {
  std::cout << e << std::endl;
}
