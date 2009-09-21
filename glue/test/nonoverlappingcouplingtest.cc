// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/common/mpihelper.hh>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
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


/** \brief trafo used for yaspgrids */
template<int dim>
class ShiftTrafo : public CoordinateTransformation<dim,dim,double>
{
  double shift;

public:
  ShiftTrafo(double x) : shift(x) {};

  virtual Dune::FieldVector<double, 2> operator()(const Dune::FieldVector<double, 2>& c) const
  {
    Dune::FieldVector<double, 2> x(c);
    x[0] += shift;
    return x;
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


template<int dim, bool par>
class MeshGenerator
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  typedef SGrid<dim,dim> GridType;

  GridType & generate()
  {
    FieldVector<int, dim> elements(2);
    FieldVector<double,dim> lower(0);
    FieldVector<double,dim> upper(1);

    if (tar)
    {
      elements = 4;
      lower[0] += 1;
      upper[0] += 1;
    }

    GridType * gridp = new GridType(elements, lower, upper);
    return *gridp;
  }

  double slice()
  {
    return 1.0;
  }

  const CoordinateTransformation<dim, dim, double> * trafo()
  {
    return 0;
  }
};


template<int dim>
class MeshGenerator<dim, true>
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  typedef YaspGrid<dim> GridType;

  GridType & generate()
  {
    FieldVector<int, dim> elements(2);
    FieldVector<double,dim> size(1);
    FieldVector<bool,dim> periodic(false);
    int overlap = 1;

    GridType * gridp = new GridType(size, elements, periodic, overlap);
    return *gridp;
  }

  double slice()
  {
    if (tar)
      return 0.0;
    return 1.0;
  }

  const CoordinateTransformation<dim, dim, double> * trafo()
  {
    double shift = 0.0;
    if (tar)
      shift = 1.0;
    return new ShiftTrafo<dim>(shift);
  }
};


template <int dim, class DomGen, class TarGen, MeshClassification::MeshType ExtractorClassification>
void testParallelCubeGrids()
{
  // ///////////////////////////////////////
  //   Make two cube grids
  // ///////////////////////////////////////

  typedef typename DomGen::GridType GridType0;
  typedef typename TarGen::GridType GridType1;

  DomGen domGen(0);
  TarGen tarGen(1);

  GridType0 & cubeGrid0 = domGen.generate();
  GridType1 & cubeGrid1 = tarGen.generate();

  // ////////////////////////////////////////
  //   Set up Traits
  // ////////////////////////////////////////

  typedef typename GridType0::LevelGridView DomGridView;
  typedef typename GridType1::LevelGridView TarGridView;

  // always test the extractor via the parallel extractor classes, even if we use a seqqential grid.
  typedef DefaultExtractionTraits<DomGridView,1,ExtractorClassification,true> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,1,ExtractorClassification,true> TarTraits;

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef GridGlue<DomTraits,TarTraits> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(cubeGrid0.levelView(0), cubeGrid1.levelView(0), &merger);

  VerticalFaceDescriptor<DomGridView> domdesc(domGen.slice());
  VerticalFaceDescriptor<TarGridView> tardesc(tarGen.slice());

  merger.setMaxDistance(0.01);

  glue.builder().setDomainTransformation(domGen.trafo());
  glue.builder().setTargetTransformation(tarGen.trafo());
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
  Dune::MPIHelper::instance(argc, argv);

  // 2d Tests
  typedef MeshGenerator<2,false>  Seq;
  typedef MeshGenerator<2,true>   Par;

  // Test two unit squares, extract boundaries using the CubeSurrfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::cube>();
  testNonMatchingCubeGrids<2,MeshClassification::cube>();
  testParallelCubeGrids<2,Seq,Seq,MeshClassification::cube>();
  testParallelCubeGrids<2,Par,Seq,MeshClassification::cube>();
  //testParallelCubeGrids<2,Seq,Par,MeshClassification::cube>();
  //testParallelCubeGrids<2,Par,Par,MeshClassification::cube>();

  // Test two unit squares, extract boundaries using the SimplexSurfaceExtractor
  // Should work, because the boundary consists of 1d simplices
  testMatchingCubeGrids<2,MeshClassification::simplex>();
  testNonMatchingCubeGrids<2,MeshClassification::simplex>();
  testParallelCubeGrids<2,Seq,Seq,MeshClassification::simplex>();
  testParallelCubeGrids<2,Par,Seq,MeshClassification::simplex>();

  // Test two unit squares, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<2,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<2,MeshClassification::hybrid>();
  testParallelCubeGrids<2,Seq,Seq,MeshClassification::hybrid>();
  testParallelCubeGrids<2,Par,Seq,MeshClassification::hybrid>();

  // 3d Tests
  typedef MeshGenerator<3,false>  Seq3d;
  typedef MeshGenerator<3,true>   Par3d;

  // Test two unit cubes, extract boundaries using the CubeSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::cube>();
  testNonMatchingCubeGrids<3,MeshClassification::cube>();
  // testParallelCubeGrids<3,Seq3d,Seq3d,MeshClassification::cube>();

  // Test two unit cubes, extract boundaries using the GeneralSurfaceExtractor
  testMatchingCubeGrids<3,MeshClassification::hybrid>();
  testNonMatchingCubeGrids<3,MeshClassification::hybrid>();
  // testParallelCubeGrids<3,Seq3d,Seq3d,MeshClassification::hybrid>();

}
catch (Exception e) {
  std::cout << e << std::endl;
}
