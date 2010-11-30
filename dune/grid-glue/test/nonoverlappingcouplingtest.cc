// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/common/mpihelper.hh>
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>
#include <dune/grid-glue/extractors/parallelextractor.hh>

#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>

using namespace Dune;


template <class GridView>
class VerticalFaceDescriptor
  : public ExtractorPredicate<GridView,1>
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

/** \brief trafo used for yaspgrids */
template<int dim, typename ctype>
class ShiftTrafo
  : public AnalyticalCoordFunction< ctype, dim, dim, ShiftTrafo<dim,ctype> >
{
  double shift;
public:
  ShiftTrafo(double x) : shift(x) {};

  virtual Dune::FieldVector<double, dim> operator()(const Dune::FieldVector<double, dim>& c) const
  {
    Dune::FieldVector<double, dim> x(c);
    x[0] += shift;
    return x;
  }

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dim> &y ) const
  {
    y = (*this)(x);
  }
};

template <int dim>
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

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

#if HAVE_PSURFACE
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}


template <int dim>
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

  VerticalFaceDescriptor<DomGridView> domdesc(1);
  VerticalFaceDescriptor<TarGridView> tardesc(1);

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

#if HAVE_PSURFACE
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#else
    #warning Not testing, because psurface backend is not available.
#endif
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

    if (tar)
    {
      elements = 4;
    }

    GridType * gridp = new GridType(size, elements, periodic, overlap);
    return *gridp;
  }

  double slice()
  {
    if (tar)
      return 0.0;
    return 1.0;
  }
};


template <int dim, class DomGen, class TarGen>
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
  typedef GeometryGrid<GridType1, ShiftTrafo<dim,double> > ShiftedGridType;
  ShiftTrafo<dim,double> trafo(tarGen.slice());   // transform dim-1 to dim
  ShiftedGridType cubeGrid1(tarGen.generate(), trafo);

  // ////////////////////////////////////////
  //   Set up Traits
  // ////////////////////////////////////////

  typedef typename GridType0::LevelGridView DomGridView;
  typedef typename ShiftedGridType::LevelGridView TarGridView;

  // always test the extractor via the parallel extractor classes, even if we use a sequential grid.

  VerticalFaceDescriptor<DomGridView> domdesc(domGen.slice());
  VerticalFaceDescriptor<TarGridView> tardesc(tarGen.slice());

  typedef ParallelExtractor< Codim1Extractor<DomGridView> > DomExtractor;
  typedef ParallelExtractor< Codim1Extractor<TarGridView> > TarExtractor;

  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////
#if HAVE_PSURFACE
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << merger.nSimplices() << " remote intersections found!" << std::endl;
  assert(merger.nSimplices() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
#else
    #warning Not testing, because psurface backend is not available.
#endif
}


int main(int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  // 2d Tests
  typedef MeshGenerator<2,false>  Seq;
  typedef MeshGenerator<2,true>   Par;

  // Test two unit squares
  std::cout << "==== 2D hybrid =============================================\n";
  testMatchingCubeGrids<2>();
  testNonMatchingCubeGrids<2>();
  // testParallelCubeGrids<2,Seq,Seq>();
  // testParallelCubeGrids<2,Par,Seq>();
  // testParallelCubeGrids<2,Seq,Par>();
  // testParallelCubeGrids<2,Par,Par>();
  std::cout << "============================================================\n";

  // 3d Tests
  typedef MeshGenerator<3,false>  Seq3d;
  typedef MeshGenerator<3,true>   Par3d;

  // Test two unit cubes
  std::cout << "==== 3D hybrid =============================================\n";
  testMatchingCubeGrids<3>();
  testNonMatchingCubeGrids<3>();
  testParallelCubeGrids<3,Seq3d,Seq3d>();
  testParallelCubeGrids<3,Par3d,Seq3d>();
  testParallelCubeGrids<3,Seq3d,Par3d>();
  testParallelCubeGrids<3,Par3d,Par3d>();
  std::cout << "============================================================\n";

}
catch (Exception e) {
  std::cout << e << std::endl;
}
