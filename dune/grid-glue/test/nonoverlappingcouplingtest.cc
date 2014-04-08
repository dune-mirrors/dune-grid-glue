// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON,2,3)
#include <dune/common/parallel/mpihelper.hh>
#else
#include <dune/common/mpihelper.hh>
#endif
#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/geometrygrid.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim1extractor.hh>

#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/gridglue.hh>

#include <dune/grid-glue/test/couplingtest.hh>
#include <dune/grid-glue/test/communicationtest.hh>

using namespace Dune;
using namespace Dune::GridGlue;

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
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    const Dune::ReferenceElement<double,dim>& refElement = Dune::ReferenceElements<double, dim>::general(eptr->type());
#else
    const Dune::GenericReferenceElement<double,dim>& refElement = Dune::GenericReferenceElements<double, dim>::general(eptr->type());
#endif

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

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dim> &y ) const
  {
    y = x;
    y[0] += shift;
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

#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
  DomExtractor domEx(cubeGrid0.levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
#endif

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  // Testing with PSurfaceMerge
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
  testCommunication(glue);
#endif

  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  ContactMergeImpl contactMerger(0.01);
  GlueType contactGlue(domEx, tarEx, &contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
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

#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
  DomExtractor domEx(cubeGrid0.levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0.levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1.levelView(0), tardesc);
#endif

#if HAVE_PSURFACE
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
  testCommunication(glue);
#endif

  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  ContactMergeImpl contactMerger(0.01);
  GlueType contactGlue(domEx, tarEx, &contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
}


template<int dim, bool par>
class MeshGenerator
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  typedef SGrid<dim,dim> GridType;

  std::shared_ptr<GridType> generate()
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

    return std::make_shared<GridType>(elements, lower, upper);
  }
};


template<int dim>
class MeshGenerator<dim, true>
{
  bool tar;
public:
  MeshGenerator(bool b) : tar(b) {};

  typedef YaspGrid<dim> HostGridType;
  typedef GeometryGrid<HostGridType, ShiftTrafo<dim,double> > GridType;

  std::shared_ptr<GridType> generate()
  {
#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
    Dune::array<int,dim> elements;
    std::fill(elements.begin(), elements.end(), 2);
    std::bitset<dim> periodic(0);
#else
    FieldVector<int, dim> elements(2);
    FieldVector<bool,dim> periodic(false);
#endif
    FieldVector<double,dim> size(1);
    int overlap = 1;
    double shift = 0.0;

    if (tar)
    {
      std::fill(elements.begin(), elements.end(), 4);
      shift = 1.0;
    }

    HostGridType * hostgridp = new HostGridType(
#if HAVE_MPI
      MPI_COMM_WORLD,
#endif // HAVE_MPI
      size, elements, periodic, overlap);
    ShiftTrafo<dim,double> * trafop = new ShiftTrafo<dim,double>(shift);
    return std::make_shared<GridType>(*hostgridp, *trafop);
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

  double slice = 1.0;

  std::shared_ptr<GridType0> cubeGrid0 = domGen.generate();
  std::shared_ptr<GridType1> cubeGrid1 = tarGen.generate();

  // ////////////////////////////////////////
  //   Set up Traits
  // ////////////////////////////////////////

  typedef typename GridType0::LevelGridView DomGridView;
  typedef typename GridType1::LevelGridView TarGridView;

  VerticalFaceDescriptor<DomGridView> domdesc(slice);
  VerticalFaceDescriptor<TarGridView> tardesc(slice);

  typedef Codim1Extractor<DomGridView> DomExtractor;
  typedef Codim1Extractor<TarGridView> TarExtractor;

#if DUNE_VERSION_NEWER(DUNE_GRID,2,3)
  DomExtractor domEx(cubeGrid0->levelGridView(0), domdesc);
  TarExtractor tarEx(cubeGrid1->levelGridView(0), tardesc);
#else
  DomExtractor domEx(cubeGrid0->levelView(0), domdesc);
  TarExtractor tarEx(cubeGrid1->levelView(0), tardesc);
#endif

  // ////////////////////////////////////////
  //   Set up coupling at their interface
  // ////////////////////////////////////////

  typedef Dune::GridGlue::GridGlue<DomExtractor,TarExtractor> GlueType;

#if HAVE_PSURFACE
  // Test using PSurfaceMerge
  typedef PSurfaceMerge<dim-1,dim,double> SurfaceMergeImpl;

  SurfaceMergeImpl merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  std::cout << "Gluing successful, " << glue.size() << " remote intersections found!" << std::endl;
  assert(glue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
  testCommunication(glue);
#endif
  // Testing with ContactMerge
  typedef ContactMerge<dim,double> ContactMergeImpl;

  ContactMergeImpl contactMerger(0.01);
  GlueType contactGlue(domEx, tarEx, &contactMerger);

  contactGlue.build();

  std::cout << "Gluing successful, " << contactGlue.size() << " remote intersections found!" << std::endl;
  assert(contactGlue.size() > 0);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(contactGlue);
  testCommunication(contactGlue);
}

#if HAVE_MPI
void eh( MPI_Comm *comm, int *err, ... )
{
  int len = 1024;
  char error_txt[len];

  MPI_Error_string(*err, error_txt, &len);
  assert(len <= 1024);
  DUNE_THROW(Dune::Exception, "MPI ERROR -- " << error_txt);
}
#endif // HAVE_MPI

int main(int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);
  Dune::dinfo.attach(std::cout);

#if HAVE_MPI
  MPI_Errhandler errhandler;
  MPI_Comm_create_errhandler(eh, &errhandler);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
#endif

  // 2d Tests
  typedef MeshGenerator<2,false>  Seq;
  typedef MeshGenerator<2,true>   Par;

  // Test two unit squares
  std::cout << "==== 2D hybrid =============================================\n";
  testMatchingCubeGrids<2>();
  std::cout << "============================================================\n";
  testNonMatchingCubeGrids<2>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Seq,Seq>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Par,Seq>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Seq,Par>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<2,Par,Par>();
  std::cout << "============================================================\n";

  // 3d Tests
  typedef MeshGenerator<3,false>  Seq3d;
  typedef MeshGenerator<3,true>   Par3d;

  // Test two unit cubes
  std::cout << "==== 3D hybrid =============================================\n";
#if ! HAVE_MPI
  testMatchingCubeGrids<3>();
  std::cout << "============================================================\n";
  testNonMatchingCubeGrids<3>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Seq3d,Seq3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Par3d,Seq3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Seq3d,Par3d>();
  std::cout << "============================================================\n";
  testParallelCubeGrids<3,Par3d,Par3d>();
  std::cout << "============================================================\n";
#endif // HAVE_MPI

  return 0;
}
catch (Exception e) {
  int i = 0; char** c = 0;
  std::cout << Dune::MPIHelper::instance(i,c).rank() << ": " << e << std::endl;
  return 1;
}
