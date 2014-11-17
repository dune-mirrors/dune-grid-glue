// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
//#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/geometrygrid.hh>

#include <dune/grid-glue/extractors/extractorpredicate.hh>
#include <dune/grid-glue/extractors/codim0extractor.hh>
#include <dune/grid-glue/adapter/gridglue.hh>
#include <dune/grid-glue/merging/mixeddimoverlappingmerge.hh>

#include <dune/grid-glue/test/couplingtest.hh>

#include <dune/grid-glue/adapter/gridgluevtkwriterroot.hh>
#include <dune/foamgrid/foamgrid.hh>
#include <dune/grid/io/file/gmshreader.hh>


using namespace Dune;
using namespace Dune::GridGlue;

/** \brief Returns always true */
template <class GridView>
class AllElementsDescriptor
  : public ExtractorPredicate<GridView,0>
{
public:
  virtual bool contains(const typename GridView::Traits::template Codim<0>::EntityPointer& element, unsigned int subentity) const
  {
    return true;
  }
};


/** \brief trafo from dim to dim+1 */
/*template<int dim, int dimw, class ctype>
class MixedDimTrafo
  : public AnalyticalCoordFunction< ctype, dim, dimw, MixedDimTrafo<dim,dimw,ctype> >
{
  //dune_static_assert(dim+1==dimw, "MixedDimTrafo assumes dim+1=dimworld");
  dune_static_assert(dim==1, "MixedDimTrafo currently assumes dim==1");
public:

  //! evaluate method for global mapping
  void evaluate ( const Dune::FieldVector<ctype, dim> &x, Dune::FieldVector<ctype, dimw> &y ) const
  {
    y[0] = 0.5;
    y[1] = 0.5;
    y[2] = x[0];

  }
};
*/

int main(int argc, char** argv)
{
  Dune::MPIHelper::instance(argc, argv);

  MixedDimOverlappingMerge<3,1,3> mixedDimOverlappingMerge;

  // /////////////////////////////////////////////////////////
  //   Make a 3d unit cube grid and a 1d grid embedded in 3d
  // /////////////////////////////////////////////////////////

  typedef SGrid<3,3> GridType3d;

  FieldVector<int, 3> elements(3);
  FieldVector<double,3> lower(0);
  FieldVector<double,3> upper(1);

  GridType3d grid0(elements, lower, upper);

  typedef FoamGrid<1,3> GridType1d;


  std::shared_ptr<FoamGrid<1, 3> > grid1( GmshReader<FoamGrid<1, 3> >::read("/temp/Natalie/DUMUX/dune-foamgrid/doc/grids/gmsh/test1d3d.msh", /*verbose*/ false, false ) );



/*typedef SGrid<1,1> GridType1d;

  FieldVector<int, 1> elements1d(1);
  FieldVector<double,1> lower1d(0.);
  FieldVector<double,1> upper1d(0.75);

  typedef GeometryGrid<GridType1d, MixedDimTrafo<1,3,double> > LiftedGridType;

  GridType1d cubeGrid1_in(elements1d, lower1d, upper1d);

  MixedDimTrafo<1,3,double> trafo;   // transform dim-1 to dim

  LiftedGridType grid1(cubeGrid1_in, trafo);
*/
  // ////////////////////////////////////////
  //   Set up an overlapping coupling
  // ////////////////////////////////////////

  typedef typename GridType3d::LeafGridView DomGridView;
  typedef typename GridType1d::LeafGridView TarGridView;

  typedef Codim0Extractor<DomGridView> DomExtractor;
  typedef Codim0Extractor<TarGridView> TarExtractor;

  AllElementsDescriptor<DomGridView> domdesc;
  AllElementsDescriptor<TarGridView> tardesc;

  DomExtractor domEx(grid0.leafGridView(), domdesc);
  TarExtractor tarEx(grid1->leafGridView(), tardesc);

  typedef ::GridGlue<DomExtractor,TarExtractor> GlueType;

  // The following code is out-commented, because the test functionality
  // doesn't actually work yet.

  MixedDimOverlappingMerge<3,1,3> merger;
  GlueType glue(domEx, tarEx, &merger);

  glue.build();

  // vtk output to test grid-glue - take care to set your path right!!
  GridGlueVtkWriterRoot glueVtk;
  glueVtk.write(glue,"/temp/Natalie/DUMUX/dune-grid-glue/dune/grid-glue/test/mixed3D");

  std::multimap<int,data2test > mymm0;
  std::multimap<int,data2test > mymm1;

  glueVtk.writeMaps(glue, mymm0, mymm1);

  // ///////////////////////////////////////////
  //   Test the coupling
  // ///////////////////////////////////////////

  testCoupling(glue);
}
