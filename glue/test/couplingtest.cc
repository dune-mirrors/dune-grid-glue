// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <iostream>

#include <dune/common/fvector.hh>
#include <dune/grid/sgrid.hh>
#include <dune/grid/common/quadraturerules.hh>

#include <dune/glue/surfaces/SurfaceDescriptor.hh>
#include <dune/glue/surfaces/GridExtractionTraits.hh>
#include <dune/glue/merging/ContactMappingSurfaceMerge.hh>
#include <dune/glue/adapter/GridGlue.hh>

using namespace Dune;

const int dim=3;

template <class GlueType>
void testCoupling(const GlueType& glue)
{

  typename GlueType::RemoteIntersectionIterator rIIt    = glue.iremotebegin();
  typename GlueType::RemoteIntersectionIterator rIEndIt = glue.iremoteend();

  for (; rIIt!=rIEndIt; ++rIIt) {

    std::cout << "Hallo Welt!" << std::endl;

    // Create a set of test points
    const QuadratureRule<double, dim-1>& quad = QuadratureRules<double, dim-1>::rule(rIIt->type(), 3);

    for (unsigned int l=0; l<quad.size(); l++) {

      Dune::FieldVector<double, dim-1> quadPos = quad[l].position();

      std::cout << rIIt->intersectionDomainGlobal().global(quadPos) << std::endl;
      std::cout << rIIt->intersectionTargetGlobal().global(quadPos) << std::endl;

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

int main(int argc, char *argv[]) try
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

  typedef GridType::LevelGridView DomGridView;
  typedef GridType::LevelGridView TarGridView;

  typedef DefaultExtractionTraits<DomGridView,false,MeshClassification::cube> DomTraits;
  typedef DefaultExtractionTraits<TarGridView,false,MeshClassification::cube> TarTraits;

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
catch (Exception e) {
  std::cout << e << std::endl;
}
