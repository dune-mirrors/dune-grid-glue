// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/grid/uggrid.hh>

#include <dune/grid-glue/common/orientedsubface.hh>

using namespace Dune;

template <class GridType>
void test(const GridType& grid) {

  const int dim = GridType::dimension;

  FieldVector<double,dim-1> localLineCorners[2];
  if (dim==2) {
    localLineCorners[0][0] = 0;
    localLineCorners[1][0] = 1;
  }

  FieldVector<double,dim-1> localTriCorners[3];
  if (dim==3) {
    localTriCorners[0][0] = 0;   localTriCorners[0][1] = 0;
    localTriCorners[1][0] = 1;   localTriCorners[1][1] = 0;
    localTriCorners[2][0] = 0;   localTriCorners[2][1] = 1;
  }

  FieldVector<double,dim-1> localQuadCorners[4];
  if (dim==3) {
    localQuadCorners[0][0] = 0;   localQuadCorners[0][1] = 0;
    localQuadCorners[1][0] = 1;   localQuadCorners[1][1] = 0;
    localQuadCorners[2][0] = 0;   localQuadCorners[2][1] = 1;
    localQuadCorners[3][0] = 1;   localQuadCorners[3][1] = 1;
  }

  double eps = 1e-6;

  typename GridType::template Codim<0>::EntityPointer element(grid.template leafbegin<0>());

  typename GridType::template Codim<0>::Entity::LevelIntersectionIterator nIt = element->ilevelbegin();
  typename GridType::template Codim<0>::Entity::LevelIntersectionIterator nEndIt = element->ilevelend();

  std::cout << "Testing: " << element->type() << std::endl;

  for (; nIt!=nEndIt; ++nIt) {

    printf("Testing subface %d\n", nIt->indexInInside());

    switch (nIt->geometry().corners()) {

    case 2 :
      for (int i=0; i<2; i++) {
        int v = orientedSubface<dim>(element->type(), nIt->indexInInside(), i);
        //std::cout << element->geometry()[v] << "  "
        //          << nIt.intersectionGlobal().global(localTriCorners[i]) << std::endl;
        assert( (element->geometry().corner(v) - nIt->geometry().global(localLineCorners[i]) ).two_norm() < eps);
      }
      break;

    case 3 :
      for (int i=0; i<3; i++) {
        int v = orientedSubface<dim>(element->type(), nIt->indexInInside(), i);
        //std::cout << element->geometry()[v] << "  "
        //          << nIt.intersectionGlobal().global(localTriCorners[i]) << std::endl;
        assert( (element->geometry().corner(v) - nIt->geometry().global(localTriCorners[i]) ).two_norm() < eps);
      }
      break;

    case 4 :
      for (int i=0; i<4; i++) {
        int v = orientedSubface<dim>(element->type(), nIt->indexInInside(), i);
        //std::cout << element->geometry()[v] << "  " << nIt.intersectionGlobal().global(localQuadCorners[i]) << std::endl;
        assert( (element->geometry().corner(v) - nIt->geometry().global(localQuadCorners[i]) ).two_norm() < eps);
      }
      break;

    default :
      DUNE_THROW(NotImplemented, "Only line, triangle, and quad faces");
    }
  }
}

void test2dElements()
{

  typedef UGGrid<2> GridType;

  // //////////////////////////////////////////
  //   Create test cube
  // //////////////////////////////////////////
  GridFactory<GridType> triangleFactory;

  //   Insert vertices
  Dune::FieldVector<double,2> pos;

  pos[0] = 0;  pos[1] = 0;  triangleFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  triangleFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  triangleFactory.insertVertex(pos);

  // Insert element
  std::vector<unsigned int> cornerIDs(3);
  for (int i=0; i<3; i++)
    cornerIDs[i] = i;

  triangleFactory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,2), cornerIDs);

  //   Finish initialization
  GridType* triangle = triangleFactory.createGrid();

  //   Test
  test(*triangle);

  //   cleanup
  delete triangle;

  // //////////////////////////////////////////
  //   Create test quadrilateral
  // //////////////////////////////////////////
  GridFactory<GridType> quadrilateralFactory;

  //   Insert vertices
  pos[0] = 0;  pos[1] = 0;     quadrilateralFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;     quadrilateralFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;     quadrilateralFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;     quadrilateralFactory.insertVertex(pos);

  // Insert element
  cornerIDs.resize(4);
  for (int i=0; i<4; i++)
    cornerIDs[i] = i;

  quadrilateralFactory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,2), cornerIDs);

  //   Finish initialization
  GridType* quadrilateral = quadrilateralFactory.createGrid();

  //   Test
  test(*quadrilateral);

  //   cleanup
  delete quadrilateral;

}

void test3dElements()
{
  typedef UGGrid<3> GridType;

  // //////////////////////////////////////////
  //   Create test tetrahedron
  // //////////////////////////////////////////
  GridFactory<GridType> tetrahedronFactory;

  //   Insert vertices
  Dune::FieldVector<double,3> pos;

  pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    tetrahedronFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    tetrahedronFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    tetrahedronFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    tetrahedronFactory.insertVertex(pos);

  // Insert element
  std::vector<unsigned int> cornerIDs(4);
  for (int i=0; i<4; i++)
    cornerIDs[i] = i;

  tetrahedronFactory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3), cornerIDs);

  //   Finish initialization
  GridType* tetrahedron =  tetrahedronFactory.createGrid();

  //   Test
  test(*tetrahedron);

  //   cleanup
  delete tetrahedron;

  // //////////////////////////////////////////
  //   Create test cube
  // //////////////////////////////////////////
  GridFactory<GridType> cubeFactory;

  //   Insert vertices
  pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;  pos[2] = 1;    cubeFactory.insertVertex(pos);

  // Insert element
  cornerIDs.resize(8);
  for (int i=0; i<8; i++)
    cornerIDs[i] = i;

  cubeFactory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3), cornerIDs);

  //   Finish initialization
  GridType* cube = cubeFactory.createGrid();

  //   Test
  test(*cube);

  //   cleanup
  delete cube;

  // //////////////////////////////////////////
  //   Create test prism
  // //////////////////////////////////////////
  GridFactory<GridType> prismFactory;

  //   Insert vertices
  pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    prismFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    prismFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    prismFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    prismFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    prismFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    prismFactory.insertVertex(pos);

  // Insert element
  cornerIDs.resize(6);
  for (int i=0; i<6; i++)
    cornerIDs[i] = i;

  prismFactory.insertElement(GeometryType(GeometryType::prism,3), cornerIDs);

  //   Finish initialization
  GridType* prism = prismFactory.createGrid();

  //   Test
  test(*prism);

  //   cleanup
  delete prism;

  // //////////////////////////////////////////
  //   Create test pyramid
  // //////////////////////////////////////////
  GridFactory<GridType> pyramidFactory;

  //   Insert vertices
  pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    pyramidFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    pyramidFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    pyramidFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    pyramidFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    pyramidFactory.insertVertex(pos);

  // Insert element
  cornerIDs.resize(5);
  for (int i=0; i<5; i++)
    cornerIDs[i] = i;

  pyramidFactory.insertElement(GeometryType(GeometryType::pyramid,3), cornerIDs);

  //   Finish initialization
  GridType* pyramid = pyramidFactory.createGrid();

  //   Test
  test(*pyramid);

  //   cleanup
  delete pyramid;

}

int main(int argc, char *argv[]) try
{
  Dune::MPIHelper::instance(argc, argv);

  test2dElements();
  test3dElements();

  return 0;
}
catch (Exception e) {
  std::cout << e << std::endl;
}
