// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <config.h>

#include <dune/common/fvector.hh>
#include <dune/grid/uggrid.hh>

#include <dune/glue/misc/geometry.hh>

using namespace Dune;

typedef UGGrid<3> GridType;

void test(const GridType& grid) {

  FieldVector<double,2> localTriCorners[3];
  localTriCorners[0][0] = 0;   localTriCorners[0][1] = 0;
  localTriCorners[1][0] = 1;   localTriCorners[1][1] = 0;
  localTriCorners[2][0] = 0;   localTriCorners[2][1] = 1;

  FieldVector<double,2> localQuadCorners[4];
  localQuadCorners[0][0] = 0;   localQuadCorners[0][1] = 0;
  localQuadCorners[1][0] = 1;   localQuadCorners[1][1] = 0;
  localQuadCorners[2][0] = 0;   localQuadCorners[2][1] = 1;
  localQuadCorners[3][0] = 1;   localQuadCorners[3][1] = 1;

  double eps = 1e-6;

  GridType::Codim<0>::EntityPointer element = grid.leafbegin<0>();

  GridType::Codim<0>::Entity::LevelIntersectionIterator nIt = element->ilevelbegin();
  GridType::Codim<0>::Entity::LevelIntersectionIterator nEndIt = element->ilevelend();

  std::cout << "Testing: " << element->type() << std::endl;

  for (; nIt!=nEndIt; ++nIt) {

    printf("Testing subface %d\n", nIt->indexInInside());

    switch (nIt->geometry().corners()) {

    case 3 :
      for (int i=0; i<3; i++) {
        int v = orientedSubface<3>(element->type(), nIt->indexInInside(), i);
        //std::cout << element->geometry()[v] << "  "
        //          << nIt.intersectionGlobal().global(localTriCorners[i]) << std::endl;
        assert( (element->geometry().corner(v) - nIt->geometry().global(localTriCorners[i]) ).two_norm() < eps);
      }
      break;

    case 4 :
      for (int i=0; i<4; i++) {
        int v = orientedSubface<3>(element->type(), nIt->indexInInside(), i);
        //std::cout << element->geometry()[v] << "  " << nIt.intersectionGlobal().global(localQuadCorners[i]) << std::endl;
        assert( (element->geometry().corner(v) - nIt->geometry().global(localQuadCorners[i]) ).two_norm() < eps);
      }
      break;

    default :
      DUNE_THROW(NotImplemented, "Only triangle and quad faces");
    }
  }
}

int main()
{

  // //////////////////////////////////////////
  //   Create test cube
  // //////////////////////////////////////////
  GridFactory<GridType> cubeFactory;

  //   Insert vertices
  Dune::FieldVector<double,3> pos;

  pos[0] = 0;  pos[1] = 0;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;  pos[2] = 0;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 0;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 0;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 0;  pos[1] = 1;  pos[2] = 1;    cubeFactory.insertVertex(pos);
  pos[0] = 1;  pos[1] = 1;  pos[2] = 1;    cubeFactory.insertVertex(pos);

  // Insert element
  std::vector<unsigned int> cornerIDs(8);
  for (int i=0; i<8; i++)
    cornerIDs[i] = i;

  cubeFactory.insertElement(Dune::GeometryType(Dune::GeometryType::cube,3), cornerIDs);

  //   Finish initialization
  std::auto_ptr<GridType> cube(cubeFactory.createGrid());

  //   Test
  test(*cube);

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
  std::auto_ptr<GridType> prism(prismFactory.createGrid());

  //   Test
  test(*prism);

  // //////////////////////////////////////////
  //   Create test cube
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
  std::auto_ptr<GridType> pyramid(pyramidFactory.createGrid());

  //   Test
  test(*pyramid);

  return 0;
}
