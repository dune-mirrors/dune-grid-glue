#include "config.h"

#include <iostream>
#include <dune/common/fvector.hh>
#include <dune/grid-glue/common/projectionhelper.hh>

typedef double Real;

/*
 * Make sure Projection::projection<2, 3, double> works independent of the
 * scale of the world geometry.
 */
bool
testScaledProjection(const Real scale)
{
  bool pass = true;

  const static int dim = 3;
  typedef Dune::FieldVector<Real, dim> Vector;
  typedef Dune::FieldVector<Real, dim-1> Local;

  Vector corner(0);

  Vector direction(0);
  direction[0] = 1;

  std::vector<Vector> target(dim);
  for (int i = 0; i < dim; ++i) {
    Vector corner(0);
    corner[0] = 1;
    corner[i] = 1;
    corner *= scale;
    target[i] = corner;
  }

  /* Project (0,0,0) on (0,0) */
  {
    Vector corner(0);

    Local image;
    bool result = Projection::projection<dim-1, dim, Real>(corner, direction, target, image);
    if (!result) {
      std::cerr << "testScaledProjection(" << scale << "): projection failed" << std::endl;
      pass = false;
    }

    /* expected image is (0,0) */
    if (!(image.infinity_norm() < 1e-6)) {
      std::cerr << "testScaledProjection(" << scale << "): image = "
                << image << " does not match expected value (0,0)" << std::endl;
      pass = false;
    }
  }

  /* Project (0,10,0), should fail */
  {
    Vector corner(0);
    corner[1] = 10 * scale;

    Local image;
    bool result = Projection::projection<dim-1, dim, Real>(corner, direction, target, image);

    if (result) {
      std::cerr << "FAIL: testScaledProjection(" << scale << "): projection succeeded for "
                << corner << ", but is expected to fail." << std::endl
                << "      image is " << image << std::endl;
      pass = false;
    }
    else {
      std::cout << "OK: projection failed, scale = " << scale << std::endl;
    }
  }

  return pass;
}

int main()
{
  bool pass = true;
  pass &= testScaledProjection(1e10);
  pass &= testScaledProjection(1);
  pass &= testScaledProjection(1e-3);
  pass &= testScaledProjection(1e-9);
  pass &= testScaledProjection(1e-100);
  return pass ? 0 : 1;
}
