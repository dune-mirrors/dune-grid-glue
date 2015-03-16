#include "config.h"

#include <iostream>

#include <dune/grid-glue/common/projectionhelper.hh>

typedef double Real;
typedef Dune::FieldVector<Real, 3> Vector;
typedef Dune::FieldVector<Real, 2> Local;

bool
testInexactInverseProjection(const std::string& name,
                             const std::vector<Vector>& corners,
                             const std::vector<Vector>& directions,
                             const Vector& target,
                             const Real overlap,
                             const bool expectedResult,
                             const Local& expectedPreimage)
{
  Local preimage;
  const bool result = Projection::ProjectionHelper<2, 3, Real>::inexactInverseProjection(corners, directions, target, preimage, overlap);

  if (result != expectedResult) {
    std::cerr << std::boolalpha
              << "ERROR: " << name << ": result from projection is " << result << " (expected " << expectedResult << ")" << std::endl;
    return false;
  }
  if (expectedResult == false)
    return true;

  if (!((preimage - expectedPreimage).two_norm() < 1e-6)) {
    std::cerr << "ERROR: " << name << ": preimage (" << preimage << ") does not match expected value (" << expectedPreimage << ")" << std::endl;
    return false;
  }

  return true;
}

int main()
{
  bool passed = true;
  std::string name;
  Vector x1, x2, x3;
  Vector n1, n2, n3;
  Vector target;
  Local preimage;

  {
    name = "test case 1";

    x1[0] = 0; x1[1] = 0; x1[2] = 0;
    x2[0] = 1; x2[1] = 0; x2[2] = 0;
    x3[0] = 0; x3[1] = 1; x2[2] = 0;

    n1[0] = 0; n1[1] = 0; n1[2] = 1;
    n2 = n1;
    n3 = n1;

    target[0] = 0.25; target[1] = 0.25; target[2] = 1;

    preimage[0] = 0.25; preimage[1] = 0.25;

    passed &= testInexactInverseProjection(name, {x1, x2, x3}, {n1, n2, n3}, target, 1e-3, true, preimage);
  }

  {
    name = "test case 2";

    x1[0] = 0; x1[1] = 0; x1[2] = 0;
    x2[0] = 1; x2[1] = 0; x2[2] = 0;
    x3[0] = 0; x3[1] = 1; x3[2] = 0;

    n1[0] = 0; n1[1] = 0; n1[2] = 1;
    n2 = n1;
    n3 = n1;

    target[0] = 0.25; target[1] = 0.25; target[2] = -1;

    passed &= testInexactInverseProjection(name, {x1, x2, x3}, {n1, n2, n3}, target, 1e-3, false, preimage);
  }

  {
    name = "test case 3";

    x1[0] = 0;    x1[1] = 0;    x1[2] = 0;
    x2[0] = 1e-2; x2[1] = 0;    x2[2] = 0;
    x3[0] = 0;    x3[1] = 1e-2; x3[2] = 0;

    /* Pythagorean triples make nice normals */
    n1[0] = 3./5.; n1[1] = 0;     n1[2] = 4./5.;
    n2[0] = 0;     n2[1] = -5./13.; n2[2] = 12./13.;
    n3[0] = 0;     n3[1] = 0;     n3[2] = 1;

    target[0] = 0.5e-2; target[1] = 0.5e-2; target[2] = -1;

    passed &= testInexactInverseProjection(name, {x1, x2, x3}, {n1, n2, n3}, target, 1e-3, false, preimage);
  }

  {
    name = "test case 4";

    x1[0] = -0.0071231350000000001; x1[1] = 0.044047700000000002; x1[2] = -0.088647350000000014;
    x2[0] = -0.0069518649999999998; x2[1] = 0.043784000000000003; x2[2] = -0.08854455;
    x3[0] = -0.0067651500000000002; x3[1] = 0.044021000000000005; x3[2] = -0.088661900000000002;
    //x1 *= 1000; x2 *= 1000; x3 *= 1000;

    n1[0] =  0.13347956722282625;   n1[1] =  0.35331475935692608; n1[2] = 0.92593298135154711;
    n2[0] =  0.051280560993381304;  n2[1] =  0.149455098537458;   n2[2] = 0.98743783479536718;
    n3[0] = -0.0093774282639142777; n3[1] = 0.35790454317655113;  n3[2] = 0.93371109119081341;

    target[0] = -0.0117477866623666; target[1] = 0.0012007769979558108; target[2] = -0.080164187395937322;
    //target *= 1000;

    passed &= testInexactInverseProjection(name, {x1, x2, x3}, {n1, n2, n3}, target, 1e-3, false, preimage);
  }

  return passed ? 0 : 1;
}
