// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    geometry.h
 *  Version:     1.0
 *  Created on:  Jan 28, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: collection of useful geometry-related functions
 *  subversion:  $Id$
 *
 */
/**
 * @file geometry.h
 * @brief header with helpful geometric computation related functions
 */

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/fixedarray.hh>
#include <dune/common/exceptions.hh>

using namespace Dune;



/**
 * @brief transforms Dune-style local coordinates in 1D and 2D to
 * barycentric coordinates.
 * Since the dimension of bar. coords is one higher than that of
 * the Dune local coords, a local coordinate for Dune's coordinate
 * system's origin is introduced which is located at corner 0 in
 * edges and triangles.
 *
 * @param
 * @return
 */
template<typename K, int dim>
inline FieldVector<K, dim-1> barycentricToReference(const FieldVector<K, dim>& bar)
{
  FieldVector<K, dim-1> result;
  if (dim == 2)
  {
    // use Dune style indexing for edge (0 - 1) where the origin (0)
    // of the local coords is at corner 0
    result[0] = bar[1];
  }
  else if (dim == 3)
  {
    // use Dune style indexing (0 - 1 - 2, ccw) where the origin (0,0)
    // of the local coords is at corner 0
    result[0] = bar[1];
    result[1] = bar[2];
  }
  else
  {
    DUNE_THROW(NotImplemented, "dimension not implemented");
  }
  return result;
}


/**
 * @brief transforms Dune-style local coordinates in 1D and 2D to
 * barycentric coordinates.
 * Since the dimension of bar. coords is one higher than that of
 * the Dune local coords, a local coordinate for Dune's coordinate
 * system's origin is introduced which is located at corner 0 in
 * edges and triangles.
 *
 * @param
 * @return
 */
template<typename K, int dim>
inline FieldVector<K, dim+1> referenceToBarycentric(const FieldVector<K, dim>& ref)
{
  FieldVector<K, dim+1> result;
  if (dim == 1)
  {
    // use Dune style indexing for edge (0 - 1) where the origin (0)
    // of the local coords is at corner 0
    result[1] = ref[0];
    result[0] = 1.0 - result[1];
  }
  else if (dim == 2)
  {
    // use Dune style indexing (0 - 1 - 2, ccw) where the origin (0,0)
    // of the local coords is at corner 0
    result[1] = ref[0];
    result[2] = ref[1];
    result[0] = 1.0 - result[2] - result[1];
  }
  else
  {
    DUNE_THROW(NotImplemented, "dimension not implemented");
  }
  return result;
}


/**
 * @brief compute the area of a 2D or 3D simplex
 * @param coords pointer to the 1st of 2 or 3 simplex corners
 * @return the "area" between the coordinates
 */
template<typename K, int dim>
inline K computeArea(const FieldVector<K, dim>& coords)
{
  if (dim == 2)
  {
    return (coords[0] - coords[1]).two_norm();
  }
  else if (dim == 3)
  {
    FieldVector<K, dim> normal, v1 = coords[0] - coords[1], v0 = coords[2] - coords[1];
    normal[0] = v0[1]*v1[2] - v0[2]*v1[1];
    normal[1] = v0[2]*v1[0] - v0[0]*v1[2];
    normal[2] = v0[0]*v1[1] - v0[1]*v1[0];
    return 0.5 * normal.two_norm();
  }
  else
  {
    DUNE_THROW(NotImplemented, "dimension not implemented");
  }
}


/**
 * @brief generic linear interpolation for convenience
 *
 * @param values a vector of values at the edge's corners,
 * order corresponding to Dune quadrilateral's indexing scheme
 * @param clocal local coordinate, starting at 1st point with 0
 * @return the interpolated result
 */
template<typename K>
inline K interpolateLinear(const FieldVector<K, 2>& values, const K clocal)
{
  return (1.0-clocal)*values[0] + clocal*values[1];
}


/**
 * @brief generic linear interpolation for convenience
 *
 * @param value0 1st point's value
 * @param value1 2nd point's value
 * @param clocal local coordinate, starting at 1st point with 0
 * @return the interpolated result
 */
template<typename K>
inline K interpolateLinear(const K value0, const K value1, const K clocal)
{
  return (1.0-clocal)*value0 + clocal*value1;
}


/**
 * @brief generic linear interpolation for convencience
 *
 * The template parameter V stands for a random accessible, default- and copy-constructible
 * container type.
 * @param values0 1st point's values
 * @param values1 2nd point's values
 * @param clocal local coordinate, starting at 1st point with 0
 * @param datadim the number of data dimensions in type V that are to be interpolated
 * (in range 0 .. datadim-1)
 * @return the interpolated result
 */
template<typename K, typename V>
inline void interpolateLinear(const V& values0, const V& values1, const K clocal, V& result, const unsigned int datadim = 1)
{
  // evaluate component-wise
  for (unsigned int j = 0; j < datadim; ++j)
    result[j] = (1.0-clocal)*values0[j] + clocal*values1[j];
}


/**
 * @brief generic bilinear interpolation on a Dune-style quadrilateral
 *
 * @param values a vector of values at the quadrilateral's corners,
 * order corresponding to Dune quadrilateral's indexing scheme
 * @param clocal a vector containing the Dune-style local coordinates,
 * i.e. at first index is point (0,0), etc.
 * @return the interpolated result
 */
template<typename K>
inline K interpolateBilinear(const FieldVector<K, 4>& values, const FieldVector<K, 2>& clocal)
{
  return (1.0-clocal[0])*((1.0-clocal[1])*values[0] + clocal[1]*values[2]) + clocal[0]*((1.0-clocal[1])*values[1] + clocal[1]*values[3]);
}


/**
 * @brief generic bilinear interpolation on a Dune-style quadrilateral
 *
 * The template parameter V stands for a random accessible, default- and copy-constructible
 * container type.
 * @param values a vector of multidimensional data associated with the quadrilateral's corners,
 * order corresponding to Dune quadrilateral's indexing scheme
 * @param clocal a vector containing the Dune-style local coordinates,
 * i.e. at first index is point (0,0), etc.
 * @param datadim the number of data dimensions in type V that are to be interpolated
 * (in range 0 .. datadim-1)
 * @return the interpolated result
 */
template<typename K, typename V>
inline void interpolateBilinear(const FieldVector<V, 4>& values, const FieldVector<K, 2>& clocal, V& result, const unsigned int datadim = 1)
{
  // evaluate component-wise
  for (unsigned int j = 0; j < datadim; ++j)
    result[j] = (1.0-clocal[0])*((1.0-clocal[1])*values[0][j] + clocal[1]*values[2][j]) + clocal[0]*((1.0-clocal[1])*values[1][j] + clocal[1]*values[3][j]);
}


/**
 * @brief generic bilinear interpolation on a Dune-style quadrilateral
 *
 * @param values a vector of values at the simplice's corners,
 * order corresponding to Dune quadrilateral's indexing scheme
 * @param clocal a vector containing the Dune-style local coordinates,
 * @return the interpolated result
 */
template<typename K>
inline K interpolateBilinear(const K value0, const K value1, const K value2, const K value3, const FieldVector<K, 2>& clocal)
{
  return (1.0-clocal[0])*((1.0-clocal[1])*value0 + clocal[1]*value2) + clocal[0]*((1.0-clocal[1])*value1 + clocal[1]*value3);
}


/**
 * @brief generic bilinear interpolation on a Dune-style quadrilateral
 *
 * The template parameter V stands for a random accessible, default- and copy-constructible
 * container type.
 * @param values0 1st point's values
 * @param values1 2nd point's values
 * @param values2 3rd point's values
 * @param values3 4th point's values
 * @param clocal a vector containing the Dune-style local coordinates,
 * @param datadim the number of data dimensions in type V that are to be interpolated
 * (in range 0 .. datadim-1)
 * @return the interpolated result
 */
template<typename K, typename V>
inline void interpolateBilinear(const V& values0, const V& values1, const V& values2, const V& values3, const FieldVector<K, 2>& clocal, V& result, const unsigned int datadim = 1)
{
  // evaluate component-wise
  for (unsigned int j = 0; j < datadim; ++j)
    result[j] = (1.0-clocal[0])*((1.0-clocal[1])*values0[j] + clocal[1]*values2[j]) + clocal[0]*((1.0-clocal[1])*values1[j] + clocal[1]*values3[j]);
}


/**
 * @brief generic evaluation of barycentric coordinates, 1D specialization
 *
 * @param values a vector of values at the simplice's corners,
 * order corresponding to the barycentric coordinates
 * @param barc a vector containing the barycentric coordinates
 * @return the interpolated result
 */
template<int dim, typename K>
inline K interpolateBarycentric(const FieldVector<K, dim>& values, const FieldVector<K, dim>& barc)
{
  K result = 0.0;
  for (typename FieldVector<K, dim>::size_type i = 0; i < dim; ++i)
    result += barc[i]*values[i];
  return result;
}


/**
 * @brief generic evaluation of barycentric coordinates
 *
 * The template parameter V stands for a random accessible, default- and copy-constructible
 * container type.
 * @param values a vector of vector-valued data at the simplice's corners,
 * order corresponding to the barycentric coordinates
 * @param barc a vector containing the barycentric coordinates
 * @param datadim the number of data dimensions in type V that are to be interpolated
 * (in range 0 .. datadim-1)
 * @return the interpolated result
 */
template<int dim, typename K, typename V>
inline void interpolateBarycentric(const array<V, dim>& values, const FieldVector<K, dim>& barc, V& result, const unsigned int datadim = 1)
{
  // initialize with zero
  for (unsigned int j = 0; j < datadim; ++j)
    result[j] = 0.0;
  // accumulate the contribution from the simplex corners
  for (typename FieldVector<K, dim>::size_type i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < datadim; ++j)
      result[j] += barc[i]*values[i][j];
}


/**
 * @brief computes barycentric coordinates for a given dim-D point and a dim-D simplex in dim-D space.
 * @param corners a vector with the (mathemat. pos. oriented) corners of the simplex
 * @param p the point (should be inside)
 * @return the barycentric coordinates (not checked, i.e. invalid if point not in simplex!)
 */
template<int dim, typename K>
FieldVector<K, dim+1> barycentric(const array<FieldVector<K, dim>, dim+1>& corners, const FieldVector<K, dim>& p)
{
  FieldVector<K, dim+1> result;
  typedef FieldVector<K, dim> CoordType;
  if (dim == 1)
  {
    // compute the ratio of the "distances"
    result[0] = (p[0] - corners[1][0]) / (corners[0][0] - corners[1][0]);

    // fill the result vector
    result[1] = 1.0 - result[0];
  }
  else if (dim == 2)
  {
    // the problem is mapped to solving a 2x2 linear system

    // build the matrix
    typename Dune::FieldMatrix<K, 2, 2> A(0.0);
    A[0][0] = corners[0][0] - corners[2][0];
    A[0][1] = corners[1][0] - corners[2][0];
    A[1][0] = corners[0][1] - corners[2][1];
    A[1][1] = corners[1][1] - corners[2][1];

    // build the rhs
    typename Dune::FieldVector<K, 2> x, b;
    b[0] = p[0] - corners[2][0];
    b[1] = p[1] - corners[2][1];

    // solve
    A.solve(x, b);

    // fill the result vector
    result[0] = x[0];
    result[1] = x[1];
    result[2] = 1.0 - result[0] - result[1];
  }
  else
  {
    DUNE_THROW(NotImplemented, "Dimension " << dim << " implementation missing.");
  }
  return result;
}


/**
 * @brief computes barycentric coordinates for in 1D
 * @param corners two values denoting the "end points" of the segment
 * @param p the "point" (should be inside)
 * @return the barycentric coordinates (not checked, i.e. invalid if point not in simplex!)
 */
template<typename K>
inline FieldVector<K, 2> barycentric(const FieldVector<K, 2>& corners, const K& p)
{
  FieldVector<K, 2> result;
  result[0] = (p - corners[1]) / (corners[0] - corners[1]);
  result[1] = 1.0 - result[0];
  return result;
}


/**
 * @brief computes barycentric coordinates for a given (dim-1)-D point and a (dim-1)-D simplex in
 * the dim-D space.
 * @param corners a vector with the (mathemat. pos. oriented) corners of the simplex
 * @param p the point (should be inside)
 * @return the barycentric coordinates (not checked, i.e. invalid if point not in simplex!)
 */
template<int dim, typename K>
FieldVector<K, dim> hyperBarycentric(const Dune::array<FieldVector<K, dim>, dim>& corners, const FieldVector<K, dim>& p)
{
  typedef FieldVector<K, dim>  CoordType;

  CoordType result(1.0);
  // In both cases projection is done to reduce the dimension of the calculation
  // by one. Consistency in the results is ensured.
  // one implementation enough, optimization will remove condition
  if (dim == 2)
  {
    // simply compare length of projected segments between the points
    // (chose the dimension that is not so "crowded")
    if (fabs(corners[0][0] - corners[1][0]) > fabs(corners[0][1] - corners[1][1]))

      result[0] = (p[0] - corners[1][0]) / (corners[0][0] - corners[1][0]);
    else
      result[0] = (p[0] - corners[1][0]) / (corners[0][0] - corners[1][0]);
    result[1] = 1.0 - result[0];
  }
  else if (dim == 3)
  {
    // the case that the point is not in the same plane as the triangle is not considered!
    // (user's responsibility to ensure!)

    // create copies of the input vectors marked with prefix "t" in name
    // (might be modified during computation)
    Dune::FieldVector<CoordType, dim> tcorners(corners);
    CoordType tp(p);

    // Rather than projecting we will be rotating the triangle until it lies
    // in a plane parallel to the xy-plane. In order to build the matrix we
    // will have to compute a rotation axis and an angle first.
    // They will be stored in these two variables (variables will be used as temporary
    // containers for other values during the computation).
    CoordType rot(0.0);             // rotation axis later, rotating ccw about it
    K cosa, sina;                   // later cosine alpha and sine alpha
    {
      CoordType v1 = corners[2] - corners[1], v2 = corners[0] - corners[1];
      CoordType n;
      n[0] = v1[1]*v2[2] - v1[2]*v2[1];
      n[1] = v1[2]*v2[0] - v1[0]*v2[2];
      n[2] = v1[0]*v2[1] - v1[1]*v2[0];

      // 1.) get projected normal's normalized normal in xy-plane ;-)
      //     (which is the rotation axis then)

      // ...therefore use cosa variable for normalization temporarily
      cosa = sqrt(n[0]*n[0] + n[1]*n[1]);

      // if cosa is now very small that means that no rotation is necessary
      // and we can just go to the part where the coordinates are actually comupted
      if (cosa < 1E-6)
        goto COMPUTE_BARYCENTRIC_COORDS;

      // otherwise just go on with the computation
      rot[0] = n[1] / cosa;
      rot[1] = -n[0] / cosa;

      // 2.) compute the angle to the z=0 plane (xy-plane)

      // get the angle between the up-vector (0,0,1)^T and the normal vector
      cosa = n[2] / n.two_norm();

      // now get the sine of alpha, too
      sina = sqrt(1 - cosa*cosa);

      // build the rotation matrix
      typename Dune::FieldMatrix<K, dim, dim> Q(1.0 - cosa);
      Q[0][0] *= rot[0]*rot[0];
      Q[0][0] += cosa;

      Q[0][1] *= rot[0]*rot[1];
      //Q[0][1] -= rot[2]*sina;   // 0

      //Q[0][2] *= rot[0]*rot[2]; // makes it 0...
      Q[0][2] = rot[1]*sina;                        // ...=> assigning here

      Q[1][0] *= rot[0]*rot[1];
      //Q[1][0] += rot[2]*sina;   // 0

      Q[1][1] = rot[1]*rot[1];
      Q[1][1] += cosa;

      //Q[1][2] = rot[1]*rot[2];  // makes it 0...
      Q[1][2] = -rot[0]*sina;                       // ...=> assigning here


      //Q[2][0] = rot[0]*rot[2];  // makes it 0...
      Q[2][0] = -rot[1]*sina;                       // ...=> assigning here

      //Q[2][1] = rot[1]*rot[2];  // makes it 0...
      Q[2][1] = rot[0]*sina;                        // ...=> assigning here

      //Q[2][2] *= rot[2]*rot[2]; // makes it 0...
      Q[2][2] = cosa;                               // ...=> assigning here

      // apply the matrix so point and triangle are now
      // lying in the xy-plane
      for (unsigned int i = 0; i < dim; ++i)
        Dune::FMatrixHelp::multAssign(Q, corners[i], tcorners[i]);
      Dune::FMatrixHelp::multAssign(Q, p, tp);
    }

COMPUTE_BARYCENTRIC_COORDS:
    // now compute the coordinates totally forgetting about the z-dimension!

    // the problem is mapped to solving a 2x2 linear system

    // build the matrix
    typename Dune::FieldMatrix<K, 2, 2> A(0.0);
    A[0][0] = tcorners[0][0] - tcorners[2][0];
    A[0][1] = tcorners[1][0] - tcorners[2][0];
    A[1][0] = tcorners[0][1] - tcorners[2][1];
    A[1][1] = tcorners[1][1] - tcorners[2][1];

    // build the rhs
    typename Dune::FieldVector<K, 2> x, b;
    b[0] = tp[0] - tcorners[2][0];
    b[1] = tp[1] - tcorners[2][1];

    // solve
    A.solve(x, b);

    // fill the result vector
    result[0] = x[0];
    result[1] = x[1];
    result[2] = 1.0 - result[0] - result[1];
  }
  else
  {
    DUNE_THROW(NotImplemented, "Dimension " << dim << " implementation missing.");
  }
  return result;
}


/**
 * @brief computes coordinate form of "plane" in 2D
 *
 * Analogous to the 3D version this computes the
 * non-normalized Hesse form of a straight line in 2D,
 * i.e. the coefficients a, b and c in
 * one valid representation of the straight line
 * a x_0 + b x_1 + c = 0
 * This form can be used to quickly determine the normal
 * vector as well as to decide on which side of the line
 * a given point is located.
 * @param p an array of coordinates holding 2 points on the line,
 * if one assumes that it's a vector pointing forward,
 * the resulting normal vector will point to the right
 * @return a vector with a, b and c
 */
template<typename K>
FieldVector<K, 3> getHyperplaneCoeffs(const array<FieldVector<K, 2>, 2>& p, bool normalize = false)
{
  FieldVector<K, 3> coeffs(0.0);
  FieldVector<K, 2> v = p[1]-p[0];
  coeffs[0] = p[1][1]-p[0][1];
  coeffs[1] = -(p[1][0]-p[0][0]);
  // plug in the point to get the offset in coeffs[2]
  coeffs[2] = -(coeffs[0]*p[0][0]+coeffs[1]*p[0][1]);
  if (normalize)
  {
    K length = 0.0;
    for (int i = 0; i < 2; ++i)
      length += coeffs[i]*coeffs[i];
    length = sqrt(length);
    for (int i = 0; i < 3; ++i)
      coeffs[i] /= length;
  }
  return coeffs;
}


/**
 * @brief computes coordinate form of plane in 3D
 *
 * Computes the non-normalized form of the Hesse form of
 * a plane in 3D, i.e. the coefficients a, b, c and d in
 * one valid representation of the plane
 * a x_0 + b x_1 + c x_2 + d = 0
 * This form can be used to quickly determine the normal
 * vector as well as to decide on which side of the plane
 * a given point is located.
 * @param p an array of coordinates holding 3 points in the plan,
 * if one assumes that it's a triangle given in cyclic order,
 * the resulting normal vector will point upwards
 * @return a vector with a, b, c and d
 */
template<typename K>
FieldVector<K, 4> getHyperplaneCoeffs(const array<FieldVector<K, 3>, 3>& p, bool normalize = false)
{
  FieldVector<K, 4> coeffs(0.0);
  // compute the cross product from (p[0]-p[1]) and (p[2]-p[1])
  FieldVector<K, 3> v2 = p[0]-p[1];
  FieldVector<K, 3> v1 = p[2]-p[1];
  coeffs[0] = v1[1]*v2[2]-v1[2]*v2[1];
  coeffs[1] = v1[2]*v2[0]-v1[0]*v2[2];
  coeffs[2] = v1[0]*v2[1]-v1[1]*v2[0];
  // plug in one point to get the offset in coeffs[3]
  coeffs[3] = -(coeffs[0]*p[0][0]+coeffs[1]*p[0][1]+coeffs[2]*p[0][2]);
  if (normalize)
  {
    K length = 0.0;
    for (int i = 0; i < 3; ++i)
      length += coeffs[i]*coeffs[i];
    length = sqrt(length);
    for (int i = 0; i < 4; ++i)
      coeffs[i] /= length;
  }
  return coeffs;
}


template<typename K>
FieldVector<K, 2> computeNormal(const array<FieldVector<K, 2>, 1>& v)
{
  FieldVector<K, 2> result(v[0][1]);
  result[1] = -v[0][0];
  return result;
}


template<typename K>
FieldVector<K, 3> computeNormal(const array<FieldVector<K, 3>, 2>& v)
{
  FieldVector<K, 3> result(0.0);
  result[0] = v[0][1]*v[1][2] - v[0][2]*v[1][1];
  result[1] = v[0][2]*v[1][0] - v[0][0]*v[1][2];
  result[2] = v[0][0]*v[1][1] - v[0][1]*v[1][0];
  return result;
}


// can not be used to do dimension independent programming (two arguments)
template<typename K>
FieldVector<K, 3> computeNormal(const FieldVector<K, 3> &v0, const FieldVector<K, 3> &v1)
{
  FieldVector<K, 3> result(0.0);
  result[0] = v0[1]*v1[2] - v0[2]*v1[1];
  result[1] = v0[2]*v1[0] - v0[0]*v1[2];
  result[2] = v0[0]*v1[1] - v0[1]*v1[0];
  return result;
}



template<typename GEO>
void printGeometry(const GEO &geo, const char *name)
{
  STDOUT(name << " :");
  for (int i = 0; i < geo.corners(); ++i)
    STDOUT(" (" << geo.corner(i) << ")");
  STDOUTLN("");
}


/**
 * @brief for a simplex geometry object (2D or 3D only!) checks if a point is contained
 *
 * The given point is in the geometry's global coordinates and it is checked
 * if the point is located inside of the simplicial shape of the geometry's domain.
 * E.g. you can use this routine to check whether a point in local codim 0 coordinates is
 * contained in one of the codim 0 entity's intersections. The intersections' checkInside
 * member fails in this case...
 * IMPORTANT: call this only for geometries that map from 1D to 2D or from 2D to 3D!!
 */
template<typename GEO, typename K>
bool checkInside(const GEO &geo, const Dune::FieldVector<K, GEO::coorddimension> &global, K tol = 1E-10)
{
  const int cdim = GEO::coorddimension;

  // idea map the point to local coordinates and then back to globals again,
  // if the result is too different (decide by given tolerance) the point has
  // not been on the face
  // NOTE: this is the most generic approach and despite of the function description
  // works for all kind of geometries!
  Dune::FieldVector<K, cdim> ownglobal = geo.global(geo.local(global));
  for (int i = 0; i < cdim; ++i)
    if (ownglobal[i] < global[i] - tol || ownglobal[i] > global[i] + tol)
      return false;
  return true;
}


template <int dim>
int orientedSubface(const Dune::GeometryType& type, int face, int vertex)
{
  const Dune::ReferenceElement<double,dim>& refElement =
    Dune::ReferenceElements<double, dim>::general(type);

  // Quadrilateral
  if (type.isQuadrilateral() && (face==0 || face==3))
    return refElement.subEntity(face,1, (vertex+1)%2, dim);

  // Pyramid, quad face:  turn 0 3 2 1 into 0 3 1 2
  if (type.isPyramid() && face==0 && vertex>=2)
    return refElement.subEntity(face,1, ((vertex-1)%2)+2, dim);

  // Prism: swap vertices 2 and 3 for all quad faces
  if (type.isPrism() && face>=1 && face<4 && vertex>=2)
    return refElement.subEntity(face,1, ((vertex-1)%2)+2, dim);

  // Hexahedron:
  // face 0: 0 2 4 6 --> 2 0 6 4
  // face 3: 2 3 6 7 --> 3 2 7 6
  if (type.isHexahedron() && (face==0 || face==3))
    return refElement.subEntity(face,1, (vertex<2) ? (vertex+1)%2 : ((vertex-1)%2)+2, dim);


  // face 4: 0 1 2 3 --> 0 2 1 3
  if (type.isHexahedron() && face==4) {
    const int renumber[4] = {0,2,1,3};
    return refElement.subEntity(face,1, renumber[vertex], dim);
  }

  return refElement.subEntity(face,1,vertex,dim);
}


/**
 * @brief computes the local coordinate of a given face's particular corner
 * if the corner was determined using orientedSubface (algorithm embedded)
 *
 * !!WARNING!!
 * This only works for codim 1 faces!
 *
 * @param gt the geometry type of the element
 * @param face the index of the element's face
 * @param corner the index of the face's corner
 * @return local coordinates for the corner in the element's local coordinates
 */
template<typename K, int dim>
Dune::FieldVector<K, dim> cornerLocalInRefElement(const Dune::GeometryType& type, int face, int vertex)
{
  unsigned int cornerlocal = 0;

  const Dune::ReferenceElement<K, dim>& refElement =
    Dune::ReferenceElements<K, dim>::general(type);

  // Quadrilateral
  if (type.isQuadrilateral() && (face==0 || face==3))
    cornerlocal = refElement.subEntity(face,1, (vertex+1)%2, dim);

  // Pyramid, quad face:  turn 0 3 2 1 into 0 3 1 2
  if (type.isPyramid() && face==0 && vertex>=2)
    cornerlocal = refElement.subEntity(face,1, ((vertex-1)%2)+2, dim);

  // Prism: swap vertices 2 and 3 for all quad faces
  if (type.isPrism() && face>=1 && face<4 && vertex>=2)
    cornerlocal = refElement.subEntity(face,1, ((vertex-1)%2)+2, dim);

  // Hexahedron:
  // face 0: 0 2 4 6 --> 2 0 6 4
  // face 3: 2 3 6 7 --> 3 2 7 6
  if (type.isHexahedron() && (face==0 || face==3))
    cornerlocal = refElement.subEntity(face,1, (vertex<2) ? (vertex+1)%2 : ((vertex-1)%2)+2, dim);


  // face 4: 0 1 2 3 --> 0 2 1 3
  if (type.isHexahedron() && face==4)
  {
    const int renumber[4] = {0,2,1,3};
    cornerlocal = refElement.subEntity(face,1, renumber[vertex], dim);
  }

  cornerlocal = refElement.subEntity(face,1,vertex,dim);

  //	STDOUTLN("cornerLocalInRefElement(" << type << ", " << face << ", " << vertex << ") => cornerlocal=" << cornerlocal << " position=(" << refElement.position(cornerlocal, dim) << ")" );

  // return the vertex' center of gravity, i.e. the vertex' location
  return refElement.position(cornerlocal, dim);
}


#endif // GEOMETRY_H_
