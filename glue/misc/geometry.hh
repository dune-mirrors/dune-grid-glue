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
 * @file
 * @brief header with helpful geometric computation related functions
 */

#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/geometrytype.hh>

#include <dune/glue/misc/orientedsubface.hh>

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
inline Dune::FieldVector<K, dim-1> barycentricToReference(const Dune::FieldVector<K, dim>& bar)
{
  Dune::FieldVector<K, dim-1> result;
  for (int i=0; i<dim-1; i++)
    result[i] = bar[i+1];

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
inline Dune::FieldVector<K, dim+1> referenceToBarycentric(const Dune::FieldVector<K, dim>& ref)
{
  Dune::FieldVector<K, dim+1> result;
  result[0] = 1.0;
  for (int i=0; i<dim; i++) {
    result[i+1] = ref[i];
    result[0] -= ref[i];
  }

  return result;
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
inline K interpolateLinear(const Dune::FieldVector<K, 2>& values, const K clocal)
{
  return (1.0-clocal)*values[0] + clocal*values[1];
}

template <typename K, int dimworld, int N>
inline Dune::FieldVector<K,dimworld> interpolateLinear(const Dune::array<Dune::FieldVector<K, dimworld>, N>& values,
                                                       const Dune::FieldVector<K, N-1>& local)
{
  Dune::FieldVector<K,dimworld> result = values[0];
  for (int i=1; i<N; i++)
    result.axpy(local[i-1], values[i]-values[0]);
  return result;
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
inline K interpolateBilinear(const Dune::FieldVector<K, 4>& values, const Dune::FieldVector<K, 2>& clocal)
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
inline void interpolateBilinear(const Dune::FieldVector<V, 4>& values, const Dune::FieldVector<K, 2>& clocal, V& result, const unsigned int datadim = 1)
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
inline K interpolateBilinear(const K value0, const K value1, const K value2, const K value3, const Dune::FieldVector<K, 2>& clocal)
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
inline void interpolateBilinear(const V& values0, const V& values1, const V& values2, const V& values3, const Dune::FieldVector<K, 2>& clocal, V& result, const unsigned int datadim = 1)
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
inline K interpolateBarycentric(const Dune::FieldVector<K, dim>& values, const Dune::FieldVector<K, dim>& barc)
{
  K result = 0.0;
  for (typename Dune::FieldVector<K, dim>::size_type i = 0; i < dim; ++i)
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
inline void interpolateBarycentric(const Dune::array<V, dim>& values, const Dune::FieldVector<K, dim>& barc, V& result, const unsigned int datadim = 1)
{
  // initialize with zero
  for (unsigned int j = 0; j < datadim; ++j)
    result[j] = 0.0;
  // accumulate the contribution from the simplex corners
  for (typename Dune::FieldVector<K, dim>::size_type i = 0; i < dim; ++i)
    for (unsigned int j = 0; j < datadim; ++j)
      result[j] += barc[i]*values[i][j];
}



template<typename K>
Dune::FieldVector<K, 2> computeNormal(const Dune::array<Dune::FieldVector<K, 2>, 1>& v)
{
  Dune::FieldVector<K, 2> result(v[0][1]);
  result[1] = -v[0][0];
  return result;
}


template<typename K>
Dune::FieldVector<K, 3> computeNormal(const Dune::array<Dune::FieldVector<K, 3>, 2>& v)
{
  Dune::FieldVector<K, 3> result(0.0);
  result[0] = v[0][1]*v[1][2] - v[0][2]*v[1][1];
  result[1] = v[0][2]*v[1][0] - v[0][0]*v[1][2];
  result[2] = v[0][0]*v[1][1] - v[0][1]*v[1][0];
  return result;
}


// can not be used to do dimension independent programming (two arguments)
template<typename K>
Dune::FieldVector<K, 3> computeNormal(const Dune::FieldVector<K, 3> &v0, const Dune::FieldVector<K, 3> &v1)
{
  Dune::FieldVector<K, 3> result(0.0);
  result[0] = v0[1]*v1[2] - v0[2]*v1[1];
  result[1] = v0[2]*v1[0] - v0[0]*v1[2];
  result[2] = v0[0]*v1[1] - v0[1]*v1[0];
  return result;
}



template<typename GEO>
void printGeometry(const GEO &geo, const char *name)
{
  std::cout << name << " :";
  for (int i = 0; i < geo.corners(); ++i)
    STDOUT(" (" << geo.corner(i) << ")");
  std::cout << std::endl;
}



/**
 * @brief computes the local coordinate of a given face's particular corner
 * if the corner was determined using orientedSubface
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
  unsigned int cornerlocal = orientedSubface<dim>(type, face, vertex);

  const Dune::GenericReferenceElement<double,dim>& refElement =
    Dune::GenericReferenceElements<double, dim>::general(type);

  // return the vertex' center of gravity, i.e. the vertex' location
  return refElement.position(cornerlocal, dim);
}


#endif // GEOMETRY_H_
