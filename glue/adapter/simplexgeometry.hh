// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SimplexGeometry.hh
 *  Version:     1.0
 *  Created on:  Feb 17, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: classes used for local and global geometries of RemoteIntersection implementations
 *  subversion:  $Id$
 *
 */
/**
 * @file SimplexGeometry.hh
 * @brief contains customized geometry implementations that can be used to
 * describe the mappings between simplex element, their faces and simplicial
 * parts inside those faces
 */

#ifndef SIMPLEXGEOMETRY_HH_
#define SIMPLEXGEOMETRY_HH_

#include <dune/common/static_assert.hh>
#include <dune/common/array.hh>
#include <dune/grid/genericgeometry/geometry.hh>
#include <dune/grid/genericgeometry/geometrytraits.hh>
#include <dune/grid/genericgeometry/cornermapping.hh>
#include <dune/grid/genericgeometry/referenceelements.hh>
#include <dune/grid/genericgeometry/misc.hh>

#include "../misc/geometry.hh"


/**
 * @class LocalSimplexGeometryTraits
 * @brief tweaked geometry traits passed to GenericGeometry::BasicGeometry
 *
 * This geometry traits class configures BasicGeometry to use only affine mappings
 * for simplicial (i.e. especially non-hybrid) grid structures.
 * The traits class is suited for the use with full-dimensional structures as well as
 * manifolds and flat hyperdimensional meshes. It chooses the correct mapping from
 * a subface simplex in the surface (or the grid mesh) to it's parent codimension 0
 * entity. Take the 3D case with (dimension, dimensionworld) for example:
 * - surfaces    (3,3): triangle -> tetrahedron
 * - manifolds   (2,3): triangle -> triangle
 * - flat meshes (2,2): triangle -> triangle
 */
template<typename G, bool is_manifold = false>
struct LocalSimplexGeometryTraits
  : public Dune::GenericGeometry::DefaultGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld - static_cast<int>(is_manifold), true>
{
  typedef typename Dune::GenericGeometry::DefaultGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld - static_cast<int>(is_manifold), true> Base;

  /// @brief we only talk about simplices here!
  static const bool hybrid = false;

  /// @brief since non-hybrid, an element type can be specified
  static const Dune::GeometryType::BasicType dunetype = Dune::GeometryType::simplex;

  /**
   * Note:
   * The coorddimension of the mapping is set to the adjusted value of dimWorld!
   * It won't work if you use G::dimensionworld, instead the construct  "G::dimensionworld - static_cast<int>(is_manifold)"
   * has to be used, which is identical with Base::dimWorld.
   * The problem only occurs when working with manifold, hence the name of the template parameter.
   */
  template<class Topology>
  struct Mapping
  {
    typedef Dune::GenericGeometry::CoordStorage<typename Base::CoordTraits, Topology, Base::dimWorld>                        CornerStorage;
    typedef Dune::GenericGeometry::CornerMapping<typename Base::CoordTraits, Topology, Base::dimWorld, CornerStorage, true>  type;
  };

  // caching configuration chosen in base class is fine
};


/**
 * @class SimplexGeometryTraits
 * @brief tweaked geometry traits passed to Dune::GenericGeometry::BasicGeometry
 *
 * This geometry traits class configures BasicGeometry to use only affine mappings
 * for simplicial (i.e. especially non-hybrid) grid structures.
 * The traits class is suited for the use with full-dimensional structures as well as
 * manifolds and flat hyperdimensional meshes. It chooses the correct mapping from
 * a subface simplex in the surface (or the grid mesh) to it's parent codimension 0
 * entity. Take the 3D case with (dimension, dimensionworld) for example:
 * - surfaces    (3,3): triangle -> tetrahedron
 * - manifolds   (2,3): triangle -> tetrahedron
 * - flat meshes (2,2): triangle -> triangle
 */
template<typename G>
struct GlobalSimplexGeometryTraits
  : public Dune::GenericGeometry::DefaultGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld, true>
{
  typedef typename Dune::GenericGeometry::DefaultGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld, true> Base;

  /// @brief we only talk about simplices here!
  static const bool hybrid = false;

  static const Dune::GeometryType::BasicType dunetype = Dune::GeometryType::simplex;

  /**
   * Note:
   * In contrast to LocalSimplexGeometryTraits it doesn't matter which one you use,
   * Base::dimWorld or G::dimensionworld. For reasons of conformance the Traits' class'
   * constant is used.
   */
  template<class Topology>
  struct Mapping
  {
    typedef Dune::GenericGeometry::CoordStorage<typename Base::CoordTraits, Topology, Base::dimWorld>                        CornerStorage;
    typedef Dune::GenericGeometry::CornerMapping<typename Base::CoordTraits, Topology, Base::dimWorld, CornerStorage, true>  type;
  };

  // caching configuration chosen in base class is fine
};





/**
 * @class SimplexGeometry
 * @brief This class is derived from BasicGeometry using tuned geometry traits.
 *
 * The used geometry traits allow for a more efficient compilation of the default
 * implementation of CornerMapping. Using only affine mappings it is specialized
 * for the case of exclusively simplicial geometries.
 */
template<int mydim, int coorddim, typename G>
class SimplexGeometry : public Dune::GenericGeometry::BasicGeometry<mydim, GlobalSimplexGeometryTraits<G> >
{
  typedef SimplexGeometry<mydim, coorddim, G> This;

  typedef Dune::GenericGeometry::BasicGeometry<mydim, GlobalSimplexGeometryTraits<G> > Base;

  enum { simplex_corners = coorddim + static_cast<int>(mydim == coorddim) };

protected:

  /// @brief normal computation module specialized for the case of full-dimensional or
  /// manifold grids
  template<bool>
  struct CommonNormalComputer
  {
    static Dune::FieldVector<typename G::ctype, coorddim> outerNormal(const This& geom, const Dune::FieldVector<typename G::ctype, mydim>& local)
    {
      Dune::FieldVector<typename G::ctype, coorddim> result = geom.corner(1);
      if (coorddim == 2)
      {
        result -= geom.corner(0);
        typename G::ctype temp = result[0];
        result[0] = result[1];
        result[1] = -temp;
        return result;
      }
      else if (coorddim == 3)
      {
        Dune::FieldVector<typename G::ctype, coorddim> v0 = geom.corner(2), v1 = geom.corner(0);
        // tris edge vectors (pos. oriented)
        v0 -= geom.corner(1);
        v1 -= geom.corner(1);
        // cross product
        result[0] = v0[1]*v1[2] - v0[2]*v1[1];
        result[1] = v0[2]*v1[0] - v0[0]*v1[2];
        result[2] = v0[0]*v1[1] - v0[1]*v1[0];
      }
      else
        DUNE_THROW(Dune::NotImplemented, "dimension not implemented!");
      return result;
    }
  };


  /// @brief for normal computation for hyperdimensional grids
  /// Note: This degenerates to determining the simplex' orientation.
  /// The normal then always points to the "added" dimension.
  template<bool>
  struct LiftingNormalComputer
  {
    static Dune::FieldVector<typename G::ctype, coorddim + 1> outerNormal(const SimplexGeometry<mydim, coorddim, G>& geom, const Dune::FieldVector<typename G::ctype, mydim>& local)
    {
      Dune::FieldVector<typename G::ctype, coorddim + 1> result(0.0);
      if (coorddim == 1)
      {
        result[1] = geom.corner(0)[0] - geom.corner(1)[0];
        return result;
      }
      else if (coorddim == 2)
      {
        Dune::FieldVector<typename G::ctype, coorddim> v0 = geom.corner(2), v1 = geom.corner(0);
        // tris edge vectors (pos. oriented)
        v0 -= geom.corner(1);
        v1 -= geom.corner(1);
        // cross product
        result[2] = v0[0]*v1[1] - v0[1]*v1[0];
        return result;
      }
      else
        DUNE_THROW(Dune::NotImplemented, "dimension not implemented");
    }
  };

  typedef Dune::GenericGeometry::ProtectedIf<mydim == coorddim, LiftingNormalComputer, CommonNormalComputer> NormalComputer;

public:

  typedef typename Base::Mapping Mapping;

  template< class CoordVector>
  SimplexGeometry(const Dune::GeometryType &type, const CoordVector &coords) : Base(type, coords)
  {}


  SimplexGeometry()
  {}


  //		SimplexGeometry& operator=(const SimplexGeometry& geom)
  //		{
  //			Base::operator=(Base(geom));
  //			return *this;
  //		}

  /**
   * @brief Setup method with a geometry type and a set of corners
   * @param type the geometry type of this subface, i.e. most likely a simplex in 1D or 2D
   * @param coordinates The corner coordinates in DUNE numbering
   */
  void setup(const Dune::GeometryType& type, const Dune::array<Dune::FieldVector<typename G::ctype, coorddim>, simplex_corners>& coordinates)
  {
    // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
    Base::operator=(Base(type, coordinates));
  }


  /**
   * @brief provided only for reasons of forward compatibility
   *
   * BasicGeometry apparently does not support the newly defined
   * corner access method BasicGeometry::corner(int i) yet.
   * @param i the index of the corner (not checked!)
   * @return a copy of the coordinate at corner
   */
  const typename Base::GlobalCoordinate corner(int i) const
  {
    return this->operator[](i);
  }


  const Dune::FieldVector<typename G::ctype, coorddim + static_cast<int>(coorddim==mydim)> outerNormal(const Dune::FieldVector<typename G::ctype, mydim>& local) const
  {
    return This::NormalComputer::outerNormal(*this, local);
  }
};



/**
 * @class LocalSimplexGeometry
 * @brief This class is derived from BasicGeometry using tuned geometry traits.
 *
 * The used geometry traits allow for a more efficient compilation of the default
 * implementation of CornerMapping. Using only affine mappings it is specialized
 * for the case of exclusively simplicial geometries.
 * Note:
 * This class mostly does the same as the SimplexGeometry class. Behaviour only differs
 * in the case of manifolds. Then the local mapping is 1D->1D resp. 2D->2D here whereas
 * the global geometry maps to the world coordinate space, i.e. 1D->2D resp. 2D->3D.
 */
template<int mydim, int coorddim, typename G>
class LocalSimplexGeometry
  : public Dune::GenericGeometry::BasicGeometry<mydim, LocalSimplexGeometryTraits<G, (coorddim < static_cast<int>(G::dimensionworld))> >
{
  typedef Dune::GenericGeometry::BasicGeometry<mydim, LocalSimplexGeometryTraits<G, (coorddim < static_cast<int>(G::dimensionworld))> > Base;

  enum { simplex_corners = coorddim + static_cast<int>(mydim == coorddim) };

public:

  typedef typename Base::Mapping Mapping;


  template< class CoordVector>
  LocalSimplexGeometry(const Dune::GeometryType &type, const CoordVector &coords) : Base(type, coords)
  {}


  LocalSimplexGeometry()
  {}


  //	SimplexGeometry& operator=(const SimplexGeometry& geom)
  //	{
  //		Base::operator=(Base(geom));
  //		return *this;
  //	}

  /**
   * @brief Setup method with a geometry type and a set of corners
   * @param type the geometry type of this subface, i.e. most likely a simplex in 1D or 2D
   * @param coordinates The corner coordinates in DUNE numbering
   */
  void setup(const Dune::GeometryType& type, const Dune::array<Dune::FieldVector<typename G::ctype, coorddim>, simplex_corners>& coordinates)
  {
    //		STDOUTLN("LocalSimplexGeometry: trying to set up");
    // set up base class
    // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
    Base::operator=(Base(type, coordinates));
    //		for (int i = 0; i < this->corners(); ++i)
    //			STDOUT(" (" << this->operator[](i) << ")");
    //		STDOUTLN("");
  }


  /**
   * @brief provided only for reasons of forward compatibility
   *
   * BasicGeometry apparently does not support the newly defined
   * corner access method BasicGeometry::corner(int i) yet.
   * @param i the index of the corner (not checked!)
   * @return a copy of the coordinate at corner
   */
  const typename Base::GlobalCoordinate corner(int i) const
  {
    return this->operator[](i);
  }

};


#endif // SIMPLEXGEOMETRY_HH_
