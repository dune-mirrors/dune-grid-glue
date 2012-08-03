// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * @author Christian Engwer
 * @author Oliver Sander
 * @author Gerrit Buse
 * @brief contains customized geometry implementations for simplices
 */

#ifndef SIMPLEXGEOMETRY_HH
#define SIMPLEXGEOMETRY_HH

#include <dune/common/static_assert.hh>
#include <dune/common/array.hh>
#include <dune/common/version.hh>
#include <dune/geometry/genericgeometry/geometry.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/cornermapping.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>

namespace Dune {
  namespace GridGlue {

    /**
     * @brief Geometry traits for simplices passed to GenericGeometry::BasicGeometry
     *
     * This geometry traits class configures BasicGeometry to use only affine mappings
     * for simplicial (i.e. especially non-hybrid) grid structures.
     */
    template<class ctype, int dim, int dimworld, bool is_manifold = false>
    struct SimplexGeometryTraits
      : public Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dimworld - static_cast<int>(is_manifold), true>
    {
      typedef typename Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dimworld - static_cast<int>(is_manifold), true> Base;

      // This traits class represents a single type only ...
      static const bool hybrid = false;
      // ... and this type is 'simplex'.
      static const unsigned int topologyId = Dune::GenericGeometry::SimplexTopology< dim >::type::id;

      template<class Topology>
      struct Mapping
      {
        typedef Dune::GenericGeometry::CoordStorage<typename Base::CoordTraits, Topology, Base::dimWorld>                        CornerStorage;
        typedef Dune::GenericGeometry::CornerMapping<typename Base::CoordTraits, Topology, Base::dimWorld, CornerStorage, true>  type;
      };

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
    class SimplexGeometry : public Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld> >
    {
      typedef Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld> > Base;

      enum { simplex_corners = mydim+1 };

    public:

      typedef typename Base::Mapping Mapping;

      template< class CoordVector>
      SimplexGeometry(const Dune::GeometryType &type, const CoordVector &coords) : Base(type, coords)
      {}


      SimplexGeometry()
      {}

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
      : public Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld, (coorddim < static_cast<int>(G::dimensionworld))> >
    {
      typedef Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<typename G::ctype, G::dimension, G::dimensionworld, (coorddim < static_cast<int>(G::dimensionworld))> > Base;

      enum { simplex_corners = mydim+1 };

    public:

      typedef typename Base::Mapping Mapping;


      template< class CoordVector>
      LocalSimplexGeometry(const Dune::GeometryType &type, const CoordVector &coords) : Base(type, coords)
      {}


      LocalSimplexGeometry()
      {}


      /**
       * @brief Setup method with a geometry type and a set of corners
       * @param type the geometry type of this subface, i.e. most likely a simplex in 1D or 2D
       * @param coordinates The corner coordinates in DUNE numbering
       */
      void setup(const Dune::GeometryType& type, const Dune::array<Dune::FieldVector<typename G::ctype, coorddim>, simplex_corners>& coordinates)
      {
        // set up base class
        // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
        Base::operator=(Base(type, coordinates));
      }

    };

  } // end namespace GridGlue
} // end namespace Dune

#endif // SIMPLEXGEOMETRY_HH
