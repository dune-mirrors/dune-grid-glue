// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * @author Christian Engwer
 * @author Oliver Sander
 * @author Gerrit Buse
 * @brief contains customized geometry implementations for simplices
 * \deprecated For use with Dune 2.2 and older.
 */

#ifndef DUNE_GRIDGLUE_COMMON_SIMPLEXGEOMETRY_HH
#define DUNE_GRIDGLUE_COMMON_SIMPLEXGEOMETRY_HH

#include <memory>

#include <dune/common/static_assert.hh>
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
    template<class ctype, int dim, int dimworld>
    struct SimplexGeometryTraits
      : public Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dimworld, true>
    {
      typedef typename Dune::GenericGeometry::DefaultGeometryTraits<ctype, dim, dimworld, true> Base;

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
    template<class ctype, int mydim, int coorddim>
    class SimplexGeometry : public Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<ctype, mydim, coorddim> >
    {
      typedef Dune::GenericGeometry::BasicGeometry<mydim, SimplexGeometryTraits<ctype, mydim, coorddim> > Base;

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
      void setup(const Dune::GeometryType& type, const std::array<Dune::FieldVector<ctype, coorddim>, simplex_corners>& coordinates)
      {
        // Yes, a strange way, but the only way, as BasicGeometry doesn't have a setup method
        Base::operator=(Base(type, coordinates));
      }

    };

  } // end namespace GridGlue
} // end namespace Dune

#endif // DUNE_GRIDGLUE_COMMON_SIMPLEXGEOMETRY_HH
