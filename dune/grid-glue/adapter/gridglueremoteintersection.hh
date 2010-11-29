// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlueRemoteIntersection.hh
 *  Version:     1.0
 *  Created on:  Mar 2, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Model of the RemoteIntersection concept provided by GridGlue.
 *  subversion:  $Id$
 *
 */
/**
 * @file GridGlueRemoteIntersection.hh
 * @brief Model of the RemoteIntersection concept provided by GridGlue.
 */

#ifndef DUNE_GRIDGLUE_REMOTEINTERSECTION_HH
#define DUNE_GRIDGLUE_REMOTEINTERSECTION_HH

#include <dune/grid-glue/adapter/simplexgeometry.hh>
#include <dune/grid-glue/adapter/gridglue.hh>

/*   I M P L E M E N T A T I O N   O F   S U B C L A S S   REMOTE INTERSECTION IMPL   */

namespace Dune {
  namespace GridGlue {

    template<typename P0, typename P1>
    class AbsoluteIntersection
    {

    private:

      typedef ::GridGlue<P0, P1> Parent;

      friend class Parent::RemoteIntersectionIterator;
      friend class Parent::Grid0IntersectionIterator;
      friend class Parent::Grid1IntersectionIterator;

    public:

      typedef typename Parent::Grid0View Grid0View;
      typedef typename Parent::Grid0Patch Grid0Patch;
      typedef typename Grid0View::Grid Grid0;

      typedef typename Parent::Grid1View Grid1View;
      typedef typename Parent::Grid1Patch Grid1Patch;
      typedef typename Grid1View::Grid Grid1;

      typedef typename Grid0View::Traits::template Codim<0>::Entity Grid0Entity;
      typedef typename Grid0View::Traits::template Codim<0>::EntityPointer Grid0EntityPointer;

      typedef typename Grid1View::Traits::template Codim<0>::Entity Grid1Entity;
      typedef typename Grid1View::Traits::template Codim<0>::EntityPointer Grid1EntityPointer;


      dune_static_assert(Grid0::dimension - Grid0Patch::codim
                         == Grid1::dimension - Grid1Patch::codim,
                         "Currently both coupling extracts need to have the same dimension!");

      /** \brief Dimension of the world space of the intersection */
      enum { coorddim = Parent::dimworld };

      /** \brief Dimension of the intersection */
      enum { mydim = Parent::Grid0View::Grid::dimension - Grid0Patch::codim };


      typedef LocalSimplexGeometry<mydim, Grid0View::dimension, Grid0> Grid0LocalGeometry;

      typedef SimplexGeometry<mydim, Grid0::dimensionworld, Grid0> Grid0Geometry;

      typedef LocalSimplexGeometry<mydim, Grid1View::dimension, Grid1> Grid1LocalGeometry;

      typedef SimplexGeometry<mydim, Grid1::dimensionworld, Grid1> Grid1Geometry;


      typedef typename Grid0::ctype ctype;


      typedef Dune::FieldVector<ctype, mydim>      LocalCoordinate;

      typedef Dune::FieldVector<ctype, coorddim>   GlobalCoordinate;

      /*   C O N S T R U C T O R S   */

      /** \brief Constructor for NULL_INTERSECTION */
      AbsoluteIntersection(const Parent* glue)
        : glue_(glue), mergeindex_((IndexType)-1), index_((IndexType)-1),
          domainindex_((IndexType)-1), targetindex_((IndexType)-1)
      {}

      /** \brief Constructor the n'th AbsoluteIntersection of a given GridGlue */
      AbsoluteIntersection(const Parent* glue, unsigned int mergeindex);

      /** \brief Copy construction from another AbsoluteIntersection */
      AbsoluteIntersection(const AbsoluteIntersection& impl)
      {
        *this = impl;
      }

      /** \brief Assignment operator */
      AbsoluteIntersection& operator=(const AbsoluteIntersection& impl)
      {
        glue_ = impl.glue_;
        mergeindex_ = impl.mergeindex_;
        index_ = impl.index_;
        domainindex_ = impl.domainindex_;
        targetindex_ = impl.targetindex_;
        domlgeom_ = impl.domlgeom_;
        domggeom_ = impl.domggeom_;
        tarlgeom_ = impl.tarlgeom_;
        targgeom_ = impl.targgeom_;
        return *this;
      }

      /*   F U N C T I O N A L I T Y   */

      /** \brief return EntityPointer to the Entity on the inside of this intersection.
              That is the Entity where we started this. */
      Grid0EntityPointer entityGrid0() const
      {
        return glue_->grid0Patch().element(glue_->merger_->template parent<0>(mergeindex_));
      }


      /** \brief return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity. */
      Grid1EntityPointer entityGrid1() const
      {
        return glue_->grid1Patch().element(glue_->merger_->template parent<1>(mergeindex_));
      }


      /** \brief Return true if intersection is conforming */
      bool conforming() const;


      /** \brief geometrical information about this intersection in local coordinates of the inside() entity.
          takes local domain intersection coords and maps them to domain parent element local coords */
      const Grid0LocalGeometry& geometryInGrid0Entity() const
      {
        return domlgeom_;
      }


      /** \brief geometrical information about this intersection in local coordinates of the outside() entity.
          takes local target intersection coords and maps them to target parent element local coords */
      const Grid1LocalGeometry& geometryInGrid1Entity() const
      {
        return tarlgeom_;
      }


      /** \brief geometrical information about this intersection in global coordinates in the domain grid.
          takes local domain intersection coords and maps them to domain grid world coords */
      const Grid0Geometry& geometryGrid0() const
      {
        return domggeom_;
      }


      /** \brief geometrical information about this intersection in global coordinates in the target grid. */
      const Grid1Geometry& geometryGrid1() const
      {
        return targgeom_;
      }

      /** \brief obtain the type of reference element for this intersection */
      Dune::GeometryType type() const
      {
        return Dune::GeometryType(Dune::GeometryType::simplex, mydim);
      }

      bool hasGrid0() const
      {
        unsigned int localindex;
        return glue_->grid0Patch().contains(glue_->merger_->template parent<0>(mergeindex_), localindex);
      }

      bool hasGrid1() const
      {
        unsigned int localindex;
        return glue_->grid1Patch().contains(glue_->merger_->template parent<1>(mergeindex_), localindex);
      }

      // Local number of codim 1 entity in the inside() Entity where intersection is contained in.
      int indexInGrid0Entity() const
      {
        return glue_->grid0Patch().indexInInside(glue_->merger_->template parent<0>(mergeindex_));
      }


      // Local number of codim 1 entity in outside() Entity where intersection is contained in.
      int indexInGrid1Entity() const
      {
        return glue_->grid1Patch().indexInInside(glue_->merger_->template parent<1>(mergeindex_));
      }


      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate outerNormalGrid0(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return domggeom_.outerNormal(local);
      }

      /** \brief Return an outer normal */
      GlobalCoordinate unitOuterNormalGrid0(const Dune::FieldVector<ctype, mydim> &local) const
      {
        Dune::FieldVector<ctype, coorddim> normal = outerNormalGrid0(local);
        normal /= normal.two_norm();
        return normal;
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate integrationOuterNormalGrid0(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return (unitOuterNormalGrid0(local) *= geometryGrid0().integrationElement(local));
      }

      GlobalCoordinate centerUnitOuterNormalGrid0() const
      {
        return domggeom_.centerUnitOuterNormal();
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate outerNormalGrid1(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return targgeom_.outerNormal(local);
      }

      /** \brief Return a unit outer normal of the target intersection */
      GlobalCoordinate unitOuterNormalGrid1(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return targgeom_.unitOuterNormal(local);
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate integrationOuterNormalGrid1(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return (unitOuterNormalGrid1(local) *= geometryGrid1().integrationElement(local));
      }

      GlobalCoordinate centerUnitOuterNormalGrid1() const
      {
        return targgeom_.centerUnitOuterNormal();
      }

#ifdef QUICKHACK_INDEX
      typedef unsigned int IndexType;

      IndexType & index()
      {
        return index_;
      }

      IndexType & domainIndex()
      {
        return domainindex_;
      }

      IndexType & targetIndex()
      {
        return targetindex_;
      }

      IndexType index() const
      {
        return index_;
      }

      IndexType globalIndex() const
      {
        return mergeindex_;
      }

      IndexType domainIndex() const
      {
        assert(domainindex_ != (IndexType)-1);
        return domainindex_;
      }

      IndexType targetIndex() const
      {
        assert(targetindex_ != (IndexType)-1);
        return targetindex_;
      }

#endif

    private:

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief the grid glue entity this is built on
      const Parent*       glue_;

      /// @brief index of this intersection after GridGlue interface
      IndexType mergeindex_;

      IndexType index_;
      IndexType domainindex_;
      IndexType targetindex_;

      Grid0LocalGeometry domlgeom_;

      Grid0Geometry domggeom_;

      Grid1LocalGeometry tarlgeom_;

      Grid1Geometry targgeom_;

    };

    template<typename P0, typename P1>
    bool AbsoluteIntersection<P0, P1>::conforming() const
    {
      std::vector<unsigned int> results;
      // first check the domain side
      bool is_conforming =
        glue_->merger_->template simplexRefined<0>(glue_->merger_->template parent<0>(mergeindex_), results) && results.size() == 1;
      results.resize(0);
      // now check the target side
      if (is_conforming)
        return glue_->merger_->template simplexRefined<1>(glue_->merger_->template parent<1>(mergeindex_), results) && results.size() == 1;
      return false;
    }


    template<typename P0, typename P1>
    AbsoluteIntersection<P0, P1>::AbsoluteIntersection(const Parent* glue, unsigned int mergeindex)
      : glue_(glue), mergeindex_(mergeindex), index_((IndexType)-1),
        domainindex_((IndexType)-1), targetindex_((IndexType)-1)
    {
      // if an invalid index is given do not proceed!
      // (happens when the parent GridGlue initializes the "end"-Intersection)
      assert (0 <= mergeindex || mergeindex < glue->index__sz);

      // initialize the local and the global geometry of the domain
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = Grid0::template Codim<0>::Geometry::mydimension;

        const int nSimplexCorners = elementdim - Parent::Grid0Patch::codim + 1;

        // coordinates within the subentity that contains the remote intersection
        Dune::array<LocalCoordinate, nSimplexCorners> corners_subEntity_local;

        for (int i = 0; i < nSimplexCorners; ++i)
          corners_subEntity_local[i] = glue->merger_->template parentLocal<0>(mergeindex, i);

        // Coordinates of the remote intersection corners wrt the element coordinate system
        Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

        // world coordinates of the remote intersection corners
        Dune::array<Dune::FieldVector<ctype, Grid0::dimensionworld>, nSimplexCorners> corners_global;

        unsigned int domainIndex = glue->merger_->template parent<0>(mergeindex);
        unsigned int unused;
        if (glue->grid0Patch().contains(domainIndex, unused))
        {
          typename Grid0Patch::Geometry domainWorldGeometry = glue->grid0Patch().geometry(domainIndex);
          typename Grid0Patch::LocalGeometry domainLocalGeometry = glue->grid0Patch().geometryLocal(domainIndex);

          for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
            corners_element_local[i] = domainLocalGeometry.global(corners_subEntity_local[i]);
            corners_global[i]        = domainWorldGeometry.global(corners_subEntity_local[i]);
          }

          // set the corners of the geometries
          domlgeom_.setup(type(), corners_element_local);
          domggeom_.setup(type(), corners_global);
        }
      }

      // do the same for the local and the global geometry of the target
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = Grid1::template Codim<0>::Geometry::mydimension;

        const int nSimplexCorners = elementdim - Parent::Grid1Patch::codim + 1;

        // coordinates within the subentity that contains the remote intersection
        Dune::array<LocalCoordinate, nSimplexCorners> corners_subEntity_local;

        for (int i = 0; i < nSimplexCorners; ++i)
          corners_subEntity_local[i] = glue->merger_->template parentLocal<1>(mergeindex, i);

        // Coordinates of the remote intersection corners wrt the element coordinate system
        Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

        // world coordinates of the remote intersection corners
        Dune::array<Dune::FieldVector<ctype, Grid1::dimensionworld>, nSimplexCorners> corners_global;

        unsigned int targetIndex = glue->merger_->template parent<1>(mergeindex);
        unsigned int unused;
        if (glue->grid1Patch().contains(targetIndex, unused))
        {
          typename Grid1Patch::Geometry targetWorldGeometry = glue->grid1Patch().geometry(targetIndex);
          typename Grid1Patch::LocalGeometry targetLocalGeometry = glue->grid1Patch().geometryLocal(targetIndex);

          for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
            corners_element_local[i] = targetLocalGeometry.global(corners_subEntity_local[i]);
            corners_global[i]        = targetWorldGeometry.global(corners_subEntity_local[i]);
          }

          // set the corners of the geometries
          tarlgeom_.setup(type(), corners_element_local);
          targgeom_.setup(type(), corners_global);
        }
      }
    }

    /*   I M P L E M E N T A T I O N   O F   S U B C L A S S   DIRECTED REMOTE INTERSECTION IMPL   */

    template<typename P0, typename P1, int inside, int outside>
    class DirectedIntersection;

    template<typename P0, typename P1>
    class DirectedIntersection<P0, P1, 0,1>
    {
    private:

      typedef Dune::GridGlue::AbsoluteIntersection<P0, P1> AbsoluteIntersection;

    public:

      typedef typename AbsoluteIntersection::Grid0View InsideGridView;

      typedef typename AbsoluteIntersection::GridView1 OutsideGridView;

      typedef typename InsideGridView::Traits::template Codim<0>::Entity InsideEntity;

      typedef typename InsideGridView::Traits::template Codim<0>::EntityPointer InsideEntityPointer;

      typedef typename OutsideGridView::Traits::template Codim<0>::Entity OutsideEntity;

      typedef typename OutsideGridView::Traits::template Codim<0>::EntityPointer OutsideEntityPointer;


      enum {
        /** \brief Dimension of the world space of the intersection */
        coorddim = AbsoluteIntersection::dimworld,
        /** \brief Dimension of the intersection */
        mydim = AbsoluteIntersection::mydim,
        /** \brief document the inside & outside domain */
        //! @{
        insidePatch = 0,
        outsidePatch = 1
                       //! @}
      };

      typedef typename AbsoluteIntersection::Grid0LocalGeometry InsideLocalGeometry;
      typedef typename AbsoluteIntersection::Grid1LocalGeometry OutsideLocalGeometry;
      typedef typename AbsoluteIntersection::Grid0Geometry Geometry;

      typedef typename AbsoluteIntersection::ctype ctype;

      typedef typename AbsoluteIntersection::LocalCoordinate LocalCoordinate;
      typedef typename AbsoluteIntersection::GlobalCoordinate GlobalCoordinate;

      // typedef unsigned int IndexType;

      /*   C O N S T R U C T O R S   */

      /** \brief Constructor for NULL_INTERSECTION */
      DirectedIntersection(const AbsoluteIntersection & i) : i_(i) {}

      /*   F U N C T I O N A L I T Y   */

      /** \brief return EntityPointer to the Entity on the inside of this intersection.
              That is the Entity where we started this. */
      InsideEntityPointer inside() const
      {
        return i_.entityGrid0();
      }

      /** \brief return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity. */
      OutsideEntityPointer outside() const
      {
        return i_.entityTarget();
      }

      /** \brief Return true if intersection is conforming */
      bool conforming() const
      {
        return i_.conforming();
      }

      /** \brief geometrical information about this intersection in local coordinates of the inside() entity.
          takes local domain intersection coords and maps them to domain parent element local coords */
      const InsideLocalGeometry& geometryInInside() const
      {
        return i_.geometryInGrid0Entity();
      }

      /** \brief geometrical information about this intersection in local coordinates of the outside() entity.
          takes local target intersection coords and maps them to target parent element local coords */
      const OutsideLocalGeometry& geometryInOutside() const
      {
        return i_.geometryInGrid1Entity();
      }

      /** \brief geometrical information about this intersection in global coordinates in the inside grid.
          takes local domain intersection coords and maps them to domain grid world coords */
      const Geometry& intersectionGrid0Global() const
      {
        return i_.geometryGrid0();
      }

      /** \brief obtain the type of reference element for this intersection */
      Dune::GeometryType type() const
      {
        return i_.type();
      }


      /** \brief return true if inside() entity exists locally */
      bool self() const
      {
        return i_.hasGrid0();
      }

      /** \brief return true if outside() entity exists locally */
      bool neighbor() const
      {
        return i_.hasGrid1();
      }

      /** \brief Local number of codim 1 entity in the inside() Entity where intersection is contained in. */
      int indexInInside() const
      {
        return i_.indexInGrid0Entity();
      }

      /** \brief Local number of codim 1 entity in outside() Entity where intersection is contained in. */
      int indexInOutside() const
      {
        return i_.indexInGrid1Entity();
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate outerNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return i_.outerNormalGrid0(local);
      }

      /** \brief Return a unit outer normal of the target intersection */
      GlobalCoordinate unitOuterNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return i_.unitOuterNormalGrid0(local);
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate integrationOuterNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return i_.integrationOuterNormalGrid0(local);
      }

      GlobalCoordinate centerUnitOuterNormal () const
      {
        return i_.centerUnitOuterNormalGrid0();
      }

#ifdef QUICKHACK_INDEX
      typedef typename AbsoluteIntersection::IndexType IndexType;

      IndexType index() const
      {
        return i_.index();
      }

      IndexType globalIndex() const
      {
        return i_.globalIndex();
      }

      IndexType domainIndex() const
      {
        return i_.domainIndex();
      }

      IndexType targetIndex() const
      {
        return i_.targetIndex();
      }

#endif

    private:

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief the grid glue entity this is built on
      const Parent*       glue_;

      /// @brief the underlying remote intersection
      const IntersectionData* i_;
    };


  } // end namesapce GridGlue
} // end namesapce Dune

#endif // DUNE_GRIDGLUE_REMOTEINTERSECTION_HH
