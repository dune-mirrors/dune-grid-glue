// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file
   @author Christian Engwer
   @brief Model of the Intersection concept provided by GridGlue.
 */

#ifndef DUNE_GRIDGLUE_REMOTEINTERSECTION_HH
#define DUNE_GRIDGLUE_REMOTEINTERSECTION_HH

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
#include <dune/geometry/affinegeometry.hh>
#else
#include <dune/grid-glue/common/simplexgeometry.hh>
#endif
#include <dune/grid-glue/adapter/gridglue.hh>

#define ONLY_SIMPLEX_INTERSECTIONS

namespace Dune {
  namespace GridGlue {

    /**
       @brief storage class for Dune::GridGlue::Intersection related data
     */
    template<typename P0, typename P1>
    struct IntersectionData
    {
      typedef ::GridGlue<P0, P1> GridGlue;

      typedef typename GridGlue::IndexType IndexType;

      /** \brief Dimension of the world space of the intersection */
      enum { coorddim = GridGlue::dimworld };

    private:
      // intermediate quantities
      static const int dim1 = GridGlue::Grid0View::Grid::dimension - GridGlue::Grid0Patch::codim;
      static const int dim2 = GridGlue::Grid1View::Grid::dimension - GridGlue::Grid1Patch::codim;

    public:
      /** \brief Dimension of the intersection */
      enum { mydim = (dim1<dim2) ? dim1 : dim2 };

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
      typedef AffineGeometry<typename GridGlue::Grid0View::ctype, mydim, GridGlue::Grid0View::dimension>
      Grid0LocalGeometry;
      typedef AffineGeometry<typename GridGlue::Grid0View::ctype, mydim, GridGlue::Grid0View::dimensionworld>
      Grid0Geometry;
      typedef AffineGeometry<typename GridGlue::Grid1View::ctype, mydim, GridGlue::Grid1View::dimension>
      Grid1LocalGeometry;
      typedef AffineGeometry<typename GridGlue::Grid1View::ctype, mydim, GridGlue::Grid1View::dimensionworld>
      Grid1Geometry;
#else
      typedef SimplexGeometry<typename GridGlue::Grid0View::ctype, mydim, GridGlue::Grid0View::dimension>
      Grid0LocalGeometry;
      typedef SimplexGeometry<typename GridGlue::Grid0View::ctype, mydim, GridGlue::Grid0View::dimensionworld>
      Grid0Geometry;
      typedef SimplexGeometry<typename GridGlue::Grid1View::ctype, mydim, GridGlue::Grid1View::dimension>
      Grid1LocalGeometry;
      typedef SimplexGeometry<typename GridGlue::Grid1View::ctype, mydim, GridGlue::Grid1View::dimensionworld>
      Grid1Geometry;
#endif

      typedef typename GridGlue::Grid0View::IndexSet::IndexType Grid0IndexType;
      typedef typename GridGlue::Grid1View::IndexSet::IndexType Grid1IndexType;

      /** \brief Constructor the n'th IntersectionData of a given GridGlue */
      IntersectionData(const GridGlue& glue, unsigned int mergeindex, unsigned int offset, bool grid0local, bool grid1local);

      /** \brief Default Constructor */
      IntersectionData() : grid0local_(false), grid1local_(false) {}

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief index of this intersection after GridGlue interface
      IndexType index_;

      bool grid0local_;              //!< true if the associated grid0 entity is local
      Grid0IndexType grid0index_;    //!< index of the associated local grid0 entity
      bool grid1local_;              //!< true if the associated grid1 entity is local
      Grid1IndexType grid1index_;    //!< index of the associated local grid1 entity

      shared_ptr<Grid0LocalGeometry>  grid0localgeom_;
      shared_ptr<Grid0Geometry>       grid0geom_;
      shared_ptr<Grid1LocalGeometry>  grid1localgeom_;
      shared_ptr<Grid1Geometry>       grid1geom_;

    };

    //! \todo move this functionality to GridGlue
    template<typename P0, typename P1>
    IntersectionData<P0, P1>::IntersectionData(const GridGlue& glue, unsigned int mergeindex, unsigned int offset,
                                               bool grid0local, bool grid1local)
      : index_(mergeindex+offset),
        grid0local_(grid0local),
        grid0index_(0),
        grid1local_(grid1local),
        grid1index_(0)
    {
      typedef typename GridGlue::ctype ctype;
      typedef Dune::FieldVector<ctype, coorddim>   GlobalCoordinate;

      // Number of corners of the intersection
      const int nSimplexCorners = mydim + 1;

      // if an invalid index is given do not proceed!
      // (happens when the parent GridGlue initializes the "end"-Intersection)
      assert (0 <= mergeindex || mergeindex < glue.index__sz);

      // initialize the local and the global geometry of the domain
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = GridGlue::Grid0View::template Codim<0>::Geometry::mydimension;

        // coordinates within the subentity that contains the remote intersection
        Dune::array<Dune::FieldVector<ctype, dim1>, nSimplexCorners> corners_subEntity_local;

        for (int i = 0; i < nSimplexCorners; ++i)
          corners_subEntity_local[i] = glue.merger_->template parentLocal<0>(mergeindex, i);

        // Coordinates of the remote intersection corners wrt the element coordinate system
        Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

        // world coordinates of the remote intersection corners
        Dune::array<Dune::FieldVector<ctype, GridGlue::Grid0View::dimensionworld>, nSimplexCorners> corners_global;

        if (grid0local)
        {
          grid0index_ = glue.merger_->template parent<0>(mergeindex);
          typename GridGlue::Grid0Patch::Geometry
          domainWorldGeometry = glue.template patch<0>().geometry(grid0index_);
          typename GridGlue::Grid0Patch::LocalGeometry
          domainLocalGeometry = glue.template patch<0>().geometryLocal(grid0index_);

          for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
            corners_element_local[i] = domainLocalGeometry.global(corners_subEntity_local[i]);
            corners_global[i]        = domainWorldGeometry.global(corners_subEntity_local[i]);
          }

          // set the corners of the geometries
#ifdef ONLY_SIMPLEX_INTERSECTIONS
          Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
#else
#error Not Implemented
#endif
          grid0localgeom_ = make_shared<Grid0LocalGeometry>(type, corners_element_local);
          grid0geom_      = make_shared<Grid0Geometry>(type, corners_global);
        }
      }

      // do the same for the local and the global geometry of the target
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = GridGlue::Grid1View::template Codim<0>::Geometry::mydimension;

        // coordinates within the subentity that contains the remote intersection
        Dune::array<Dune::FieldVector<ctype, dim2>, nSimplexCorners> corners_subEntity_local;

        for (int i = 0; i < nSimplexCorners; ++i)
          corners_subEntity_local[i] = glue.merger_->template parentLocal<1>(mergeindex, i);

        // Coordinates of the remote intersection corners wrt the element coordinate system
        Dune::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

        // world coordinates of the remote intersection corners
        Dune::array<Dune::FieldVector<ctype, GridGlue::Grid1View::dimensionworld>, nSimplexCorners> corners_global;

        if (grid1local)
        {
          grid1index_ = glue.merger_->template parent<1>(mergeindex);
          typename GridGlue::Grid1Patch::Geometry
          targetWorldGeometry = glue.template patch<1>().geometry(grid1index_);
          typename GridGlue::Grid1Patch::LocalGeometry
          targetLocalGeometry = glue.template patch<1>().geometryLocal(grid1index_);

          for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
            corners_element_local[i] = targetLocalGeometry.global(corners_subEntity_local[i]);
            corners_global[i]        = targetWorldGeometry.global(corners_subEntity_local[i]);
          }

          // set the corners of the geometries
#ifdef ONLY_SIMPLEX_INTERSECTIONS
          Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
#else
#error Not Implemented
#endif
          grid1localgeom_ = make_shared<Grid1LocalGeometry>(type, corners_element_local);
          grid1geom_      = make_shared<Grid1Geometry>(type, corners_global);
        }
      }
    }

    /**
       @brief
       @todo doc me
     */
    template<typename P0, typename P1, int P>
    struct IntersectionDataView;

    template<typename P0, typename P1>
    struct IntersectionDataView<P0, P1, 0>
    {
      typedef const typename IntersectionData<P0,P1>::Grid0LocalGeometry LocalGeometry;
      typedef const typename IntersectionData<P0,P1>::Grid0Geometry Geometry;
      typedef const typename IntersectionData<P0,P1>::Grid0IndexType IndexType;
      static LocalGeometry& localGeometry(const IntersectionData<P0,P1> & i)
      {
        return *i.grid0localgeom_;
      }
      static Geometry& geometry(const IntersectionData<P0,P1> & i)
      {
        return *i.grid0geom_;
      }
      static bool local(const IntersectionData<P0,P1> & i)
      {
        return i.grid0local_;
      }
      static IndexType index(const IntersectionData<P0,P1> & i)
      {
        return i.grid0index_;
      }
    };

    template<typename P0, typename P1>
    struct IntersectionDataView<P0, P1, 1>
    {
      typedef const typename IntersectionData<P0,P1>::Grid1LocalGeometry LocalGeometry;
      typedef const typename IntersectionData<P0,P1>::Grid1Geometry Geometry;
      typedef const typename IntersectionData<P0,P1>::Grid1IndexType IndexType;
      static LocalGeometry& localGeometry(const IntersectionData<P0,P1> & i)
      {
        return *i.grid1localgeom_;
      }
      static Geometry& geometry(const IntersectionData<P0,P1> & i)
      {
        return *i.grid1geom_;
      }
      static IndexType local(const IntersectionData<P0,P1> & i)
      {
        return i.grid1local_;
      }
      static IndexType index(const IntersectionData<P0,P1> & i)
      {
        return i.grid1index_;
      }
    };

    /**
       @brief
       @todo doc me
     */
    template<typename P0, typename P1, int inside, int outside>
    struct IntersectionTraits;

    template<typename P0, typename P1>
    struct IntersectionTraits<P0,P1,0,1>
    {
      typedef ::GridGlue<P0, P1> GridGlue;
      typedef Dune::GridGlue::IntersectionData<P0,P1> IntersectionData;

      typedef typename GridGlue::Grid0View InsideGridView;
      typedef typename GridGlue::Grid1View OutsideGridView;

      typedef const typename IntersectionData::Grid0LocalGeometry InsideLocalGeometry;
      typedef const typename IntersectionData::Grid1LocalGeometry OutsideLocalGeometry;
      typedef const typename IntersectionData::Grid0Geometry Geometry;
      typedef const typename IntersectionData::Grid1Geometry OutsideGeometry;

      enum {
        coorddim = IntersectionData::coorddim,
        mydim = IntersectionData::mydim,
        insidePatch = 0,
        outsidePatch = 1
      };

      typedef typename GridGlue::ctype ctype;
      typedef Dune::FieldVector<ctype, mydim>      LocalCoordinate;
      typedef Dune::FieldVector<ctype, coorddim>   GlobalCoordinate;
    };

    template<typename P0, typename P1>
    struct IntersectionTraits<P0,P1,1,0>
    {
      typedef ::GridGlue<P0, P1> GridGlue;
      typedef Dune::GridGlue::IntersectionData<P0,P1> IntersectionData;

      typedef typename GridGlue::Grid1View InsideGridView;
      typedef typename GridGlue::Grid0View OutsideGridView;

      typedef const typename IntersectionData::Grid1LocalGeometry InsideLocalGeometry;
      typedef const typename IntersectionData::Grid0LocalGeometry OutsideLocalGeometry;
      typedef const typename IntersectionData::Grid1Geometry Geometry;
      typedef const typename IntersectionData::Grid0Geometry OutsideGeometry;
      typedef const typename IntersectionData::Grid1IndexType InsideIndexType;
      typedef const typename IntersectionData::Grid0IndexType OutsideIndexType;

      enum {
        coorddim = IntersectionData::coorddim,
        mydim = IntersectionData::mydim,
        insidePatch = 1,
        outsidePatch = 0
      };

      typedef typename GridGlue::ctype ctype;
      typedef Dune::FieldVector<ctype, mydim>      LocalCoordinate;
      typedef Dune::FieldVector<ctype, coorddim>   GlobalCoordinate;
    };

    /** @brief The intersection of two entities of the two patches of a GridGlue
     */
    template<typename P0, typename P1, int I, int O>
    class Intersection
    {

    public:

      typedef IntersectionTraits<P0,P1,I,O> Traits;

      typedef typename Traits::GridGlue GridGlue;
      typedef typename Traits::IntersectionData IntersectionData;


      typedef typename Traits::InsideGridView InsideGridView;
      typedef typename Traits::InsideLocalGeometry InsideLocalGeometry;

      typedef typename Traits::OutsideGridView OutsideGridView;
      typedef typename Traits::OutsideLocalGeometry OutsideLocalGeometry;
      typedef typename Traits::OutsideGeometry OutsideGeometry;

      typedef typename Traits::Geometry Geometry;
      typedef typename Traits::ctype ctype;

      typedef typename InsideGridView::Traits::template Codim<0>::Entity InsideEntity;
      typedef typename InsideGridView::Traits::template Codim<0>::EntityPointer InsideEntityPointer;

      typedef typename OutsideGridView::Traits::template Codim<0>::Entity OutsideEntity;
      typedef typename OutsideGridView::Traits::template Codim<0>::EntityPointer OutsideEntityPointer;

      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;

      enum {
        /** \brief Dimension of the world space of the intersection */
        coorddim = Traits::coorddim,
        /** \brief Dimension of the intersection */
        mydim = Traits::mydim,
        /** \brief document the inside & outside domain */
        //! @{
        insidePatch = Traits::insidePatch,
        outsidePatch = Traits::outsidePatch
                       //! @}
      };

      // typedef unsigned int IndexType;

      /*   C O N S T R U C T O R S   */

      /** \brief Constructor for a given Dataset */
      Intersection(const GridGlue* glue, const IntersectionData*  i) :
        glue_(glue), i_(i) {}

      /*   F U N C T I O N A L I T Y   */

      /** \brief return EntityPointer to the Entity on the inside of this intersection.
              That is the Entity where we started this. */
      InsideEntityPointer inside() const
      {
        assert(self());
        return glue_->template patch<I>().element(
                 IntersectionDataView<P0,P1,I>::index(*i_));
      }

      /** \brief return EntityPointer to the Entity on the outside of this intersection. That is the neighboring Entity. */
      OutsideEntityPointer outside() const
      {
        assert(neighbor());
        return glue_->template patch<O>().element(
                 IntersectionDataView<P0,P1,O>::index(*i_));
      }

      /** \brief Return true if intersection is conforming */
      bool conforming() const
      {
        assert(("not implemented", false));
      }

      /** \brief geometrical information about this intersection in local coordinates of the inside() entity.
          takes local domain intersection coords and maps them to domain parent element local coords */
      const InsideLocalGeometry& geometryInInside() const
      {
        return IntersectionDataView<P0,P1,I>::localGeometry(*i_);
      }

      /** \brief geometrical information about this intersection in local coordinates of the outside() entity.
          takes local target intersection coords and maps them to target parent element local coords */
      const OutsideLocalGeometry& geometryInOutside() const
      {
        return IntersectionDataView<P0,P1,O>::localGeometry(*i_);
      }

      /** \brief geometrical information about this intersection in global coordinates of the inside grid.
          takes local domain intersection coords and maps them to domain grid world coords */
      const Geometry& geometry() const
      {
        return IntersectionDataView<P0,P1,I>::geometry(*i_);
      }

      /** \brief geometrical information about this intersection in global coordinates of the outside grid.
          takes local domain intersection coords and maps them to target grid world coords */
      const OutsideGeometry& geometryOutside() const // DUNE_DEPRECATED
      {
        return IntersectionDataView<P0,P1,O>::geometry(*i_);
      }

      /** \brief obtain the type of reference element for this intersection */
      Dune::GeometryType type() const
      {
        #ifdef ONLY_SIMPLEX_INTERSECTIONS
        static const Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
        return type;
        #else
        #error Not Implemented
        #endif
      }


      /** \brief return true if inside() entity exists locally */
      bool self() const
      {
        return IntersectionDataView<P0,P1,I>::local(*i_);
      }

      /** \brief return true if outside() entity exists locally */
      bool neighbor() const
      {
        return IntersectionDataView<P0,P1,O>::local(*i_);
      }

      /** \brief Local number of codim 1 entity in the inside() Entity where intersection is contained in. */
      int indexInInside() const
      {
        assert(self());
        return glue_->template patch<I>().indexInInside(
                 IntersectionDataView<P0,P1,I>::index(*i_));
      }

      /** \brief Local number of codim 1 entity in outside() Entity where intersection is contained in. */
      int indexInOutside() const
      {
        assert(neighbor());
        return glue_->template patch<O>().indexInInside(
                 IntersectionDataView<P0,P1,O>::index(*i_));
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate outerNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        Dune::FieldVector<ctype, coorddim> normal;

        // Codimension with respect to the world(!)
        int codimension = coorddim - mydim;

        if (codimension == 0)

          DUNE_THROW(Dune::Exception, "There is no normal vector to a full-dimensional intersection");

        else if (codimension == 1) {

          /** \todo Implement the general n-ary cross product here */
          FieldMatrix<ctype, mydim,coorddim> jacobianTransposed = geometry().jacobianTransposed(local);
          if (mydim==1) {
            normal[0] = - jacobianTransposed[0][1];
            normal[1] =   jacobianTransposed[0][0];
          } else if (mydim==2) {
            normal[0] =   (jacobianTransposed[0][1] * jacobianTransposed[1][2] - jacobianTransposed[0][2] * jacobianTransposed[1][1]);
            normal[1] = - (jacobianTransposed[0][0] * jacobianTransposed[1][2] - jacobianTransposed[0][2] * jacobianTransposed[1][0]);
            normal[2] =   (jacobianTransposed[0][0] * jacobianTransposed[1][1] - jacobianTransposed[0][1] * jacobianTransposed[1][0]);
          } else
            DUNE_THROW(Dune::NotImplemented, "Remote intersections don't implement the 'outerNormal' method for " << mydim << "-dimensional intersections yet");

        } else
          DUNE_THROW(Dune::NotImplemented, "Remote intersections don't implement the 'outerNormal' method for intersections with codim >= 2 yet");

        return normal;
      }

      /** \brief Return a unit outer normal of the target intersection */
      GlobalCoordinate unitOuterNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        Dune::FieldVector<ctype, coorddim> normal = outerNormal(local);
        normal /= normal.two_norm();
        return normal;
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate integrationOuterNormal(const Dune::FieldVector<ctype, mydim> &local) const
      {
        return (unitOuterNormal(local) *= geometry().integrationElement(local));
      }

      /** \brief Unit outer normal at the center of the intersection
       *
       * Used for some grids that do not implement element geometries
       */
      GlobalCoordinate centerUnitOuterNormal () const
      {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
        return unitOuterNormal(ReferenceElements<ctype,mydim>::general(type()).position(0,0));
#else
        return unitOuterNormal(GenericReferenceElements<ctype,mydim>::general(type()).position(0,0));
#endif
      }

#ifdef QUICKHACK_INDEX
      typedef typename GridGlue::IndexType IndexType;

      IndexType index() const
      {
        return i_->index_;
      }

#endif

    private:

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief the grid glue entity this is built on
      const GridGlue*       glue_;

      /// @brief the underlying remote intersection
      const IntersectionData* i_;
    };


  } // end namesapce GridGlue
} // end namesapce Dune

#endif // DUNE_GRIDGLUE_REMOTEINTERSECTION_HH
