// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file
   @author Christian Engwer
   @brief Model of the Intersection concept provided by GridGlue.
 */

#ifndef DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH
#define DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH

#include <memory>

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
#include <dune/geometry/affinegeometry.hh>
#else
#include <dune/grid-glue/common/simplexgeometry.hh>
#endif
#include <dune/grid-glue/gridglue.hh>

#define ONLY_SIMPLEX_INTERSECTIONS

namespace Dune {
  namespace GridGlue {

    // forward declaration
    template<typename P0, typename P1>
    class IntersectionIndexSet;

    /**
       @brief storage class for Dune::GridGlue::Intersection related data
     */
    template<typename P0, typename P1>
    class IntersectionData
    {
    public:
      typedef ::Dune::GridGlue::GridGlue<P0, P1> GridGlue;

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
      std::vector<Grid0IndexType> grid0indices_;    //!< indices of the associated local grid0 entity
      bool grid1local_;              //!< true if the associated grid1 entity is local
      std::vector<Grid1IndexType> grid1indices_;    //!< indices of the associated local grid1 entity

      /**
       * Embedding of intersection into local grid0 entity coordinates.
       */
      std::vector<shared_ptr<Grid0LocalGeometry> >  grid0localgeom_;
      /**
       * Global intersection geometry on grid0 side.
       *
       * This is the same as gₚ∘iₚ for any embedding iₚ into a grid0
       * entity as stored in grid0localgeom_ and that entities global
       * geometry gₚ.
       */
      shared_ptr<Grid0Geometry>       grid0geom_;
      /**
       * Embedding of intersection into local grid1 entity coordinates.
       */
      std::vector<shared_ptr<Grid1LocalGeometry> >  grid1localgeom_;
      /**
       * Global intersection geometry on grid1 side.
       *
       * This is the same as gₚ∘iₚ for any embedding iₚ into a grid1
       * entity as stored in grid1localgeom_ and that entities global
       * geometry gₚ.
       */
      shared_ptr<Grid1Geometry>       grid1geom_;

    };

    //! \todo move this functionality to GridGlue
    template<typename P0, typename P1>
    IntersectionData<P0, P1>::IntersectionData(const GridGlue& glue, unsigned int mergeindex, unsigned int offset,
                                               bool grid0local, bool grid1local)
      : index_(mergeindex+offset),
        grid0local_(grid0local),
        grid1local_(grid1local)
    {
      unsigned int n_grid0Parents = glue.merger_->template parents<0>(mergeindex);
      unsigned int n_grid1Parents = glue.merger_->template parents<1>(mergeindex);

      assert (0 <= n_grid0Parents || 0 <= n_grid1Parents);

      // init containers
      grid0indices_.resize(n_grid0Parents);
      grid0localgeom_.resize(n_grid0Parents);

      grid1indices_.resize(n_grid1Parents);
      grid1localgeom_.resize(n_grid1Parents);

      // default values
      grid0indices_[0] = 0;
      grid1indices_[0] = 0;

      typedef typename GridGlue::ctype ctype;

      // Number of corners of the intersection
      const int nSimplexCorners = mydim + 1;

      // if an invalid index is given do not proceed!
      // (happens when the parent GridGlue initializes the "end"-Intersection)
      assert (0 <= mergeindex || mergeindex < glue.index__sz);

      // initialize the local and the global geometries of grid0
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = GridGlue::Grid0View::template Codim<0>::Geometry::mydimension;

        // coordinates within the subentity that contains the remote intersection
        std::array<Dune::FieldVector<ctype, dim1>, nSimplexCorners> corners_subEntity_local;

        for (unsigned int par = 0; par < n_grid0Parents; ++par) {
            for (int i = 0; i < nSimplexCorners; ++i)
              corners_subEntity_local[i] = glue.merger_->template parentLocal<0>(mergeindex, i, par);

            // Coordinates of the remote intersection corners wrt the element coordinate system
            std::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

            if (grid0local)
            {
              grid0indices_[par] = glue.merger_->template parent<0>(mergeindex,par);

              typename GridGlue::Grid0Patch::LocalGeometry
              grid0LocalGeometry = glue.template patch<0>().geometryLocal(grid0indices_[par]);
              typename GridGlue::Grid0Patch::Geometry grid0WorldGeometry1 =  glue.template patch<0>().geometry(grid0indices_[par]);
              for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
                corners_element_local[i] = grid0LocalGeometry.global(corners_subEntity_local[i]);
              }

              // set the corners of the local geometry
    #ifdef ONLY_SIMPLEX_INTERSECTIONS
              Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
    #else
    #error Not Implemented
    #endif
              grid0localgeom_[par] = make_shared<Grid0LocalGeometry>(type, corners_element_local);

              // Add world geometry only for 0th parent
              if (par == 0) {
                typename GridGlue::Grid0Patch::Geometry
                grid0WorldGeometry = glue.template patch<0>().geometry(grid0indices_[par]);

                // world coordinates of the remote intersection corners
                std::array<Dune::FieldVector<ctype, GridGlue::Grid0View::dimensionworld>, nSimplexCorners> corners_global;

                for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
                  corners_global[i]        = grid0WorldGeometry.global(corners_subEntity_local[i]);
                }

                grid0geom_      = make_shared<Grid0Geometry>(type, corners_global);
              }
            }
        }
      }

      // do the same for the local and the global geometry of grid1
      {
        // compute the coordinates of the subface's corners in codim 0 entity local coordinates
        const int elementdim = GridGlue::Grid1View::template Codim<0>::Geometry::mydimension;

        // coordinates within the subentity that contains the remote intersection
        std::array<Dune::FieldVector<ctype, dim2>, nSimplexCorners> corners_subEntity_local;

        for (unsigned int par = 0; par < n_grid1Parents; ++par) {

            for (int i = 0; i < nSimplexCorners; ++i)
              corners_subEntity_local[i] = glue.merger_->template parentLocal<1>(mergeindex, i, par);

            // Coordinates of the remote intersection corners wrt the element coordinate system
            std::array<Dune::FieldVector<ctype, elementdim>, nSimplexCorners> corners_element_local;

            if (grid1local)
            {
              grid1indices_[par] = glue.merger_->template parent<1>(mergeindex, par);

              typename GridGlue::Grid1Patch::LocalGeometry
              grid1LocalGeometry = glue.template patch<1>().geometryLocal(grid1indices_[par]);

              for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
                corners_element_local[i] = grid1LocalGeometry.global(corners_subEntity_local[i]);
              }

              // set the corners of the geometries
    #ifdef ONLY_SIMPLEX_INTERSECTIONS
              Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
    #else
    #error Not Implemented
    #endif
              grid1localgeom_[par] = make_shared<Grid1LocalGeometry>(type, corners_element_local);

              // Add world geomety only for 0th parent
              if (par == 0) {
                typename GridGlue::Grid1Patch::Geometry
                grid1WorldGeometry = glue.template patch<1>().geometry(grid1indices_[par]);

                // world coordinates of the remote intersection corners
                std::array<Dune::FieldVector<ctype, GridGlue::Grid1View::dimensionworld>, nSimplexCorners> corners_global;

                for (std::size_t i=0; i<corners_subEntity_local.size(); i++) {
                  corners_global[i]        = grid1WorldGeometry.global(corners_subEntity_local[i]);
                }

                grid1geom_      = make_shared<Grid1Geometry>(type, corners_global);
              }
            }
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
      static LocalGeometry& localGeometry(const IntersectionData<P0,P1> & i, unsigned int parentId = 0)
      {
          return *i.grid0localgeom_[parentId];
      }
      static Geometry& geometry(const IntersectionData<P0,P1> & i)
      {
          return *i.grid0geom_;
      }
      static bool local(const IntersectionData<P0,P1> & i)
      {
          return i.grid0local_;
      }
      static IndexType index(const IntersectionData<P0,P1> & i, unsigned int parentId = 0)
      {
          return i.grid0indices_[parentId];
      }
      static IndexType parents(const IntersectionData<P0,P1> & i)
      {
          return i.grid0indices_.size();
      }
    };

    template<typename P0, typename P1>
    struct IntersectionDataView<P0, P1, 1>
    {
      typedef const typename IntersectionData<P0,P1>::Grid1LocalGeometry LocalGeometry;
      typedef const typename IntersectionData<P0,P1>::Grid1Geometry Geometry;
      typedef const typename IntersectionData<P0,P1>::Grid1IndexType IndexType;
      static LocalGeometry& localGeometry(const IntersectionData<P0,P1> & i, unsigned int parentId = 0)
      {
          return *i.grid1localgeom_[parentId];
      }
      static Geometry& geometry(const IntersectionData<P0,P1> & i)
      {
          return *i.grid1geom_;
      }
      static IndexType local(const IntersectionData<P0,P1> & i)
      {
          return i.grid1local_;
      }
      static IndexType index(const IntersectionData<P0,P1> & i, unsigned int parentId = 0)
      {
          return i.grid1indices_[parentId];
      }
      static IndexType parents(const IntersectionData<P0,P1> & i)
      {
          return i.grid1indices_.size();
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
      typedef ::Dune::GridGlue::GridGlue<P0, P1> GridGlue;
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
      typedef ::Dune::GridGlue::GridGlue<P0, P1> GridGlue;
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
        /** \brief document the inside & outside patch */
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

      /** \brief Return entity on the inside of this intersection.
       */
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4) || DOXYGEN
      InsideEntity
#else
      InsideEntityPointer
#endif
      inside(unsigned int parentId = 0) const
      {
        assert(self());
        return glue_->template patch<I>().element(
                 IntersectionDataView<P0,P1,I>::index(*i_, parentId));
      }

      /** \brief Return entity on the outside of this intersection.
       */
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4) || DOXYGEN
      OutsideEntity
#else
      OutsideEntityPointer
#endif
      outside(unsigned int parentId = 0) const
      {
        assert(neighbor());
        return glue_->template patch<O>().element(
                 IntersectionDataView<P0,P1,O>::index(*i_, parentId));
      }

      /** \brief Return true if intersection is conforming */
      bool conforming() const
      {
        throw Dune::NotImplemented();
      }

      /** \brief Geometric information about this intersection in local coordinates of the inside() entity.
       */
      const InsideLocalGeometry& geometryInInside(unsigned int parentId = 0) const
      {
        return IntersectionDataView<P0,P1,I>::localGeometry(*i_, parentId);
      }

      /** \brief Geometric information about this intersection in local coordinates of the outside() entity.
       */
      const OutsideLocalGeometry& geometryInOutside(unsigned int parentId = 0) const
      {
        return IntersectionDataView<P0,P1,O>::localGeometry(*i_, parentId);
      }

      /** \brief Geometric information about this intersection as part of the inside grid.
       */
      const Geometry& geometry() const
      {
        return IntersectionDataView<P0,P1,I>::geometry(*i_);
      }

      /** \brief Geometric information about this intersection as part of the outside grid.
       */
      const OutsideGeometry& geometryOutside() const // DUNE_DEPRECATED
      {
        return IntersectionDataView<P0,P1,O>::geometry(*i_);
      }

      /** \brief Type of reference element for this intersection */
      Dune::GeometryType type() const
      {
        #ifdef ONLY_SIMPLEX_INTERSECTIONS
        static const Dune::GeometryType type(Dune::GeometryType::simplex, mydim);
        return type;
        #else
        #error Not Implemented
        #endif
      }


      /** \brief For parallel computations: Return true if inside() entity exists locally */
      bool self() const
      {
        return IntersectionDataView<P0,P1,I>::local(*i_);
      }

      /** \brief Return number of embeddings into local grid0 (grid1) entities. */
      size_t neighbor(unsigned int g = 0) const
      {
          if (g == 0 && IntersectionDataView<P0,P1,O>::local(*i_)) {
            return IntersectionDataView<P0,P1,O>::parents(*i_);
          } else if (g == 1  && IntersectionDataView<P0,P1,I>::local(*i_)) {
            return IntersectionDataView<P0,P1,I>::parents(*i_);
          }
          return 0;
      }

      /** \brief Local number of codim 1 entity in the inside() Entity where intersection is contained in. */
      int indexInInside(unsigned int parentId = 0) const
      {
        assert(self());
        return glue_->template patch<I>().indexInInside(
                 IntersectionDataView<P0,P1,I>::index(*i_, parentId));
      }

      /** \brief Local number of codim 1 entity in outside() Entity where intersection is contained in. */
      int indexInOutside(unsigned int parentId = 0) const
      {
        assert(neighbor());
        return glue_->template patch<O>().indexInInside(
                 IntersectionDataView<P0,P1,O>::index(*i_, parentId));
      }

      /** \brief Return an outer normal (length not necessarily 1) */
      GlobalCoordinate outerNormal(const LocalCoordinate &local) const
      {
        GlobalCoordinate normal;

        // Codimension with respect to the world(!)
        int codimension = coorddim - mydim;

        if (codimension == 0)

          DUNE_THROW(Dune::Exception, "There is no normal vector to a full-dimensional intersection");

        else if (codimension == 1) {

          /** \todo Implement the general n-ary cross product here */
          const auto jacobianTransposed = geometry().jacobianTransposed(local);
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

      /** \brief Return a unit outer normal */
      GlobalCoordinate unitOuterNormal(const LocalCoordinate &local) const
      {
        Dune::FieldVector<ctype, coorddim> normal = outerNormal(local);
        normal /= normal.two_norm();
        return normal;
      }

      /** \brief Return an outer normal with the length of the integration element */
      GlobalCoordinate integrationOuterNormal(const LocalCoordinate &local) const
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

      Intersection<P0,P1,O,I> flip() const
      {
        return Intersection<P0,P1,O,I>(glue_,i_);
      }

#ifdef QUICKHACK_INDEX
      typedef typename GridGlue::IndexType IndexType;

      IndexType index() const
      {
        return i_->index_;
      }

#endif

    private:

      friend class IntersectionIndexSet<P0,P1>;

      /*   M E M B E R   V A R  I A B L E S   */

      /// @brief the grid glue entity this is built on
      const GridGlue*       glue_;

      /// @brief the underlying remote intersection
      const IntersectionData* i_;
    };


  } // end namespace GridGlue
} // end namespace Dune

#endif // DUNE_GRIDGLUE_ADAPTER_INTERSECTION_HH
