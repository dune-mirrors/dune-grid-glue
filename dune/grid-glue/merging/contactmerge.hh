// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
  @file
  @brief Merge two grid boundary surfaces that may be a positive distance apart
 */

#ifndef DUNE_GRID_GLUE_MERGING_CONTACTMERGE_HH
#define DUNE_GRID_GLUE_MERGING_CONTACTMERGE_HH


#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <limits>

#include <dune/common/fvector.hh>
#include <dune/common/function.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/bitsetvector.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/standardmerge.hh>
#include <dune/grid-glue/gridglue.hh>
#include <dune/grid-glue/extractors/extractorpredicate.hh>

/** \brief Merge two codimension-1 surfaces that may be a positive distance apart

  \tparam dimworld  Dimension of the world coordinates.
  \tparam T Type used for coordinates
 */
template<int dimworld, typename T = double>
class ContactMerge
: public StandardMerge<T,dimworld-1,dimworld-1,dimworld>
{
    enum {dim = dimworld-1};

    static_assert( dim==1 || dim==2,
            "ContactMerge yet only handles the cases dim==1 and dim==2!");

    typedef StandardMerge<T,dim,dim,dimworld> Base;
public:

    /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

    /// @brief the numeric type used in this interface
    typedef T ctype;

    /// @brief the coordinate type used in this interface
    typedef Dune::FieldVector<T, dimworld>  WorldCoords;

    /// @brief the coordinate type used in this interface
    typedef Dune::FieldVector<T, dim>  LocalCoords;

    ContactMerge() {}

    ContactMerge(const T allowedOverlap,
            const Dune::VirtualFunction<WorldCoords,WorldCoords>* domainDirections = NULL,
            const Dune::VirtualFunction<WorldCoords,WorldCoords>* targetDirections = NULL)
        : domainDirections_(domainDirections), targetDirections_(targetDirections),
          overlap_(allowedOverlap)
    {}

    /**
     * @brief Set surface direction functions
     *
     * The matching of the geometries offers the possibility to specify a function for
     * the exact evaluation of domain surface normals. If no such function is specified
     * (default) normals are interpolated.
     * @param value the new function (or NULL to unset the function)
     */
    inline
    void setSurfaceDirections(const Dune::VirtualFunction<WorldCoords,WorldCoords>* domainDirections,
                              const Dune::VirtualFunction<WorldCoords,WorldCoords>* targetDirections)
    {
            domainDirections_ = domainDirections;
            targetDirections_ = targetDirections;
            this->valid = false;
    }

    //! Set the allowed overlap of the surfaces.
    void setOverlap(T overlap)
    {
        overlap_ = overlap;
    }

private:
    /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
     */
    const Dune::VirtualFunction<WorldCoords,WorldCoords>* domainDirections_;
    std::vector<WorldCoords> nodalDomainDirections_;

    /** \brief Vector field on the target surface which prescribes a 'forward'
      direction.

      We use the normals of the target side to increase projection
      robustness.  If these cannot be computed from the surface directly
      (e.g. because it is not properly oriented), they can be given
      explicitly through the targetDirections field.
     */
    const Dune::VirtualFunction<WorldCoords,WorldCoords>* targetDirections_;
    std::vector<WorldCoords> nodalTargetDirections_;

    //! Allow some overlap, i.e. also look in the negative projection directions
    T overlap_;

    /** \brief Compute the intersection between two overlapping elements
     *
     *   The result is a set of simplices.
     */
    void computeIntersection(const Dune::GeometryType& grid1ElementType,
            const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
            unsigned int grid1Index,
            std::bitset<(1<<dim)>& neighborIntersects1,
            const Dune::GeometryType& grid2ElementType,
            const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
            unsigned int grid2Index,
            std::bitset<(1<<dim)>& neighborIntersects2);

    /**
      * @copydoc StandardMerge<T,grid1Dim,grid2Dim,dimworld>::build
     */
    void build(const std::vector<Dune::FieldVector<T,dimworld> >& grid1Coords,
        const std::vector<unsigned int>& grid1Elements,
        const std::vector<Dune::GeometryType>& grid1ElementTypes,
        const std::vector<Dune::FieldVector<T,dimworld> >& grid2Coords,
        const std::vector<unsigned int>& grid2Elements,
        const std::vector<Dune::GeometryType>& grid2ElementTypes)
    {
        std::cout<<"ContactMerge building grid!\n";
        // setup the nodal direction vectors
        setupNodalDirections(grid1Coords, grid1Elements, grid1ElementTypes,
                grid2Coords, grid2Elements, grid2ElementTypes);

        Base::build(grid1Coords, grid1Elements, grid1ElementTypes,
                   grid2Coords, grid2Elements, grid2ElementTypes);

    }

private:

    //! Compute local coordinates of a corner
    static LocalCoords localCornerCoords(int i, const Dune::GeometryType& gt)
    {
        const Dune::ReferenceElement<T,dim>& ref = Dune::ReferenceElements<T,dim>::general(gt);
        return ref.position(i,dim);
    }

    //! Check if local coordinates correspond to a vertex and if yes returns its local index
    static int isCorner(const Dune::GeometryType& gt, const LocalCoords& local)
    {
        const Dune::ReferenceElement<T,dim>& ref = Dune::ReferenceElements<T,dim>::general(gt);

        for (int i=0; i<ref.size(dim); i++)
            if ((ref.position(i,dim)-local).two_norm()<1e-8)
                return i;
        return -1;
    }

    //! Transform local to barycentric coordinates
    static std::vector<T> localToBarycentric(const Dune::GeometryType& gt, const LocalCoords& local)
    {
        std::vector<T> bar(Dune::ReferenceElements<T,dim>::general(gt).size(dim));
        if (bar.size()<4) {
            bar[0] = 1.0;
            for (int i=0; i<dim; i++) {
                bar[i+1] = local[i];
                bar[0] -= local[i];
            }
        } else {
            // quadrilateral case
            for (int i=0; i<4; i++) {
                bar[i] = 1;

                for (int j=0; j<dim; j++)
                    // if j-th bit of i is set multiply with in[j], else with 1-in[j]
                    bar[i] *= (i & (1<<j)) ? local[j] :  1-local[j];
            }
        }
        return bar;
    }

    //! Compute global coordinates of a local point using linear interpolation of the corners
    static WorldCoords interpolate(const std::vector<WorldCoords>& p,
            const Dune::GeometryType& gt, const LocalCoords& local)
    {
        // Compute barycentric coordinates
        std::vector<T> bar = localToBarycentric(gt, local);

        WorldCoords global(0);
        for (size_t i=0; i<p.size(); i++)
            global.axpy(bar[i],p[i]);

        return global;
    }

    //! Order the corners of the intersection polytope in cyclic order
    void computeCyclicOrder(const std::vector<Dune::array<LocalCoords,2> >& polytopeCorners,
            const LocalCoords& center, std::vector<int>& ordering) const;

    //! Setup the direction vectors containing the directions for each vertex
    void setupNodalDirections(const std::vector<WorldCoords>& coords1,
            const std::vector<unsigned int>& elements1,
            const std::vector<Dune::GeometryType>& elementTypes1,
            const std::vector<WorldCoords>& coords2,
            const std::vector<unsigned int>& elements2,
            const std::vector<Dune::GeometryType>& elementTypes2);

    //! If no direction field was specified compute the outer normal field
    void computeOuterNormalField(const std::vector<WorldCoords>& coords,
            const std::vector<unsigned int>& elements,
            const std::vector<Dune::GeometryType>& elementTypes,
            std::vector<WorldCoords>& normals);

    //! Remove all multiples
    void removeDoubles(std::vector<Dune::array<LocalCoords,2> >& polytopeCorners);
};

#include "contactmerge.cc"

#endif // DUNE_GRID_GLUE_MERGING_CONTACTMERGE_HH
