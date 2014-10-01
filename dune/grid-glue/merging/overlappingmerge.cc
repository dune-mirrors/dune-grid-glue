// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRIDGLUE_OVERLAPPINGMERGE_CC
#define DUNE_GRIDGLUE_OVERLAPPINGMERGE_CC
//#include <algorithm>

namespace Dune {
namespace GridGlue {

template<int dim1, int dim2, int dimworld, typename T>
bool OverlappingMerge<dim1,dim2,dimworld, T>::inPlane(std::vector<FieldVector<T,dimworld> >& points) {

    T eps = 1e-8;

    assert(dim1 == 3 && dim2 == 3 && dimworld == 3);
    assert(points.size() == 4);

    FieldVector<T,dimworld> v1 = points[1]-points[0];
    FieldVector<T,dimworld> v2 = points[2]-points[0];
    FieldVector<T,dimworld> v3 = points[3]-points[0];

    FieldVector<T,dimworld> v1xv2;
    v1xv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
    v1xv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
    v1xv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

    return (std::abs(v3.dot(v1xv2)) < eps);
}

template<int dim1, int dim2, int dimworld, typename T>
void OverlappingMerge<dim1,dim2,dimworld, T>::computeIntersections(const Dune::GeometryType& grid1ElementType,
                                               const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                               std::bitset<(1<<dim1)>& neighborIntersects1,
                                               unsigned int grid1Index,
                                               const Dune::GeometryType& grid2ElementType,
                                               const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                               std::bitset<(1<<dim2)>& neighborIntersects2,
                                               unsigned int grid2Index,
                                               std::vector<RemoteSimplicialIntersection>& intersections)
{
    this->counter++;
    intersections.clear();

    typedef SimplexMethod<dimworld,dim1,dim2,T> CM;

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
    const Dune::ReferenceElement<T,dim1>& refElement1 = Dune::ReferenceElements<T,dim1>::general(grid1ElementType);
    const Dune::ReferenceElement<T,dim2>& refElement2 = Dune::ReferenceElements<T,dim2>::general(grid2ElementType);
#else
    const Dune::GenericReferenceElement<T,dim1>& refElement1 = Dune::GenericReferenceElements<T,dim1>::general(grid1ElementType);
    const Dune::GenericReferenceElement<T,dim2>& refElement2 = Dune::GenericReferenceElements<T,dim2>::general(grid2ElementType);
#endif

    // A few consistency checks
    assert((unsigned int)(refElement1.size(dim1)) == grid1ElementCorners.size());
    assert((unsigned int)(refElement2.size(dim2)) == grid2ElementCorners.size());

    // Make generic geometries representing the grid1- and grid2 element.
    // this eases computation of local coordinates.
    typedef CachedMultiLinearGeometry<T,dim1,dimworld> Geometry1;
    typedef CachedMultiLinearGeometry<T,dim2,dimworld> Geometry2;

    Geometry1 grid1Geometry(grid1ElementType, grid1ElementCorners);
    Geometry2 grid2Geometry(grid2ElementType, grid2ElementCorners);

    // get size_type for all the vectors we are using
    typedef typename std::vector<Empty>::size_type size_type;

    const int dimis = dim1 < dim2 ? dim1 : dim2;
    const size_type n_intersectionnodes = dimis+1;
    size_type i;

    std::vector<FieldVector<T,dimworld> >  P(0);
    std::vector<std::vector<int> >          H,SX(1<<dim1),SY(1<<dim2);
    FieldVector<T,dimworld>                centroid;
    std::vector<FieldVector<T,dim1> > g1local(n_intersectionnodes);
    std::vector<FieldVector<T,dim2> > g2local(n_intersectionnodes);

    // compute the intersection nodes
    bool b = IntersectionComputation<CM>::computeIntersection(grid1ElementCorners,grid2ElementCorners,SX,SY,P);

    for (size_type i = 0; i < neighborIntersects1.size(); ++i) {
        if (i < SX.size())
            neighborIntersects1[i] = (SX[i].size() > 0);
        else
            neighborIntersects1[i] = false;
    }
    for (size_type i = 0; i < neighborIntersects2.size(); ++i) {
        if (i < SY.size())
            neighborIntersects2[i] = (SY[i].size() > 0);
        else
            neighborIntersects2[i] = false;
    }

    // P is an simplex of dimension dimis
    if (P.size() == n_intersectionnodes) {

        for (i = 0; i < n_intersectionnodes; ++i) {
            g1local[i] = grid1Geometry.local(P[i]);
            g2local[i] = grid2Geometry.local(P[i]);
        }

        bool isinplane = false;
        if (dimis == 3)
            isinplane = inPlane(P);

        intersections.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));
        for (i = 0; i < n_intersectionnodes; ++i) {
            intersections.back().grid1Local_[0][i] = g1local[i];
            intersections.back().grid2Local_[0][i] = g2local[i];
        }

    } else if (P.size() > n_intersectionnodes) {  // P is a union of simplices of dimension dimis

        assert(dimis != 1);
        std::vector<FieldVector<T,dimworld> > global(n_intersectionnodes);

        // Compute the centroid
        centroid=0;
        for (size_type i=0; i < P.size(); i++)
            centroid += P[i] ;
        centroid /= static_cast<T>(P.size()) ;

        // order the points and get intersection face indices
        H.clear() ;
        IntersectionComputation<CM>::template orderPoints<dimis,dimworld>(centroid,SX,SY,P,H);

        //  Loop over all intersection elements
        for (size_type i=0; i < H.size(); i++) {
            int hs = H[i].size(); // number of nodes of the intersection

            // if the intersection element is not degenerated
            if (hs==dimis) {

                // create the intersection geometry
                for ( size_type j=0 ; j < dimis; ++j) {
                    global[j]= P[H[i][j]]; // get the intersection face
                }

                // intersection face + centroid = new element
                global[dimis]=centroid;

                // create local representation of the intersection
                for (size_type j = 0; j < n_intersectionnodes; ++j) {
                    g1local[j] = grid1Geometry.local(global[j]);
                    g2local[j] = grid2Geometry.local(global[j]);
                }

                // do not insert degenerated elements
                bool isinplane = false;
                if (dimis == 3)
                    isinplane = inPlane(global);

                intersections.push_back(RemoteSimplicialIntersection(grid1Index,grid2Index));
                for (size_type j = 0; j < n_intersectionnodes; ++j) {
                    intersections.back().grid1Local_[0][j] = g1local[j];
                    intersections.back().grid2Local_[0][j] = g2local[j];
                }
            }
        }
    }
}

} /* namespace Dune::GridGlue */
} /* namespace Dune */

#endif // DUNE_GRIDGLUE_OVERLAPPINGMERGE_CC
