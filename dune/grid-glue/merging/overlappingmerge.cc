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
int OverlappingMerge<dim1,dim2,dimworld, T>::intersectionIndex(int grid1Index, int grid2Index,
                                                                  std::vector<FieldVector<T,dim1> > & g1Local,
                                                                  std::vector<FieldVector<T,dim2> > & g2Local) {

    int n_intersections = this->intersections_.size();

    if (dim1 != dim2)  {
        int i,j,k,l;
        int isDim = std::min(dim1,dim2);
        int n_parents1, n_parents2;
        bool b1,b2;
        T d;
        T eps = 1e-10;

        for (i=n_intersections-1; i >=  0; --i) {
            n_parents1 = this->intersections_[i].grid1Entities_.size();
            n_parents2 = this->intersections_[i].grid2Entities_.size();

            if (n_parents1 > 0 || n_parents2 > 0) {
                for (j = 0; j < n_parents1; ++j) {
                    b2 = true;

                    // an existing intersection candidate must be contained in at least one grid elem already
                    if (this->intersections_[i].grid1Entities_[j] == grid1Index) {
                        // start comparing local grid 1 intersection geometries to local grid 1 nodes
                        for (k = 0; k < isDim+1; ++k) {
                            b1 = false;
                            FieldVector<T,dim1>  v = this->intersections_[i].grid1Local_[j][k];
                            for (l = 0; l < isDim+1; ++l) {
                                FieldVector<T,dim1> w = g1Local[l];
                                d = (v-w).infinity_norm();
                                b1 = b1 || (d < eps);
                                if (d < eps)
                                    break;
                            }

                            b2 = b2 && b1;

                            if (!b2)
                                break;

                        }

                        if (b2){
                            return i;
                        }
                    }
                }
                for (j = 0; j < n_parents2;++j) {
                    b2 = true;

                    // an existing intersection candidate must be contained in at least one grid elem already
                    if (this->intersections_[i].grid2Entities_[j] == grid2Index) {

                        // start comparing local grid 1 intersection geometries to local grid 1 nodes
                        for (k = 0; k < isDim+1; ++k) {
                            b1 = false;
                            FieldVector<T,dim2> v = this->intersections_[i].grid2Local_[j][k];
                            for (l = 0; l < isDim+1; ++l) {
                                FieldVector<T,dim2> w = g2Local[l];
                                d = (v-w).infinity_norm();
                                b1 = b1 || (d < eps);
                                if (d < eps)
                                    break;
                            }

                            b2 = b2 && b1;
                            if (!b2)
                                break;
                        }

                        if (b2)  {
                            return i;
                        }
                    }
                }
            }
        }
    }

    return n_intersections;
}

template<int dim1, int dim2, int dimworld, typename T>
void OverlappingMerge<dim1,dim2,dimworld, T>::computeIntersection(const Dune::GeometryType& grid1ElementType,
                                                                     const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                                                     unsigned int grid1Index,
                                                                     std::bitset<(1<<dim1)>& neighborIntersects1,
                                                                     const Dune::GeometryType& grid2ElementType,
                                                                     const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                                                     unsigned int grid2Index,
                                                                     std::bitset<(1<<dim2)>& neighborIntersects2)
{
    this->counter++;

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
    size_type i, isindex;

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
        isindex = intersectionIndex(grid1Index,grid2Index,g1local,g2local);  // check whether the intersection is contained in the intersections already
        bool isinplane = false;
        if (dimis == 3)
            isinplane = inPlane(P);

        if (isindex >= this->intersections_.size() && !isinplane) { // a new intersection is found
            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));
            for (i = 0; i < n_intersectionnodes; ++i) {
                this->intersections_.back().grid1Local_[0][i] = g1local[i];
                this->intersections_.back().grid2Local_[0][i] = g2local[i];
            }

        } else if (isindex < this->intersections_.size() && !isinplane) {// intersections equal

            array<FieldVector<T,dim1>, n_intersectionnodes > g1localarr;
            array<FieldVector<T,dim2>, n_intersectionnodes > g2localarr;

            for (i = 0; i < n_intersectionnodes; ++i) {
                g1localarr[i] = g1local[i];
                g2localarr[i] = g2local[i];
            }

            // not efficient, we push back local coordinates of at least one element twice
            this->intersections_[isindex].grid1Local_.push_back(g1localarr);
            this->intersections_[isindex].grid2Local_.push_back(g2localarr);
            this->intersections_[isindex].grid1Entities_.push_back(grid1Index);
            this->intersections_[isindex].grid2Entities_.push_back(grid2Index);
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

                // check whether the intersection is contained in the intersections already
                isindex = intersectionIndex(grid1Index,grid2Index,g1local,g2local);

                if (isindex >= this->intersections_.size()  && !isinplane) { // a new intersection is found
                    this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));
                    for (size_type j = 0; j < n_intersectionnodes; ++j) {
                        this->intersections_.back().grid1Local_[0][j] = g1local[j];
                        this->intersections_.back().grid2Local_[0][j] = g2local[j];
                    }

                } else if (isindex < this->intersections_.size() && !isinplane) {// intersections equal

                    array<FieldVector<T,dim1>, n_intersectionnodes > g1localarr;
                    array<FieldVector<T,dim2>, n_intersectionnodes > g2localarr;

                    for (size_type j = 0; j < n_intersectionnodes; ++j) {
                        g1localarr[j] = g1local[j];
                        g2localarr[j] = g2local[j];
                    }

                    // not efficient, we push back local coordinates of at least one element twice
                    this->intersections_[isindex].grid1Local_.push_back(g1localarr);
                    this->intersections_[isindex].grid2Local_.push_back(g2localarr);
                    this->intersections_[isindex].grid1Entities_.push_back(grid1Index);
                    this->intersections_[isindex].grid2Entities_.push_back(grid2Index);
                }
            }
        }
    }
}

} /* namespace Dune::GridGlue */
} /* namespace Dune */

#endif // DUNE_GRIDGLUE_OVERLAPPINGMERGE_HH
