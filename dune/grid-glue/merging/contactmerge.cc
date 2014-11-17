// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include <dune/grid-glue/common/projectionhelper.hh>

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeIntersection(const Dune::GeometryType& grid1ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                                   unsigned int grid1Index,
                                   std::bitset<(1<<dim)>& neighborIntersects1,
                                   const Dune::GeometryType& grid2ElementType,
                                   const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                                   unsigned int grid2Index,
                                   std::bitset<(1<<dim)>& neighborIntersects2)
{
    std::vector<Dune::array<LocalCoords,2> > polytopeCorners;

    // Initialize
    neighborIntersects1.reset();
    neighborIntersects2.reset();

    const int nCorners1 = grid1ElementCorners.size();
    const int nCorners2 = grid2ElementCorners.size();

    // The grid1 projection directions
    std::vector<WorldCoords> directions1(nCorners1);
    for (size_t i=0; i<directions1.size(); i++)
        directions1[i] = nodalDomainDirections_[this->grid1ElementCorners_[grid1Index][i]];

    // The grid2 projection directions
    // These are only needed in a check for validity of the projection and they should be chosen
    // to be some outer pointing directions, e.g. outer normals
    std::vector<WorldCoords> directions2(nCorners1);
    for (size_t i=0; i<directions1.size(); i++)
        directions2[i] = nodalTargetDirections_[this->grid2ElementCorners_[grid2Index][i]];

    /////////////////////////////////////////////////////
    //  Compute all corners of the intersection polytope
    /////////////////////////////////////////////////////


    // If we hit a corner of the grid1 element we store its index
    std::vector<int> hitCorners(nCorners2,-1);
    // If a corner was hit already we don't have to project it
    std::vector<bool> skip(nCorners1);
    std::vector<bool> proj(nCorners2);

    // local coordinates
    LocalCoords localCoords;
    T eps = 1e-4;

    // Compute which points of the grid1 element lie inside the grid2 element
    for (size_t i=0; i<grid2ElementCorners.size(); i++)
        if (Projection::template inverseProjection<dim,dimworld,T>(grid1ElementCorners, directions1,
                                                                grid2ElementCorners[i], localCoords, overlap_)) {

            // There are case where the intersection is not valid
            WorldCoords base = interpolate(grid1ElementCorners,grid1ElementType,localCoords);
            WorldCoords baseNormal = interpolate(directions1,grid1ElementType,localCoords);

            WorldCoords segment = grid2ElementCorners[i] - base;
            T distance = segment.two_norm();

            // If both criterions are invalid
            if ((segment*directions2[i]) > eps
                 && (segment*baseNormal) < -eps)
            {
                // then the projection is still valid if the overlap is small
                if (distance > overlap_)
                    continue;

            // if only one is invalid then the projection is invalid
            } else if ((segment*directions2[i]) > eps
                    || (segment*baseNormal) < -eps)
                continue;

            // corner could be projected
            proj[i] = true;

            Dune::array<LocalCoords,2> corner;
            corner[1] = localCornerCoords(i,grid2ElementType);
            corner[0] = localCoords;

            polytopeCorners.push_back(corner);

            hitCorners[i] = isCorner(grid1ElementType,localCoords);

            if (hitCorners[i] != -1)
                skip[hitCorners[i]]=true;
        }

    // Compute which points of the grid2 element lie inside the grid1 element
    for (size_t i=0; i<grid1ElementCorners.size(); i++) {

        if (skip[i])
            continue;

        if (Projection::template projection<dim,dimworld,T>(grid1ElementCorners[i], directions1[i],
                                                                grid2ElementCorners, localCoords, overlap_)) {
            Dune::array<LocalCoords,2> corner;
            corner[0] = localCornerCoords(i,grid1ElementType);
            corner[1] = localCoords;

            // corner could be projected
            skip[i] = true;

            polytopeCorners.push_back(corner);
        }
    }

    // check which neighbors might also intersect
    const Dune::ReferenceElement<T,dim>& ref2 = Dune::ReferenceElements<T,dim>::general(grid2ElementType);
    for (int i=0; i<ref2.size(1); i++) {

        // if all face corners hit the the other element then
        // the neighbor might also intersect

        bool intersects(true);
        for (int k=0; k<ref2.size(i,1,dim); k++)
            intersects &= proj[ref2.subEntity(i,1,k,dim)];

        if (intersects)
            neighborIntersects2[i] = true;
    }

    const Dune::ReferenceElement<T,dim>& ref1 = Dune::ReferenceElements<T,dim>::general(grid1ElementType);
    for (int i=0; i<ref1.size(1); i++) {

        // if all face corners hit the the other element then
        // the neighbor might also intersect

        bool intersects(true);
        for (int k=0; k<ref1.size(i,1,dim); k++)
            intersects &= skip[ref1.subEntity(i,1,k,dim)];

        if (intersects)
            neighborIntersects1[i] = true;
    }

    // Compute the edge intersections
    Projection::template addEdgeIntersections<dim,dimworld,T>(grid1ElementCorners,grid2ElementCorners,
                                                              directions1, grid1ElementType,
                                                              grid2ElementType, polytopeCorners, hitCorners,
                                                              neighborIntersects1, neighborIntersects2, overlap_);

    // remove possible doubles
    removeDoubles(polytopeCorners);

    // Compute an interior point of the polytope
    int nPolyCorners = polytopeCorners.size();

    // If the polytope is degenerated then there is no intersection
    if (nPolyCorners<dimworld)
        return;

    // If the polytope is a simplex return it
    if (nPolyCorners==dim+1) {

     //   std::cout<<"Add intersection: 1\n";
        typename Base::RemoteSimplicialIntersection intersect;
        intersect.grid1Entities_[0] = grid1Index;
        intersect.grid2Entities_[0] = grid2Index;

        for (int j=0;j<dim+1; j++) {
            intersect.grid1Local_[0][j]=polytopeCorners[j][0];
            intersect.grid2Local_[0][j]=polytopeCorners[j][1];
        }
        this->intersections_.push_back(intersect);

        return;
    }

    // At this point we must have dimworld>=3

    ///////////////////////////////////////////////////////////////////////////////
    //  Compute a point in the middle of the polytope and order all corners cyclic
    //////////////////////////////////////////////////////////////////////////////

    Dune::array<LocalCoords,2> center;
    center[0] = 0; center[1] = 0;
    for (int i=0; i<nPolyCorners; i++) {
        center[0].axpy(1.0/nPolyCorners,polytopeCorners[i][0]);
        center[1].axpy(1.0/nPolyCorners,polytopeCorners[i][1]);
    }

    // Order cyclic
    std::vector<int> ordering;
    computeCyclicOrder(polytopeCorners,center[0],ordering);

    //////////////////////////////////////
    // Add intersections
    ////////////////////////////////

      //  std::cout<<"Add intersection: "<<polytopeCorners.size()<<"\n";
    for (size_t i=0; i<polytopeCorners.size(); i++) {

        typename Base::RemoteSimplicialIntersection intersect;
        intersect.grid1Entities_[0] = grid1Index;
        intersect.grid2Entities_[0] = grid2Index;

        for (int j=0;j<dim; j++) {
            intersect.grid1Local_[0][j]=polytopeCorners[ordering[(i+j)%nPolyCorners]][0];
            intersect.grid2Local_[0][j]=polytopeCorners[ordering[(i+j)%nPolyCorners]][1];
        }

        // last corner is the center for all intersections
        intersect.grid1Local_[0][dim]=center[0];
        intersect.grid2Local_[0][dim]=center[1];

        this->intersections_.push_back(intersect);
    }
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeCyclicOrder(const std::vector<Dune::array<LocalCoords,2> >& polytopeCorners,
                        const LocalCoords& center, std::vector<int>& ordering) const
{
    ordering.resize(polytopeCorners.size());

    for (size_t k=0; k<ordering.size(); k++)
        ordering[k] = k;

    //TODO Do I have to order triangles to get some correct orientation?
    if (polytopeCorners.size()<=3)
        return;

    // compute angles inside the polygon plane w.r.t to this axis
    LocalCoords  edge0 = polytopeCorners[0][0] - center;

    // Compute a vector that is perpendicular to the edge but lies in the polytope plane
    // So we have a unique ordering
    LocalCoords  edge1 = polytopeCorners[1][0] - center;
    LocalCoords normal0 = edge1;
    normal0.axpy(-(edge0*edge1),edge0);

    std::vector<T> angles(polytopeCorners.size());

    for (size_t i=0; i<polytopeCorners.size(); i++) {

        LocalCoords  edge = polytopeCorners[i][0] - center;

        T x(edge*edge0);
        T y(edge*normal0);

        angles[i] = std::atan2(y, x);
        if (angles[i]<0)
            angles[i] += 2*M_PI;
    }

    // bubblesort

    for (int i=polytopeCorners.size(); i>1; i--){
        bool swapped = false;

        for (int j=0; j<i-1; j++){

            if (angles[j] > angles[j+1]){
                swapped = true;
                std::swap(angles[j], angles[j+1]);
                std::swap(ordering[j], ordering[j+1]);
            }
        }

        if (!swapped)
            break;
    }
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::setupNodalDirections(const std::vector<WorldCoords>& coords1,
                              const std::vector<unsigned int>& elements1,
                              const std::vector<Dune::GeometryType>& elementTypes1,
                              const std::vector<WorldCoords>& coords2,
                              const std::vector<unsigned int>& elements2,
                              const std::vector<Dune::GeometryType>& elementTypes2)
{
    if (domainDirections_) {

        // Sample the provided analytical contact direction field
        nodalDomainDirections_.resize(coords1.size());
        for (size_t i=0; i<coords1.size(); i++)
            domainDirections_->evaluate(coords1[i], nodalDomainDirections_[i]);
    } else
        computeOuterNormalField(coords1,elements1,elementTypes1, nodalDomainDirections_);

    if (targetDirections_) {

        // Sample the provided analytical target direction field
        nodalTargetDirections_.resize(coords2.size());
        for (size_t i=0; i<coords2.size(); i++)
            targetDirections_->evaluate(coords2[i], nodalTargetDirections_[i]);
    } else
        computeOuterNormalField(coords2,elements2,elementTypes2, nodalTargetDirections_);
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::computeOuterNormalField(const std::vector<WorldCoords>& coords,
                                const std::vector<unsigned int>& elements,
                                const std::vector<Dune::GeometryType>& elementTypes,
                                std::vector<WorldCoords>& normals)
{
    normals.assign(coords.size(),WorldCoords(0));


    int offset = 0;

    for (size_t i=0; i<elementTypes.size(); i++) {

        int nCorners = Dune::ReferenceElements<T,dim>::general(elementTypes[i]).size(dim);

        // For segments 1, for triangles or quadrilaterals take the first 2
        std::vector<WorldCoords> edges(dim);
        for (int j=1; j<=dim; j++)
            edges[j-1] = coords[elements[offset + j]] - coords[elements[offset]];

        WorldCoords elementNormal;

        if (dim==1) {
            elementNormal[0] = edges[0][1]; elementNormal[1] = -edges[0][0];
        } else
            elementNormal = Projection::crossProduct(edges[0],edges[1]);

        elementNormal /= elementNormal.two_norm();

        for (int j=0; j<nCorners;j++)
            normals[elements[offset + j]] += elementNormal;

        offset += nCorners;
    }

    for (size_t i=0; i<coords.size(); i++)
        normals[i] /= normals[i].two_norm();
}

template<int dimworld, typename T>
void ContactMerge<dimworld, T>::removeDoubles(std::vector<Dune::array<LocalCoords,2> >& polytopeCorners)
{

    size_t counter(1);
    for (size_t i=1; i<polytopeCorners.size(); i++) {
        bool contained = false;
        for (size_t j=0; j<counter; j++)
            if ( (polytopeCorners[j][0]-polytopeCorners[i][0]).two_norm()<1e-10) {
                assert((polytopeCorners[j][1]-polytopeCorners[i][1]).two_norm()<1e-10);
                contained = true;
                break;
            }

        if (!contained) {
            if (counter < i)
                polytopeCorners[counter] = polytopeCorners[i];
            counter++;
        }
    }
    polytopeCorners.resize(counter);
}
