// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER_HH
#define DUNE_GRIDGLUE_COMMON_PROJECTIONHELPER_HH

#include <bitset>
#include <memory>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

/**
 *  \brief This namespace contains helper functions for the projection of two triangular surface on each other.
 */
namespace Projection
{
    //! Compute the cross product of two vectors
    template <class T, int dim>
        static Dune::FieldVector<T,dim> crossProduct(const Dune::FieldVector<T,dim>& a,
                const Dune::FieldVector<T,dim>& b)
        {
            if (dim!=3)
              DUNE_THROW(Dune::NotImplemented, "crossProduct does not work for dimension " << dim);

            Dune::FieldVector<T,dim> c;
            c[0] = a[1]*b[2] - a[2]*b[1];
            c[1] = a[2]*b[0] - a[0]*b[2];
            c[2] = a[0]*b[1] - a[1]*b[0];

            return c;
        }

    //! Helper class that provides static methods to compute the projection and inverse projection of a point along some given directions
    template <int dim, int dimworld, class T=double>
        class ProjectionHelper
        {
            typedef Dune::FieldVector<T,dimworld> WorldCoords;
            typedef Dune::FieldVector<T,dim> LocalCoords;

            public:

            /**  \brief Compute the inverse projection of a point onto some surface element where the projection is done across directions which are associated to corners of an surface element.
             *
             *  \param corners The coordinates of the corners.
             *  \param directions The directions along which the projection is done.
             *  \param target The point whose inverse projection is computed.
             *  \param preImage The pre-image of the target point in local coordinates of the surface element.
             *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
             *
             *  \return Returns true if the computed pre-image is within the convex hull of the corner points.
             */
            static bool inverseProjection(const std::vector<WorldCoords>& corners,
                    const std::vector<WorldCoords>& directions,
                    const WorldCoords& target, LocalCoords& preImage, const T overlap=1e-1) {

                DUNE_THROW(Dune::NotImplemented, "inverseProjection is not implemented for dimworld=="<<dimworld
                        <<" and dim=="<<dim);
                return false;
            }

            /**  \brief Compute the projection of a point along a given direction into the convex hull of some target points.
             *
             *  \param corner The coordinates of the point that is projected.
             *  \param direction The direction along which an intersection with the target surface element is searched.
             *  \param targetCorners The corner coordinates of the target surface element.
             *  \param image The projected corner in local coordinates of the target surface element.
             *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
             *
             *  \return Returns true if the computed image is within the convex hull of the target corner points.
             */
            static bool projection(const WorldCoords& corner, const WorldCoords& direction,
                    const std::vector<WorldCoords>& targetCorners,
                    LocalCoords& image, const T overlap=1e-1) {

                DUNE_THROW(Dune::NotImplemented, "projection is not implemented for dimworld=="<<dimworld
                        <<" and dim=="<<dim);
                return false;

            }
        };

    template <class T>
        class ProjectionHelper<2,3,T>
        {
            typedef Dune::FieldVector<T,3> WorldCoords;
            typedef Dune::FieldVector<T,2> LocalCoords;

            public:
            static bool projection(const WorldCoords& corner, const WorldCoords& direction,
                    const std::vector<WorldCoords>& targetCorners, LocalCoords& image,
                    const T overlap=1e-1)
            {

                if (targetCorners.size() != 3 && targetCorners.size() != 4)
                    DUNE_THROW(Dune::Exception, "projection is not implemented for polygone with "<<targetCorners.size()<<" vertices");

                // split up quadrilateral into triangles
                if (targetCorners.size() == 4) {

                    // lower triangle
                    std::vector<WorldCoords> triCorners(3);
                    for (int i=0; i<3; i++)
                        triCorners[i] = targetCorners[i];

                    if (projection(corner,direction,triCorners,image,overlap))
                        return true;

                    // upper triangle
                    triCorners[0] = targetCorners[3];

                    if (projection(corner,direction,triCorners,image,overlap)) {
                        // Transform local coordinates of the upper triangle to quadrilateral ones
                        T a = image[0];
                        image[0] = 1-image[1];
                        image[1] = 1-a;
                        return true;
                    }
                }


                //////////////////////////////////////////////
                //
                //  Solve the linear system:
                //
                //  (1-x-y)*t_0 + x*t_1 + y*t_2 + -z*direction = corner
                //
                //  with  (x,y) local coordinates in the target element
                //           z  the 'direction' scaling
                ////////////////////////////////////////////////////


                T eps = 1e-6;

                image = 0;

                Dune::FieldMatrix<T,3, 3> mat(0);
                std::array<WorldCoords, 2> directions;
                std::array<T, 2> scales;
                for (unsigned i = 0; i < 2; ++i) {
                    directions[i] = targetCorners[i+1] - targetCorners[0];
                    scales[i] = directions[i].infinity_norm();
                    directions[i] /= scales[i];
                }
                for (int i=0; i<3; i++) {
                    mat[i][0] = directions[0][i];
                    mat[i][1] = directions[1][i];
                    mat[i][2] = -direction[i];
                }

                // Solve the linear system
                WorldCoords rhs = corner - targetCorners[0];
                WorldCoords x(0);
                mat.solve(x,rhs);

                for (unsigned i = 0; i < 2; ++i)
                  x[i] /= scales[i];

                // only allow a certain overlap (we solved for '-z')
                if (x[2]<-overlap)
                    return false;

                if (x[0]<-eps || x[1]<-eps || (x[0] + x[1]>1+eps) )
                    return false;

                image[0] = x[0];
                image[1] = x[1];

                return true;
            }


            static bool inverseProjection(const std::vector<WorldCoords>& corners,
                    const std::vector<WorldCoords>& directions,
                    const WorldCoords& target, LocalCoords& preImage,
                    const T overlap = 1e-1)
            {
                const T eps = 1e-6;

                if (corners.size() != 3 && corners.size() != 4)
                    DUNE_THROW(Dune::NotImplemented, "inverseProjection is not implemented for elements with "<<corners.size()<<" corners!");

                // Split quadrilateral into triangles
                if (corners.size() == 4) {

                    // lower triangle
                    std::vector<WorldCoords> triCorners(3),triDirections(3);
                    for (int i=0; i<3; i++) {
                        triCorners[i] = corners[i];
                        triDirections[i] = directions[i];
                    }

                    if (inverseProjection(triCorners,triDirections,target,preImage,overlap))
                        return true;

                    // upper triangle
                    triCorners[0] = corners[3];
                    triDirections[0] = directions[3];

                    if (inverseProjection(triCorners,triDirections,target,preImage,overlap)) {
                        // Transform local coordinates to quadrilateral ones
                        T a=preImage[0];
                        preImage[0] = 1-preImage[1];
                        preImage[1] = 1-a;

                        return true;
                    }
                }

                // try to solve a cubic equation for the distance parameter, then compute the barycentric coordinates from it

                // the barycentric coordinates and the distance of the projected point
                WorldCoords x(3);

                // cubic coefficient
                WorldCoords n02 = directions[0] - directions[2];
                WorldCoords n12 = directions[1] - directions[2];
                WorldCoords n02n12 = crossProduct(n02,n12);

                T cubic = (directions[2]*n02n12);

                // quadratic coefficient

                WorldCoords p02 = corners[0] - corners[2];
                WorldCoords p12 = corners[1] - corners[2];
                WorldCoords p2q = corners[2] -target;
                WorldCoords p02n12 = crossProduct(p02,n12);
                WorldCoords n02p12 = crossProduct(n02,p12);

                T quadratic = (directions[2]*p02n12)+ (directions[2]*n02p12)+ (p2q*n02n12);

                // constant coefficient
                WorldCoords p02p12 = crossProduct(p02,p12);
                T constant = p2q*p02p12;

                // linear coefficient
                T linear = (directions[2]*p02p12) + (p2q*n02p12) + (p2q*p02n12);

                // save all zeros we find
                std::vector<T> zeros;

                if (std::fabs(cubic) <1e-10 && std::fabs(quadratic)<1e-10 && std::fabs(linear)<1e-10) {
                    return false;
                } else if (std::fabs(cubic) <1e-10 && std::fabs(quadratic)<1e-10) {

                    // problem is linear
                    zeros.push_back(-constant/linear);

                } else if(std::fabs(cubic)<1e-10) {

                    // problem is quadratic
                    T p = linear/quadratic;
                    T q = constant/quadratic;

                    T sqt = 0.25*p*p -q;

                    // no real solution
                    if (sqt<-1e-10)
                        return false;

                    zeros.push_back(-0.5*p + std::sqrt(sqt));
                    zeros.push_back(-0.5*p -std::sqrt(sqt));

                } else {

                    // problem is cubic
                    quadratic /= cubic;
                    linear /= cubic;
                    constant /= cubic;

                    // Transform to reduced form z^3 + p*z + q = 0 where x = z-quadratic/3
                    T p= linear - quadratic*quadratic/3;
                    T q=quadratic*(2*quadratic*quadratic/27 - linear/3) + constant;

                    // use Cardano's method to solve the problem
                    T D = 0.25*q*q + std::pow(p,3)/27;

                    if (D>1e-10) {
                        // one real zero

                        // be careful when computing the cubic roots
                        T nu = -q/2+std::sqrt(D);
                        T zer = std::pow(std::fabs(nu),1.0/3.0) * ((nu<-1e-10) ? -1 : 1);

                        nu = -q/2-std::sqrt(D);
                        zer += std::pow(std::fabs(nu),1.0/3.0) * ((nu<-1e-10) ? -1 : 1);

                        zeros.push_back(zer-quadratic/3);

                    } else if (D<-1e-10) {

                        // three real zeros, using trigonometric functions to compute them
                        T a = std::sqrt(-4*p/3);
                        T b = std::acos(-0.5*q*std::sqrt(-27/(std::pow(p,3))));

                        for (int i=0;i<3; i++)
                            zeros.push_back(std::pow(-1,i+1)*a*std::cos((b+(1-i)*M_PI)/3) -quadratic/3);


                    } else {
                        // one single and one double zero

                        if (std::fabs(q)<1e-10) {
                            zeros.push_back(-quadratic/3);

                            if (p<-1e-10)
                                zeros.push_back(std::sqrt(-p)-quadratic/3);

                        } else if (std::fabs(p)<1e-10) { // is this case correct?

                            T nu = std::pow(std::fabs(q),1.0/3.0) * ((q<-eps) ? -1 : 1);
                            zeros.push_back(nu-quadratic/3);

                        } else {
                            zeros.push_back(3*q/p - quadratic/3);
                            zeros.push_back(-1.5*q/p - quadratic/3);
                        }
                    }
                }

                int index = -1;
                WorldCoords r;

                for (size_t i=0;i<zeros.size();i++) {

                    T nu=zeros[i];
                    // only look in the direction of the outer normals
                    if (nu<-overlap) // allowed overlap
                        continue;

                    if (index != -1)
                        if (nu > zeros[index]) // is this one really closer ?
                            continue;

                    r[2] = nu;

                    // the computation of the other components might lead to nan or inf
                    // if this happens use a different equation to compute them
                    WorldCoords e = p2q; e.axpy(nu,directions[2]);
                    WorldCoords f = p12; f.axpy(nu,n12);
                    WorldCoords g = p02; g.axpy(nu,n02);
                    WorldCoords c = crossProduct(e,g);
                    WorldCoords d = crossProduct(g,f);

                    // computation of the other components is unstable
                    for (int j=0;j<3; j++) {

                        r[1] = c[j]/d[j];

                        if (isnan(r[1]) || isinf(r[1]))
                            continue;

                        r[0] = -(e[(j+1)%3]+r[1]*f[(j+1)%3])/g[(j+1)%3];

                        // Compute the residual
                        WorldCoords residual = p2q;
                        residual.axpy(r[0],p02);
                        residual.axpy(r[1],p12);
                        residual.axpy(r[2]*r[0],n02);
                        residual.axpy(r[2]*r[1],n12);
                        residual.axpy(r[2],directions[2]);


                        if (!(isnan(r[0]) || isinf(r[0])) && (residual.two_norm()<1e-5))
                            break;

                        r[0] = -(e[(j+2)%3]+r[1]*f[(j+2)%3])/g[(j+2)%3];

                        // Compute the residual
                        residual = p2q;
                        residual.axpy(r[0],p02);
                        residual.axpy(r[1],p12);
                        residual.axpy(r[2]*r[0],n02);
                        residual.axpy(r[2]*r[1],n12);
                        residual.axpy(r[2],directions[2]);

                        if (!(isnan(r[0]) || isinf(r[0])) && (residual.two_norm()<1e-5))
                            break;

                    }

                    if (r[0] > -eps && r[1]> -eps && (r[0]+r[1] < 1+eps)) {
                        index = i;
                        x = r;
                    }
                }

                // Check if we found a feasible zero
                WorldCoords residual = p2q;
                residual.axpy(x[0],p02);
                residual.axpy(x[1],p12);
                residual.axpy(x[2]*x[0],n02);
                residual.axpy(x[2]*x[1],n12);
                residual.axpy(x[2],directions[2]);

                if (residual.two_norm()<1e-6) {
                    if (index >= 0) {
                        preImage[0] = x[1];
                        preImage[1] = 1-x[0]-x[1];

                        return true;
                    }
                    return false;
                }


                //std::cout<<"Direct solution failed, use Newton method\n";
                // In some cases the direct solution of the cubic equation is unstable, in this case we use
                // Newton's method to compute at least one solution

                // Some problems have two solutions and the Newton converges to the wrong one
                return inexactInverseProjection(corners,directions, target, preImage,overlap);
            }

            /**
             *  \brief Compute the inverse projection using Newton's method.
             *         This is more stable but as the problem to solve is cubic we only compute one solution
             *         which might not be the one we are looking for, i.e. the one in the surface element.
             */
            static bool inexactInverseProjection(const std::vector<WorldCoords>& corners,
                    const std::vector<WorldCoords>& directions,
                    const WorldCoords& target, LocalCoords& preImage, const T overlap=1e-1)
            {
                assert(corners.size() == 3);
                assert(directions.size() == 3);

                // feasible initial Newton iterate
                const int nCorners = 3;
                Dune::FieldVector<T,nCorners> x(1.0/((T) nCorners));
                {
                    WorldCoords d(0);
                    for (std::size_t i = 0; i < corners.size(); ++i)
                        d += corners[i];
                    d *= 1./3;
                    d -= target;
                    x[2] = d.two_norm();
                }

                for (int i=0; i<30; i++) {

                    // compute Newton correction
                    WorldCoords Fxk =  target -corners[2];

                    for (size_t i=0; i<corners.size()-1; i++) {
                        Fxk.axpy(x[i],corners[2]-corners[i]);
                        Fxk.axpy(x[2]*x[i],directions[2]-directions[i]);
                    }
                    Fxk.axpy(-x[2],directions[2]);

                    Dune::FieldMatrix<T,nCorners,nCorners> FPrimexk(0);

                    FPrimexk[0] = corners[0] - corners[2];
                    FPrimexk[0].axpy(x[2],directions[0]-directions[2]);
                    FPrimexk[1] = corners[1] - corners[2];
                    FPrimexk[1].axpy(x[2],directions[1]-directions[2]);
                    FPrimexk[2] = directions[2];
                    FPrimexk[2].axpy(x[0],directions[0]-directions[2]);
                    FPrimexk[2].axpy(x[1],directions[1]-directions[2]);

                    try {
                        FPrimexk.invert();
                    }
                    catch (const Dune::FMatrixError&) {
                        return false;
                    }

                    WorldCoords newtonCorrection; // = (-1) * FPrimexk.inverse() * Fxk;

                    FPrimexk.mtv(Fxk, newtonCorrection);

                    x += newtonCorrection;
                }


                if (x[0]>-1e-6 && x[1]>-1e-6 && (x[0]+x[1] <1+1e-6)) {

                    WorldCoords residual = corners[2]-target;
                    residual.axpy(x[0],corners[0]-corners[2]);
                    residual.axpy(x[1],corners[1]-corners[2]);
                    residual.axpy(x[2]*x[0],directions[0]-directions[2]);
                    residual.axpy(x[2]*x[1],directions[1]-directions[2]);
                    residual.axpy(x[2],directions[2]);

                    // Newton did not converge
                    if (residual.two_norm()>1e-6 || x[2]<-overlap)
                        return false;

                    preImage[0] = x[1];
                    preImage[1] = 1-x[0]-x[1];

                    return true;
                }

                return false;
            }
        };

    template <class T>
        class ProjectionHelper<1,2,T>
        {

            typedef Dune::FieldVector<T,2> WorldCoords;
            typedef Dune::FieldVector<T,1> LocalCoords;

            public:
            static bool projection(const WorldCoords& corner, const WorldCoords& direction,
                    const std::vector<WorldCoords>& targetCorners, LocalCoords& image, const T overlap=1e-1)
            {
                T eps = 1e-8;
                // we solve the equation basePoint + x_0 * normal = a + x_1 * (b-a)
                image = 0;

                Dune::FieldMatrix<T,2,2> mat;
                mat[0][0] = direction[0];
                mat[1][0] = direction[1];
                mat[0][1] = targetCorners[0][0]-targetCorners[1][0];
                mat[1][1] = targetCorners[0][1]-targetCorners[1][1];

                WorldCoords rhs = targetCorners[0] - corner;
                WorldCoords x(0);

                // Solve the system.  If it is singular the normal and the segment
                // are parallel and there is no intersection

                T detinv = mat[0][0]*mat[1][1]-mat[0][1]*mat[1][0];
                if (std::abs(detinv)<1e-14)
                    return false;

                detinv = 1/detinv;

                x[0] = detinv*(mat[1][1]*rhs[0]-mat[0][1]*rhs[1]);
                x[1] = detinv*(mat[0][0]*rhs[1]-mat[1][0]*rhs[0]);

                // x[0] is the distance, x[1] is the intersection point
                // in local coordinates on the segment
                if (x[1]<-eps || x[1] > 1+eps)
                    return false;

                // allow some overlap
                if (x[0]<-overlap)
                    return false;

                image[0] = x[1];

                return true;
            }

            static bool inverseProjection(const std::vector<WorldCoords>& corners,
                    const std::vector<WorldCoords>& directions,
                    const WorldCoords& target, LocalCoords& preImage,
                    const T overlap=1e-1)
            {
                T eps = 1e-8;

                preImage = 0;
                WorldCoords p10 = corners[1]-corners[0];
                WorldCoords v10 = directions[1]-directions[0];
                WorldCoords qp0 = target-corners[0];

                // Compute coefficients of the polynomial by eliminating the distance variable
                T a = p10[1]*v10[0] - p10[0]*v10[1];

                T b = -qp0[1]*v10[0] +p10[1]*directions[0][0]
                    + qp0[0]*v10[1] -p10[0]*directions[0][1];

                T c = -qp0[1]*directions[0][0] +qp0[0]*directions[0][1];

                // Is the quadratic formula degenerated to a linear one?
                if (std::abs(a) < 1e-10) {
                    T x = -c/b;

                    if (x >= -eps && x <= 1+eps) {

                        //check if projection is along the positive directions
                        T dist = (qp0[0]-x*p10[0])/(directions[0][0]+x*v10[0]);
                        if (dist <-overlap)
                            return false;

                        preImage[0] = x;
                        return true;
                    }
                    return false;
                }

                // The abc-formula
                T mu_0 = (-b + std::sqrt(b*b - 4*a*c))/(2*a);
                T mu_1 = (-b - std::sqrt(b*b - 4*a*c))/(2*a);

                if (mu_0 >= -eps && mu_0 <= 1+eps) {

                    T dist = (qp0[0]-mu_0*p10[0])/(directions[0][0]+mu_0*v10[0]);
                    if (dist <-overlap)
                        return false;

                    preImage[0] = mu_0;
                    return true;

                } else if (mu_1 >= -eps && mu_1 <= 1+eps) {

                    T dist = (qp0[0]-mu_1*p10[0])/(directions[0][0]+mu_1*v10[0]);
                    if (dist <-overlap)
                        return false;

                    preImage[0] = mu_1;
                    return true;
                }
                return false;

            }
        };
    /**
     *  \brief Compute the projection of a point along a given direction into the convex hull of some target points.
     *
     *  \param corner The coordinates of the point that is projected.
     *  \param direction The direction along which an intersection with the target surface element is searched.
     *  \param targetCorners The corner coordinates of the target surface element.
     *  \param image The projected corner in local coordinates of the target surface element.
     *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
     *
     *  \return Returns true if the computed image is within the convex hull of the target corner points.
     */

    template <int dim, int dimworld, class T>
        static bool projection(const Dune::FieldVector<T,dimworld>& corner, const Dune::FieldVector<T,dimworld>& direction,
                const std::vector<Dune::FieldVector<T,dimworld> >& targetCorners, Dune::FieldVector<T,dim>& image,
                const T overlap=1e-1)
        {
            return ProjectionHelper<dim,dimworld,T>::projection(corner,direction,
                    targetCorners,image,overlap);
        }

    /**
     *  \brief Compute the inverse projection of a point onto some surface element where the projection is done across directions which are associated to corners of an surface element.
     *
     *  \param corners The coordinates of the corners.
     *  \param directions The directions along which the projection is done.
     *  \param target The point whose inverse projection is computed.
     *  \param preImage The pre-image of the target point in local coordinates of the surface element.
     *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
     *
     *  \return Returns true if the computed pre-image is within the convex hull of the corner points.
     */
    template <int dim, int dimworld, class T>
        static bool inverseProjection(const std::vector<Dune::FieldVector<T,dimworld> >& corners,
                const std::vector<Dune::FieldVector<T,dimworld> >& directions,
                const Dune::FieldVector<T,dimworld>& target, Dune::FieldVector<T,dim>& preImage,
                const T overlap=1e-1)
        {
            return ProjectionHelper<dim,dimworld,T>::inverseProjection(corners,directions,
                    target, preImage, overlap);
        }
    /** \brief Helper class that provides static methods for the computation of the intersection of surface element edges projected onto each other. */
    template <int dim, int dimworld, class T>
        class EdgeIntersectionHelper
        {

            typedef Dune::FieldVector<T,dimworld> WorldCoords;
            typedef Dune::FieldVector<T,dim> LocalCoords;

            public:
            /**
             * \brief Compute the projection along given directions of surface element edges onto target edges.
             *
             *  \param corners1 The coordinates of the surface element corners whose edges are projected.
             *  \param corners2 The coordinates of the surface element corners on which is projected.
             *  \param directions1 The directions along which the projection is done.
             *  \param gt1 The geometry type of the projected surface element.
             *  \param gt2 The geometry type of the target surface element.
             *  \param polygonCorners If intersection points are found their local coordinates are added to this vector.
             *  \param hitCorners Vector containing information on which surface element corners are projected on each other.
            *  \param neighborIntersects1 If two edges intersect then the corresponding surface element neighbors also intersect. This information for the projected surface element is stored in the bitfield.
            *  \param neighborIntersects2 If two edges intersect then the corresponding surface element neighbors also intersect. This information for the surface element on which is projected is stored in the bitfield.
             *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
             */
            static void addEdgeIntersections(const std::vector<WorldCoords>& corners1,
                    const std::vector<WorldCoords>& corners2,
                    const std::vector<WorldCoords>& directions1,
                    const Dune::GeometryType& gt1, Dune::GeometryType& gt2,
                    std::vector<std::array<LocalCoords,2> >& polygonCorners,
                    const std::vector<int>& hitCorners, std::bitset<(1<<dim)>& neighborIntersects1,
                    std::bitset<(1<<dim)>& neighborIntersects2, const T overlap = 1e-1)
            {
                DUNE_THROW(Dune::NotImplemented, "addEdgeIntersections is not implemented for dimworld=="<<dimworld
                        <<" and dim=="<<dim);
            }
        };

    template <typename T>
        class EdgeIntersectionHelper<1, 2, T>
        {

            typedef Dune::FieldVector<T,2> WorldCoords;
            typedef Dune::FieldVector<T,1> LocalCoords;

            public:
            //! For 1D surfaces there are no edges.
            static void addEdgeIntersections(const std::vector<WorldCoords>& corners1,
                    const std::vector<WorldCoords>& corners2,
                    const std::vector<WorldCoords>& directions1,
                    const Dune::GeometryType& gt1, const Dune::GeometryType& gt2,
                    std::vector<std::array<LocalCoords,2> >& polygonCorners,
                    const std::vector<int>& hitCorners, std::bitset<2>& neighborIntersects1,
                    std::bitset<2>& neighborIntersects2, const T overlap=1e-1)
            {}

        };

    template <typename T>
        class EdgeIntersectionHelper<2, 3, T>
        {

            typedef Dune::FieldVector<T,3> WorldCoords;
            typedef Dune::FieldVector<T,2> LocalCoords;

            public:
            static void addEdgeIntersections(const std::vector<WorldCoords>& corners1,
                    const std::vector<WorldCoords>& corners2,
                    const std::vector<WorldCoords>& directions1,
                    const Dune::GeometryType& gt1, const Dune::GeometryType& gt2,
                    std::vector<std::array<LocalCoords,2> >& polygonCorners,
                    const std::vector<int>& hitCorners, std::bitset<4>& neighborIntersects1,
                    std::bitset<4>& neighborIntersects2, const T overlap=1e-1)
            {
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
                const Dune::ReferenceElement<T,2>& ref1 = Dune::ReferenceElements<T,2>::general(gt1);
                const Dune::ReferenceElement<T,2>& ref2 = Dune::ReferenceElements<T,2>::general(gt2);
#else
                const Dune::GenericReferenceElement<T,2>& ref1 = Dune::GenericReferenceElements<T,2>::general(gt1);
                const Dune::GenericReferenceElement<T,2>& ref2 = Dune::GenericReferenceElements<T,2>::general(gt2);
#endif

                //loop  over edges
                std::array<Dune::FieldVector<T,1>,2> intersection;
                for (int i=0; i<ref1.size(1); i++) {

                    std::vector<int> edgeCorners1(2);
                    edgeCorners1[0] = ref1.subEntity(i,1,0,2);
                    edgeCorners1[1] = ref1.subEntity(i,1,1,2);

                    int nIntersects(0);
                    for (int j=0; j<ref2.size(1); j++) {

                        std::vector<int> edgeCorners2(2);
                        edgeCorners2[0] = ref2.subEntity(j,1,0,2);
                        edgeCorners2[1] = ref2.subEntity(j,1,1,2);

                        // if any edge endpoints hit each other we don't have to compute the intersections
                        // either there is none or we already have computed it before
                        if (hitCorners[edgeCorners2[0]]==edgeCorners1[0] || hitCorners[edgeCorners2[0]]==edgeCorners1[1]
                                || hitCorners[edgeCorners2[1]]==edgeCorners1[0] || hitCorners[edgeCorners2[1]]==edgeCorners1[1])
                            continue;

                        if (edgeIntersection(corners1[edgeCorners1[0]],corners1[edgeCorners1[1]],
                                    corners2[edgeCorners2[0]],corners2[edgeCorners2[1]],
                                    directions1[edgeCorners1[0]],directions1[edgeCorners1[1]],
                                    intersection, overlap) )
                        {
                            nIntersects++;

                            // compute the local coordinates
                            std::array<LocalCoords,2> corner;
                            corner[0] = ref1.template geometry<1>(i).global(intersection[0]);
                            corner[1] = ref2.template geometry<1>(j).global(intersection[1]);
                            polygonCorners.push_back(corner);

                            // if the grid2 edge intersected then the neighbor will also intersect
                            neighborIntersects2[j] = true;
                        }
                    }

                    // if grid1 edge intersected an other edge then the neighbor will also intersect
                    if (nIntersects>0)
                        neighborIntersects1[i] = true;
                }
            }

            private:
            /**
             * \brief Compute the projection along given directions of a surface element edge onto another edge.
             *
             *  \param corner1 The coordinates of the first edge corner whose edge is projected.
             *  \param corner2 The coordinates of the second edge corner whose edge is projected.
             *  \param target1 The coordinates of the first edge corner on whose edge is projected.
             *  \param target2 The coordinates of the second edge corner on whose edge is projected.
             *  \param direction1 The direction corresponding to the first edge corner along which is projected.
             *  \param direction2 The direction corresponding to the second edge corner along which is projected.
             *  \param intersection If an intersection points is found its local coordinates stored in this vector.
             *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
             *
             *  \returns Returns true if an intersection point is found.
             */
            static bool edgeIntersection(const WorldCoords& corner1, const WorldCoords& corner2,
                    const WorldCoords& target1, const WorldCoords& target2,
                    const WorldCoords& direction1, const WorldCoords& direction2,
                    std::array<Dune::FieldVector<T,1>,2>& intersection, const T overlap = 1e-1)
            {

                T eps = 1e-6;
                // solve a quadratic scalar equation for the distance parameter eta, then compute the barycentric coordinates from it

                WorldCoords n21 = direction2 - direction1;
                WorldCoords p21 = corner2 - corner1;
                WorldCoords q21 = target2 - target1;
                WorldCoords q21n21 = crossProduct(q21,n21);
                WorldCoords q21p21 = crossProduct(q21,p21);
                WorldCoords p1q1 = corner1 -target1;

                // quadratic coefficient
                T quadratic = direction1*q21n21;

                // linear coefficient
                T linear = (direction1*q21p21) + (p1q1*q21n21);

                // constant coefficient
                T constant = p1q1*q21p21;

                // save all zeros we find
                std::vector<T> zeros;

                if (std::fabs(quadratic)<1e-10 && std::fabs(linear)<1e-10) {
                    return false;
                } else if (std::fabs(quadratic)<1e-10) {

                    // problem is linear
                    zeros.push_back(-constant/linear);

                } else {

                    // problem is quadratic
                    T p = linear/quadratic;
                    T q = constant/quadratic;

                    T sqt = 0.25*p*p -q;

                    // no real solution
                    if (sqt<-1e-10)
                        return false;

                    zeros.push_back(-0.5*p + std::sqrt(sqt));
                    zeros.push_back(-0.5*p -std::sqrt(sqt));

                }

                int index = -1;
                WorldCoords r,x(-1);

                for (size_t i=0;i<zeros.size();i++) {

                    T eta=zeros[i];

                    // only look in the direction of the outer normals
                    if (eta<-overlap)
                        continue;

                    r[2] = eta;

                    // the computation of the other components might lead to nan or inf
                    // if this happens use a different equation to compute them
                    WorldCoords dummy = p1q1; dummy.axpy(eta,direction1);
                    WorldCoords c =crossProduct(dummy,q21);
                    dummy = p21; dummy.axpy(eta,n21);
                    WorldCoords d =crossProduct(q21,dummy);

                    for (int j=0;j<3; j++) {

                        r[0] = c[j]/d[j];
                        if (isnan(r[0]) || isinf(r[0]))
                            continue;

                        r[1] = (p1q1[(j+1)%3]+eta*direction1[(j+1)%3] + r[0]*(p21[(j+1)%3]+eta*n21[(j+1)%3]))/q21[(j+1)%3];

                        // computation of the other components can be instable
                        WorldCoords residual = p1q1;
                        residual.axpy(r[0],p21); residual.axpy(r[2],direction1);
                        residual.axpy(r[2]*r[0],n21); residual.axpy(-r[1],q21);

                        if (!(isnan(r[1]) || isinf(r[1])) && residual.two_norm()<1e-5)
                            break;

                        residual.axpy(r[1],q21);
                        r[1] = (p1q1[(j+2)%3]+eta*direction1[(j+2)%3] + r[0]*(p21[(j+2)%3]+eta*n21[(j+2)%3]))/q21[(j+2)%3];

                        residual.axpy(r[1],q21);
                        // computation of the other components can be instable
                        if (!(isnan(r[1]) || isinf(r[1])) && residual.two_norm()<1e-5)
                            break;

                    }
                    if (r[0] >= -eps && r[1]>= -eps && (r[0]<=1+eps)  && (r[1] <= 1+eps)) {
                        index = i;
                        x=r;
                    }
                }

                // TODO CLEAN THIS UP
                WorldCoords residual = p1q1;
                residual.axpy(x[0],p21); residual.axpy(x[2],direction1);
                residual.axpy(x[2]*x[0],n21); residual.axpy(-x[1],q21);

                if (residual.two_norm()<eps)
                {
                    if (index >= 0 && x[0] >= -eps && x[1]>= -eps && (x[0]<=1+eps)  && (x[1] <= 1+eps)) {

                        intersection[0] = x[0];
                        intersection[1] = x[1];
                        return true;
                    }
                    return false;
                }

                // if the direct compuation failed, use a Newton method to compute at least one zero

                // Fix some initial value
                // sometimes it only works when the initial value is an intersection...
                x[0] = x[1] = 0.5;
                x[2] = 0.5;
                WorldCoords newtonCorrection(0);

                for (int i=0; i<30; i++) {

                    // compute Newton correction

                    WorldCoords Fxk = target1 -corner1;
                    Fxk.axpy(-x[0], p21);
                    Fxk.axpy(-x[2],direction1);
                    Fxk.axpy(-x[2]*x[0],n21);
                    Fxk.axpy(x[1],q21);

                    Dune::FieldMatrix<T,3,3> FPrimexk;
                    FPrimexk[0] = p21;
                    FPrimexk[0].axpy(x[2],n21);
                    FPrimexk[1] = target1-target2;
                    FPrimexk[2] = direction1;
                    FPrimexk[2].axpy(x[0],n21);
                    if (FPrimexk.determinant()<1e-10)
                        return false;

                    FPrimexk.invert();

                    FPrimexk.mtv(Fxk, newtonCorrection);

                    x += newtonCorrection;
                }

                residual = p1q1;
                residual.axpy(x[0],p21); residual.axpy(x[2],direction1);
                residual.axpy(x[2]*x[0],n21); residual.axpy(-x[1],q21);

                if (residual.two_norm()<=eps) {

                    if (x[0]>=-eps && x[0]<=(1+eps) && x[1]>=-eps && x[1]<=(1+eps)) {
                        if (x[2]<-overlap)
                            return false;

                        intersection[0] = x[0];
                        intersection[1] = x[1];

                        return true;
                    }
                    return false;

                }

                //std::cout<<"Newton did not converge either!\n";

                return false;
            }
        };

    /**
     * \brief Compute the projection along given directions of surface element edges onto target edges.
     *
     *  \param corners1 The coordinates of the surface element corners whose edges are projected.
     *  \param corners2 The coordinates of the surface element corners on which is projected.
     *  \param directions1 The directions along which the projection is done.
     *  \param gt1 The geometry type of the projected surface element.
     *  \param gt2 The geometry type of the target surface element.
     *  \param polygonCorners If intersection points are found their local coordinates are added to this vector.
     *  \param hitCorners Vector containing information on which surface element corners are projected on each other.
     *  \param neighborIntersects1 If two edges intersect then the corresponding surface element neighbors also intersect. This information for the projected surface element is stored in the bitfield.
     *  \param neighborIntersects2 If two edges intersect then the corresponding surface element neighbors also intersect. This information for the surface element on which is projected is stored in the bitfield.
     *  \param overlap The amount of overlap that is allowed, i.e. projection among the opposite direction is valid if the scaling is smaller than overlap.
     */
    template <int dim, int dimworld, class T=double>
        static void addEdgeIntersections(const std::vector<Dune::FieldVector<T,dimworld> >& corners1,
                const std::vector<Dune::FieldVector<T,dimworld> >& corners2,
                const std::vector<Dune::FieldVector<T,dimworld> >& directions1,
                const Dune::GeometryType& gt1, const Dune::GeometryType& gt2,
                std::vector<std::array<Dune::FieldVector<T,dim>,2> >& polygonCorners,
                const std::vector<int>& hitCorners, std::bitset<(1<<dim)>& neighborIntersects1,
                std::bitset<(1<<dim)>& neighborIntersects2, const T overlap = 1e-1)
        {
            EdgeIntersectionHelper<dim,dimworld,T>::addEdgeIntersections(corners1,corners2,
                    directions1, gt1, gt2, polygonCorners,
                    hitCorners, neighborIntersects1,
                    neighborIntersects2, overlap);
        }

} // end namespace
#endif
