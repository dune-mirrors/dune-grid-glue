// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_GLUE_MIXED_DIM_OVERLAPPING_MERGE_HH
#define DUNE_GRID_GLUE_MIXED_DIM_OVERLAPPING_MERGE_HH

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/common/version.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/multilineargeometry.hh>

#include <dune/grid/common/grid.hh>

#include <dune/grid-glue/merging/standardmerge.hh>

namespace Dune {

  namespace GridGlue {

/** \brief Computing overlapping grid intersections for grids of different dimensions

   \tparam dim1 Grid dimension of grid 1
   \tparam dim2 Grid dimension of grid 2
   \tparam dimworld World dimension
   \tparam T Type used for coordinates
 */
template<int dim1, int dim2, int dimworld, typename T = double>
class MixedDimOverlappingMerge
  : public StandardMerge<T,dim1,dim2,dimworld>
{

public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  //typedef Dune::FieldVector<T, dim>  LocalCoords;

  MixedDimOverlappingMerge()
  {}

private:

  typedef typename StandardMerge<T,dim1,dim2,dimworld>::RemoteSimplicialIntersection RemoteSimplicialIntersection;

  /** \brief Compute the intersection between two overlapping elements

     The result is a set of simplices.

     \param grid1ElementType Type of the first element to be intersected
     \param grid1ElementCorners World coordinates of the corners of the first element

     \param grid2ElementType Type of the second element to be intersected
     \param grid2ElementCorners World coordinates of the corners of the second element

   */
  void computeIntersection(const Dune::GeometryType& grid1ElementType,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid1ElementCorners,
                           unsigned int grid1Index,
                           std::bitset<(1<<dim1)>& neighborIntersects1,
                           const Dune::GeometryType& grid2ElementType,
                           const std::vector<Dune::FieldVector<T,dimworld> >& grid2ElementCorners,
                           unsigned int grid2Index,
                           std::bitset<(1<<dim2)>& neighborIntersects2)
  {
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
    typedef Dune::CachedMultiLinearGeometry<T, dim1, dimworld> Geometry1;
    typedef Dune::CachedMultiLinearGeometry<T, dim2, dimworld> Geometry2;

    Geometry1 grid1Geometry(grid1ElementType, grid1ElementCorners);
    Geometry2 grid2Geometry(grid2ElementType, grid2ElementCorners);

    if (dim1==1 and dim2==2 and dimworld==2)
    {
      assert(grid1ElementCorners.size() == 2);  // line segment

      // Compute intersections of segment with all edges of the 2d element
      std::vector<Dune::FieldVector<T,dimworld> > intersections;
      // Used to check if itersection in an element was already found
      bool sameElement = false;

      for (int i=0; i<refElement2.size(1); i++)
      {
        Dune::array<Dune::FieldVector<T,dimworld>, 2> edge;
        edge[0] = grid2ElementCorners[refElement2.subEntity(i,dim2-1,0,dim2)];
        edge[1] = grid2ElementCorners[refElement2.subEntity(i,dim2-1,1,dim2)];

        segmentSegmentIntersection2D(grid1ElementCorners,edge,intersections, sameElement);
      }

      // Compute whether the segment endpoints are contained in the triangle
      std::bitset<2> inTriangle("00");

      for (size_t i=0; i<2; i++)
#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
        inTriangle[i] = ReferenceElements<T,dim2>::general(grid2ElementType).checkInside(grid2Geometry.local(grid1ElementCorners[i]));
#else
        inTriangle[i] = Dune::GenericReferenceElements<T,dim2>::general(grid2ElementType).checkInside(grid2Geometry.local(grid1ElementCorners[i]));
#endif

      // Everything is easy if both segment endpoints are contained in the 2d element
      if (inTriangle[0] and inTriangle[1])
      {
        this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0][0] = 0;
        this->intersections_.back().grid1Local_[0][1] = 1;

        // Compute local coordinates in the grid2 element
        this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(grid1ElementCorners[0]);
        this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(grid1ElementCorners[1]);

        return;
      }

      switch (intersections.size()) {

          case 1:
          {
            if (inTriangle[0])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0][0] = 0;
              this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[0]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(grid1ElementCorners[0]);
              this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[0]);

            } else if (inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
              this->intersections_.back().grid1Local_[0][1] = 1;

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
              this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(grid1ElementCorners[1]);

            }

            return;
          }


          case 2:
          case 3:
          {
            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
            this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[1]);

            // Compute local coordinates in the grid2 element
            this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
            this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[1]);

            return;
          }

      }

    } else if (dim1==2 and dim2==1 and dimworld==2)
    {
      assert(grid2ElementCorners.size() == 2);  // line segment

      // Compute intersections of segment with all edges of the 2d element
      std::vector<Dune::FieldVector<T,dimworld> > intersections;
      bool sameElement = false;
      for (int i=0; i<refElement1.size(1); i++)
      {
        Dune::array<Dune::FieldVector<T,dimworld>, 2> edge;
        edge[0] = grid1ElementCorners[refElement1.subEntity(i,dim1-1,0,dim1)];
        edge[1] = grid1ElementCorners[refElement1.subEntity(i,dim1-1,1,dim1)];
        segmentSegmentIntersection2D(grid2ElementCorners,edge,intersections, sameElement);
      }

      // Compute whether the segment endpoints are contained in the triangle
      std::bitset<2> inTriangle("00");

      for (size_t i=0; i<2; i++)

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
        inTriangle[i] = ReferenceElements<T,dim1>::general(grid1ElementType).checkInside(grid1Geometry.local(grid2ElementCorners[i]));
#else
        inTriangle[i] = Dune::GenericReferenceElements<T,dim1>::general(grid1ElementType).checkInside(grid1Geometry.local(grid2ElementCorners[i]));
#endif

      // Everything is easy if both segment endpoints are contained in the 2d element
      if (inTriangle[0] and inTriangle[1])
      {
        this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(grid2ElementCorners[0]);
        this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(grid2ElementCorners[1]);

        // Compute local coordinates in the grid2 element
        this->intersections_.back().grid2Local_[0][0] = 0;
        this->intersections_.back().grid2Local_[0][1] = 1;

        return;
      }

      switch (intersections.size()) {

          case 1:
          {
            if (inTriangle[0])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(grid2ElementCorners[0]);
              this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[0]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0][0] = 0;
              this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[0]);

            } else if (inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
              this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(grid2ElementCorners[1]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
              this->intersections_.back().grid2Local_[0][1] = 1;

            }

            return;
          }

          case 2:
          case 3:
          {
            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
            this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[1]);

            // Compute local coordinates in the grid2 element
            this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
            this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[1]);

            return;
          }
      }
    }
    else if (dim1==3 and dim2==1 and dimworld==3)
      {
        assert(grid2ElementCorners.size() == 2);  // line segment

      // Compute intersections of segment with all edges of the 2d element
      std::vector<Dune::FieldVector<T,dimworld> > intersections;
      bool sameElement = false;
      double eps = 10e-9;

      for (int i=0; i<refElement1.size(1); i++)
        {
          Dune::array<Dune::FieldVector<T,dimworld>, 4> face;
          face[0] = grid1ElementCorners[refElement1.subEntity(i,dim1-2,0,dim1)];
          face[1] = grid1ElementCorners[refElement1.subEntity(i,dim1-2,1,dim1)];
          face[2] = grid1ElementCorners[refElement1.subEntity(i,dim1-2,2,dim1)];
          face[3] = grid1ElementCorners[refElement1.subEntity(i,dim1-2,3,dim1)];

          segmentPlaneIntersection3D(grid2ElementCorners,face,intersections,sameElement);
        }

      // Compute whether the segment endpoints are contained in the triangle
      std::bitset<2> inCube("00");

      for (size_t i=0; i<2; i++)
        inCube[i] = Dune::ReferenceElements<T,dim1>::general(grid1ElementType).checkInside(grid1Geometry.local(grid2ElementCorners[i]));

      // Everything is easy if both segment endpoints are contained in the 2d element
      if (inCube[0] and inCube[1])
      {
        this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(grid2ElementCorners[0]);
        this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(grid2ElementCorners[1]);

        // Compute local coordinates in the grid2 element
        this->intersections_.back().grid2Local_[0][0] = 0;
        this->intersections_.back().grid2Local_[0][1] = 1;

        return;
      }

      switch (intersections.size()) {

      case 1:
        {
          if( (std::abs(intersections[0][0]) > eps) &&
              (std::abs(intersections[0][1]) > eps) &&
              (std::abs(intersections[0][2]) > eps) ) {

            if (inCube[0])
              {
                this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(grid2ElementCorners[0]);
              this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[0]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0][0] = 0;
              this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[0]);
              }
            else if (inCube[1])
              {
                this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

                // Compute local coordinates in the grid1 element
                this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
                this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(grid2ElementCorners[1]);

                // Compute local coordinates in the grid2 element
                this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
                this->intersections_.back().grid2Local_[0][1] = 1;
              }
          }

          return;
        }

      case 2:
      case 3:
        {
          if( (std::abs(intersections[0][0]) > eps) &&
              (std::abs(intersections[0][1]) > eps) &&
              (std::abs(intersections[0][2]) > eps) &&
              (std::abs(intersections[1][0]) > eps) &&
              (std::abs(intersections[1][1]) > eps) &&
              (std::abs(intersections[1][2]) > eps) ) {

            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(intersections[0]);
            this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(intersections[1]);

            // Compute local coordinates in the grid2 element
            this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(intersections[0]);
            this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(intersections[1]);

            return;
          }
        }
      }
      }
    else if (((dim1==3 and dim2 == 2) || (dim1==2 and dim2 == 3)) and dimworld == 3) {


        typedef typename std::vector<Dune::Empty>::size_type size_type;

        assert(grid1ElementCorners.size() == 4);  // tetrahedron
        assert(grid2ElementCorners.size() == 3);  // triangle

        std::vector<Dune::FieldVector<T,dimworld> >  P;              // boundary Intersection points

        if (dim1 == 3 && dim2==2)
            intersections3D(grid1ElementCorners ,grid2ElementCorners, P);
        else
            intersections3D(grid2ElementCorners ,grid1ElementCorners, P);

        if (P.size() == 3) {

            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(P[0]);
            this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(P[1]);
            this->intersections_.back().grid1Local_[0][2] = grid1Geometry.local(P[2]);

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(P[0]);
            this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(P[1]);
            this->intersections_.back().grid2Local_[0][2] = grid2Geometry.local(P[2]);

        } else if (P.size() > 3) {
            Dune::FieldVector<T,dimworld> centroid(0.0);
            std::vector<int> no;

            for (size_type i=0; i < P.size(); i++)
                centroid += P[i] ;
            centroid /= P.size();

            orderPoints3D(centroid,no,P);

            for (size_type i=0; i < P.size() ; i++) {

                this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

                // Compute local coordinates in the grid1 element
                this->intersections_.back().grid1Local_[0][0] = grid1Geometry.local(P[no[i]]);
                this->intersections_.back().grid1Local_[0][1] = grid1Geometry.local(P[no[(i+1)%(P.size())]]);
                this->intersections_.back().grid1Local_[0][2] = grid1Geometry.local(centroid);

                // Compute local coordinates in the grid2 element
                this->intersections_.back().grid2Local_[0][0] = grid2Geometry.local(P[no[i]]);
                this->intersections_.back().grid2Local_[0][1] = grid2Geometry.local(P[no[(i+1)%(P.size())]]);
                this->intersections_.back().grid2Local_[0][2] = grid2Geometry.local(centroid);

            }
        }




    } else
      std::cout << "MixedDimOverlappingMerge for dim1 = " << dim1 << ", dim2 = "
                                                          << dim2 << ", dimworld = "
                                                          << dimworld << ": Not implemented yet!" << std::endl;
  }

  static void segmentSegmentIntersection2D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                            const Dune::array<Dune::FieldVector<T,dimworld>, 2>& edge,
                                            std::vector<Dune::FieldVector<T,dimworld> > & P ,
                                            bool& sameElement);


  static void segmentPlaneIntersection3D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                          const Dune::array<Dune::FieldVector<T,dimworld>, 4>& face,
                                          std::vector<Dune::FieldVector<T,dimworld> > & P ,
                                          bool& sameElement);

  static void intersections3D( const std::vector<Dune::FieldVector<T,dimworld> >   X,
                                               const std::vector<Dune::FieldVector<T,dimworld> >   Y,
                                               std::vector<Dune::FieldVector<T,dimworld> > & P );

  static bool segmentSegmentIntersection3D( const Dune::FieldVector<T,dimworld>    X0,
                                      const Dune::FieldVector<T,dimworld>    X1,
                                      const Dune::FieldVector<T,dimworld>    Y0,
                                      const Dune::FieldVector<T,dimworld>    Y1,
                                      Dune::FieldVector<T,dimworld>   & p);

  static bool triangleLineIntersection3D( const Dune::FieldVector<T,dimworld>    X0,
                                          const Dune::FieldVector<T,dimworld>    X1,
                                          const Dune::FieldVector<T,dimworld>    Y0,
                                          const Dune::FieldVector<T,dimworld>    Y1,
                                          const Dune::FieldVector<T,dimworld>    Y2,
                                          Dune::FieldVector<T,dimworld>   & p);

  static bool pointOnSegment3D(const Dune::FieldVector<T,dimworld> X,
                               const Dune::FieldVector<T,dimworld> Y0,
                               const Dune::FieldVector<T,dimworld> Y1);

  static bool pointInTriangle3D(const Dune::FieldVector<T,dimworld> X,
                                const Dune::FieldVector<T,dimworld> Y1,
                                const Dune::FieldVector<T,dimworld> Y2,
                                const Dune::FieldVector<T,dimworld> Y3);

  static bool pointInTetrahedra3D( const Dune::FieldVector<T,dimworld> X,
                                                      const std::vector<Dune::FieldVector<T,dimworld> >   Y);

  static int insertPoint3D( const Dune::FieldVector<T,dimworld>  p, std::vector<Dune::FieldVector<T,dimworld> > &  P);

  static void orderPoints3D(const Dune::FieldVector<T,dimworld>     centroid,
                            std::vector<int> &                      id,
                            std::vector<Dune::FieldVector<T,dimworld> >  & P);

};

template<int dim1, int dim2, int dimworld, typename T>
void MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::segmentSegmentIntersection2D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                                    const Dune::array<Dune::FieldVector<T,dimworld>, 2>& Y,
                                                    std::vector<Dune::FieldVector<T,dimworld> > & P,
                                                    bool& sameElement )
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  Dune::FieldVector<T,dimworld>  p,r ;
  Dune::FieldMatrix<T,dimworld,dimworld>  A ;

  int i=0;
  int j=0;

  Dune::FieldVector<T,dimworld> B = Y[j] - X[i] ;

  int I = (i+1)%3 ;
  int J = (j+1)%3 ;

  A[0][0] =  X[I][0] - X[i][0] ;  A[1][0] =  X[I][1] - X[i][1] ;
  A[0][1] =  Y[j][0] - Y[J][0] ;  A[1][1] =  Y[j][1] - Y[J][1] ;

  if (A.determinant()!=0) {

    A.solve(r,B) ;

    if ((r[0]>=0)&&(r[0]<=1)&&(r[1]>=0)&&(r[1]<=1)) {
      p = X[I] - X[i] ;
      p *= r[0] ;
      p += X[i] ;

      for (int k = 0; k < P.size(); k++){
        if( p == P[k] && sameElement )
          return;
      }
      P.push_back(p);
      sameElement = true;
    }
  }
}

template<int dim1, int dim2, int dimworld, typename T>
void MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::segmentPlaneIntersection3D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                                    const Dune::array<Dune::FieldVector<T,dimworld>, 4>& Y,
                                                    std::vector<Dune::FieldVector<T,dimworld> > & P ,
                                                    bool& sameElement)
{
  // X: Segment coords
  // P: intersection coords

  /* from http://en.wikipedia.org/wiki/Line-plane_intersection */

  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  Dune::FieldVector<T,dimworld>  p, q, r,  v1, v2, v3;
  Dune::FieldMatrix<T,dimworld,dimworld> A ;

  int i=0;
  double eps = 1e-12;

  Dune::FieldVector<T,dimworld> B = X[i+1] - Y[i] ;

  A [0][0] = X[i+1][0] - X[i][0] ;  A [0][1] = Y[i+1][0] -  Y[i][0] ;  A [0][2] = Y[i+2][0] -  Y[i][0] ;
  A [1][0] = X[i+1][1] - X[i][1] ;  A [1][1] = Y[i+1][1] -  Y[i][1] ;  A [1][2] = Y[i+2][1] -  Y[i][1] ;
  A [2][0] = X[i+1][2] - X[i][2] ;  A [2][1] = Y[i+1][2] -  Y[i][2] ;  A [2][2] = Y[i+2][2] -  Y[i][2] ;

  if (A.determinant()!=0) {
    A.solve(r,B) ;

    if ( ((r[0])>= 0.0)&&
         ((r[0])<= 1.0)&&
         ((r[1])>= 0.0)&&
         ((r[1])<= 1.0)&&
         ((r[2])>= 0.0)&&
         ((r[2])<= 1.0) )
      {

      p = X[i] - X[i+1] ;
      p *= r[0] ;
      p += X[i+1] ;

      for (int k = 0; k < P.size(); k++){

        v3[0] = std::abs(P[k][0] - p[0]);
        v3[1] = std::abs(P[k][1] - p[1]);
        v3[2] = std::abs(P[k][2] - p[2]);

        // skip if this intersection was already found in this element
        if( ((v3[0] < eps)&& (v3[1] < eps) && (v3[2] < eps) ) && sameElement )
          return;
      }

      if( ( std::abs(p[0]) < eps )||  (std::abs(p[1]) < eps) || ( std::abs(p[2]) < eps) )
        return;

      for (int m = 0; m < 2; m++)
        if(  (std::abs(p[0]- X[0+m][0]) < eps ) &&
             (std::abs(p[1]- X[0+m][1]) < eps ) &&
             (std::abs(p[2]- X[0+m][2]) < eps ) )
          return;

      P.push_back(p);
      sameElement = true;
    }
  }
}

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
//
//                               3D subroutines !
//
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

//  INTERSECTIONS computes edge intersections of two triangles for the the given tetrahedron X and triangle Y.
//  (point coordinates are stored column-wise, in counter clock order) the points P where their edges intersect.

template<int dim1, int dim2, int dimworld, typename T>
void MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::intersections3D( const std::vector<Dune::FieldVector<T,dimworld> >   X,
                                                                         const std::vector<Dune::FieldVector<T,dimworld> >   Y,
                                                                         std::vector<Dune::FieldVector<T,dimworld> > & P )
{

    assert(X.size() == 4);  // simplex
    assert(Y.size() == 3);  // triangle

    // get size_type for all the vectors we are using
    typedef typename std::vector<Dune::Empty>::size_type size_type;

    int l1[6],l2[6] ;
    Dune::FieldVector<T,dimworld>  p ;
    int k ;

    l1[0]= 0 ; l1[1]= 1 ; l1[2]= 2 ; l1[3]= 3 ; l1[4]= 0 ; l1[5]= 1 ;   //  enumeration of lines
    l2[0]= 1 ; l2[1]= 2 ; l2[2]= 3 ; l2[3]= 0 ; l2[4]= 2 ; l2[5]= 3 ;

    // find the points of Y in X
    for ( size_type i=0; i<3; ++i)
    {
        if (pointInTetrahedra3D(Y[i],X))
            k=insertPoint3D(Y[i],P) ;
    }

    // find the points of X in Y (simplex face aligns with hypersurface grid element
    for ( size_type i=0; i<4; ++i)
    {
        if (pointInTriangle3D(X[i],Y[0],Y[1],Y[2])) {
            k=insertPoint3D(X[i],P);

            std::cout << "point " << X[i]  << " in triangle " << Y[0] << ", " << Y[1] << ", " << Y[2] << std::endl;

        }
    }

    // find points of X on edge of Y
    for ( size_type i=0; i<4; ++i) for ( size_type j=0; j<3; ++j) {
        if (pointOnSegment3D(X[i],Y[i],Y[(i+1)%3]))
            k=insertPoint3D(X[i],P);
    }


    // find points of Y on edge of X
    for ( size_type i=0; i<3; ++i) for ( size_type j=0; j<6; ++j) {
        if (pointOnSegment3D(Y[i],X[l1[j]],X[l2[j]]))
            k=insertPoint3D(Y[i],P);
    }

    // find intersection points on element segments
    for ( size_type i=0; i<6; ++i) for ( size_type j=0; j<3; ++j)
    {
        if (segmentSegmentIntersection3D(X[l1[i]],X[l2[i]],Y[j],Y[(j+1)%3],p))
             k=insertPoint3D(p,P);

    }

    //  find the intersection of edges of X vs triangle Y

    for ( size_type i=0; i<6; ++i)
    {
        p = 0 ;
        if (triangleLineIntersection3D(X[l1[i]],X[l2[i]],Y[0],Y[1],Y[2],p))
            k=insertPoint3D(p,P) ;
    }

    //  find the intersection of edges of Y vs triangles of X

    for ( size_type i=0; i<3; ++i) for ( size_type j=0; j<4; ++j)
    {
        p = 0 ;
        if (triangleLineIntersection3D(Y[i],Y[(i+1)%3],X[j],X[(j+1)%4],X[(j+2)%4],p))
            k=insertPoint3D(p,P) ;
    }

}


// SEGMENTSEGMENTINTERSECTOPM check if segment with nodes X0 and X1 intersects with segment
// with nodes Y0 and Y1 and compute the intersection point.
template<int dim1, int dim2, int dimworld, typename T>
bool MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::segmentSegmentIntersection3D( const Dune::FieldVector<T,dimworld>    X0,
                                    const Dune::FieldVector<T,dimworld>    X1,
                                    const Dune::FieldVector<T,dimworld>    Y0,
                                    const Dune::FieldVector<T,dimworld>    Y1,
                                    Dune::FieldVector<T,dimworld>   & p) {

    Dune::FieldVector<T,dimworld>  dX, dY, dZ, cXY, cYZ;

    dX = X1-X0;
    dY = Y1-Y0;
    dZ = Y0-X0;

    cXY[0] = dX[1]* dY[2] - dX[2]* dY[1];
    cXY[1] = dX[2]* dY[0] - dX[0]* dY[2];
    cXY[2] = dX[0]* dY[1] - dX[1]* dY[0];

    if (dZ.dot(cXY)!= 0)
        return false;

    cYZ[0] = dY[1]* dZ[2] - dY[2]* dZ[1];
    cYZ[1] = dY[2]* dZ[0] - dY[0]* dZ[2];
    cYZ[2] = dY[0]* dZ[1] - dY[1]* dZ[0];

    T s = cYZ.dot(cXY) / cXY.two_norm();

    if (s >= 0 && s <= 1) {
        p = dX;
        p*= s;
        p+= X0;
        return true;
    }

    return false;
}


//   TRIANGLELINEINTERSECTION intersection of a line and a triangle TriangeLineIntersection(X,Y,P);
//   computes for a given line X in 3d,  and a triangle Y in 3d, the point p of intersection in the
//   triangle, and otherwise returns p=[];

template<int dim1, int dim2, int dimworld, typename T>
bool MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::triangleLineIntersection3D( const Dune::FieldVector<T,dimworld>    X0,
                                                                                    const Dune::FieldVector<T,dimworld>    X1,
                                                                                    const Dune::FieldVector<T,dimworld>    Y0,
                                                                                    const Dune::FieldVector<T,dimworld>    Y1,
                                                                                    const Dune::FieldVector<T,dimworld>    Y2,
                                                                                    Dune::FieldVector<T,dimworld>   & p)
{
    Dune::FieldVector<T,dimworld>      B,r ;
    Dune::FieldMatrix<T,dimworld,dimworld>  A ;
    bool found = false ;

    B = Y0 - X0 ;

    A[0][0] =  X1[0] - X0[0] ;  A[0][1] =  Y0[0] - Y1[0] ;  A[0][2] =  Y0[0] - Y2[0] ;
    A[1][0] =  X1[1] - X0[1] ;  A[1][1] =  Y0[1] - Y1[1] ;  A[1][2] =  Y0[1] - Y2[1] ;
    A[2][0] =  X1[2] - X0[2] ;  A[2][1] =  Y0[2] - Y1[2] ;  A[2][2] =  Y0[2] - Y2[2] ;

    if (A.determinant()!=0) {

        A.solve(r,B) ;

        if ((r[0]>=0)&&(r[0]<=1)&&(r[1]>=0)&&(r[1]<=1)&&(r[2]>=0)&&((r[1]+r[2])<=1)) {
            p =  X1 - X0 ;
            p *= r[0] ;
            p += X0 ;
            found = true ;
        }

    }
    return found ;
}

// POINTONSEGMENT check if the point X is contained in the line segment with nodes Y0 and Y1
template<int dim1, int dim2, int dimworld, typename T>
bool MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::pointOnSegment3D(const Dune::FieldVector<T,dimworld> X,
                              const Dune::FieldVector<T,dimworld> Y0,
                              const Dune::FieldVector<T,dimworld> Y1) {

    T eps = 1e-10;

    Dune::FieldVector<T,dimworld> u,v;

    u = Y0 - X;
    v = Y1 - Y0;

    T d = sqrt((u.dot(u) * v.dot(v) - (u.dot(v)*u.dot(v))) / (v.dot(v)));

    if (fabs(d) < eps) {
        Dune::FieldVector<T,dimworld> w = X - Y0;
        T t = w.dot(v) / v.dot(v);

        return (t >=0 && t <= 1);
    }

    return false;
}


// POINTINTRIANGLE check if the point X is contained in the triangle with nodes (Y0,Y1,Y2).
template<int dim1, int dim2, int dimworld, typename T>
bool MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::pointInTriangle3D(const Dune::FieldVector<T,dimworld> X,
                                                                          const Dune::FieldVector<T,dimworld> Y0,
                                                                          const Dune::FieldVector<T,dimworld> Y1,
                                                                          const Dune::FieldVector<T,dimworld> Y2) {


    double eps= 1.e-10 ;     // tolerance for relative error

    Dune::FieldVector<T,dimworld> v0,v1,v2,r;

    v0 = Y1 - Y0;
    v1 = Y2 - Y0;
    v2 = X - Y0;

    T s,t,d;

    d =  ((v0.dot(v0))*(v1.dot(v1)) - (v0.dot(v1))*(v0.dot(v1)));

    s = ((v1.dot(v1))*(v0.dot(v2)) - (v0.dot(v1))*(v1.dot(v2))) / d;
    t = ((v0.dot(v0))*(v1.dot(v2)) - (v0.dot(v1))*(v0.dot(v2))) / d;

    v0*=s;
    v1*=t;
    r = Y0 + v0 + v1 - X;

    return (s >= 0 && t >= 0 && (s+t)<= 1 && r.infinity_norm() < eps);

}

//   POINTINTETRAHEDRA check if the point X is contained in the tetrahedra Y.

template<int dim1, int dim2, int dimworld, typename T>
bool MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::pointInTetrahedra3D( const Dune::FieldVector<T,dimworld> X,
                                                    const std::vector<Dune::FieldVector<T,dimworld> >   Y)
{
  Dune::FieldMatrix<T,dimworld+1,dimworld+1>  D,DD ;
  T D0,D1,D2,D3,D4 ;

  D[0][0] =  Y[0][0] ;  D[0][1] =  Y[1][0] ;  D[0][2] =  Y[2][0] ;  D[0][3] =  Y[3][0] ;
  D[1][0] =  Y[0][1] ;  D[1][1] =  Y[1][1] ;  D[1][2] =  Y[2][1] ;  D[1][3] =  Y[3][1] ;
  D[2][0] =  Y[0][2] ;  D[2][1] =  Y[1][2] ;  D[2][2] =  Y[2][2] ;  D[2][3] =  Y[3][2] ;
  D[3][0] =        1 ;  D[3][1] =        1 ;  D[3][2] =        1 ;  D[3][3] =        1 ;

  D0 = D.determinant() ;

  DD = D ; DD[0][0] = X[0] ; DD[1][0] = X[1] ; DD[2][0] = X[2] ; D1 = DD.determinant() ;
  DD = D ; DD[0][1] = X[0] ; DD[1][1] = X[1] ; DD[2][1] = X[2] ; D2 = DD.determinant() ;
  DD = D ; DD[0][2] = X[0] ; DD[1][2] = X[1] ; DD[2][2] = X[2] ; D3 = DD.determinant() ;
  DD = D ; DD[0][3] = X[0] ; DD[1][3] = X[1] ; DD[2][3] = X[2] ; D4 = DD.determinant() ;

  return ((D0*D1>=0)&&(D0*D2>=0)&&(D0*D3>=0)&&(D0*D4>=0)) ;
}

//  INSERTPOINT inserts an intersection point p into the list of intersection points P. If the point p is already
//  contained in P, the point is not inserted, and k is the index of the point found in the list. Otherwise,
//  the new point is inserted, and k is its index in the list.

template<int dim1, int dim2, int dimworld, typename T>
int MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::insertPoint3D( const Dune::FieldVector<T,dimworld>  p, std::vector<Dune::FieldVector<T,dimworld> > &  P)
{
    double eps= 1.e-10 ;     // tolerance for identical nodes
    int k=0 ;

    if (P.size()>0) {

        while ((k<P.size())&&((p - P[k]).infinity_norm()>eps))
            k++ ;

        if (k>=P.size())
            P.push_back(p) ;        //  new node is not contained in P

    }
    else
        P.push_back(p);

    return k ;
}

//   ORDERPOINTS order points counterclockwise [p,id]=OrderPoints(p,centroid); orders the points in a plane in 3d
//   stored columnwise in the matrix p counterclockwise looking from the side opposite to the point centroid
//   outside the plane. There must be more than two points. p contains the reorderd points, and id
//   contains the index reordering.

template<int dim1, int dim2, int dimworld, typename T>
void MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::orderPoints3D(const Dune::FieldVector<T,dimworld>     centroid,
                                                                      std::vector<int> &                      id,
                                                                      std::vector<Dune::FieldVector<T,dimworld> >  & P)
{
    // get size_type for all the vectors we are using
    typedef typename std::vector<Dune::Empty>::size_type size_type;

    Dune::FieldVector<T,dimworld> c,d1,d2,dr,dn,cross,d ;
    std::vector<T> ai ;

    d1 = P[1] - P[0] ;    // two reference vectors
    d2 = P[2] - P[0] ;

    cross[0] = d1[1]*d2[2] - d1[2]*d2[1] ;    // cross product
    cross[1] = d1[2]*d2[0] - d1[0]*d2[2] ;
    cross[2] = d1[0]*d2[1] - d1[1]*d2[0] ;

    if (((centroid - P[0])*cross)<0)   // good orientation ?
    {
        dr = d1 ;
        dr /= dr.two_norm()  ;       // 'x-axis' unit vector
        dn = dr ;
        dn *= -(d2*dr) ;
        dn += d2 ;
        dn /= dn.two_norm()  ;       // 'y-axis' unit vector
    }
    else
    {
        dr = d2 ;
        dr /= dr.two_norm()  ;       // 'y-axis' unit vector
        dn = dr ;
        dn *= -(d1*dr) ;
        dn += d1 ;
        dn /= dn.two_norm()  ;        // 'x-axis' unit vector
    }

    // definition of angles, using projection on the local reference, ie by scalarly multipliying by dr and dn resp.
    for ( size_type j=1 ; j < P.size() ; j++)
    {
        ai.push_back(atan2((P[j]-P[0])*dn,(P[j]-P[0])*dr)) ;
        id.push_back(j) ;
    }

    // sort according to increasing angles
    for ( size_type j=1; j < ai.size(); j++)
        for ( size_type i=0; i < j; i++)
            if (ai[j]<ai[i]) {
                std::swap<T>(ai[i],ai[j]) ;
                std::swap<int>(id[i],id[j]) ;
            }

    id.insert(id.begin(),0) ;
}


}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRID_GLUE_MIXED_DIM_OVERLAPPING_MERGE_HH
