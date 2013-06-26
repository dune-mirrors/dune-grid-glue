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
      std::vector<Dune::FieldVector<T,2> > intersections;
      for (int i=0; i<refElement2.size(1); i++)
      {
        Dune::array<Dune::FieldVector<T,2>, 2> edge;
        edge[0] = grid2ElementCorners[refElement2.subEntity(i,1,0,2)];
        edge[1] = grid2ElementCorners[refElement2.subEntity(i,1,1,2)];
        segmentSegmentIntersection2D(grid1ElementCorners,edge,intersections);
      }

      // Compute whether the segment endpoints are contained in the triangle
      std::bitset<2> inTriangle("00");

      for (size_t i=0; i<2; i++)
        inTriangle[i] = Dune::GenericReferenceElements<T,dim2>::general(grid2ElementType).checkInside(grid2Geometry.local(grid1ElementCorners[i]));

      switch (intersections.size()) {

          case 0:
            if (inTriangle[0] and inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = 0;
              this->intersections_.back().grid1Local_[1] = 1;

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = grid2Geometry.local(grid1ElementCorners[0]);
              this->intersections_.back().grid2Local_[1] = grid2Geometry.local(grid1ElementCorners[1]);
            }
            return;

          case 1:
          {
            if (inTriangle[0])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = 0;
              this->intersections_.back().grid1Local_[1] = grid1Geometry.local(intersections[0]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = grid2Geometry.local(grid1ElementCorners[0]);
              this->intersections_.back().grid2Local_[1] = grid2Geometry.local(intersections[0]);

            } else if (inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = grid1Geometry.local(intersections[0]);
              this->intersections_.back().grid1Local_[1] = 1;

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = grid2Geometry.local(intersections[0]);
              this->intersections_.back().grid2Local_[1] = grid2Geometry.local(grid1ElementCorners[1]);

            }

            return;


          }


          case 2:
          case 3:
          {
            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0] = grid1Geometry.local(intersections[0]);
            this->intersections_.back().grid1Local_[1] = grid1Geometry.local(intersections[1]);

            // Compute local coordinates in the grid2 element
            this->intersections_.back().grid2Local_[0] = grid2Geometry.local(intersections[0]);
            this->intersections_.back().grid2Local_[1] = grid2Geometry.local(intersections[1]);

            return;


          }

      }

    } else if (dim1==2 and dim2==1 and dimworld==2)
    {
      assert(grid2ElementCorners.size() == 2);  // line segment

      // Compute intersections of segment with all edges of the 2d element
      std::vector<Dune::FieldVector<T,2> > intersections;
      for (int i=0; i<refElement1.size(1); i++)
      {
        Dune::array<Dune::FieldVector<T,2>, 2> edge;
        edge[0] = grid1ElementCorners[refElement1.subEntity(i,1,0,2)];
        edge[1] = grid1ElementCorners[refElement1.subEntity(i,1,1,2)];
        segmentSegmentIntersection2D(grid2ElementCorners,edge,intersections);
      }

      // Compute whether the segment endpoints are contained in the triangle
      std::bitset<2> inTriangle("00");

      for (size_t i=0; i<2; i++)
        inTriangle[i] = Dune::GenericReferenceElements<T,dim1>::general(grid1ElementType).checkInside(grid1Geometry.local(grid2ElementCorners[i]));

      switch (intersections.size()) {

          case 0:
            if (inTriangle[0] and inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = grid1Geometry.local(grid2ElementCorners[0]);
              this->intersections_.back().grid1Local_[1] = grid1Geometry.local(grid2ElementCorners[1]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = 0;
              this->intersections_.back().grid2Local_[1] = 1;
            }
            return;

          case 1:
          {
            if (inTriangle[0])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = grid1Geometry.local(grid2ElementCorners[0]);
              this->intersections_.back().grid1Local_[1] = grid1Geometry.local(intersections[0]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = 0;
              this->intersections_.back().grid2Local_[1] = grid2Geometry.local(intersections[0]);

            } else if (inTriangle[1])
            {
              this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = grid1Geometry.local(intersections[0]);
              this->intersections_.back().grid1Local_[1] = grid1Geometry.local(grid2ElementCorners[1]);

              // Compute local coordinates in the grid2 element
              this->intersections_.back().grid2Local_[0] = grid2Geometry.local(intersections[0]);
              this->intersections_.back().grid2Local_[1] = 1;

            }

            return;


          }


          case 2:
          case 3:
          {
            this->intersections_.push_back(RemoteSimplicialIntersection(grid1Index, grid2Index));

            // Compute local coordinates in the grid1 element
            this->intersections_.back().grid1Local_[0] = grid1Geometry.local(intersections[0]);
            this->intersections_.back().grid1Local_[1] = grid1Geometry.local(intersections[1]);

            // Compute local coordinates in the grid2 element
            this->intersections_.back().grid2Local_[0] = grid2Geometry.local(intersections[0]);
            this->intersections_.back().grid2Local_[1] = grid2Geometry.local(intersections[1]);

            return;


          }

      }

    } else
      std::cout << "MixedDimOverlappingMerge for dim1 = " << dim1 << ", dim2 = "
                                                          << dim2 << ", dimworld = "
                                                          << dimworld << ": Not implemented yet!" << std::endl;
  }

  static void segmentSegmentIntersection2D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                              const Dune::array<Dune::FieldVector<T,dimworld>, 2>& edge,
                                              std::vector<Dune::FieldVector<T,dimworld> > & P );

};

template<int dim1, int dim2, int dimworld, typename T>
void MixedDimOverlappingMerge<dim1, dim2, dimworld, T>::segmentSegmentIntersection2D( const std::vector<Dune::FieldVector<T,dimworld> >& X,
                                                    const Dune::array<Dune::FieldVector<T,dimworld>, 2>& Y,
                                                    std::vector<Dune::FieldVector<T,dimworld> > & P )
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
          P.push_back(p);
        }
      }

}


}  // namespace Dune

#endif // DUNE_GRID_GLUE_MIXED_DIM_OVERLAPPING_MERGE_HH
