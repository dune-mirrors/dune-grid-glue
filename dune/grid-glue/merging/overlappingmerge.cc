// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include <algorithm>
#include <dune/common/typetraits.hh>

template<int dim, typename T>
void OverlappingMerge<dim, T>::
computeIntersection(const Dune::GeometryType& grid1ElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& grid1ElementCorners,
                    unsigned int grid1Index,
                    std::bitset<(1<<dim)>& neighborIntersects1,
                    const Dune::GeometryType& grid2ElementType,
                    const std::vector<Dune::FieldVector<T,dim> >& grid2ElementCorners,
                    unsigned int grid2Index,
                    std::bitset<(1<<dim)>& neighborIntersects2)
{
  this->counter++;

  // A few consistency checks
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid1ElementType).size(dim)) == grid1ElementCorners.size());
  assert((unsigned int)(Dune::GenericReferenceElements<T,dim>::general(grid2ElementType).size(dim)) == grid2ElementCorners.size());

  // Make generic geometries representing the grid1- and grid2 element.
  // this eases computation of local coordinates.
  typedef Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<T,dim,dim> > Geometry;

  Geometry grid1Geometry(grid1ElementType, grid1ElementCorners);
  Geometry grid2Geometry(grid2ElementType, grid2ElementCorners);

  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  switch (dim) {
  case 1 : {

    // Check consistent orientation
    // \todo Reverse the orientation if this check fails
    assert(grid1ElementCorners[0][0] <= grid1ElementCorners[1][0]);
    assert(grid2ElementCorners[0][0] <= grid2ElementCorners[1][0]);

    T lowerBound = std::max(grid1ElementCorners[0][0], grid2ElementCorners[0][0]);
    T upperBound = std::min(grid1ElementCorners[1][0], grid2ElementCorners[1][0]);

    if (lowerBound <= upperBound) {      // Intersection is non-empty

      this->intersections_.push_back(RemoteSimplicialIntersection());

      // Compute local coordinates in the grid1 element
      this->intersections_.back().grid1Local_[0] = grid1Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().grid1Local_[1] = grid1Geometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Compute local coordinates in the grid2 element
      this->intersections_.back().grid2Local_[0] = grid2Geometry.local(Dune::FieldVector<T,dim>(lowerBound));
      this->intersections_.back().grid2Local_[1] = grid2Geometry.local(Dune::FieldVector<T,dim>(upperBound));

      // Set indices
      this->intersections_.back().grid1Entity_ = grid1Index;
      this->intersections_.back().grid2Entity_ = grid2Index;

      //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

    }
    break;
  }
  case 2 : {

    std::vector<Dune::FieldVector<T,dim> >  P;

    // find the intersections of any segement of the two triangles.
    edgeIntersections2D(grid1ElementCorners,grid2ElementCorners,P);

    // add the points of grid1 in grid2 and of grid2 in grid1.
    pointsofXinY2D(grid1ElementCorners,grid2ElementCorners,P);
    pointsofXinY2D(grid2ElementCorners,grid1ElementCorners,P);

    // sort points counter clock wise and removes duplicates
    if (P.size()>=3)
      sortAndRemoveDoubles2D(P);

    //      TO check if the previous function made a good job, uncomment the next lines.
    //
    //        if (P.size()>=3) {   for ( size_type i=0 ; i < 3 ; ++i) std::cout << " grid1 " << grid1ElementCorners[i] << std::endl ;
    //                             for ( size_type i=0 ; i < 3 ; ++i) std::cout << " grid2 " << grid2ElementCorners[i] << std::endl ;
    //                             std::cout << " Size E " << P.size() << std::endl ;
    //                             for ( size_type i=0 ; i < P.size() ; ++i) std::cout << " P " << P[i] << std::endl ;  }

    if (P.size()>=3)
      for ( size_type i=0 ; i < P.size() - 2 ; ++i) {

        this->intersections_.push_back(RemoteSimplicialIntersection());

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0] = grid1Geometry.local(P[0]);
        this->intersections_.back().grid1Local_[1] = grid1Geometry.local(P[i+1]);
        this->intersections_.back().grid1Local_[2] = grid1Geometry.local(P[i+2]);

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid2Local_[0] = grid2Geometry.local(P[0]);
        this->intersections_.back().grid2Local_[1] = grid2Geometry.local(P[i+1]);
        this->intersections_.back().grid2Local_[2] = grid2Geometry.local(P[i+2]);

        // Set indices
        this->intersections_.back().grid1Entity_ = grid1Index;
        this->intersections_.back().grid2Entity_ = grid2Index;

      }

    break;
  }
  case 3 : {

    std::vector<Dune::FieldVector<T,dim> >  P;
    std::vector<std::vector<int> >          H,SX(4),SY(4);
    Dune::FieldVector<T,dim>                centroid;

    // Compute intersections ( Create SX, SY and P )
    intersections3D(grid1ElementCorners,grid2ElementCorners,SX,SY,P) ;

    if (P.size()>=4) {

      if (P.size()==4) {         // if the intersection is one tetrahedron, no need to go further

        this->intersections_.push_back(RemoteSimplicialIntersection());

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid1Local_[0] = grid1Geometry.local(P[0]);
        this->intersections_.back().grid1Local_[1] = grid1Geometry.local(P[1]);
        this->intersections_.back().grid1Local_[2] = grid1Geometry.local(P[2]);
        this->intersections_.back().grid1Local_[3] = grid1Geometry.local(P[3]);

        // Compute local coordinates in the grid1 element
        this->intersections_.back().grid2Local_[0] = grid2Geometry.local(P[0]);
        this->intersections_.back().grid2Local_[1] = grid2Geometry.local(P[1]);
        this->intersections_.back().grid2Local_[2] = grid2Geometry.local(P[2]);
        this->intersections_.back().grid2Local_[3] = grid2Geometry.local(P[3]);

        // Set indices
        this->intersections_.back().grid1Entity_ = grid1Index;
        this->intersections_.back().grid2Entity_ = grid2Index;
      }
      else
      {
        // Compute the centroid
        centroid=0;
        for (size_type i=0; i < P.size(); i++)
          centroid += P[i] ;
        centroid /= P.size() ;

        // Sorte each faces ( Create H )
        H.clear() ;
        sorting3D(centroid,SX,SY,P,H) ;

        /*      TO check if the previous routines made a good job, uncomment the next lines.

                std::cout << " ------------------------------------------------------------ " << std::endl ;

                   for ( size_type i=0 ; i < 4 ; ++i) std::cout << " grid1 " << grid1ElementCorners[i] << std::endl ;
                   for ( size_type i=0 ; i < 4 ; ++i) std::cout << " grid2 " << grid2ElementCorners[i] << std::endl ;
                   std::cout << " +++ Size P " << P.size() << std::endl ;
                   for ( size_type i=0 ; i < P.size() ; ++i) std::cout << " P " << P[i] << std::endl ;

                   std::cout << " +++ SX " << std::endl ;
                   for ( size_type i=0; i < SX.size() ; i++)
                        { for ( size_type j=0 ; j < SX[i].size() ; ++j)  {   std::cout  << SX[i][j] << "  ;  "  ; }  std::cout << std::endl ; }

                   std::cout << " +++ SY " << std::endl ;
                   for ( size_type i=0; i < SY.size() ; i++)
                        { for ( size_type j=0 ; j < SY[i].size() ; ++j)  {   std::cout  << SY[i][j] << "  ;  "  ; }  std::cout << std::endl ; }

                   std::cout << " +++ H " << std::endl ;
                   for ( size_type i=0; i < H.size() ; i++)
                        { for ( size_type j=0 ; j < H[i].size() ; ++j)   {   std::cout  << H[i][j] << "  ;  "  ; }  std::cout << std::endl ; }
         */

        //  Loop over all facets
        for (size_type i=0; i < H.size() ; i++) {
          //  Loop over all triangles of facets if the face is not degenerated
          if (H[i].size()>=3)
            for ( size_type j=0 ; j < H[i].size() - 2 ; ++j) {

              // Output the tetrahedron (anchor, next, nextNext, centroid)
              this->intersections_.push_back(RemoteSimplicialIntersection());

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid1Local_[0] = grid1Geometry.local(P[H[i][0]]);
              this->intersections_.back().grid1Local_[1] = grid1Geometry.local(P[H[i][j+1]]);
              this->intersections_.back().grid1Local_[2] = grid1Geometry.local(P[H[i][j+2]]);
              this->intersections_.back().grid1Local_[3] = grid1Geometry.local(centroid);

              // Compute local coordinates in the grid1 element
              this->intersections_.back().grid2Local_[0] = grid2Geometry.local(P[H[i][0]]);
              this->intersections_.back().grid2Local_[1] = grid2Geometry.local(P[H[i][j+1]]);
              this->intersections_.back().grid2Local_[2] = grid2Geometry.local(P[H[i][j+2]]);
              this->intersections_.back().grid2Local_[3] = grid2Geometry.local(centroid);

              // Set indices
              this->intersections_.back().grid1Entity_ = grid1Index;
              this->intersections_.back().grid2Entity_ = grid2Index;

            }
        }
      }
    }

    break;
  }
  default :
    DUNE_THROW(Dune::NotImplemented, "OverlappingMerge is not implemented for dim==" << dim << "!");

  }

}

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
//
//                               2 D   sub  routines !
//
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------


//  EDGEINTERSECTIONS computes edge intersections of two triangles for the two given triangles X and Y
//  (point coordinates are stored column-wise, in counter clock order) the points P where their edges intersect.

template<int dim, typename T>
void OverlappingMerge<dim, T>::edgeIntersections2D( const std::vector<Dune::FieldVector<T,dim> >   X,
                                                    const std::vector<Dune::FieldVector<T,dim> >   Y,
                                                    std::vector<Dune::FieldVector<T,dim> > & P )
{

  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  int I,J ;
  Dune::FieldVector<T,dim>  p,B,r ;
  Dune::FieldMatrix<T,dim,dim>  A ;

  for ( size_type i=0; i<3; ++i) for ( size_type j=0; j<3; ++j)
    {

      B = Y[j] - X[i] ;

      I = (i+1)%3 ;
      J = (j+1)%3 ;

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

}

// PointsofXinY finds corners of one triangle within another one, computes for the two given triangles X
// and Y (point coordinates are stored column-wise, in counter clock  order) the corners P of X which lie in the interior of Y.

template<int dim, typename T>
void OverlappingMerge<dim, T>::pointsofXinY2D( const std::vector<Dune::FieldVector<T,dim> >   X,
                                               const std::vector<Dune::FieldVector<T,dim> >   Y,
                                               std::vector<Dune::FieldVector<T,dim> > & P )
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  Dune::FieldMatrix<double,dim,dim> A ;
  Dune::FieldVector<double,dim>  v0 ,v1 ,v2, U, B  ;

  v0 = Y[1] - Y[0] ;
  v1 = Y[2] - Y[0] ;

  A[0][0] = v0*v0 ;
  A[0][1] = v0*v1 ;
  A[1][0] = v1*v0 ;
  A[1][1] = v1*v1 ;

  for ( size_type i=0; i<3; ++i) {
    v2 = X[i] - Y[0];
    B[0] = v0*v2 ;
    B[1] = v1*v2 ;
    A.solve(U,B) ;
    if ((U[0]>=0)&&(U[1]>=0)&&(U[0]+U[1]<=1))
      P.push_back(X[i]);
  }

}

//  SortAndRemoveDoubles: orders polygon corners in P counter clock wise and removes duplicates

template<int dim, typename T>
void OverlappingMerge<dim, T>::sortAndRemoveDoubles2D( std::vector<Dune::FieldVector<T,dim> > & P )
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  Dune::FieldVector<double,dim> c ;
  std::vector< double > ai ;

  // build barycentre c of all these points
  c = 0 ;
  for ( size_type i=0; i < P.size(); i++)
    c += P[i] ;
  c /= P.size() ;

  // definition of angles
  for ( size_type i=0; i < P.size(); i++)
    ai.push_back(atan2(P[i][1]-c[1],P[i][0]-c[0]));

  // sort according to increasing angles
  for ( size_type j=1; j < ai.size(); j++)
    for ( size_type i=0; i < j; i++) if (ai[j]<ai[i]) {
        std::swap<double>(ai[i],ai[j]);
        std::swap<Dune::FieldVector<double,dim> >(P[i],P[j]);
      }

  //  Remove exact Doubles
  P.erase(std::unique(P.begin(), P.end()), P.end());

  // If for some reasons, one wants to eliminate also very close neighbouring (up to a distance eps)
  // then just coment the previous line, ad uncomment the next ones.

  //        double eps=1.e-10 ;
  //        std::vector<Dune::FieldVector<double,dim> >   Q = P ;
  //        P.clear() ;
  //        P.push_back(Q[0]) ;
  //        for ( size_type j=1; j < Q.size(); j++) if ((P[P.size()-1] - Q[j]).infinity_norm()>eps) P.push_back(Q[j]) ;

}

// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------
//
//                               3D subroutines !
//
// -------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------------------

//  INTERSECTIONS computes edge intersections of two triangles for the two given triangles X and Y
//  (point coordinates are stored column-wise, in counter clock order) the points P where their edges intersect.

template<int dim, typename T>
void OverlappingMerge<dim, T>::intersections3D( const std::vector<Dune::FieldVector<T,dim> >   X,
                                                const std::vector<Dune::FieldVector<T,dim> >   Y,
                                                std::vector<std::vector<int> >         & SX,
                                                std::vector<std::vector<int> >         & SY,
                                                std::vector<Dune::FieldVector<T,dim> > & P )
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  int l1[6],l2[6],s1[6],s2[6],ni[4][3] ;
  Dune::FieldVector<T,dim>  p ;
  int k ;

  l1[0]= 0 ; l1[1]= 1 ; l1[2]= 2 ; l1[3]= 3 ; l1[4]= 0 ; l1[5]= 1 ;   //  enumeration of lines
  l2[0]= 1 ; l2[1]= 2 ; l2[2]= 3 ; l2[3]= 0 ; l2[4]= 2 ; l2[5]= 3 ;

  s1[0]= 0 ; s1[1]= 0 ; s1[2]= 1 ; s1[3]= 2 ; s1[4]= 0 ; s1[5]= 1 ;   // faces touching each line
  s2[0]= 3 ; s2[1]= 1 ; s2[2]= 2 ; s2[3]= 3 ; s2[4]= 2 ; s2[5]= 3 ;

  ni[0][0]= 0 ;   ni[0][1]= 2 ;   ni[0][2]= 3 ;   // faces touching each node
  ni[1][0]= 0 ;   ni[1][1]= 1 ;   ni[1][2]= 3 ;
  ni[2][0]= 0 ;   ni[2][1]= 1 ;   ni[2][2]= 2 ;
  ni[3][0]= 1 ;   ni[3][1]= 2 ;   ni[3][2]= 3 ;

  //  find the intersection of edges of X vs triangles of Y

  for ( size_type i=0; i<6; ++i) for ( size_type j=0; j<4; ++j)
    {
      p = 0 ;
      if (triangleLineIntersection3D(X[l1[i]],X[l2[i]],Y[j],Y[(j+1)%4],Y[(j+2)%4],p))
      {
        k=insertPoint3D(p,P) ;
        SY[j].push_back(k) ;                 // remember to which surfaces
        SX[s1[i]].push_back(k) ;
        SX[s2[i]].push_back(k) ;
      }
    }

  //  find the intersection of edges of Y vs triangles of X

  for ( size_type i=0; i<6; ++i) for ( size_type j=0; j<4; ++j)
    {
      p = 0 ;
      if (triangleLineIntersection3D(Y[l1[i]],Y[l2[i]],X[j],X[(j+1)%4],X[(j+2)%4],p))
      {
        k=insertPoint3D(p,P) ;
        SX[j].push_back(k) ;                // remember to which surfaces
        SY[s1[i]].push_back(k) ;
        SY[s2[i]].push_back(k) ;
      }
    }

  // find the points of X in Y

  for ( size_type i=0; i<4; ++i)
  {
    if (pointInTetrahedra3D(X[i],Y))
    {
      k=insertPoint3D(X[i],P) ;
      SX[ni[i][0]].push_back(k) ;             // remember to which surfaces
      SX[ni[i][1]].push_back(k) ;
      SX[ni[i][2]].push_back(k) ;
    }
  }

  // find the points of Y in X

  for ( size_type i=0; i<4; ++i)
  {
    if (pointInTetrahedra3D(Y[i],X))
    {
      k=insertPoint3D(Y[i],P) ;
      SY[ni[i][0]].push_back(k) ;               // remember to which surfaces
      SY[ni[i][1]].push_back(k) ;
      SY[ni[i][2]].push_back(k) ;
    }
  }

}

// SORTING ROUTINE

template<int dim, typename T>
void OverlappingMerge<dim, T>::sorting3D( const Dune::FieldVector<T,dim>    centroid,
                                          const std::vector<std::vector<int> >           SX,
                                          const std::vector<std::vector<int> >           SY,
                                          const std::vector<Dune::FieldVector<T,dim> >   P,
                                          std::vector<std::vector<int> >               & H)
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  // sorting
  int m ;
  std::vector<int> no,id,temp ;
  std::vector<Dune::FieldVector<T,dim> > p ;

  if (P.size()>3)
  {
    for ( size_type i=0; i<4; ++i)
    {
      if (SX[i].size()>0)         // loop on faces of X
      {
        no = SX[i] ;
        removeDuplicates(no) ;
        m = no.size() ;
        if ((m>2) && newFace3D(no,H))               // don't compute degenerate polygons and check if face is new
        {
          for ( size_type l=0; l<m; ++l)
            p.push_back(P[no[l]]) ;
          orderPoints3D(centroid,id,p) ;                                   // order points counter-clock-wise
          for ( size_type l=0; l<m; ++l)
            temp.push_back(no[id[l]]) ;
          H.push_back(temp) ;
          temp.clear();
          p.clear();
          id.clear();                 // clean
        }
        no.clear() ;             // clean
      }
      if (SY[i].size()>0)         // loop on faces of Y
      {
        no = SY[i] ;
        removeDuplicates(no) ;
        m = no.size() ;
        if ((m>2) && newFace3D(no,H))               // don't compute degenerate polygons  and check if face is new
        {
          m = no.size() ;
          for ( size_type l=0; l<m; ++l)
            p.push_back(P[no[l]]) ;
          orderPoints3D(centroid,id,p) ;                                   // order points counter-clock-wise
          for ( size_type l=0; l<m; ++l)
            temp.push_back(no[id[l]]) ;
          H.push_back(temp) ;
          temp.clear();
          p.clear();
          id.clear();                 // clean
        }
        no.clear() ;             // clean
      }
    }
  }

}

//   TRIANGLELINEINTERSECTION intersection of a line and a triangle TriangeLineIntersection(X,Y,P);
//   computes for a given line X in 3d,  and a triangle Y in 3d, the point p of intersection in the
//   triangle, and otherwise returns p=[];

template<int dim, typename T>
bool OverlappingMerge<dim, T>::triangleLineIntersection3D( const Dune::FieldVector<T,dim>    X0,
                                                           const Dune::FieldVector<T,dim>    X1,
                                                           const Dune::FieldVector<T,dim>    Y0,
                                                           const Dune::FieldVector<T,dim>    Y1,
                                                           const Dune::FieldVector<T,dim>    Y2,
                                                           Dune::FieldVector<T,dim>   & p)
{
  Dune::FieldVector<T,dim>      B,r ;
  Dune::FieldMatrix<T,dim,dim>  A ;
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

//   POINTINTETRAHEDRA check if the point X is contained in the tetrahedra Y.

template<int dim, typename T>
bool OverlappingMerge<dim, T>::pointInTetrahedra3D( const Dune::FieldVector<T,dim>                 X,
                                                    const std::vector<Dune::FieldVector<T,dim> >   Y)
{
  Dune::FieldMatrix<T,dim+1,dim+1>  D,DD ;
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

template<int dim, typename T>
int OverlappingMerge<dim, T>::insertPoint3D( const Dune::FieldVector<T,dim>  p, std::vector<Dune::FieldVector<T,dim> > &  P)
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

// REMOVEDUPLICATES removes duplicate entries from the vector p.

template<int dim, typename T>
void OverlappingMerge<dim, T>::removeDuplicates( std::vector<int> & p)
{
  sort(p.begin(),p.end());
  std::vector<int>::iterator it = std::unique(p.begin(),p.end());
  p.erase(it,p.end());
}

//   ORDERPOINTS order points counterclockwise [p,id]=OrderPoints(p,centroid); orders the points in a plane in 3d
//   stored columnwise in the matrix p counterclockwise looking from the side opposite to the point centroid
//   outside the plane. There must be more than two points. p contains the reorderd points, and id
//   contains the index reordering.

template<int dim, typename T>
void OverlappingMerge<dim, T>::orderPoints3D(const Dune::FieldVector<T,dim>     centroid,
                                             std::vector<int> &                      id,
                                             std::vector<Dune::FieldVector<T,dim> >  & P)
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  Dune::FieldVector<T,dim> c,d1,d2,dr,dn,cross,d ;
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

// NEWFACE checks if index set is contained already in H  b=NewFace(H,id) checks if a permutation
// of the vector id is contained in a row of the matrix H.

template<int dim, typename T>
bool OverlappingMerge<dim, T>::newFace3D(const std::vector<int> id, const std::vector<std::vector<int> > H)
{
  // get size_type for all the vectors we are using
  typedef typename std::vector<Dune::Empty>::size_type size_type;

  int n = H.size() ;
  int m = id.size() ;
  std::vector<int> A ;
  std::vector<int> B = id ;
  sort(B.begin(),B.end()) ;
  int i = 0 ;
  bool b = true ;
  double tp ;

  while ( b && (i<n) )
  {
    if ((H[i].size())>=m)
    {
      A=H[i] ;
      sort(A.begin(),A.end());
      tp = 0 ;
      for ( size_type j=0 ; j < m; j++)
        tp += std::fabs(A[j]-B[j]) ;
      b = (tp>0) ;
    }

    i += 1 ;
  }

  return b ;
}
