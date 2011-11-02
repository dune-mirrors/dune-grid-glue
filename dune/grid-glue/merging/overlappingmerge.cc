// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:


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

    std::vector<Dune::FieldVector<T,dim> >  P ;
    P.clear();

    // find the intersections of any segement of the two triangles.
    EdgeIntersections(grid1ElementCorners,grid2ElementCorners,P) ;

    // add the points of grid1 in grid2 and of grid2 in grid1.
    PointsofXinY(grid1ElementCorners,grid2ElementCorners,P) ;
    PointsofXinY(grid2ElementCorners,grid1ElementCorners,P) ;

    // sort points counter clock wise and removes duplicates
    if (P.size()>=3) SortAndRemoveDoubles(P) ;

    //      TO check if the previous function made a good job, uncomment the next lines.
    //
    //        if (P.size()>=3) {   for ( int i=0 ; i < 3 ; ++i) std::cout << " grid1 " << grid1ElementCorners[i] << std::endl ;
    //                             for ( int i=0 ; i < 3 ; ++i) std::cout << " grid2 " << grid2ElementCorners[i] << std::endl ;
    //                             std::cout << " Size E " << P.size() << std::endl ;
    //                             for ( int i=0 ; i < P.size() ; ++i) std::cout << " P " << P[i] << std::endl ;  }

    if (P.size()>=3) for ( int i=0 ; i < P.size() - 2 ; ++i) {

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

        //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;
      }

    break;
  }

  default :
    DUNE_THROW(Dune::NotImplemented, "GanderMerge is not implemented for dim==" << dim << "!");

  }

}


//  EDGEINTERSECTIONS computes edge intersections of two triangles for the two given triangles X and Y
//  (point coordinates are stored column-wise, in counter clock order) the points P where their edges intersect.

template<int dim, typename T>
void OverlappingMerge<dim, T>::EdgeIntersections( const std::vector<Dune::FieldVector<T,dim> >   X,
                                                  const std::vector<Dune::FieldVector<T,dim> >   Y,
                                                  std::vector<Dune::FieldVector<T,dim> > & P )
{

  int I,J ;
  T b0,b1,A00,A01,A10,A11,D,r0,r1 ;
  Dune::FieldVector<T,dim>  CO ;

  for (int i=0; i<3; ++i) for (int j=0; j<3; ++j)
    {

      b0 = Y[j][0] - X[i][0] ;
      b1 = Y[j][1] - X[i][1] ;

      if (i==2) I = 0 ;else I = i+1 ;
      if (j==2) J = 0 ;else J = j+1 ;

      A00 =  X[I][0] - X[i][0] ;
      A10 =  X[I][1] - X[i][1] ;
      A01 =  Y[j][0] - Y[J][0] ;
      A11 =  Y[j][1] - Y[J][1] ;

      D = A00*A11 - A01*A10 ;

      if (D!=0) {
        r0 = ( b0*A11 - b1*A01 ) / D ;
        r1 = ( A00*b1 - A10*b0 ) / D ;

        if ((r0>=0)&&(r0<=1)&&(r1>=0)&&(r1<=1)) {
          CO[0] = X[i][0] + r0*(X[I][0] - X[i][0]);
          CO[1] = X[i][1] + r0*(X[I][1] - X[i][1]) ;
          P.push_back(CO);
        }
      }
    }

}

// PointsofXinY finds corners of one triangle within another one, computes for the two given triangles X
// and Y (point coordinates are stored column-wise, in counter clock  order) the corners P of X which lie in the interior of Y.

template<int dim, typename T>
void OverlappingMerge<dim, T>::PointsofXinY( const std::vector<Dune::FieldVector<T,dim> >   X,
                                             const std::vector<Dune::FieldVector<T,dim> >   Y,
                                             std::vector<Dune::FieldVector<T,dim> > & P )
{

  Dune::FieldVector<double,dim>  v0 ,v1 ,v2  ;
  double d00,d01,d11,d02,d12,u,v,id ;

  v0 = Y[1] - Y[0] ;
  v1 = Y[2] - Y[0] ;

  d00 = v0*v0 ;
  d01 = v0*v1 ;
  d11 = v1*v1 ;

  id=1./(d00*d11-d01*d01);

  for (int i=0; i<3; ++i) {
    v2 = X[i] - Y[0];
    d02 = v0*v2 ;
    d12 = v1*v2 ;
    u=(d11*d02-d01*d12)*id ;
    v=(d00*d12-d01*d02)*id ;
    if ((u>=0)&&(v>=0)&&(u+v<=1)) P.push_back(X[i]);
  }

}

//  SortAndRemoveDoubles: orders polygon corners in P counter clock wise and removes duplicates

template<int dim, typename T>
void OverlappingMerge<dim, T>::SortAndRemoveDoubles( std::vector<Dune::FieldVector<T,dim> > & P )
{
  std::vector<Dune::FieldVector<double,dim> >   Q ;
  Dune::FieldVector<double,dim>  c,d ;
  double eps=1.e-10 ;

  std::vector< double > air,ail ;
  std::vector< int >    iit,iir,iil,iib ;

  // build barycentre c of all these points
  c = 0 ;
  for (int im=0; im < P.size(); im++) c += P[im] ;
  c /= P.size() ;

  // sort the points respectively into 4 sets : strict right half plan (r for right), half top line x=0 (t for top)
  //                                            strict left half plan (l for left),  half bottom line x=0 (b for bottom)
  // air & ail are the rates of vectors \vec{cP} for each point P, to be used next to sort the cloud of points in each sets.

  for (int im=0; im < P.size(); im++) {
    d = P[im] - c ;
    if      (d[0]>0)     { air.push_back(d[1]/d[0]) ; iir.push_back(im) ; }
    else if (d[0]<0)     { ail.push_back(d[1]/d[0]) ; iil.push_back(im) ; }
    else {  if (d[1]>=0) {                            iit.push_back(im) ; }
            else         {                            iib.push_back(im) ; }   }
  }

  // sort the cloud of points in each half plans.

  for (int j=1; j < air.size(); j++)
    for (int i=0; i < j; i++) if (air[j]<air[i]) { std::swap<double>(air[i],air[j]) ;
                                                   std::swap<int>   (iir[i],iir[j]) ; }

  for (int j=1; j < ail.size(); j++)
    for (int i=0; i < j; i++) if (ail[j]<ail[i]) { std::swap<double>(ail[i],ail[j]) ;
                                                   std::swap<int>   (iil[i],iil[j]) ; }

  // filling the new set of points, respectively with sets r,t,l and b.
  // Only one point can be taken from t and b, since any other would be an indesirable double.

  for (int j=0; j < iir.size(); j++) Q.push_back(P[iir[j]]) ;
  if  (iit.size()>0) Q.push_back(P[iit[0]]) ;
  for (int j=0; j < iil.size(); j++) Q.push_back(P[iil[j]]) ;
  if  (iib.size()>0) Q.push_back(P[iib[0]]) ;

  //  Remove Doubles, we might use the comand  Q.erase(std::unique(Q.begin(), Q.end()), Q.end());
  //  however it works only to eliminate exact doubles,
  //  and we also want to eliminate very close neighbouring (in a distance eps).

  P.clear() ;
  P.push_back(Q[0]) ;
  for (int j=1; j < Q.size(); j++)  { if ((P[P.size()-1] - Q[j]).infinity_norm()>eps) P.push_back(Q[j]) ;}

}
