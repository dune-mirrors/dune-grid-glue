// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CGAL_EXTERN
#include "config.h"
#endif

#include <dune/grid-glue/merging/cgalmerge.hh>
#include <dune/grid-glue/merging/cgalmergeimp.hh>

#ifdef HAVE_CGAL  // without CGAL we can still handle 1d problems
// 2d
//#include "bso_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <list>


#include <CGAL/Gmpz.h>
#include <CGAL/Extended_homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#endif

#ifdef HAVE_CGAL

#ifdef CGAL_USE_GMP

// GMP is installed. Use the GMP rational number-type.
  #include <CGAL/Gmpq.h>

typedef CGAL::Gmpq Number_type;

#else

// GMP is not installed. Use CGAL's exact rational number-type.
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif
#endif

/* IMPLEMENTATION */

template<int dim, typename T>
void CGALMerge<dim, T>::
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

  // /////////////////////////////////////////////////////////////////////////////////////
  //   Compute the intersection between the two elements.  The 1d case is implemented
  //   by hand; 2d and 3d use CGAL.
  // /////////////////////////////////////////////////////////////////////////////////////

  switch (dim) {
  case 1 : {

    CGALMergeImp<dim,T>::compute1dIntersection(grid1Geometry, grid1ElementCorners, grid1Index,
                                               grid2Geometry, grid2ElementCorners, grid2Index,
                                               this->intersections_
                                               );

    break;
  }
#ifdef HAVE_CGAL
  case 2 : {

    CGALMergeImp<dim,T>::compute2dIntersection(grid1ElementType, grid1Geometry, grid1ElementCorners, grid1Index,
                                               grid2ElementType, grid2Geometry, grid2ElementCorners, grid2Index,
                                               this->intersections_
                                               );

    break;
  }
  case 3 : {

    typedef CGAL::Cartesian<Number_type>   Kernel;
    typedef typename Kernel::Point_3 Point_3;
    typedef CGAL::Polyhedron_3<Kernel>     Polyhedron_3;
    typedef CGAL::Nef_polyhedron_3<Kernel> Nef_Polyhedron_3;

    // Construct the two input polyhedra.
    Polyhedron_3 P;

    if (grid1ElementType.isSimplex()) {
      P.make_tetrahedron(Point_3 (grid1ElementCorners[0][0], grid1ElementCorners[0][1], grid1ElementCorners[0][2]),
                         Point_3 (grid1ElementCorners[1][0], grid1ElementCorners[1][1], grid1ElementCorners[1][2]),
                         Point_3 (grid1ElementCorners[2][0], grid1ElementCorners[2][1], grid1ElementCorners[2][2]),
                         Point_3 (grid1ElementCorners[3][0], grid1ElementCorners[3][1], grid1ElementCorners[3][2]));
    } if (grid1ElementType.isCube()) {

      CGALMergeImp<dim,T>::makeHexahedron(P, grid1ElementCorners);

    } else
      DUNE_THROW(Dune::GridError, "Type " << grid1ElementType << " not supported by CGALMerge yet");
    //std::cout << "P = " << P << std::endl;

    // Turn polyhedron into a Nef polyhedron
    Nef_Polyhedron_3 NP(P);
    Polyhedron_3 Q;

    if (grid2ElementType.isSimplex()) {
      Q.make_tetrahedron(Point_3 (grid2ElementCorners[0][0], grid2ElementCorners[0][1], grid2ElementCorners[0][2]),
                         Point_3 (grid2ElementCorners[1][0], grid2ElementCorners[1][1], grid2ElementCorners[1][2]),
                         Point_3 (grid2ElementCorners[2][0], grid2ElementCorners[2][1], grid2ElementCorners[2][2]),
                         Point_3 (grid2ElementCorners[3][0], grid2ElementCorners[3][1], grid2ElementCorners[3][2]));
    } if (grid2ElementType.isCube()) {

      CGALMergeImp<dim,T>::makeHexahedron(Q, grid2ElementCorners);

    } else
      DUNE_THROW(Dune::GridError, "Type " << grid2ElementType << " not supported by CGALMerge yet");

    Nef_Polyhedron_3 NQ(Q);
    //std::cout << "Q = " << Q << std::endl;

    //////////////////////////////////////////////////////////
    // Compute the intersection of P and Q.
    //////////////////////////////////////////////////////////
    Nef_Polyhedron_3 intersection = NP * NQ;

    if (intersection.is_empty())
      return;

    Polyhedron_3 intersectionP;
    if(intersection.is_simple()) {
      intersection.convert_to_polyhedron(intersectionP);
      //std::cout << intersectionP;
    } else
      std::cerr << "N1 is not a 2-manifold." << std::endl;

    //std::cout << "Intersection has " << intersection.number_of_vertices() << " vertices\n";

    //////////////////////////////////////////////////////////////////////
    //  Compute for each face of each element whether it will also intersect
    //////////////////////////////////////////////////////////////////////

    CGALMergeImp<dim,T>::computeNeighborIntersections(grid1ElementType, grid1ElementCorners, NQ, intersection, neighborIntersects1);

    CGALMergeImp<dim,T>::computeNeighborIntersections(grid2ElementType, grid2ElementCorners, NP, intersection, neighborIntersects2);

    //////////////////////////////////////////////////////////////////////
    //  Triangulate the intersection polyhedron.
    //  For this we first compute the centroid to have a point that
    //  is certainly within the intersection.  (this of course presupposes
    //  that the intersection is convex).  Then we connect the centroid
    //  to all facets.
    //////////////////////////////////////////////////////////////////////

    assert(intersectionP.is_closed());

    // If the intersection is a tetrahedron we take a shortcut
    if (intersectionP.is_tetrahedron(intersectionP.halfedges_begin())) {

      // Make Dune types from CGAL types

      // The following check is there, because we are using CGAL::to_double().
      // If necessary the code can be generalized somewhat.
      dune_static_assert((Dune::is_same<T,double>::value), "T must be 'double'");

      Dune::array<Dune::FieldVector<T,dim>, 4> duneVertexPos;

      typename Polyhedron_3::Point_iterator vIt = intersectionP.points_begin();

      for (int i=0; i<4; ++i, ++vIt)
        for (int j=0; j<3; j++)
          duneVertexPos[i][j] = CGAL::to_double((*vIt)[j]);

      assert(vIt==intersectionP.points_end());

      // ///////////////////////////////////////////////////
      // Output the tetrahedron
      // ///////////////////////////////////////////////////

      this->intersections_.push_back(RemoteSimplicialIntersection<T,dim,dim,dim>());

      // Compute local coordinates in the grid1 and grid2 elements
      for (int i=0; i<4; i++) {
        this->intersections_.back().grid1Local_[i] = grid1Geometry.local(duneVertexPos[i]);
        this->intersections_.back().grid2Local_[i] = grid2Geometry.local(duneVertexPos[i]);
      }

      // Set indices
      this->intersections_.back().grid1Entity_ = grid1Index;
      this->intersections_.back().grid2Entity_ = grid2Index;

      //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

    } else {

      //////////////////////////////////////////////////////////////////
      //   Compute the centroid
      //////////////////////////////////////////////////////////////////

      Dune::FieldVector<T,dim> centroid(0);

      for (typename Polyhedron_3::Point_iterator vIt = intersectionP.points_begin();
           vIt != intersectionP.points_end();
           ++vIt) {
        for (int i=0; i<dim; i++)
          centroid[i] += CGAL::to_double((*vIt)[i]);
      }

      centroid /= intersectionP.size_of_vertices();

      //////////////////////////////////////////////////////////////////////////
      //  Loop over all facets, triangulate the facet and create
      //  an intersection from each such triangle together with the centroid
      //////////////////////////////////////////////////////////////////////////

      for (typename Polyhedron_3::Facet_iterator fIt = intersectionP.facets_begin();
           fIt != intersectionP.facets_end();
           ++fIt) {

        ////////////////////////////////////////////////////////////
        // Get the vertices of this facet in circular order
        ////////////////////////////////////////////////////////////
        unsigned int nVertices = fIt->facet_degree();

        std::vector<Dune::FieldVector<T,dim> > facetVertices(nVertices);

        typename Polyhedron_3::Facet::Halfedge_around_facet_circulator fVCirc = fIt->facet_begin();
        assert (fVCirc != 0);          // an empty circulator would be an error
        int i=0;           // array iterator

        do {

          for (int j=0; j<dim; j++)
            facetVertices[i][j] = CGAL::to_double(fVCirc->vertex()->point()[j]);

          ++i;

        } while (++fVCirc != fIt->facet_begin());

        /////////////////////////////////////////////////////////////////////////
        //  Triangulate the facet and enter an intersection for each triangle
        /////////////////////////////////////////////////////////////////////////

        typename Polyhedron_3::Facet::Halfedge_around_facet_circulator anchor = fIt->facet_begin();

        typename Polyhedron_3::Facet::Halfedge_around_facet_circulator next = anchor;
        ++next;

        typename Polyhedron_3::Facet::Halfedge_around_facet_circulator nextNext = next;
        ++nextNext;

        do {

          // Make Dune types from CGAL types

          // The following check is there, because we are using CGAL::to_double().
          // If necessary the code can be generalized somewhat.
          dune_static_assert((Dune::is_same<T,double>::value), "T must be 'double'");

          Dune::FieldVector<T,dim> anchorFieldVector;
          Dune::FieldVector<T,dim> nextFieldVector;
          Dune::FieldVector<T,dim> nextNextFieldVector;

          for (int i=0; i<dim; i++) {
            anchorFieldVector[i]   = CGAL::to_double(anchor->vertex()->point()[i]);
            nextFieldVector[i]     = CGAL::to_double(next->vertex()->point()[i]);
            nextNextFieldVector[i] = CGAL::to_double(nextNext->vertex()->point()[i]);
          }

          // ////////////////////////////////////////////////////////////
          // Output the tetrahedron (anchor, next, nextNext, centroid)
          // ////////////////////////////////////////////////////////////

          this->intersections_.push_back(RemoteSimplicialIntersection<T,dim,dim,dim>());

          // Compute local coordinates in the grid1 element
          this->intersections_.back().grid1Local_[0] = grid1Geometry.local(anchorFieldVector);
          this->intersections_.back().grid1Local_[1] = grid1Geometry.local(nextFieldVector);
          this->intersections_.back().grid1Local_[2] = grid1Geometry.local(nextNextFieldVector);
          this->intersections_.back().grid1Local_[3] = grid1Geometry.local(centroid);

          // Compute local coordinates in the grid1 element
          this->intersections_.back().grid2Local_[0] = grid2Geometry.local(anchorFieldVector);
          this->intersections_.back().grid2Local_[1] = grid2Geometry.local(nextFieldVector);
          this->intersections_.back().grid2Local_[2] = grid2Geometry.local(nextNextFieldVector);
          this->intersections_.back().grid2Local_[3] = grid2Geometry.local(centroid);

          // Set indices
          this->intersections_.back().grid1Entity_ = grid1Index;
          this->intersections_.back().grid2Entity_ = grid2Index;

          //std::cout << "Intersection between elements " << grid1Index << " and " << grid2Index << std::endl;

          // move to the next triangle
          ++next;
          ++nextNext;

        } while (nextNext != fIt->facet_begin());

      }

    }

  }
  break;
#endif

  default :
    DUNE_THROW(Dune::NotImplemented, "CGALMerge is not implemented for dim==" << dim << "!");

  }

}

template<int dim, typename T>
inline bool CGALMerge<dim, T>::grid1SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid1 entities
  return true;
}


template<int dim, typename T>
inline bool CGALMerge<dim, T>::grid2SimplexMatched(unsigned int idx) const
{
  // naive: we assume that there is a partner for all grid2 entities
  return true;
}


// Explicit instantiation
#ifdef CGAL_EXTERN
#define DECL extern
#else
#define DECL
#endif
#define CGAL_INSTANTIATION(D,T)                                          \
  DECL template void CGALMerge<D, T>::computeIntersection(const Dune::GeometryType& grid1ElementType, \
                                                          const std::vector<Dune::FieldVector<T,D> >& grid1ElementCorners, \
                                                          unsigned int grid1Index, \
                                                          std::bitset<(1<<D)>& neighborIntersects1, \
                                                          const Dune::GeometryType& grid2ElementType, \
                                                          const std::vector<Dune::FieldVector<T,D> >& grid2ElementCorners, \
                                                          unsigned int grid2Index, \
                                                          std::bitset<(1<<D)>& neighborIntersects2); \
\
  DECL template bool CGALMerge<D, T>::grid1SimplexMatched(unsigned int idx) const; \
\
  DECL template bool CGALMerge<D, T>::grid2SimplexMatched(unsigned int idx) const

CGAL_INSTANTIATION(1, double);
CGAL_INSTANTIATION(2, double);
CGAL_INSTANTIATION(3, double);
// CGAL_INSTANTIATION(1, float);
// CGAL_INSTANTIATION(2, float);
// CGAL_INSTANTIATION(3, float);

#undef DECL
