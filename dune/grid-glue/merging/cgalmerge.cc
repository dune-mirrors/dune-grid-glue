// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef CGAL_EXTERN
#include "config.h"
#endif

#include <dune/grid-glue/merging/cgalmerge.hh>
#include <dune/grid-glue/merging/cgalmergeimp.hh>

#if HAVE_CGAL  // without CGAL we can still handle 1d problems
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

#ifdef CGAL_USE_GMP

// GMP is installed. Use the GMP rational number-type.
  #include <CGAL/Gmpq.h>

typedef CGAL::Gmpq Number_type;

#else

// GMP is not installed. Use CGAL's exact rational number-type.
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

typedef CGAL::Quotient<CGAL::MP_Float>                Number_type;

#endif // CGAL_USE_GMP
#else
typedef double Number_type;
#endif // HAVE_CGAL

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

    CGALMergeImp<dim,T,Number_type>::compute1dIntersection(grid1Geometry, grid1ElementCorners, grid1Index,
                                                           grid2Geometry, grid2ElementCorners, grid2Index,
                                                           this->intersections_
                                                           );

    break;
  }
#if HAVE_CGAL
  case 2 : {

    CGALMergeImp<dim,T,Number_type>::compute2dIntersection(grid1ElementType, grid1Geometry, grid1ElementCorners, grid1Index,
                                                           grid2ElementType, grid2Geometry, grid2ElementCorners, grid2Index,
                                                           this->intersections_
                                                           );

    break;
  }
  case 3 : {

    CGALMergeImp<dim,T,Number_type>::compute3dIntersection(grid1ElementType, grid1Geometry, grid1ElementCorners, grid1Index, neighborIntersects1,
                                                           grid2ElementType, grid2Geometry, grid2ElementCorners, grid2Index, neighborIntersects2,
                                                           this->intersections_
                                                           );

    break;
  }
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
#ifndef CGAL_EXTERN
template class CGALMerge<1,double>;
template class CGALMerge<2,double>;
template class CGALMerge<3,double>;
#endif
