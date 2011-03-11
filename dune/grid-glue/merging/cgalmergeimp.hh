// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * \brief Implementation of the Merger concept using the CGAL library
 */

#ifndef CGAL_MERGE_IMP_HH
#define CGAL_MERGE_IMP_HH

#include <iostream>
#include <vector>
#include <algorithm>
#include <bitset>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/grid/common/genericreferenceelements.hh>
#include <dune/grid/common/grid.hh>

#include <dune/grid/genericgeometry/geometry.hh>

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


/** \brief Implementation of the Merger concept using the CGAL library
 *
 * This is a separate implementation class to keep all CGAL stuff out of the cgalmerge.hh header.

   \tparam dim Grid dimension of the coupling grids.  The world dimension is assumed to be the same.
   \tparam T Type used for coordinates
 */
template<int dim, typename T = double>
class CGALMergeImp
{

public:
  typedef CGAL::Cartesian<Number_type>   Kernel;
  typedef Kernel::Point_3 Point_3;
  typedef CGAL::Polyhedron_3<Kernel>     Polyhedron_3;
  typedef CGAL::Nef_polyhedron_3<Kernel> Nef_Polyhedron_3;

  static void makeHexahedron(Polyhedron_3& P,
                             const std::vector<Dune::FieldVector<T,dim> >& c);


  static void computeNeighborIntersections(const Dune::GeometryType& elementType,
                                           const std::vector<Dune::FieldVector<T,dim> >& elementCorners,
                                           const Nef_Polyhedron_3& NQ,
                                           const Nef_Polyhedron_3& intersection,
                                           std::bitset<(1<<dim)>& neighborIntersects
                                           );



private:

};

#ifdef CGAL_EXTRA_TYPES
#define CGAL_EXTERN
#include "cgalmergeimp.cc"
#undef CGAL_EXTERN
#endif

#endif // CGAL_MERGE_IMP_HH
