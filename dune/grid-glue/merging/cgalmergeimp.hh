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

#include <dune/grid-glue/merging/remotesimplicialintersection.hh>

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
#endif // HAVE_CGAL


/** \brief Implementation of the Merger concept using the CGAL library
 *
 * This is a separate implementation class to keep all CGAL stuff out of the cgalmerge.hh header.

   \tparam dim Grid dimension of the coupling grids.  The world dimension is assumed to be the same.
   \tparam Dune_number_type Type used by Dune for coordinates
   \tparam CGAL_number_type Type used by CGAL for coordinates
 */
template<int dim, class Dune_number_type, class CGAL_number_type>
class CGALMergeImp
{

public:
#if HAVE_CGAL
  typedef CGAL::Cartesian<CGAL_number_type>   Kernel;

  typedef typename Kernel::Point_2 Point_2;
  typedef CGAL::Polygon_2<Kernel>                    Polygon_2;
  typedef CGAL::Polygon_with_holes_2<Kernel>         Polygon_with_holes_2;

  typedef typename Kernel::Point_3 Point_3;
  typedef CGAL::Polyhedron_3<Kernel>     Polyhedron_3;
  typedef CGAL::Nef_polyhedron_3<Kernel> Nef_Polyhedron_3;

  static void makeHexahedron(Polyhedron_3& P,
                             const std::vector<Dune::FieldVector<Dune_number_type,dim> >& c);


  static void computeNeighborIntersections(const Dune::GeometryType& elementType,
                                           const std::vector<Dune::FieldVector<Dune_number_type,dim> >& elementCorners,
                                           const Nef_Polyhedron_3& NQ,
                                           const Nef_Polyhedron_3& intersection,
                                           std::bitset<(1<<dim)>& neighborIntersects
                                           );
#endif // HAVE_CGAL

  static void compute1dIntersection(const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                    unsigned int grid1Index,
                                    const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                    unsigned int grid2Index,
                                    std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                    );

  static void compute2dIntersection(const Dune::GeometryType& grid1ElementType,
                                    const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                    unsigned int grid1Index,
                                    const Dune::GeometryType& grid2ElementType,
                                    const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                    unsigned int grid2Index,
                                    std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                    );

  static void compute3dIntersection(const Dune::GeometryType& grid1ElementType,
                                    const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid1Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid1ElementCorners,
                                    unsigned int grid1Index,
                                    std::bitset<(1<<dim)>& neighborIntersects1,
                                    const Dune::GeometryType& grid2ElementType,
                                    const Dune::GenericGeometry::BasicGeometry<dim, Dune::GenericGeometry::DefaultGeometryTraits<Dune_number_type,dim,dim> >& grid2Geometry,
                                    const std::vector<Dune::FieldVector<Dune_number_type,dim> >& grid2ElementCorners,
                                    unsigned int grid2Index,
                                    std::bitset<(1<<dim)>& neighborIntersects2,
                                    std::vector<RemoteSimplicialIntersection<Dune_number_type,dim,dim,dim> >& intersections
                                    );

private:

};

#ifdef CGAL_EXTRA_TYPES
#define CGAL_EXTERN
#include "cgalmergeimp.cc"
#undef CGAL_EXTERN
#endif

#endif // CGAL_MERGE_IMP_HH
