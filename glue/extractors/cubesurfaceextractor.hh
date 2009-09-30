// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    CubeSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Apr 25, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for cube surface grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief grid extractor implementation for hypercube surface grids
 */

#ifndef CUBESURFACEEXTRACTOR_HH_
#define CUBESURFACEEXTRACTOR_HH_


#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>
#include "surfacedescriptor.hh"


/**
 * @class CubeSurfaceExtractor
 * @brief grid extractor implementation for cube surface grids
 *
 * Provides methods that build topology information for given grids.
 * The template parameters
 * @li GV the grid class type
 */
template<typename GV>
class CubeSurfaceExtractor
  : public Codim1Extractor<GV>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    dimworld = GV::dimensionworld
  };

  enum
  {
    dim = GV::dimension
  };

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dim
  };

  enum
  {
    cube_corners = 1 << (dim-1)
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                           Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  static const Dune::PartitionIteratorType PType = Dune::All_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim1Extractor<GV>::FaceInfo FaceInfo;

public:

  using Codim1Extractor<GV>::codim;

  typedef Dune::GenericGeometry::BasicGeometry<dim-codim, Dune::GenericGeometry::DefaultGeometryTraits<ctype,dim-codim,dimworld> > Geometry;

  typedef Dune::GenericGeometry::BasicGeometry<dim-codim, Dune::GenericGeometry::DefaultGeometryTraits<ctype,dim-codim,dim> > LocalGeometry;

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  CubeSurfaceExtractor(const GV& gv)
    : Codim1Extractor<GV>(gv)
  {
    std::cout << "This is CubeSurfaceExtractor on a <"
              << GV::dimension << "," << GV::dimensionworld
              << "> grid!" << std::endl;
  }

  /*  F U N C T I O N A L I T Y  */


  /**
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a mapper or set object
   * @c m specifying an index set with codimension 0 entities near or on the boundary.
   * It is assumed that only one geometric shape exists on the boundary.
   * The template parameter n then denotes the number of corners per boundary element
   * (e.g. n==3 for triangles in a 3D grid of tetrahedra).
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are cubes.
   * @param descr a descriptor that "selects" the faces to add to the surface
   */
  void update(const FaceDescriptor<GV>& descr);

  /** \brief Get world geometry of the extracted face */
  Geometry geometry(unsigned int index) const;

  /** \brief Get geometry of the extracted face in element coordinates */
  LocalGeometry geometryLocal(unsigned int index) const;

}; // end of class CubeSurfaceExtractor



template<typename GV>
void CubeSurfaceExtractor<GV>::update(const FaceDescriptor<GV>& descr)
{
  // Currently the extractor splits (if necessary) quadrilaterals into triangles (because psurface)
  // can only handle triangles). Since this is done manually we cannot handle grids
  // of dimension higher than 3.
  dune_static_assert(dim<=3, "CubeSurfaceExtractor only implemented for 1d, 2d, 3d");

  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.

  // free everything there is in this object
  this->clear();

  // several counter for consecutive indexing are needed
  size_t simplex_index = 0;
  size_t vertex_index = 0;

  // a temporary container where newly acquired face
  // information can be stored at first
  std::deque<FaceInfo> temp_faces;

  // iterate over all codim 0 elemets on the grid
  for (ElementIter elit = this->_gv.template begin<0, PType>();
       elit != this->_gv.template end<0, PType>(); ++elit)
  {
    // count number of selected faces for this element
    size_t simplex_count = 0;

    IndexType eindex = this->indexSet().template index<0>(*elit);

    // iterate over all intersections of codim 1 and test if the
    // boundary intersections are to be added to the surface
    for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
    {
      // only look at boundary faces
      if (is->boundary() && descr.contains(elit, is->indexInInside())) {

        // Make sure the face is a cube
        if (!is->type().isCube())
          DUNE_THROW(Dune::GridError, "found non-cube boundary entity: " << is->type());

        unsigned int face_index = is->indexInInside();

        // now we only have to care about the 3D case, i.e. the quadrilateral
        // face has to be divided into two triangles
        unsigned int vertex_indices[4];
        unsigned int vertex_numbers[4];

        // get the vertex pointers for the quadrilateral's corner vertices
        // and try for each of them whether it is already inserted or not
        for (int i = 0; i < cube_corners; ++i)
        {
          // get the number of the vertex in the parent element
          vertex_numbers[i] = orientedSubface<dim>(elit->type(), face_index, i);

          // get the vertex pointer and the index from the index set
          VertexPtr vptr(elit->template subEntity<dim>(vertex_numbers[i]));
          IndexType vindex = this->indexSet().template index<dim>(*vptr);

          // if the vertex is not yet inserted in the vertex info map
          // it is a new one -> it will be inserted now!
          typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
          if (vimit == this->_vtxInfo.end())
          {
            // insert into the map
            this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, vptr);
            // remember this vertex' index
            vertex_indices[i] = vertex_index;
            // increase the current index
            vertex_index++;
          }
          else
          {
            // only remember the vertex' index
            vertex_indices[i] = vimit->second->idx;
          }
        }

        // add a new face to the temporary collection
        simplex_count ++;
        temp_faces.push_back(
          FaceInfo(simplex_index, eindex, face_index, FaceInfo::first));
        for (int i=0; i<dim; i++) {
          temp_faces.back().corners[i].idx = vertex_indices[i];
          // remember the vertices' numbers in parent element's vertices
          temp_faces.back().corners[i].num = vertex_numbers[i];
        }

        // now introduce the second triangle subdividing the quadrilateral
        // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
        // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1
        if (dim==3) {
          // add a new face to the temporary collection for the second tri
          simplex_count++;
          temp_faces.push_back(
            FaceInfo(simplex_index+1, eindex, face_index, FaceInfo::second));
          temp_faces.back().corners[0].idx = vertex_indices[3];
          temp_faces.back().corners[1].idx = vertex_indices[2];
          temp_faces.back().corners[2].idx = vertex_indices[1];
          // remember the vertices' numbers in parent element's vertices
          temp_faces.back().corners[0].num = vertex_numbers[3];
          temp_faces.back().corners[1].num = vertex_numbers[2];
          temp_faces.back().corners[2].num = vertex_numbers[1];
        }
      }        // end if is contained

    }             // end loop over intersections

    // add element entry with the offset of the first simplex and the number of simplices
    if (simplex_count != 0)
      this->_elmtInfo[eindex] = new typename Codim1Extractor<GV>::ElementInfo(simplex_index, elit, simplex_count);
    simplex_index += simplex_count;

  }       // end loop over elements

  // allocate the array for the face specific information...
  this->_faces.resize(simplex_index);
  // ...and fill in the data from the temporary containers
  copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());

  // now first write the array with the coordinates...
  this->_coords.resize(this->_vtxInfo.size());
  typename Codim1Extractor<GV>::VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    typename Codim1Extractor<GV>::CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }

}


/** \brief Get World geometry of the extracted face */
template<typename GV>
typename CubeSurfaceExtractor<GV>::Geometry CubeSurfaceExtractor<GV>::geometry(unsigned int index) const
{
  std::vector<Coords> corners(simplex_corners);
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;

  return Geometry(Dune::GeometryType(Dune::GeometryType::simplex,dim-codim), corners);
}


/** \brief Get Geometry of the extracted face in element coordinates */
template<typename GV>
typename CubeSurfaceExtractor<GV>::LocalGeometry CubeSurfaceExtractor<GV>::geometryLocal(unsigned int index) const
{
  std::vector<Coords> corners(simplex_corners);

  unsigned int num_in_self = this->indexInInside(index);
  Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->type();

  // computing the locals is straight forward for flat rectangles,
  // we only need the triangle's corners in element coordinate
  if (this->_faces[index].category == 1) {
    // the triangle's corners are (0 1 2) using face indices
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimworld>(gt, num_in_self, i);

  } else {
    // the triangle's corners are (3 2 1) using face indices
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimworld>(gt, num_in_self, 3-i);
  }

  return LocalGeometry(Dune::GeometryType(Dune::GeometryType::simplex,dim-codim), corners);
}

#endif // CUBESURFACEEXTRACTOR_HH_
