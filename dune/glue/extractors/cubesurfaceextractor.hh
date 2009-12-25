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
#include <dune/grid/common/genericreferenceelements.hh>
#include "../misc/geometry.hh"
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

  using Codim1Extractor<GV>::codim;
  using Codim1Extractor<GV>::dim;
  using Codim1Extractor<GV>::dimworld;
  using Codim1Extractor<GV>::simplex_corners;
  using Codim1Extractor<GV>::cube_corners;

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  static const Dune::PartitionIteratorType PType = Dune::All_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim1Extractor<GV>::SubEntityInfo SubEntityInfo;
  typedef typename Codim1Extractor<GV>::ElementInfo ElementInfo;
  typedef typename Codim1Extractor<GV>::CoordinateInfo CoordinateInfo;
  typedef typename Codim1Extractor<GV>::VertexInfo VertexInfo;
  typedef typename Codim1Extractor<GV>::VertexInfoMap VertexInfoMap;

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
  std::deque<SubEntityInfo> temp_faces;

  // iterate over all codim 0 elemets on the grid
  for (ElementIter elit = this->_gv.template begin<0, PType>();
       elit != this->_gv.template end<0, PType>(); ++elit)
  {
    // count number of selected faces for this element
    size_t simplex_count = 0;

    // element index
    IndexType eindex = this->indexSet().template index<0>(*elit);

    // make sure each face is handle only once
    std::set<unsigned int> faces;

    // cell reference element
    const Dune::GenericReferenceElement<ctype, dim> & refElem =
      Dune::GenericReferenceElements<ctype, dim>::general(elit->type());

    // iterate over all intersections of codim 1 and test if the
    // boundary intersections are to be added to the surface
    for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
    {
      unsigned int face_index = is->indexInInside();

      // only look at boundary faces
      if (is->boundary() &&
          faces.find(face_index) == faces.end() &&
          descr.contains(elit, face_index)) {

        // Make sure the face is a cube
        if (refElem.size(face_index,1,dim) != cube_corners)
          DUNE_THROW(Dune::GridError, "found non-cube boundary entity (" << refElem.size(face_index,1,dim) << " corners)");

        unsigned int vertex_indices[cube_corners];         // index in global vector
        unsigned int vertex_numbers[cube_corners];         // index in parent entity

        // get the vertex pointers for the quadrilateral's corner vertices
        // and try for each of them whether it is already inserted or not
        for (int i = 0; i < cube_corners; ++i)
        {
          // get the number of the vertex in the parent element
          // vertex_numbers[i] = refElem.subEntity(face_index,1,i,dim);
          vertex_numbers[i] = orientedSubface<dim>(elit->type(), face_index, i);

          // get the vertex pointer and the index from the index set
          VertexPtr vptr(elit->template subEntity<dim>(vertex_numbers[i]));
          IndexType vindex = this->indexSet().template index<dim>(*vptr);

          // if the vertex is not yet inserted in the vertex info map
          // it is a new one -> it will be inserted now!
          typename VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
          if (vimit == this->_vtxInfo.end())
          {
            // insert into the map
            this->_vtxInfo[vindex] = new VertexInfo(vertex_index, vptr);
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

        // check normal vector & flip face if necessary
#if 0
        {
          static Dune::FieldVector<ctype, dim-1> c(0.5);
          static Dune::array<Dune::FieldVector<ctype, dim>, dim-1> A;
          for (int i=0; i<dim-1; i++)
            A[i] = elit->geometry().corner(vertex_numbers[1<<i]) -
                   // here we use some expert knowledge      ^^^^
                   // the ref element and about psurface
                   elit->geometry().corner(vertex_numbers[0]);
          Dune::FieldVector<ctype, dim> Nis = is->outerNormal(c);
          Dune::FieldVector<ctype, dim> Nface = computeNormal(A);
          assert(Nis * Nface != 0.0);
          if (Nis * Nface < 0.0)
          {
            // do the flip
            for (int i = 0; i < cube_corners; i+=2)
            {
              // swap i and i+1
              int x = vertex_indices[i];
              vertex_indices[i] = vertex_indices[i+1];
              vertex_indices[i+1] = x;

              x = vertex_numbers[i];
              vertex_numbers[i] = vertex_numbers[i+1];
              vertex_numbers[i+1] = x;
            }
          }
        }
#endif

        // add a new face to the temporary collection
        simplex_count ++;
        temp_faces.push_back(SubEntityInfo(eindex, face_index));
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
          temp_faces.push_back(SubEntityInfo(eindex, face_index));
          temp_faces.back().corners[0].idx = vertex_indices[3];
          temp_faces.back().corners[1].idx = vertex_indices[2];
          temp_faces.back().corners[2].idx = vertex_indices[1];
          // remember the vertices' numbers in parent element's vertices
          temp_faces.back().corners[0].num = vertex_numbers[3];
          temp_faces.back().corners[1].num = vertex_numbers[2];
          temp_faces.back().corners[2].num = vertex_numbers[1];
        }
      }        // end if is contained

      // remember we already checked this face
      faces.insert(face_index);

    }             // end loop over intersections

    // add element entry with the offset of the first simplex and the number of simplices
    if (simplex_count != 0)
      this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit, simplex_count);
    simplex_index += simplex_count;

  }       // end loop over elements

  // allocate the array for the face specific information...
  this->subEntities_.resize(simplex_index);
  // ...and fill in the data from the temporary containers
  copy(temp_faces.begin(), temp_faces.end(), this->subEntities_.begin());

  // now first write the array with the coordinates...
  this->_coords.resize(this->_vtxInfo.size());
  typename VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }

}

#endif // CUBESURFACEEXTRACTOR_HH_
