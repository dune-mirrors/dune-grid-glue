// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SimplicialSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Jan 15, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for simplicial surface grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief grid extractor implementation for simplicial surface grids
 */

#ifndef SIMPLICIALSURFACEEXTRACTOR_HH
#define SIMPLICIALSURFACEEXTRACTOR_HH


#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/genericgeometry/geometry.hh>
#include "surfacedescriptor.hh"

#include <dune/glue/extractors/codim1extractor.hh>



/**
 * @brief provides static methods for grid surface extraction
 *
 * Provides methods that build topology information for given grids.
 * \tparam GV the grid view  type
 */
template<typename GV>
class SimplicialSurfaceExtractor
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
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim1Extractor<GV>::FaceInfo FaceInfo;

public:

  using Codim1Extractor<GV>::codim;

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  SimplicialSurfaceExtractor(const GV& gv)
    : Codim1Extractor<GV>(gv)
  {
    std::cout << "This is SimplicialSurfaceExtractor on a <" << GV::dimension
              << "," << GV::dimensionworld << "> grid!" << std::endl;
  }

  /*  F U N C T I O N A L I T Y  */


  /**
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a decider function
   * that returns a boolean value for each codim 0 entity's every face.
   * Note that the function is only called for faces that are part of the domain boundary.
   * Inner faces cannot be part of the surface.
   * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If we the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are simplices.
   * @param descr a descriptor that "selects" the faces to add to the surface
   */
  void update(const FaceDescriptor<GV>& descr);

}; // end of class SimplicialSurfaceExtractor


template<typename GV>
void SimplicialSurfaceExtractor<GV>::update(const FaceDescriptor<GV>& descr)
{
  // free everything there is in this object
  this->clear();

  // In this first pass iterate over all entities of codim 0.
  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.

  {
    // several counter for consecutive indexing are needed
    int simplex_index = 0;
    int vertex_index = 0;
    IndexType eindex = 0;             // supress warning

    // needed later for insertion into a std::set which only
    // works with const references

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<FaceInfo> temp_faces;

    // iterate over all codim 0 elemets on the grid
    for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
    {

      // remember the indices of the faces that shall become
      // part of the surface
      std::set<int> boundary_faces;

      // iterate over all intersections of codim 1 and test if the
      // boundary intersections are to be added to the surface
      for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
      {
        // only look at boundary faces
        if (is->boundary() && descr.contains(elit, is->indexInInside())) {

          // Make sure the face is a simplex
          if (!is->type().isSimplex())
            DUNE_THROW(Dune::GridError, "found non-simplicial boundary entity: " << is->type());

          boundary_faces.insert(is->indexInInside());
        }

      }

      // if some face is part of the surface add it!
      if (boundary_faces.size() != 0)
      {
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new typename Codim1Extractor<GV>::ElementInfo(simplex_index, elit, boundary_faces.size());

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (typename std::set<int>::const_iterator sit = boundary_faces.begin(); sit != boundary_faces.end(); ++sit)
        {
          // add a new face to the temporary collection
          temp_faces.push_back(FaceInfo(eindex, *sit));

          // try for each of the faces vertices whether it is already inserted or not
          for (int i = 0; i < simplex_corners; ++i)
          {
            // get the number of the vertex in the parent element
            int vertex_number = orientedSubface<dim>(elit->type(), *sit, i);

            // get the vertex pointer and the index from the index set
            VertexPtr vptr(elit->template subEntity<dim>(vertex_number));
            IndexType vindex = this->indexSet().template index<dim>(*vptr);

            // remember the vertex' number in parent element's vertices
            temp_faces.back().corners[i].num = vertex_number;

            // if the vertex is not yet inserted in the vertex info map
            // it is a new one -> it will be inserted now!
            typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
            if (vimit == this->_vtxInfo.end())
            {
              // insert into the map
              this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, vptr);
              // remember the vertex as a corner of the current face in temp_faces
              temp_faces.back().corners[i].idx = vertex_index;
              // increase the current index
              vertex_index++;
            }
            else
            {
              // only insert the index into the simplices array
              temp_faces.back().corners[i].idx = vimit->second->idx;
            }
          }

          // now increase the current face index
          simplex_index++;

        }                         // end loop over found surface parts
      }
    }             // end loop over elements

    // allocate the array for the face specific information...
    this->_faces.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());
  }


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

#endif // SIMPLICIALSURFACEEXTRACTOR_HH_
