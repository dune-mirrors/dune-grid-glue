// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SimplicialMeshExtractor.hh
 *  Version:     1.0
 *  Created on:  Feb 19, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for "flat" simplicial grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief
 */

#ifndef SIMPLICIALMESHEXTRACTOR_HH_
#define SIMPLICIALMESHEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>

#include "surfacedescriptor.hh"
#include <dune/glue/extractors/codim0extractor.hh>


/**
 * @brief grid extractor implementation for simplicial grids
 *
 * Provides methods that build topology information for given grids.
 * Note that these methods only operate on the grid.
 *
 * \tparam GV the grid view class type
 */
template<typename GV>
class SimplicialMeshExtractor
  : public Codim0Extractor<GV>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  using Codim0Extractor<GV>::codim;
  using Codim0Extractor<GV>::dim;
  using Codim0Extractor<GV>::dimworld;
  using Codim0Extractor<GV>::simplex_corners;
  using Codim0Extractor<GV>::cube_corners;

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                                Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  static const Dune::PartitionIteratorType PType = Dune::All_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim0Extractor<GV>::FaceInfo FaceInfo;
  typedef typename Codim0Extractor<GV>::ElementInfo ElementInfo;
  typedef typename Codim0Extractor<GV>::VertexInfo VertexInfo;
  typedef typename Codim0Extractor<GV>::CoordinateInfo CoordinateInfo;
  typedef typename Codim0Extractor<GV>::VertexInfoMap VertexInfoMap;

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  SimplicialMeshExtractor(const GV& gv)
    : Codim0Extractor<GV>(gv)
  {
    std::cout << "This is SimplicialMeshExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid!" << std::endl;
  }



  /*  F U N C T I O N A L I T Y  */

  /**
   */
  void update(const ElementDescriptor<GV>& descr);

}; // end of class SimplicialMeshExtractor


template<typename GV>
void SimplicialMeshExtractor<GV>::update(const ElementDescriptor<GV>& descr)
{
  // In this first pass iterate over all entities of codim 0.
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
  for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
  {
    // check if there are unwanted geometric shapes
    // if one appears => exit with error
    if (!elit->type().isSimplex())
      DUNE_THROW(Dune::GridError, "Expected simplicial grid but found a " << elit->type());

    IndexType eindex = this->indexSet().template index<0>(*elit);

    // only do sth. if this element is "interesting"
    // implicit cast is done automatically
    if (descr.contains(elit))
    {
      // add an entry to the element info map, the index will be set properly later
      this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit, 1);

      unsigned int vertex_indices[simplex_corners];       // index in global vector
      unsigned int vertex_numbers[simplex_corners];       // index in parent entity

      // try for each of the faces vertices whether it is already inserted or not
      for (int i = 0; i < simplex_corners; ++i)
      {
        vertex_numbers[i] = i;

        // get the vertex pointer and the index from the index set
        // Note that the orientation is always the same for all simplices,
        // i.e. CCW which is 0,1 in 2D and 0,1,2 in 3D
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

      // flip cell if necessary
      if (this->positiveNormalDirection())
      {
        assert(false);
      }

      // add a new face to the temporary collection
      temp_faces.push_back(FaceInfo(eindex,0));
      simplex_index++;
      for (int i=0; i<simplex_corners; i++) {
        temp_faces.back().corners[i].idx = vertex_indices[i];
        // remember the vertices' numbers in parent element's vertices
        temp_faces.back().corners[i].num = vertex_numbers[i];
      }

    }
  }   // end loop over elements

  // allocate the array for the face specific information...
  this->_faces.resize(simplex_index);
  // ...and fill in the data from the temporary containers
  copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());

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

#if 0
  // now add the vertices' parent faces in the _vertexFaces map.
  // therefore iterate over all indices in the _index array...
  {
    std::vector<unsigned int> refcount(this->_coords.size(), 0);

    // for each coordinate count the references in the _indices array
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
      for (unsigned int j = 0; j < simplex_corners; ++j)
        refcount[this->_faces[i].corners[j]]++;

    // allocate the right amount of storage for each vertex's references
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
    {
      // allocate an array and initialize its first element with its length
      refcount[i] = 0;                   // used as "pointer" in next loop
    }

    // add the references
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
    {
      for (unsigned int j = 0; j < simplex_corners; ++j)
      {
        unsigned int ref = this->_faces[i].corners[j];
        refcount[ref]++;
      }
    }
  }
#endif

}

#endif // SIMPLICIALMESHEXTRACTOR_HH_
