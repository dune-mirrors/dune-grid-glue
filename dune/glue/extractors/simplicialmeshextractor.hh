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

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                                Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  static const Dune::PartitionIteratorType PType = Dune::All_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim0Extractor<GV>::SubEntityInfo SubEntityInfo;
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
  // Get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the corners.

  // free everything there is in this object
  this->clear();

  // several counter for consecutive indexing are needed
  size_t element_index = 0;
  size_t vertex_index = 0;

  // a temporary container where newly acquired face
  // information can be stored at first
  std::deque<SubEntityInfo> temp_faces;

  // iterate over all codim 0 elements on the grid
  for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
  {

    IndexType eindex = this->indexSet().template index<0>(*elit);

    // only do sth. if this element is "interesting"
    // implicit cast is done automatically
    if (descr.contains(elit))
    {
      // add an entry to the element info map, the index will be set properly later
      this->_elmtInfo[eindex] = new ElementInfo(element_index, elit, 1);

      int numCorners = elit->template count<dim>();
      unsigned int vertex_indices[numCorners];       // index in global vector
      unsigned int vertex_numbers[numCorners];       // index in parent entity

      // try for each of the faces vertices whether it is already inserted or not
      for (int i = 0; i < numCorners; ++i)
      {
        vertex_numbers[i] = i;

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

      // add a new face to the temporary collection
      temp_faces.push_back(SubEntityInfo(eindex,0,elit->type()));
      element_index++;
      for (int i=0; i<numCorners; i++) {
        temp_faces.back().corners[i].idx = vertex_indices[i];
        // remember the vertices' numbers in parent element's vertices
        temp_faces.back().corners[i].num = vertex_numbers[i];
      }

    }
  }   // end loop over elements

  // allocate the array for the face specific information...
  this->subEntities_.resize(element_index);
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

#endif // SIMPLICIALMESHEXTRACTOR_HH_
