// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    codim0extractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: base class for grid extractors extracting surface grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Mesh grid extractor base class
 */

#ifndef DUNE_CODIM_0_EXTRACTOR_HH
#define DUNE_CODIM_0_EXTRACTOR_HH

#include <deque>

#include <dune/grid/common/mcmgmapper.hh>

#include "extractor.hh"
#include "surfacedescriptor.hh"

template<typename GV>
class Codim0Extractor : public Extractor<GV,0>
{

protected:

  /** \brief Layout class needed to create a consecutive index for all element types together */
  template <int dim>
  struct P0Layout
  {
    /// Return true if gt has the same dimension as the grid
    bool contains (Dune::GeometryType gt) const
    {
      return (gt.dim()==dim);
    }
  };

public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */
  using Extractor<GV,0>::codim;
  using Extractor<GV,0>::dim;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;

  static const Dune::PartitionIteratorType PType = Dune::All_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet::IndexType IndexType;

  // import typedefs from base class
  typedef typename Codim0Extractor<GV>::SubEntityInfo SubEntityInfo;
  typedef typename Codim0Extractor<GV>::ElementInfo ElementInfo;
  typedef typename Codim0Extractor<GV>::VertexInfo VertexInfo;
  typedef typename Codim0Extractor<GV>::CoordinateInfo CoordinateInfo;
  typedef typename Codim0Extractor<GV>::VertexInfoMap VertexInfoMap;

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  Codim0Extractor(const GV& gv)
    :  Extractor<GV,0>(gv), _positiveNormalDirection(false)
  {
    std::cout << "This is Codim0Extractor on a <"
              << GV::dimension << "," << GV::dimensionworld << "> grid!" << std::endl;
  }

  /**
   */
  void update(const ElementDescriptor<GV>& descr);

  bool & positiveNormalDirection() { return _positiveNormalDirection; }
  const bool & positiveNormalDirection() const { return _positiveNormalDirection; }

protected:
  bool _positiveNormalDirection;
};


template<typename GV>
void Codim0Extractor<GV>::update(const ElementDescriptor<GV>& descr)
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

  // Make a set of consecutive indices for the elements, irrespective of the element types
  Dune::MultipleCodimMultipleGeomTypeMapper<GV, P0Layout > elementMapper(this->_gv);

  // iterate over all codim 0 elements on the grid
  for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
  {

    IndexType eindex = elementMapper.map(*elit);

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

      // flip cell if necessary
      if (_positiveNormalDirection)
      {
        // choose element type
        if (dim > elit->type().isCube())
        {
          for (int i = 0; i < (1<<dim); i+=2)
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
        if (elit->type().isSimplex())
        {
          switch ((int)dim) {
          case 1 :
          case 2 :
          {
            int x = vertex_indices[0];
            vertex_indices[0] = vertex_indices[1];
            vertex_indices[1] = x;

            x = vertex_numbers[0];
            vertex_numbers[0] = vertex_numbers[1];
            vertex_numbers[1] = x;

            break;
          }
          default :
            assert(false);
          }
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
  this->_subEntities.resize(element_index);
  // ...and fill in the data from the temporary containers
  copy(temp_faces.begin(), temp_faces.end(), this->_subEntities.begin());

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

#endif // DUNE_CODIM_0_EXTRACTOR_HH_
