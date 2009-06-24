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


#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>




/**
 * @brief provides static methods for codim 0 grid extraction
 *
 * Provides methods that build topology information for given grids.

   \tparam GV the grid view type
 */
template<typename GV>
class Codim0Extractor
{
  /** \todo This should rather be protected */
public:

  enum {dimworld = GV::dimensionworld};
  enum {dim      = GV::dimension};

  /** \brief This class extracts codim-1 stuff (surfaces) */
  enum {codim    = 0};

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;


  /************************** PRIVATE SUBCLASSES **********************/

  struct CoordinateInfo
  {
    CoordinateInfo()
    {}

    CoordinateInfo(unsigned int index_, IndexType vtxindex_)
      :       self(vtxindex_), index(index_), num_faces(0), faces(NULL)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType self;

    /// @brief the coordinate
    Coords coord;

    /// @brief the index of this coordinate (in internal storage scheme) // NEEDED??
    unsigned int index : 28;

    /// @brief the number of extracted faces with this coord as corner,
    /// major purpose is holding the length of the array "faces"
    unsigned int num_faces : 4;

    /// @brief holds the indices of the faces of which the coordinate is a corner
    unsigned int* faces;
  };

  /**
   * @class EntityInfo
   * @brief simple struct holding a vertex pointer and an index
   */
  struct VertexInfo
  {
    VertexInfo(unsigned int idx_, VertexPtr& p_) : idx(idx_), p(p_)
    {}
    unsigned int idx;
    VertexPtr p;
  };

  /**
   * @class ElementInfo
   * @brief simple struct holding an entity pointer and an index
   */
  struct ElementInfo
  {
    ElementInfo(unsigned int idx_, ElementPtr& p_) : idx(idx_), p(p_)
    {}

    /// @brief the index of this element's first face in the internal list of extracted faces
    unsigned int idx;

    /// @brief the entity pointer to the element
    ElementPtr p;
  };

  typedef std::map<IndexType, ElementInfo* >  ElementInfoMap;
  typedef std::map<IndexType, VertexInfo* >   VertexInfoMap;

};



#endif // DUNE_CODIM_0_EXTRACTOR_HH_
