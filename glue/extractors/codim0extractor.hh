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

  /** \brief This class extracts codim-0 stuff (elements) */
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
      :       self(vtxindex_), index(index_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType self;

    /// @brief the coordinate
    Coords coord;

    /// @brief the index of this coordinate (in internal storage scheme) // NEEDED??
    unsigned int index : 28;

  };

  /**
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

  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /// @brief the grid object to extract the surface from
  const GV&                                              _gv;


  /*        Geometrical and Topological Information                */

  /// @brief all information about the corner vertices of the extracted
  std::vector<CoordinateInfo>   _coords;

  /// @brief a map enabling faster access to vertices and coordinates
  ///
  /// Maps a vertex' index (from index set) to an object holding the locally
  /// associated index of the vertex' coordinate in _coords and an entity
  /// pointer to the codim<dim> entity.
  VertexInfoMap _vtxInfo;

  /// @brief a map enabling faster access to elements and faces
  ///
  /// Maps an element's index (from index set) to an object holding the locally
  /// associated index of its first face in _indices (if there are more they are
  /// positioned consecutively) and an entity pointer to the codim<0> entity.
  ElementInfoMap _elmtInfo;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  Codim0Extractor(const GV& gv) :
    _gv(gv)
  {}

  ~Codim0Extractor()
  {
    // only the objects that have been allocated manually have to be
    // deallocated manually again
    // free all the manually allocated memory
    for (typename VertexInfoMap::iterator it = _vtxInfo.begin(); it != _vtxInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    for (typename ElementInfoMap::iterator it = _elmtInfo.begin(); it != _elmtInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
  }

  void clear()
  {
    // this is an inofficial way on how to free the memory allocated
    // by a std::vector
    std::vector<typename Codim0Extractor<GV>::CoordinateInfo> dummy;
    this->_coords.swap(dummy);

    // first free all manually allocated vertex/element info items...
    for (typename VertexInfoMap::iterator it = _vtxInfo.begin(); it != _vtxInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    for (typename ElementInfoMap::iterator it = _elmtInfo.begin(); it != _elmtInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;

    // ...then clear the maps themselves, too
    _vtxInfo.clear();
    _elmtInfo.clear();
  }

  bool contains (unsigned int global, unsigned int & local) const
  {
    local = global;
    return true;
  }

  /*  G E T T E R S  */

  /**
   * @brief getter for the coordinates array
   * @param coords a vector that will be resized (!) and filled with the coordinates
   */
  void getCoords(std::vector<Dune::FieldVector<ctype, dimworld> >& coords) const
  {
    coords.resize(this->_coords.size());
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
      coords[i] = this->_coords[i].coord;
  }

  /**
   * @brief getter for the count of coordinates
   * @return the count
   */
  unsigned int nCoords() const
  {
    return this->_coords.size();
  }


  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return this->_gv.indexSet();
  }

};



#endif // DUNE_CODIM_0_EXTRACTOR_HH_
