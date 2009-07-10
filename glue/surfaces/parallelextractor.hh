// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    parallelextractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: extractor for parallel grids, uses a local extractor
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief parallel grid extractor
 */

#ifndef DUNE_PARALLEL_EXTRACTOR_HH
#define DUNE_PARALLEL_EXTRACTOR_HH

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
 * @brief provides static methods for grid surface extraction
 *
 * Provides methods that build topology information for given grids.

   \tparam GV the grid view type
 */
template<typename LX>
class ParallelExtractor
{
  /** \todo This should rather be protected */
public:

  typedef LX LocalExtractor;

  // retrieve typedefs etc. from LocalExtractor
  enum { dimworld = LX::dimworld };
  enum { dim      = LX::dim };
  enum { codim    = LX::codim };

  typedef LX::GridView GridView;
  typedef LX::Coords Coords;
  typedef LX::GridView GridView;
  typedef LX::GridView GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;
  typedef Dune::array<unsigned int, simplex_corners>               SimplexTopology;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

private:

  LX & _lx;
  std::map<unsigned int, unsigned int>  _global2local;
  // global data
  std::vector<Coords> _coords;
  std::vector<SimplexTopology> _cells;
  std::vector<GlobalID> _globalId;

  // methods
  bool trymap (unsigned int global, unsigned int & local) const
  {
    std::map<unsigned int, unsigned int>::iterrator i = _global2local.find(global);
    if (i != data.end())
    {
      local = i->second;
      return true;
    }
    return false;
  }

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  ParallelExtractor(const LX& lx)
    :  _lx(lx)
  {}

  /** \brief Destructor frees allocated memory */
  ~ParallelExtractor()
  {
    clear();
  }

  /*  F U N C T I O N A L I T Y  */

  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear()
  {
    map.clear();
    assert(false);     // actually it would be a bit strange to call clear of the parallel extractor.
  }

  // TODO: doku, get Descriptor Type from LX
  template<typename Descriptor>
  void update(const Descriptor& descr)
  {
    _lx.update(descr);
    _lx.getCoords(_coords);
    _lx.getFaces(_faces);

    // TODO: communicate coordinates
    // TODO: communicate elements
    // TODO: set parallel info for elements
  };

  /*  G E T T E R S  */

  /**
   * @brief getter for the coordinates array
   * @param coords a vector that will be resized (!) and filled with the coordinates,
   * note that the single components are written consecutively
   */
  void getCoords(std::vector<Dune::FieldVector<ctype, dimworld> >& coords) const
  {
    coords.resize(_coords.size());
    for (unsigned int i = 0; i < _coords.size(); ++i)
      coords[i] = _coords[i];
  }


  /**
   * @brief getter for the count of coordinates
   * @return the count
   */
  unsigned int nCoords() const
  {
    return _coords.size();
  }


  /**
   * @brief getter for the indices array
   * It is strongly recommended not to modify its contents.
   * Deallocation is done in this class.
   * @return the _indices array
   */
  void getFaces(std::vector<SimplexTopology>& faces) const
  {
    faces.resize(this->_faces.size());
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
      for (unsigned int j = 0; j < simplex_corners; ++j)
        faces[i][j] = this->_faces[i][j];
  }


  /**
   * @brief gets index of first face as well as the total number of faces that
   * were extracted from this element
   * @param e the element
   * @param first will contain the first index if found, else -1
   * @param count will contain the number of faces if found, else 0
   * @return success
   */
  bool faceIndices(const Element& e, int& first, int& count) const
  {
    typename Codim1Extractor<GV>::ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->template index<0>(e));
    if (it == this->_elmtInfo.end())
    {
      first = -1;
      count = 0;
      return false;
    }
    // the iterator is valid, fill the out params
    first = it->second->idx;
    count = it->second->faces;
    return true;
  }


  /**
   * @brief gets the number face in the parent element
   * @param index the index of the face
   * @return if failed -1, else the index
   */
  int numberInSelf(unsigned int index) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.numberInSelf(l_index);
  }


  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return _lx.indexSet();
  }


  /**
   * @brief gets the parent element for a given face index,
   * throws an exception if index not valid
   * @param index the index of the face
   * @return a reference to the element's stored pointer
   */
  const ElementPtr& element(unsigned int index) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.element(l_index);
  }


  /**
   * @brief gets the vertex for a given coordinate index
   * throws an exception if index not valid
   * @param index the index of the coordinate
   * @return a reference to the vertex' stored pointer
   */
  const VertexPtr& vertex(unsigned int index) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.vertex(l_index);
  }


  /**
   * @brief for given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  void globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.globalCoords(l_index, bcoords, wcoords);
  }


  /**
   * @brief for given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  void localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.localCoords(l_index, bcoords, ecoords);
  }


  /**
   * @brief for given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.localAndGlobalCoords(l_index, bcoords, ecoords, wcoords);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  template<typename CoordContainer>
  void globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.globalCoords(l_index, bcoords, wcoords, size);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  template<typename CoordContainer>
  void localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, int size) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.localCoords(l_index, bcoords, ecoords, size);
  }


  /**
   * @brief for several given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   * @return
   */
  template<typename CoordContainer>
  void localAndGlobalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, CoordContainer &wcoords, int size) const
  {
    int l_index = 0;
    bool have = trymap(index, l_index);
    assert(have);
    return _lx.localAndGlobalCoords(l_index, bcoords, ecoords, wcoords, size);
  }

};

#endif // DUNE_CODIM_1_EXTRACTOR_HH_
