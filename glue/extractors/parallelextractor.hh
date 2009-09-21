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

#include <dune/common/collectivecommunication.hh>
#include <stdlib.h>

template <typename COMM>
class CommHelper
{
public:
  CommHelper (const COMM & c) {}

  template<typename BUF>
  void allgather(const BUF & in, BUF & out, int)
  {
    out = in;
  }
};

#if HAVE_MPI
template <>
class CommHelper< Dune::CollectiveCommunication<MPI_Comm> >
{
  typedef Dune::CollectiveCommunication<MPI_Comm> COMM;
  const COMM & comm;
public:
  CommHelper (const COMM & c) : comm(c) {}

  template<typename BUF>
  void allgather(BUF & in, BUF & out, int outsz)
  {
    typedef typename BUF::value_type V;
    MPI_Allgather(&in[0], in.size(),
                  Dune::Generic_MPI_Datatype<V>::get(),
                  &out[0], outsz,
                  Dune::Generic_MPI_Datatype<V>::get(), comm);
  }
};
#endif

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
  enum { dimw     = LX::dimworld };
  enum { codim    = LX::codim };
  enum { simplex_corners = LX::simplex_corners };

  typedef typename LX::GridView GV;
  typedef typename LX::GridView GridView;
  typedef typename LX::Coords Coords;

  typedef typename GV::Grid Grid;
  typedef typename GV::Grid::ctype ctype;
  typedef typename LX::SimplexTopology SimplexTopology;
  typedef typename LX::FaceInfo FaceInfo;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  typedef typename GV::Grid::GlobalIdSet::IdType GlobalId;

  // types for the global setup
  struct GlobalCoordInfo
  {
    Coords c;
    GlobalId i;
    bool valid;

    GlobalCoordInfo() : valid(false) {}

    bool operator < (const GlobalCoordInfo & other) const
    {
      return valid && (!other.valid || i < other.i);
      if (!valid) return false;
      if (!other.valid) return true;
      return i < other.i;
    }

    bool operator == (const GlobalCoordInfo & other) const
    {
      return (! valid && ! other.valid) || (valid && other.valid && i == other.i);
    }
  };
  typedef Dune::array<GlobalId, simplex_corners> GlobalSimplexTopology;
  struct GlobalFaceInfo
  {
    GlobalSimplexTopology v;
    GlobalId i;
    bool valid;

    GlobalFaceInfo() : valid(false) {}

    bool operator < (const GlobalFaceInfo & other) const
    {
      if (!valid) return false;
      if (!other.valid) return true;
      return i < other.i;
    }

    bool operator == (const GlobalFaceInfo & other) const
    {
      return (! valid && ! other.valid) || (valid && other.valid && i == other.i);
    }
  };

private:

  const GV& _gv;
  LX _lx;
  std::vector<unsigned int>  _local2global;
  std::vector<unsigned int>  _global2local;
  // global data
  std::vector<Coords> _coords;
  std::vector<SimplexTopology> _faces;
  std::vector<GlobalId> _faceId;

public:

  // methods
  bool contains (unsigned int global, unsigned int & local) const
  {
    local = _global2local[global];
    return (local != (unsigned int)(-1));
  }

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  ParallelExtractor(const GV& gv)
    :  _gv(gv), _lx(gv)
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
    _lx.clear();
    _global2local.clear();
  }

  // TODO: doku, get Descriptor Type from LX
  template<typename Descriptor>
  void update(const Descriptor& descr)
  {
    // setup local extractor
    _lx.update(descr);

    // obtain local data
    std::vector<GlobalCoordInfo> localCoordInfos;
    std::vector<GlobalFaceInfo> localFaceInfos;

    // get vertex data
    {
      std::vector<Coords> coords;
      _lx.getCoords(coords);
      localCoordInfos.resize(coords.size());
      for (unsigned int i=0; i<coords.size(); i++)
      {
        localCoordInfos[i].c = coords[i];
        localCoordInfos[i].i = _gv.grid().globalIdSet().id(* _lx.vertex(i));
        localCoordInfos[i].valid = true;
      }
    }

    // get face data
    {
      std::vector<SimplexTopology> faces;
      _lx.getFaces(faces);
      localFaceInfos.resize(faces.size());
      for (unsigned int i=0; i<faces.size(); i++)
      {
        for (int v=0; v<simplex_corners; v++)
        {
          localFaceInfos[i].v[v] = localCoordInfos[faces[i][v]].i;
        }
        localFaceInfos[i].i = _gv.grid().globalIdSet().id(* _lx.element(i));
        localFaceInfos[i].valid = true;
      }
    }

    // merge parallel data
    std::vector<GlobalCoordInfo> globalCoordInfos;
    std::vector<GlobalFaceInfo>  globalFaceInfos;

    // communicate coordinates
    // TODO: use more efficient communication
    typedef typename Grid::CollectiveCommunication Comm;
    const Comm & comm = _gv.grid().comm();
    CommHelper<Comm> commHelper(comm);

    size_t globalCoordLocalSize = localCoordInfos.size();
    globalCoordLocalSize = comm.max(globalCoordLocalSize) + 2;
    size_t globalCoordSize = globalCoordLocalSize * comm.size();;
    globalCoordInfos.resize(globalCoordSize);
    commHelper.allgather(localCoordInfos, globalCoordInfos, globalCoordLocalSize);

    // communicate faces
    // TODO: use more efficient communication
    size_t globalFaceLocalSize = localFaceInfos.size();
    globalFaceLocalSize = comm.max(globalFaceLocalSize);
    size_t globalFaceSize = globalFaceLocalSize * comm.size();
    globalFaceInfos.resize(globalFaceSize);
    commHelper.allgather(localFaceInfos, globalFaceInfos, globalFaceLocalSize);

    // sort and shrink vectors
    {
      std::sort(globalCoordInfos.begin(), globalCoordInfos.end());
      typename std::vector<GlobalCoordInfo>::iterator where =
        std::unique(globalCoordInfos.begin(), globalCoordInfos.end());
      --where;
      while(where->valid == false) --where;
      ++where;
      globalCoordInfos.erase(where, globalCoordInfos.end());
    }
    assert(globalCoordInfos.front().valid == true);
    globalCoordSize = globalCoordInfos.size();

    std::sort(globalFaceInfos.begin(), globalFaceInfos.end());
    {
      typename std::vector<GlobalFaceInfo>::iterator where =
        std::unique(globalFaceInfos.begin(), globalFaceInfos.end());
      --where;
      while(where->valid == false) --where;
      ++where;
      globalFaceInfos.erase(where, globalFaceInfos.end());
    }
    globalFaceSize = globalFaceInfos.size();

    // setup parallel coords and faces
    {
      std::map<GlobalId, unsigned int> coordIndex;       // quick access to vertex index, give an ID
      _coords.resize(globalCoordSize);
      _faces.resize(globalFaceSize);
      _faceId.resize(globalFaceSize);

      for (size_t c=0; c<globalCoordSize; ++c)
      {
        coordIndex[globalCoordInfos[c].i] = c;
        _coords[c] = globalCoordInfos[c].c;
      }

      for (size_t f=0; f<globalFaceSize; ++f)
      {
        for (size_t c=0; c<simplex_corners; ++c)
          _faces[f][c] = coordIndex[globalFaceInfos[f].v[c]];
        _faceId[f] = globalFaceInfos[f].i;
      }
    }

    // setup local/global mappings
    {
      // create a temporary map for quick lookup
      std::map<GlobalId, unsigned int> globalIndex;
      for (size_t f=0; f<globalFaceInfos.size(); ++f)
      {
        globalIndex[globalFaceInfos[f].i] = f;
      }
      // "copy" map to _local2global
      _local2global.resize(localFaceInfos.size());
      for (unsigned int i = 0; i<localFaceInfos.size(); i++)
      {
        _local2global[i] = globalIndex[localFaceInfos[i].i];
      }
    }
    {
      // create a temporary map for quick lookup
      std::map<GlobalId, unsigned int> localIndex;
      for (size_t f=0; f<localFaceInfos.size(); ++f)
      {
        localIndex[localFaceInfos[f].i] = f;
      }
      // "copy" map to _global2local
      // not all entries are contained in the map, if not the entry is "-1"
      _global2local.resize(globalFaceInfos.size());
      for (unsigned int i = 0; i<globalFaceInfos.size(); i++)
      {
        typename std::map<GlobalId, unsigned int>::iterator where =
          localIndex.find(globalFaceInfos[i].i);
        if (where != localIndex.end())
          _global2local[i] = where->second;
        else
          _global2local[i] = (unsigned int)-1;
      }
    }
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
    if (_lx.faceIndices(e, first, count))
    {
      // convert local index to global index
      first = _local2global[first];
      return true;
    }
    return false;
  }


  /**
   * @brief gets the number face in the parent element
   * @param index the index of the face
   * @return if failed -1, else the index
   */
  int indexInInside(unsigned int index) const
  {
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.indexInInside(l_index);
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
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
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
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.vertex(l_index);
  }

#if 0

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
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
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
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
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
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.localAndGlobalCoords(l_index, bcoords, ecoords, wcoords);
  }

#endif

  /**
   * @brief for several given barycentric coords in a simplex compute world coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param wcoords to be filled with world coordinates
   */
  template<typename BCoordContainer, typename CoordContainer>
  void globalCoords(unsigned int index, const BCoordContainer &bcoords, CoordContainer &wcoords) const
  {
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.globalCoords(l_index, bcoords, wcoords);
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
  template<typename BCoordContainer, typename CoordContainer>
  void localCoords(unsigned int index, const BCoordContainer &bcoords, CoordContainer &ecoords) const
  {
    int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.localCoords(l_index, bcoords, ecoords);
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
  template<typename BCoordContainer, typename ECoordContainer, typename WCoordContainer>
  void localAndGlobalCoords(unsigned int index, const BCoordContainer &bcoords, ECoordContainer &ecoords, WCoordContainer &wcoords) const
  {
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return _lx.localAndGlobalCoords(l_index, bcoords, ecoords, wcoords);
  }

};

#endif // DUNE_CODIM_1_EXTRACTOR_HH_
