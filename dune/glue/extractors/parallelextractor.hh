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

#if HAVE_MPI
#include <dune/common/mpicollectivecommunication.hh>
#endif

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

  typedef typename LX::GridView GV;
  typedef typename LX::GridView GridView;
  typedef typename LX::Coords Coords;

  typedef typename GV::Grid Grid;
  typedef typename GV::Grid::ctype ctype;
  typedef typename LX::VertexVector VertexVector;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

  typedef typename GV::Grid::GlobalIdSet::IdType GlobalId;
  typedef std::pair<GlobalId, size_t>                              GlobalFaceId;

  typedef typename LX::Geometry Geometry;
  typedef typename LX::LocalGeometry LocalGeometry;

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
  // typedef std::vector<GlobalId> GlobalEntityTopology;
  struct GlobalFaceInfo
  {
    // GlobalEntityTopology v;
    GlobalId v[(1<<(dim-codim))];
    char corners;
    GlobalFaceId i;
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

  const GV& gv_;
  LX lx_;
  // map between local and global _face_ index
  std::vector<unsigned int>  local2global_;
  std::vector<unsigned int>  global2local_;
  // global data
  std::vector<Coords> coords_;
  std::vector<VertexVector> faces_;
  std::vector<Dune::GeometryType> geometryTypes_;

  Dune::GeometryType guessGeomType(int corners) const
  {
    Dune::GeometryType type;
    switch (dim-codim)
    {
    case 0 :
      type.makeVertex();
      break;
    case 1 :
      type.makeLine();
      break;
    case 2 :
      switch (corners) {
      case 3 :
        type.makeTriangle();
        break;
      case 4 :
        type.makeQuadrilateral();
        break;
      default :
        DUNE_THROW(Dune::Exception,
                   "impossible intersection (dim=" <<
                   dim-codim << ", " << corners << " corners)");
      }
      break;
    case 3 :
      switch (corners) {
      case 4 :
        type.makeTetrahedron();
        break;
      case 5 :
        type.makePyramid();
        break;
      case 6 :
        type.makePrism();
        break;
      case 8 :
        type.makeHexahedron();
        break;
      default :
        DUNE_THROW(Dune::Exception,
                   "impossible intersection (dim=" <<
                   dim-codim << ", " << corners << " corners)");
      }
      break;
    default :
      DUNE_THROW(Dune::NotImplemented,
                 "Only intersections up to dim=3 are supported");
    }
    return type;
  }

public:

  // methods
  bool contains (unsigned int global, unsigned int & local) const
  {
    local = global2local_[global];
    // assert(local == global);
    return (local != (unsigned int)(-1));
  }

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  ParallelExtractor(const GV& gv)
    :  gv_(gv), lx_(gv)
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
    lx_.clear();
    global2local_.clear();
  }

  // TODO: doku, get Descriptor Type from LX
  template<typename Descriptor>
  void update(const Descriptor& descr)
  {
    // setup local extractor
    lx_.update(descr);

    // obtain local data
    std::vector<GlobalCoordInfo> localCoordInfos;
    std::vector<GlobalFaceInfo> localFaceInfos;
    std::vector<Dune::GeometryType> localGeometryTypes;

    // get vertex data
    {
      std::vector<Coords> coords;
      lx_.getCoords(coords);
      localCoordInfos.resize(coords.size());
      for (size_t i=0; i<coords.size(); i++)
      {
        localCoordInfos[i].c = coords[i];
        localCoordInfos[i].i = gv_.grid().globalIdSet().id(* lx_.vertex(i));
        localCoordInfos[i].valid = true;
      }
    }

    // get face data
    {
      std::vector<VertexVector> faces;
      lx_.getFaces(faces);
      localFaceInfos.resize(faces.size());
      size_t Xsub = 0;
      for (size_t i=0; i<faces.size(); i++)
      {
        size_t corners = faces[i].size();
        // localFaceInfos[i].v.resize(corners);
        localFaceInfos[i].corners = corners;
        for (size_t v=0; v<corners; v++)
        {
          localFaceInfos[i].v[v] = localCoordInfos[faces[i][v]].i;
        }
        localFaceInfos[i].i.first = gv_.grid().globalIdSet().id(* lx_.element(i));
        if (i>0 && localFaceInfos[i].i.first != localFaceInfos[i-1].i.first)
          Xsub = 0;
        localFaceInfos[i].i.second = Xsub++;
        localFaceInfos[i].valid = true;
      }
    }

    // merge parallel data
    std::vector<GlobalCoordInfo> globalCoordInfos;
    std::vector<GlobalFaceInfo>  globalFaceInfos;
    std::vector<Dune::GeometryType> globalGeometryTypes;

    // communicate coordinates
    // TODO: use more efficient communication
    typedef typename Grid::CollectiveCommunication Comm;
    const Comm & comm = gv_.grid().comm();
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

    // sort and shrink coordinate vector
    std::sort(globalCoordInfos.begin(), globalCoordInfos.end());
    {
      typename std::vector<GlobalCoordInfo>::iterator where =
        std::unique(globalCoordInfos.begin(), globalCoordInfos.end());
      --where;
      while(where->valid == false) --where;
      ++where;
      globalCoordInfos.erase(where, globalCoordInfos.end());
    }
    assert(globalCoordInfos.front().valid == true);
    globalCoordSize = globalCoordInfos.size();

    // sort and shrink face vector
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

    // sort and shrink geometry types vector
    std::sort(globalGeometryTypes.begin(), globalGeometryTypes.end());
    {
      typename std::vector<Dune::GeometryType>::iterator where =
        std::unique(globalGeometryTypes.begin(), globalGeometryTypes.end());
      globalGeometryTypes.erase(where, globalGeometryTypes.end());
    }
    geometryTypes_ = globalGeometryTypes;

    // setup parallel coords and faces
    {
      std::map<GlobalId, unsigned int> coordIndex;       // quick access to vertex index, give an ID
      coords_.resize(globalCoordSize);
      faces_.resize(globalFaceSize);

      for (size_t c=0; c<globalCoordSize; ++c)
      {
        coordIndex[globalCoordInfos[c].i] = c;
        coords_[c] = globalCoordInfos[c].c;
      }

      for (size_t f=0; f<globalFaceSize; ++f)
      {
        // size_t corners = globalFaceInfos[f].v.size();
        size_t corners = globalFaceInfos[f].corners;
        faces_[f].resize(corners);
        for (size_t c=0; c<corners; ++c)
          faces_[f][c] = coordIndex[globalFaceInfos[f].v[c]];
      }
    }

    // setup local/global mappings
    {
      // create a temporary map for quick lookup
      std::map<GlobalFaceId, unsigned int> globalIndex;
      for (size_t f=0; f<globalFaceInfos.size(); ++f)
      {
        globalIndex[globalFaceInfos[f].i] = f;
      }
      // "copy" map to _local2global
      local2global_.resize(localFaceInfos.size());
      for (size_t i = 0; i<localFaceInfos.size(); i++)
      {
        local2global_[i] = globalIndex[localFaceInfos[i].i];
      }
    }
    {
      // create a temporary map for quick lookup
      std::map<GlobalFaceId, unsigned int> localIndex;
      for (size_t f=0; f<localFaceInfos.size(); ++f)
      {
        localIndex[localFaceInfos[f].i] = f;
      }
      // "copy" map to _global2local
      // not all entries are contained in the map, if not the entry is "-1"
      global2local_.resize(globalFaceInfos.size());
      for (size_t i = 0; i<globalFaceInfos.size(); i++)
      {
        typename std::map<GlobalFaceId, unsigned int>::iterator where =
          localIndex.find(globalFaceInfos[i].i);
        if (where != localIndex.end())
          global2local_[i] = where->second;
        else
          global2local_[i] = (unsigned int)-1;
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
    coords.resize(coords_.size());
    std::copy(coords_.begin(), coords_.end(), coords.begin());
  }


  /**
   * @brief getter for the count of coordinates
   * @return the count
   */
  unsigned int nCoords() const
  {
    return coords_.size();
  }

  /**
   * @brief Get the list of geometry types
   */
  void getGeometryTypes(std::vector<Dune::GeometryType>& geometryTypes) const
  {
    geometryTypes.resize(faces_.size());
    for (size_t i=0; i<faces_.size(); i++)
      geometryTypes[i] = guessGeomType(faces_[i].size());
  }

  /**
   * @brief getter for the indices array
   * It is strongly recommended not to modify its contents.
   * Deallocation is done in this class.
   * @return the _indices array
   */
  void getFaces(std::vector<VertexVector>& faces) const
  {
    faces.resize(this->faces_.size());
    for (size_t i = 0; i < this->faces_.size(); ++i)
      faces[i] = this->faces_[i];
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
    if (lx_.faceIndices(e, first, count))
    {
      // convert local index to global index
      first = local2global_[first];
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
    return lx_.indexInInside(l_index);
  }


  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return lx_.indexSet();
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
    return lx_.element(l_index);
  }

  /** \brief Get World geometry of the extracted face */
  Geometry geometry(unsigned int index) const
  {
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return lx_.geometry(l_index);
  }

  /** \brief Get Geometry of the extracted face in element coordinates */
  Geometry geometryLocal(unsigned int index) const
  {
    unsigned int l_index = 0;
    bool have = contains(index, l_index);
    assert(have);
    return lx_.geometryLocal(l_index);
  }

  // export codim0 method
  bool & positiveNormalDirection()
  {
    return lx_.positiveNormalDirection();
  }
};

#endif // DUNE_CODIM_1_EXTRACTOR_HH_
