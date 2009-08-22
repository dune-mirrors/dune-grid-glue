// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SimplicialSurfaceExtractor.hh
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
 * @brief Surface grid extractor base class
 */

#ifndef DUNE_CODIM_1_EXTRACTOR_HH
#define DUNE_CODIM_1_EXTRACTOR_HH


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
template<typename GV>
class Codim1Extractor
{
  /** \todo This should rather be protected */
public:

  enum {dimworld = GV::dimensionworld};
  enum {dim      = GV::dimension};

  /** \brief This class extracts codim-1 stuff (surfaces) */
  enum {codim    = 1};

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


  /************************** PRIVATE SUBCLASSES **********************/

  /**
   * @brief Helpful struct holding one index for the coordinate (vertex)
   * to which it is associated and the element's corner index;
   */
  struct CornerInfo
  {
    unsigned int idx : 28;   ///< index of the vertex
    unsigned int num : 4;   ///< element corner
  };

  struct CoordinateInfo
  {
    CoordinateInfo()
    {}

    CoordinateInfo(unsigned int index_, IndexType vtxindex_)
      :   vtxindex(vtxindex_), index(index_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType vtxindex;

    /// @brief the coordinate
    Coords coord;

    /// @brief the index of this coordinate (in internal storage scheme) // NEEDED??
    unsigned int index;
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
    ElementInfo(unsigned int idx_, ElementPtr& p_, unsigned int f_) : idx(idx_), faces(f_), p(p_)
    {}

    /// @brief the index of this element's first face in the internal list of extracted faces
    unsigned int idx : 28;

    /// @brief the number of extracted faces for this element
    unsigned int faces : 4;

    /// @brief the entity pointer to the element
    ElementPtr p;
  };


  /**
   * @class FaceInfo
   * @brief simple struct holding some packed information about a codimension 1 entity
   * to its parent codim 0 entity
   */
  struct FaceInfo
  {
    FaceInfo()
    {}

    FaceInfo(unsigned int index_, IndexType parent_, unsigned int num_in_parent_, unsigned int category_)
      :       parent(parent_), index(index_), num_in_parent(num_in_parent_), category(category_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType parent;

    /// @brief the index of this face (in internal storage scheme) // NEEDED??
    unsigned int index : 27;

    /// @brief the number of the face in the parent element
    unsigned int num_in_parent : 3;

    /// @brief flag telling whether this is the first of two triangles
    /// refining the associated face

    /// 0: triangle face
    /// 1: first of two tris refining a quadrilateral
    /// 2: second of two tris refining a quadrilateral
    unsigned int category : 2;

    enum {
      triangle = 0,
      first = 1,
      second = 2
    };

    /// @brief the corner indices plus the numbers of the vertices in the parent element
    CornerInfo corners[dim];     // sim = numer of vertices in a simplex
  };


  typedef std::map<IndexType, ElementInfo* >  ElementInfoMap;
  typedef std::map<IndexType, VertexInfo* >   VertexInfoMap;


  /************************** MEMBER VARIABLES ************************/

  /// @brief the grid object to extract the surface from
  const GV&                                       _gv;

  /*        Geometrical and Topological Information                */

  /// @brief all information about the corner vertices of the extracted
  std::vector<CoordinateInfo>   _coords;

  /// @brief all information about the extracted faces
  std::vector<FaceInfo>         _faces;

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
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  Codim1Extractor(const GV& gv)
    :  _gv(gv)
  {}

  /** \brief Destructor frees allocated memory */
  ~Codim1Extractor();

  /*  F U N C T I O N A L I T Y  */

  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear()
  {
    // this is an inofficial way on how to free the memory allocated
    // by a std::vector
    {
      std::vector<typename Codim1Extractor<GV>::CoordinateInfo> dummy;
      this->_coords.swap(dummy);
    }
    {
      std::vector<FaceInfo> dummy;
      this->_faces.swap(dummy);
    }

    // first free all manually allocated vertex/element info items...
    for (typename Codim1Extractor<GV>::VertexInfoMap::iterator it = this->_vtxInfo.begin();
         it != this->_vtxInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    for (typename Codim1Extractor<GV>::ElementInfoMap::iterator it = this->_elmtInfo.begin();
         it != this->_elmtInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    // ...then clear the maps themselves, too
    this->_vtxInfo.clear();
    this->_elmtInfo.clear();
  }


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
      coords[i] = _coords[i].coord;
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
        faces[i][j] = this->_faces[i].corners[j].idx;
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
    typename Codim1Extractor<GV>::ElementInfoMap::const_iterator it =
      this->_elmtInfo.find(this->indexSet().template index<0>(e));
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
  int indexInInside(unsigned int index) const
  {
    return index < this->_faces.size() ? this->_faces[index].num_in_parent : -1;
  }

  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return _gv.indexSet();
  }


  bool contains (unsigned int global, unsigned int & local) const
  {
    local = global;
    return true;
  }

#if 0
  /**
   * @brief getter for the index of an entity of codim cc
   * @return the index specified by the grid's index set
   */
  template<int cc>
  IndexType index(const typename GV::Traits::template Codim<cc>::Entity& e) const
  {
    return indexSet().template index<cc>(e);
  }


  /**
   * @brief gets index of coordinate in _coords associated with given vertex
   * @return the index if possible, -1 else
   */
  int coordinateIndex(const Vertex& v) const
  {
    typename VertexInfoMap::const_iterator it = _vtxInfo.find(index<dim>(v));
    return (it == _vtxInfo.end()) ? -1 : it->second->idx;
  }


  /**
   * @brief gets for each vertex corner of given face (by index) the number of
   * the vertex in parent element's ordering
   * @param index the face's index
   * @param corner the index of the corner
   * @return -1 <=> index invalid or array not filled, else index
   */
  int numCornerInParent(unsigned int index, unsigned int corner) const
  {
    return (index >= this->_faces.size() || corner >= simplex_corners) ?
           -1 : this->_faces[index].corners[corner].num;
  }
#endif

  /**
   * @brief gets the parent element for a given face index,
   * throws an exception if index not valid
   * @param index the index of the face
   * @return a reference to the element's stored pointer
   */
  const ElementPtr& element(unsigned int index) const
  {
    if (index >= this->_faces.size())
      DUNE_THROW(Dune::GridError, "invalid face index");
    return (this->_elmtInfo.find(this->_faces[index].parent))->second->p;
  }

#if 1
  /**
   * @brief gets the vertex for a given coordinate index
   * throws an exception if index not valid
   * @param index the index of the coordinate
   * @return a reference to the vertex' stored pointer
   */
  const VertexPtr& vertex(unsigned int index) const
  {
    if (index >= this->_coords.size())
      DUNE_THROW(Dune::GridError, "invalid coordinate index");
    return (this->_vtxInfo.find(this->_coords[index].vtxindex))->second->p;
  }
#endif
};


template<typename GV>
Codim1Extractor<GV>::~Codim1Extractor()
{
  clear();
}




#endif // DUNE_CODIM_1_EXTRACTOR_HH_
