// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    extractor.hh
 *  Version:     1.0
 *  Created on:  Oct 05, 2009
 *  Author:      Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: base class for all grid extractors
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief extractor base class
 */

#ifndef DUNE_EXTRACTOR_HH
#define DUNE_EXTRACTOR_HH

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

/**
 * @brief provides static methods for grid extraction
 *
 * \tparam GV the grid view type
 * \tparam cd codimension of the extracted entities
 */
template<typename GV, int cd>
class Extractor
{
public:

  enum {dimworld = GV::dimensionworld};
  enum {dim      = GV::dimension};
  enum {codim    = cd};

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dim-codim+1
  };

  enum
  {
    cube_corners = 1 << (dim-codim)
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

  // transformations
  typedef Dune::GenericGeometry::BasicGeometry<dim-codim, Dune::GenericGeometry::DefaultGeometryTraits<ctype,dim-codim,dimworld> > Geometry;
  typedef Dune::GenericGeometry::BasicGeometry<dim-codim, Dune::GenericGeometry::DefaultGeometryTraits<ctype,dim-codim,dim> > LocalGeometry;

protected:
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

    FaceInfo(IndexType parent_, unsigned int num_in_parent_)
      :   parent(parent_), num_in_parent(num_in_parent_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType parent;

    /// @brief the number of the face in the parent element
    unsigned int num_in_parent : 3;

    /// @brief the corner indices plus the numbers of the vertices in the parent element
    CornerInfo corners[simplex_corners];     // sim = numer of vertices in a simplex
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
  Extractor(const GV& gv)
    :  _gv(gv)
  {}

  /** \brief Destructor frees allocated memory */
  ~Extractor();

  /*  F U N C T I O N A L I T Y  */

  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear()
  {
    // this is an inofficial way on how to free the memory allocated
    // by a std::vector
    {
      std::vector<CoordinateInfo> dummy;
      this->_coords.swap(dummy);
    }
    {
      std::vector<FaceInfo> dummy;
      this->_faces.swap(dummy);
    }

    // first free all manually allocated vertex/element info items...
    for (typename VertexInfoMap::iterator it = this->_vtxInfo.begin();
         it != this->_vtxInfo.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    for (typename ElementInfoMap::iterator it = this->_elmtInfo.begin();
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
    typename ElementInfoMap::const_iterator it =
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

  /** \brief Get world geometry of the extracted face */
  Geometry geometry(unsigned int index) const;

  /** \brief Get geometry of the extracted face in element coordinates */
  LocalGeometry geometryLocal(unsigned int index) const;

};


template<typename GV, int cd>
Extractor<GV,cd>::~Extractor()
{
  clear();
}


/** \brief Get World geometry of the extracted face */
template<typename GV, int cd>
typename Extractor<GV,cd>::Geometry Extractor<GV,cd>::geometry(unsigned int index) const
{
  std::vector<Coords> corners(simplex_corners);
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;

  return Geometry(Dune::GeometryType(Dune::GeometryType::simplex,dim-codim), corners);
}


/** \brief Get Geometry of the extracted face in element coordinates */
template<typename GV, int cd>
typename Extractor<GV,cd>::LocalGeometry Extractor<GV,cd>::geometryLocal(unsigned int index) const
{
  std::vector<Coords> corners(simplex_corners);

  // get face info
  const FaceInfo & face = this->_faces[index];
  Dune::GeometryType facetype(Dune::GeometryType::simplex, dim-codim);
  // get reference element
  Dune::GeometryType celltype(Dune::GeometryType::cube, dim);
  const Dune::GenericReferenceElement<ctype, dim> & re =
    Dune::GenericReferenceElements<ctype, dim>::general(celltype);

  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = re.position(face.corners[i].num,dim);

  return LocalGeometry(facetype, corners);
}

#endif // DUNE_EXTRACTOR_HH
