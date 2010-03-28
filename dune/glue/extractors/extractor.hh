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
#include <map>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/genericgeometry/geometry.hh>

/**
 * @brief Provides codimension-independent methods for grid extraction
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

  enum
  {
    cube_corners = 1 << (dim-codim)
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;
  typedef Dune::FieldVector<ctype, dim>                            LocalCoords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  typedef std::vector<unsigned int>                                VertexVector;

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
      : vtxindex(vtxindex_), index(index_)
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
   * @brief Holds some information about an element's subentity involved in a coupling
   */
  struct SubEntityInfo
  {
    SubEntityInfo()
    {
      geometryType_.makeSimplex(dim-codim);
    }

    SubEntityInfo(IndexType parent_, unsigned int num_in_parent_,
                  const Dune::GeometryType& geometryType)
      : parent(parent_), num_in_parent(num_in_parent_), geometryType_(geometryType)
    {}

    unsigned int nCorners() const
    {
      return Dune::GenericReferenceElements<ctype, dim-codim>::general(geometryType_).size(dim-codim);
    }

    /// @brief the index of the parent element (from index set)
    IndexType parent;

    /// @brief the number of the face in the parent element
    unsigned int num_in_parent : 3;

    /** \brief The GeometryType of the subentity */
    Dune::GeometryType geometryType_;

    /** @brief the corner indices plus the numbers of the vertices in the parent element

        This array has the length cube_corners, because currently that is an upper boun d
        for the number of corners of an element.  If more general element types appear we
        need to change this.
     */
    CornerInfo corners[cube_corners];     // sim = numer of vertices in a simplex
  };


  typedef std::map<IndexType, ElementInfo* >  ElementInfoMap;
  typedef std::map<IndexType, VertexInfo* >   VertexInfoMap;


  /************************** MEMBER VARIABLES ************************/

  /// @brief the grid object to extract the surface from
  const GV&                       gv_;

  /*        Geometrical and Topological Information                */

  /// @brief all information about the corner vertices of the extracted
  std::vector<CoordinateInfo>   coords_;

  /// @brief all information about the extracted subEntities
  std::vector<SubEntityInfo> subEntities_;

  /// @brief a map enabling faster access to vertices and coordinates
  ///
  /// Maps a vertex' index (from index set) to an object holding the locally
  /// associated index of the vertex' coordinate in coords_ and an entity
  /// pointer to the codim<dim> entity.
  VertexInfoMap vtxInfo_;

  /// @brief a map enabling faster access to elements and faces
  ///
  /// Maps an element's index (from index set) to an object holding the locally
  /// associated index of its first face in _indices (if there are more they are
  /// positioned consecutively) and an entity pointer to the codim<0> entity.
  ElementInfoMap elmtInfo_;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   */
  Extractor(const GV& gv)
    :  gv_(gv)
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
      this->coords_.swap(dummy);
    }
    {
      std::vector<SubEntityInfo> dummy;
      this->subEntities_.swap(dummy);
    }

    // first free all manually allocated vertex/element info items...
    for (typename VertexInfoMap::iterator it = this->vtxInfo_.begin();
         it != this->vtxInfo_.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    for (typename ElementInfoMap::iterator it = this->elmtInfo_.begin();
         it != this->elmtInfo_.end(); ++it)
      if (it->second != NULL)
        delete it->second;
    // ...then clear the maps themselves, too
    this->vtxInfo_.clear();
    this->elmtInfo_.clear();
  }


  /*  G E T T E R S  */

  /**
   * @brief getter for the coordinates array
   * @param coords a vector that will be resized (!) and filled with the coordinates,
   * note that the single components are written consecutively
   */
  void getCoords(std::vector<Dune::FieldVector<ctype, dimworld> >& coords) const
  {
    coords.resize(coords_.size());
    for (unsigned int i = 0; i < coords_.size(); ++i)
      coords[i] = coords_[i].coord;
  }


  /**
   * @brief getter for the count of coordinates
   * @return the count
   */
  unsigned int nCoords() const
  {
    return coords_.size();
  }

  /** \brief Get the list of geometry types */
  void getGeometryTypes(std::vector<Dune::GeometryType>& geometryTypes) const
  {
    geometryTypes.resize(subEntities_.size());
    for (size_t i=0; i<subEntities_.size(); i++)
      geometryTypes[i] = subEntities_[i].geometryType_;
  }


  /**
   * @brief Get the corners of the extracted subentities
   */
  void getFaces(std::vector<VertexVector>& faces) const
  {
    faces.resize(this->subEntities_.size());
    for (unsigned int i = 0; i < this->subEntities_.size(); ++i) {
      faces[i].resize(subEntities_[i].nCorners());
      for (unsigned int j = 0; j < subEntities_[i].nCorners(); ++j)
        faces[i][j] = this->subEntities_[i].corners[j].idx;
    }
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
      this->elmtInfo_.find(this->indexSet().template index<0>(e));
    if (it == this->elmtInfo_.end())
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
    return index < this->subEntities_.size() ? this->subEntities_[index].num_in_parent : -1;
  }

  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return gv_.indexSet();
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
    if (index >= this->subEntities_.size())
      DUNE_THROW(Dune::GridError, "invalid face index");
    return (this->elmtInfo_.find(this->subEntities_[index].parent))->second->p;
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
    if (index >= this->coords_.size())
      DUNE_THROW(Dune::GridError, "invalid coordinate index");
    return (this->vtxInfo_.find(this->coords_[index].vtxindex))->second->p;
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
  std::vector<Coords> corners(subEntities_[index].nCorners());
  for (unsigned int i = 0; i < subEntities_[index].nCorners(); ++i)
    corners[i] = this->coords_[this->subEntities_[index].corners[i].idx].coord;

  return Geometry(subEntities_[index].geometryType_, corners);
}


/** \brief Get Geometry of the extracted face in element coordinates */
template<typename GV, int cd>
typename Extractor<GV,cd>::LocalGeometry Extractor<GV,cd>::geometryLocal(unsigned int index) const
{
  std::vector<LocalCoords> corners(subEntities_[index].nCorners());

  // get face info
  const SubEntityInfo & face = this->subEntities_[index];
  Dune::GeometryType facetype = subEntities_[index].geometryType_;

  // get reference element
  Dune::GeometryType celltype = elmtInfo_.find(face.parent)->second->p->type();
  const Dune::GenericReferenceElement<ctype, dim> & re =
    Dune::GenericReferenceElements<ctype, dim>::general(celltype);

  for (unsigned int i = 0; i < subEntities_[index].nCorners(); ++i)
    corners[i] = re.position(face.corners[i].num,dim);

  return LocalGeometry(facetype, corners);
}

#endif // DUNE_EXTRACTOR_HH
