// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    CubeManifoldExtractor.hh
 *  Version:     1.0
 *  Created on:  Apr 27, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for cube manifolds
 *  subversion:  $Id$
 *
 */
/**
 * @file CubeManifoldExtractor.hh
 * @brief
 */

#ifndef CUBEMANIFOLDEXTRACTOR_HH_
#define CUBEMANIFOLDEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/common/geometry.hh>

#include "surfacedescriptor.hh"
#include "simplicialmanifoldextractor.hh"


/**
 * @class CubeManifoldExtractor
 * @brief grid extractor implementation for cube manifolds
 *
 *
 **!!! WARNING !!!
 **!!! As of now manifold grids are not supported in Dune.
 **!!! This code is supposed to do the job when support is finally
 **!!! implemented but there is no guarantee.
 **!!! Not a single test hast been done!!
 * Provides methods that build topology information for given grids.
 * Note that these methods only operate on the grid.
 * The template parameters
 * @li GV the grid class type
 */
template<typename GV, bool rect = false, int dimG = GV::dimension>
class CubeManifoldExtractor
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */


  enum
  {
    dim = GV::dimension
  };

  enum
  {
    dimw = GV::dimensionworld
  };

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dimw
  };

  enum
  {
    cube_corners = 1 << dim
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimw>                           Coords;
  typedef Dune::array<unsigned int, simplex_corners>               SimplexTopology;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<dim>::Iterator VertexIter;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;

private:

  /************************** PRIVATE SUBCLASSES **********************/

  /**
   * @class FaceInfo
   * @brief simple struct holding some packed information about this codim 0 entity
   */
  struct FaceInfo
  {
    FaceInfo()
    {}

    FaceInfo(unsigned int index_, IndexType elmtindex_, unsigned int first_)
      :       self(elmtindex_), index(index_), first(first_)
    {}

    /// @brief the index of this element (from index set)
    IndexType self;

    /// @brief the index of this face (in internal storage scheme) // NEEDED??
    unsigned int index : 31;

    /// @brief flag telling whether this is the first of two triangles
    /// refining the associated face
    unsigned int first : 1;

    /// @brief the corner indices plus the numbers of the vertices in this element
    unsigned int corners[simplex_corners];
  };


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


  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /// @brief the grid object to extract the surface from
  const GV&                                              _gv;


  /*        Geometrical and Topological Information                */

  /// @brief all information about the extracted faces
  std::vector<FaceInfo>         _faces;

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


  /// @brief geometry type of the surface patches
  const Dune::GeometryType _codim0element;


public:


  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  CubeManifoldExtractor(const GV& gv) :
    _gv(gv), _codim0element(Dune::GeometryType::cube, dim)
  {
    STDOUTLN("This is CubeManifoldExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid working in " << dimw << " space expecting faces of type " << _codim0element << "!");
  }


  /**
   * @brief default destructor, frees memory
   */
  ~CubeManifoldExtractor();


  /*  F U N C T I O N A L I T Y  */

  /**
   * ASSUMPTION:
   * dim == dimw
   *
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a mapper or set object
   * @c m specifying an index set with codimension 0 entities near or on the boundary.
   * It is assumed that only one geometric shape exists on the boundary.
   * The template parameter n then denotes the number of corners per boundary element
   * (e.g. n==3 for triangles in a 3D grid of tetrahedra).
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If we the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are simplices.
   */
  void update(const ElementDescriptor<GV>& descr);


  /**
   * @brief delete everything build up so far and free the memory
   */
  void clear();


  /*  S E T T E R S  */


  /*  G E T T E R S  */

  /**
   * @brief getter for the coordinates array
   * @param coords a vector that will be resized (!) and filled with the coordinates,
   * note that the single components are written consecutively
   */
  void getCoords(std::vector<Dune::FieldVector<ctype, dimw> >& coords) const
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
        faces[i][j] = this->_faces[i].corners[j];
  }


  /**
   * @brief getter for internally used index set (grid's index set)
   * @return the index set
   */
  const IndexSet& indexSet() const
  {
    return this->_gv.indexSet();
  }


  /**
   * @brief getter for the index of an entity of codim cc
   * @return the index specified by the grid's index set
   */
  template<int cc>
  IndexType index(const typename GV::Traits::template Codim<cc>::Entity& e) const
  {
    return this->indexSet().template index<cc>(e);
  }


  /**
   * @brief gets index of coordinate in _coords associated with given vertex
   * @return the index if possible, -1 else
   */
  int coordinateIndex(const Vertex& v) const
  {
    typename VertexInfoMap::const_iterator it = this->_vtxInfo.find(this->index<dim>(v));
    return it == this->_vtxInfo.end() ? -1 : it->second->idx;
  }


  /**
   * @brief gets index of first face as well as the total number of faces that
   * were extracted from this element
   * @param e the element
   * @param first will contain the first index if found, else -1
   * @param count will be 1 if found, else 0.
   * Note: parameter is kept only to have the same interface as the surface extractor
   * @return success
   */
  bool faceIndices(const Element& e, int& first, int& count) const
  {
    typename ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->index<0>(e));
    if (it == this->_elmtInfo.end())
    {
      first = -1;
      count = 0;
      return false;
    }
    // the iterator is valid, fill the out params
    first = it->second->idx;
    count = 2;
    return true;
  }


  /**
   * @brief gets the number of the face in the element (which is the face)
   *
   * Note: this method is only kept to have the same interface as the surface extractor
   * @param index the index of the element face
   * @return if failed -1, else 0 (consistent with Dune speaking of entity<0>(0) )
   */
  int numberInSelf(unsigned int index) const
  {
    return index < this->_faces.size() ? 0 : -1;
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
  void globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const;


  /**
   * @brief for given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  void localCoords(unsigned int index, const Coords &bcoords, Dune::FieldVector<ctype, dim> &ecoords) const;


  /**
   * @brief for given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const Coords &bcoords, Dune::FieldVector<ctype, dim> &ecoords, Coords &wcoords) const;


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
  void globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const;


  /**
   * @brief for several given barycentric coords in a simplex compute element coordinates
   *
   * If both are to be computed, element and world coordinates, then use the
   * combined method for efficiency!
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   */
  template<typename CoordContainerB, typename CoordContainerE>
  void localCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, int size) const;


  /**
   * @brief for several given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   * @return
   */
  template<typename CoordContainerB, typename CoordContainerE>
  void localAndGlobalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, CoordContainerB &wcoords, int size) const;


  /**
   * @brief gets for each vertex corner of given face (by index) the number of
   * the vertex in parent element's ordering
   * @param index the face's index
   * @param corner the index of the corner
   * @return -1 <=> index invalid or array not filled, else index
   */
  int numCornerInParent(unsigned int index, unsigned int corner) const
  {
    return (index >= this->_faces.size() || corner >= simplex_corners) ? -1 :
           ((this->_faces[index].first == 0) ? 3 - corner : corner);
  }


  /**
   * @brief gets the the element for a given face/element index,
   * throws an exception if index not valid
   * @param index the index of the face
   * @return a reference to the element's stored pointer
   */
  const ElementPtr& element(unsigned int index) const
  {
    if (index >= this->_faces.size())
      DUNE_THROW(Dune::GridError, "invalid face index");
    return (this->_elmtInfo.find(this->_faces[index].self))->second->p;
  }


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
    return (this->_vtxInfo.find(this->_coords[index].self))->second->p;
  }


  /**
   * @brief gets the indices of all faces with the given coordinate as corner
   * @param index the index of the coordinate
   * @param faces array if given index was legal
   * @param count length of array if successful
   * @return TRUE <=> if successful
   * DO NOT MODIFY THE ARRAY'S CONTENT!
   */
  bool parentFaces(unsigned int index,  unsigned int const*& faces, unsigned int& count) const
  {
    if (index >= this->_coords.size())
      return false;
    // index valid
    faces = this->_coords[index].faces;
    count = this->_coords[index].num_faces;
    return true;
  }

}; // end of class CubeManifoldExtractor



template<typename GV, bool rect, int dimG>
CubeManifoldExtractor<GV, rect, dimG>::~CubeManifoldExtractor()
{
  // only the objects that have been allocated manually have to be
  // deallocated manually again
  // free all the manually allocated memory
  for (unsigned int i = 0; i < this->_coords.size(); ++i)
    if (this->_coords[i].faces != NULL)
      delete this->_coords[i].faces;
  for (typename VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}



template<typename GV, bool rect, int dimG>
void CubeManifoldExtractor<GV, rect, dimG>::clear()
{
  // this is an inofficial way on how to free the memory allocated
  // by a std::vector
  {
    // free all the manually allocated memory
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
      if (this->_coords[i].faces != NULL)
        delete this->_coords[i].faces;
    std::vector<CoordinateInfo> dummy;
    this->_coords.swap(dummy);
  }
  {
    std::vector<FaceInfo> dummy;
    this->_faces.swap(dummy);
  }

  // first free all manually allocated vertex/element info items...
  for (typename VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  // ...then clear the maps themselves, too
  this->_vtxInfo.clear();
  this->_elmtInfo.clear();
}



template<typename GV, bool rect, int dimG>
void CubeManifoldExtractor<GV, rect, dimG>::update(const ElementDescriptor<GV>& descr)
{
  // free everything there is in this object
  this->clear();

  // In this first pass iterate over all entities of codim 0.
  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.

  {
    // several counter for consecutive indexing are needed
    int simplex_index = 0;
    int vertex_indices[cube_corners];
    int vertex_index = 0;
    IndexType eindex = 0;             // supress warning

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<FaceInfo> temp_faces;

    // iterate over all codim 0 elemets on the grid
    for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
    {
      // check if there are unwanted geometric shapes
      // if one appears => exit with error
      if (elit->type() != this->_codim0element)
        DUNE_THROW(Dune::GridError, "expected cube grid but found non-cube entity of codimension 0: " << elit->type());

      // only do sth. if this element is "interesting"
      // implicit cast is done automatically
      if (descr.contains(elit))
      {
        // add an entry to the element info map, the index will be set properly later
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit);

        // try for each of the faces vertices whether it is already inserted or not
        for (int i = 0; i < cube_corners; ++i)
        {
          // TODO Think about a way to ensure that orientation of all faces is consistent!
          // get the vertex pointer and the index from the index set
          // Note that the orientation is always the same for all simplices,
          // i.e. CCW which is 0,1 in 2D and 0,1,2 in 3D
          VertexPtr vptr(elit->template entity<dim>(i));
          IndexType vindex = this->index<dim>(*vptr);

          // if the vertex is not yet inserted in the vertex info map
          // it is a new one -> it will be inserted now!
          typename VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
          if (vimit == this->_vtxInfo.end())
          {
            // insert into the map
            this->_vtxInfo[vindex] = new VertexInfo(vertex_index, vptr);
            // remember the vertex as a corner of the current face in temp_faces
            vertex_indices[i] = vertex_index;
            // increase the current index
            vertex_index++;
          }
          else
          {
            // only insert the index into the simplices array
            vertex_indices[i] = vimit->second->idx;
          }
        }


        // add a new face to the temporary collection for the 1st triangle
        temp_faces.push_back(FaceInfo(simplex_index++, eindex, 1));
        temp_faces.back().corners[0] = vertex_indices[0];
        temp_faces.back().corners[1] = vertex_indices[1];
        temp_faces.back().corners[2] = vertex_indices[2];

        // add a new face to the temporary collection for the 2nd triangle
        temp_faces.push_back(FaceInfo(simplex_index++, eindex, 0));
        temp_faces.back().corners[0] = vertex_indices[3];
        temp_faces.back().corners[1] = vertex_indices[2];
        temp_faces.back().corners[2] = vertex_indices[1];
      }
    }             // end loop over elements

    // allocate the array for the face specific information...
    this->_faces.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());
  }


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
    current->self = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }


  // now add the vertices' parent faces in the _vertexFaces map.
  // therefore iterate over all indices in the _index array...
  {
    std::vector<unsigned int> refcount(this->_coords.size(), 0);

    // for each coordinate count the references in the _indices array
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
      for (unsigned int j = 0; j < simplex_corners; ++j)
        refcount[this->_faces[i].corners[j]]++;

    // allocate the right amount of storage for each vertex's references
    for (unsigned int i = 0; i < this->_coords.size(); ++i)
    {
      // allocate an array and initialize its first element with its length
      this->_coords[i].num_faces = refcount[i];
      this->_coords[i].faces = new unsigned int[refcount[i]];
      refcount[i] = 0;                   // used as "pointer" in next loop
    }

    // add the references
    for (unsigned int i = 0; i < this->_faces.size(); ++i)
    {
      for (unsigned int j = 0; j < simplex_corners; ++j)
      {
        unsigned int ref = this->_faces[i].corners[j];
        this->_coords[ref].faces[refcount[ref]] = i;
        refcount[ref]++;
      }
    }
  }


  //	const char prefix[] = "CubeManifoldExtractor: ";
  //
  //	STDOUTLN(prefix << "Extracted Coordinates (size=" << this->_coords.size() << ")");
  //	for (unsigned int i = 0; i < this->_coords.size(); ++i)
  //	{
  ////		if (i % 100 == 0)
  //		{
  //			STDOUT(prefix << "self=" << this->_coords[i].self << " index=" << this->_coords[i].index
  //					<< " coord=(" << this->_coords[i].coord << ") num_faces=" << this->_coords[i].num_faces
  //					<< " faces={");
  //			for (unsigned int j = 0; j < this->_coords[i].num_faces; ++j)
  //				STDOUT(" " << this->_coords[i].faces[j]);
  //			STDOUTLN("}");
  //		}
  //	}
  //
  //	STDOUTLN("\n" << prefix << "Extracted faces (size=" << this->_faces.size() << ")");
  //	for (unsigned int i = 0; i < this->_faces.size(); ++i)
  //	{
  ////		if (i % 100 == 0)
  //		{
  //			STDOUT(prefix << "self=" << this->_faces[i].self << " index=" << this->_faces[i].index << " corners={");
  //			for (unsigned int j = 0; j < simplex_corners; ++j)
  //				STDOUT("(" << this->_faces[i].corners[j] << ")");
  //			STDOUTLN("}");
  //		}
  //
  //	}

}


template<typename GV, bool rect, int dimG>
inline void CubeManifoldExtractor<GV, rect, dimG>::globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
{
  // only interpolate barycentric in the given triangle => for flat quads this is exact!
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i]].coord;
  interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dim> >(corners, bcoords, wcoords, dim);
}


template<typename GV, bool rect, int dimG>
inline void CubeManifoldExtractor<GV, rect, dimG>::localCoords(unsigned int index, const Coords &bcoords, Dune::FieldVector<ctype, dim> &ecoords) const
{
  Coords wcoords;
  this->localAndGlobalCoords(index, bcoords, ecoords, wcoords);
}


template<typename GV, bool rect, int dimG>
inline void CubeManifoldExtractor<GV, rect, dimG>::localAndGlobalCoords(unsigned int index, const Coords &bcoords, Dune::FieldVector<ctype, dim> &ecoords, Coords &wcoords) const
{
  this->globalCoords(index, bcoords, wcoords);
  ecoords = this->_elmtInfo.find(this->_faces[index].self)->second->p->geometry().local(wcoords);
}


template<typename GV, bool rect, int dimG>
template<typename CoordContainer>
void CubeManifoldExtractor<GV, rect, dimG>::globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i]].coord;
  for (int i = 0; i < size; ++i)
    interpolateBarycentric<simplex_corners, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], wcoords[i], dimw);
}


template<typename GV, bool rect, int dimG>
template<typename CoordContainerB, typename CoordContainerE>
void CubeManifoldExtractor<GV, rect, dimG>::localCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, int size) const
{
  CoordContainerE wcoords;
  this->localAndGlobalCoords(index, bcoords, ecoords, wcoords, size);
}


template<typename GV, bool rect, int dimG>
template<typename CoordContainerB, typename CoordContainerE>
void CubeManifoldExtractor<GV, rect, dimG>::localAndGlobalCoords(unsigned int index, const CoordContainerB &bcoords, CoordContainerE &ecoords, CoordContainerB &wcoords, int size) const
{
  this->globalCoords(index, bcoords, wcoords, size);
  ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].self)->second->p;
  for (int i = 0; i < size; ++i)
    ecoords[i] = eptr->geometry().local(wcoords[i]);
}


/*   S P E C I A L I Z A T I O N   F O R   2 D   G R I D S   */

template<typename GV, bool rect>
class CubeManifoldExtractor<GV, rect, 1> : public SimplicialManifoldExtractor<GV>
{
private:

  typedef SimplicialManifoldExtractor<GV>  Base;


public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    cube_corners = 1 << (Base::dim-1)
  };


  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  CubeManifoldExtractor(const GV& gv) :
    Base(gv, Dune::GeometryType(Dune::GeometryType::cube, Base::dim))
  {
    STDOUTLN("This is CubeManifoldExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid working in " << Base::dimw << " space expecting faces of type " << Dune::GeometryType(Dune::GeometryType::cube, Base::dim) << "!");
  }
};


#endif // CUBEMANIFOLDEXTRACTOR_HH_
