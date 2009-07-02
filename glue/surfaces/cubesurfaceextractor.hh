// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    CubeSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Apr 25, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for cube surface grids
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief grid extractor implementation for hypercube surface grids
 */

#ifndef CUBESURFACEEXTRACTOR_HH_
#define CUBESURFACEEXTRACTOR_HH_


#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <dune/common/static_assert.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>
#include "surfacedescriptor.hh"


/**
 * @class CubeSurfaceExtractor
 * @brief grid extractor implementation for cube surface grids
 *
 * Provides methods that build topology information for given grids.
 * The template parameters
 * @li GV the grid class type
 */
template<typename GV, bool rect = false>
class CubeSurfaceExtractor
  : public Codim1Extractor<GV>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    dimw = GV::dimensionworld
  };

  enum
  {
    dim = GV::dimension
  };

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
  typedef Dune::FieldVector<ctype, dimw>                           Coords;
  typedef Dune::array<unsigned int, simplex_corners>               SimplexTopology;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<dim>::Iterator VertexIter;

  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;
  typedef typename GV::Traits::template Codim<0>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // index sets and index types
  typedef typename GV::IndexSet IndexSet;
  typedef typename IndexSet::IndexType IndexType;


private:



  /************************** PRIVATE SUBCLASSES **********************/

  /**
   * @class FaceInfo
   * @brief simple struct holding some packed information about a codimension 1 entity
   * to its parent codim 0 entity
   */
  struct FaceInfo
  {
    FaceInfo()
    {}

    FaceInfo(unsigned int index_, IndexType parent_, unsigned int num_in_parent_, unsigned int first_)
      :       parent(parent_), index(index_), num_in_parent(num_in_parent_), first(first_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType parent;

    /// @brief the index of this face (in internal storage scheme) // NEEDED??
    unsigned int index : 27;

    /// @brief the number of the face in the parent element
    unsigned int num_in_parent : 4;

    /// @brief flag telling whether this is the first of two triangles
    /// refining the associated face
    unsigned int first : 1;

    /// @brief the corner indices plus the numbers of the vertices in the parent element
    typename Codim1Extractor<GV>::CornerInfo corners[simplex_corners];
  };


  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /// @brief the grid object to extract the surface from
  const GV&                                           _gv;


  /*        Geometrical and Topological Information                */

  /// @brief all information about the extracted faces
  std::vector<FaceInfo>         _faces;

  /// @brief all information about the corner vertices of the extracted
  std::vector<typename Codim1Extractor<GV>::CoordinateInfo>   _coords;

  /// @brief a map enabling faster access to vertices and coordinates
  ///
  /// Maps a vertex' index (from index set) to an object holding the locally
  /// associated index of the vertex' coordinate in _coords and an entity
  /// pointer to the codim<dim> entity.
  typename Codim1Extractor<GV>::VertexInfoMap _vtxInfo;

  /// @brief a map enabling faster access to elements and faces
  ///
  /// Maps an element's index (from index set) to an object holding the locally
  /// associated index of its first face in _indices (if there are more they are
  /// positioned consecutively) and an entity pointer to the codim<0> entity.
  typename Codim1Extractor<GV>::ElementInfoMap _elmtInfo;


  /// @brief geometry type of the surface patches
  const Dune::GeometryType _codim0element;


public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  CubeSurfaceExtractor(const GV& gv) :
    _gv(gv), _codim0element(Dune::GeometryType::cube, dim)
  {
    std::cout << "This is CubeSurfaceExtractor on a <"
              << GV::dimension << "," << GV::dimensionworld
              << "> grid working in " << dimw
              << " space expecting faces of type " << _codim0element << "!" << std::endl;
  }


  /**
   * @brief default destructor, frees memory
   */
  ~CubeSurfaceExtractor();


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
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are cubes.
   * @param descr a descriptor that "selects" the faces to add to the surface
   */
  void update(const FaceDescriptor<GV>& descr);


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
        faces[i][j] = this->_faces[i].corners[j].idx;
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
    typename Codim1Extractor<GV>::VertexInfoMap::const_iterator it = this->_vtxInfo.find(this->index<dim>(v));
    return it == this->_vtxInfo.end() ? -1 : it->second->idx;
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
    typename Codim1Extractor<GV>::ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->index<0>(e));
    if (it == this->_elmtInfo.end())
    {
      first = -1;
      count = 0;
      return false;
    }
    // the iterator is valid, fill the out params
    first = it->second->idx;
    count = it->second->num;
    return true;
  }


  /**
   * @brief gets the number face in the parent element
   * @param index the index of the face
   * @return if failed -1, else the index
   */
  int numberInSelf(unsigned int index) const
  {
    return index < this->_faces.size() ? this->_faces[index].num_in_parent : -1;
  }


  //	/**
  //	 * @brief transforms local triangle coordinates to quadrilateral coordinates
  //	 *
  //	 * Necessary because quadrilaterals on the surface are internally subdivided
  //	 * to triangles which are ordinarily accessed with the usual local coordinates.
  //	 * But if the geometry of the face or the element is about to be used the
  //	 * triangle coordinates first have to be transformed to quadrilateral coordinates.
  //	 * @param index the index of the face
  //	 * @param coords barycentric coords
  //	 * @return @c local face coords
  //	 */
  //	Dune::FieldVector<dimw-1> toFaceCoords(unsigned int index, const Coords &coords) const
  //	{
  //		// nothing to do for the first triangle
  //		if (this->_faces[index].first == 1)
  //			return coords;
  //
  //		// it is a "2nd" face => corners (3 2 1) i.e. ((1 1) (0 1) (1 0))
  //		Coords result;
  //		ctype first_coord = 1.0 - result[0] - result[1];
  //		// point 1 = (1 0)  ->  (0 1)
  //		// point 2 = (0 1)  ->  (1 0)
  //		// point 3 = (0 0)  ->  (1 1)
  //		result[0] = 1.0-coords[0];
  //		result[1] = 1.0-coords[1];
  //		return result;
  //	}


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
  void localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const;


  /**
   * @brief for given barycentric coords in a simplex compute element and world coordinates
   *
   * @param index the index of the simplex
   * @param bcoords the barycentric coordinates
   * @param ecoords to be filled with element coordinates
   * @param wcoords to be filled with world coordinates
   */
  void localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const;


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
  template<typename CoordContainer>
  void localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, int size) const;


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
  void localAndGlobalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, CoordContainer &wcoords, int size) const;


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

}; // end of class CubeSurfaceExtractor



template<typename GV, bool rect>
CubeSurfaceExtractor<GV, rect>::~CubeSurfaceExtractor()
{
  // only the objects that have been allocated manually have to be
  // deallocated manually again
  for (typename Codim1Extractor<GV>::VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename Codim1Extractor<GV>::ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}


template<typename GV, bool rect>
void CubeSurfaceExtractor<GV, rect>::clear()
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
  for (typename Codim1Extractor<GV>::VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename Codim1Extractor<GV>::ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  // ...then clear the maps themselves, too
  this->_vtxInfo.clear();
  this->_elmtInfo.clear();
}


template<typename GV, bool rect>
void CubeSurfaceExtractor<GV, rect>::update(const FaceDescriptor<GV>& descr)
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
    int vertex_index = 0;
    IndexType eindex = 0;             // supress warning

    // needed later for insertion into a std::set which only
    // works with const references

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

      // remember the indices of the faces that shall become
      // part of the surface
      std::set<int> boundary_faces;

      // iterate over all intersections of codim 1 and test if the
      // boundary intersections are to be added to the surface
      for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
      {
        // only look at boundary faces
        if (is->boundary() && descr.contains(elit, is->indexInInside()))
          boundary_faces.insert(is->indexInInside());
      }

      // if some face is part of the surface add it!
      if (boundary_faces.size() != 0)
      {
        dune_static_assert(dim == 2 || dim == 3, "CubeSurfaceExtractor works only for 2D and 3D");
        static const int factor = (dim == 2 ? 1 : 2);
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new typename Codim1Extractor<GV>::ElementInfo(simplex_index, elit, factor * boundary_faces.size());

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (typename std::set<int>::const_iterator sit = boundary_faces.begin(); sit != boundary_faces.end(); ++sit)
        {
          // now we only have to care about the 3D case, i.e. the quadrilateral
          // face has to be divided into two triangles

          VertexPtr* vptrs[4];
          unsigned int vertex_indices[4];
          unsigned int vertex_numbers[4];

          // get the vertex pointers for the quadrilateral's corner vertices
          // and try for each of them whether it is already inserted or not
          for (int i = 0; i < cube_corners; ++i)
          {
            // get the number of the vertex in the parent element
            vertex_numbers[i] = orientedSubface<dim>(this->_codim0element, *sit, i);

            // get the vertex pointer and the index from the index set
            vptrs[i] = new typename Codim1Extractor<GV>::VertexPtr(elit->template subEntity<dim>(vertex_numbers[i]));
            IndexType vindex = this->index<dim>(*(*vptrs[i]));

            // if the vertex is not yet inserted in the vertex info map
            // it is a new one -> it will be inserted now!
            typename Codim1Extractor<GV>::VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
            if (vimit == this->_vtxInfo.end())
            {
              // insert into the map
              this->_vtxInfo[vindex] = new typename Codim1Extractor<GV>::VertexInfo(vertex_index, *vptrs[i]);
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

          // Currently the extractor splits quadrilaterals into triangles (because psurface)
          // can only handle triangles). Since this is done manually we cannot handle grids
          // of dimension higher than 3.
          dune_static_assert(dim<=3, "CubeSurfaceExtractor only implemented for 1d, 2d, 3d");

          // add a new face to the temporary collection
          temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, 1));
          for (int i=0; i<dim; i++) {
            temp_faces.back().corners[i].idx = vertex_indices[i];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[i].num = vertex_numbers[i];
          }

          // now introduce the second triangle subdividing the quadrilateral
          // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
          // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1
          if (dim==3) {
            // add a new face to the temporary collection for the second tri
            temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, 0));
            temp_faces.back().corners[0].idx = vertex_indices[3];
            temp_faces.back().corners[1].idx = vertex_indices[2];
            temp_faces.back().corners[2].idx = vertex_indices[1];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[3];
            temp_faces.back().corners[1].num = vertex_numbers[2];
            temp_faces.back().corners[2].num = vertex_numbers[1];
          }

        }                         // end loop over found surface parts
      }
    }             // end loop over elements

    // allocate the array for the face specific information...
    this->_faces.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->_faces.begin());
  }


  // now first write the array with the coordinates...
  this->_coords.resize(this->_vtxInfo.size());
  typename Codim1Extractor<GV>::VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    typename Codim1Extractor<GV>::CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }


  //	const char prefix[] = "CubeSurfaceExtractor: ";
  //
  //	STDOUTLN(prefix << "Extracted Coordinates (size=" << this->_coords.size() << ")");
  //	for (unsigned int i = 0; i < this->_coords.size(); ++i)
  //	{
  ////		if (i % 100 == 0)
  //		{
  //			STDOUT(prefix << "vtxindex=" << this->_coords[i].vtxindex << " index=" << this->_coords[i].index
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
  //			STDOUT(prefix << "parent=" << this->_faces[i].parent << " index=" << this->_faces[i].index
  //					<< " nip=" << this->_faces[i].num_in_parent << " corners={");
  //			for (unsigned int j = 0; j < simplex_corners; ++j)
  //				STDOUT("(" << this->_faces[i].corners[j].idx << ", " << this->_faces[i].corners[j].num << ")");
  //			STDOUTLN("}");
  //		}
  //
  //	}
}


template<typename GV, bool rect>
inline void CubeSurfaceExtractor<GV, rect>::globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
{
  // only interpolate barycentric in the given triangle => for flat quads this is exact!
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, wcoords, dimw);
}


template<typename GV, bool rect>
inline void CubeSurfaceExtractor<GV, rect>::localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const
{
  if (rect)
  {
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->indexInInside(index);
    // computing the locals is straight forward for flat rectangles,
    // we only need the triangle's corners in element coordinate
    if (this->_faces[index].first)
    {
      // the triangle's corners are (0 1 2) using face indices
      for (int i = 0; i < simplex_corners; ++i)
        corners[i] = cornerLocalInRefElement<ctype, dimw>(this->_codim0element, num_in_self, i);
    }
    else
    {
      // the triangle's corners are (3 2 1) using face indices
      for (int i = 0; i < simplex_corners; ++i)
        corners[i] = cornerLocalInRefElement<ctype, dimw>(this->_codim0element, num_in_self, 3-i);
    }
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, ecoords, dimw);
  }
  else
  {
    Coords wcoords;
    this->localAndGlobalCoords(index, bcoords, ecoords, wcoords);
  }
}


template<typename GV, bool rect>
inline void CubeSurfaceExtractor<GV, rect>::localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const
{
  this->globalCoords(index, bcoords, wcoords);
  // for rectangles avoid using world coordinates
  if (rect)
    this->localCoords(index, bcoords, ecoords);
  else
    ecoords = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().local(wcoords);
}


template<typename GV, bool rect>
template<typename CoordContainer>
void CubeSurfaceExtractor<GV, rect>::globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  for (int i = 0; i < size; ++i)
    interpolateBarycentric<simplex_corners, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], wcoords[i], dimw);
}


template<typename GV, bool rect>
template<typename CoordContainer>
void CubeSurfaceExtractor<GV, rect>::localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, int size) const
{
  if (rect)
  {
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->numberInSelf(index);
    // computing the locals is straight forward for flat rectangles,
    // we only need the triangle's corners in element coordinate
    if (this->_faces[index].first)
    {
      // the triangle's corners are (0 1 2) using face indices
      for (int i = 0; i < simplex_corners; ++i)
        corners[i] = cornerLocalInRefElement<ctype, dimw>(this->_codim0element, num_in_self, i);
    }
    else
    {
      // the triangle's corners are (3 2 1) using face indices
      for (int i = 0; i < simplex_corners; ++i)
        corners[i] = cornerLocalInRefElement<ctype, dimw>(this->_codim0element, num_in_self, 3-i);
    }
    for (int i = 0; i < size; ++i)
      interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], ecoords[i], dimw);
  }
  else
  {
    CoordContainer wcoords;
    this->localAndGlobalCoords(index, bcoords, ecoords, wcoords, size);
  }
}


template<typename GV, bool rect>
template<typename CoordContainer>
void CubeSurfaceExtractor<GV, rect>::localAndGlobalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, CoordContainer &wcoords, int size) const
{
  this->globalCoords(index, bcoords, wcoords, size);
  // for rectangles avoid using world coordinates
  if (rect)
    this->localCoords(index, bcoords, ecoords, size);
  else
  {
    ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].parent)->second->p;
    for (int i = 0; i < size; ++i)
      ecoords[i] = eptr->geometry().local(wcoords[i]);
  }
}

#endif // CUBESURFACEEXTRACTOR_HH_
