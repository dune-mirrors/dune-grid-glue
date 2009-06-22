// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GeneralSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Feb 3, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for hybrid surface grids (tris and/or quads on surface)
 *  subversion:  $Id$
 *
 */
/**
 * @file GeneralSurfaceExtractor.hh
 * @brief grid extractor implementation for hybrid surface grids (tris and/or quads on surface)
 */

#ifndef GENERALSURFACEEXTRACTOR_HH_
#define GENERALSURFACEEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/referenceelements.hh>
#include <dune/grid/common/geometry.hh>
#include "surfacedescriptor.hh"
#include "generaledgeextractor.hh"


/**
 * @class GeneralSurfaceExtractor
 * @brief grid extractor implementation for arbitrary surface grids (tris and/or quads on surface)
 *
 * Provides methods that build topology information for given grids.
 * The template parameters
 * @li GV the grid class type
 */
template<typename GV, int dimG = GV::dimension>
class GeneralSurfaceExtractor
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
   * @class CornerInfo
   * @brief
   * Helpful struct holding one index for the coordinate (vertex)
   * to which it is associated and the element's corner index;
   */
  struct CornerInfo
  {
    unsigned int idx : 28;           ///< index of the vertex
    unsigned int num : 4;           ///< element corner
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
    /// refining the associated face or a even a triangle face
    ///
    /// 0: triangle face
    /// 1: first of two tris refining a quadrilateral
    /// 2: second of two tris refining a quadrilateral
    unsigned int category : 2;

    /// @brief the corner indices plus the numbers of the vertices in the parent element
    CornerInfo corners[simplex_corners];
  };


  struct CoordinateInfo
  {
    CoordinateInfo()
    {}

    CoordinateInfo(unsigned int index_, IndexType vtxindex_)
      :       vtxindex(vtxindex_), index(index_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType vtxindex;

    /// @brief the coordinate
    Coords coord;

    /// @brief the index of this coordinate (in internal storage scheme) // NEEDED??
    unsigned int index;
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
    ElementInfo(unsigned int idx_, ElementPtr& p_, unsigned int num_) : idx(idx_), num(num_), p(p_)
    {}

    /// @brief the index of this element's first face in the internal list of extracted faces
    unsigned int idx : 28;

    /// @brief the number of extracted faces for this element
    unsigned int num : 4;

    /// @brief the entity pointer to the element
    ElementPtr p;
  };


  typedef std::map<IndexType, ElementInfo* >  ElementInfoMap;
  typedef std::map<IndexType, VertexInfo* >   VertexInfoMap;


  /************************** MEMBER VARIABLES ************************/

  // these values are filled on surface extraction and can be
  // asked by the corresponding getters

  /// @brief the grid object to extract the surface from
  const GV&                                           _gv;


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


public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief except from the GridView initializes all member variables with null values
   * @param gv the grid view object to work with
   */
  GeneralSurfaceExtractor(const GV& gv) :
    _gv(gv)
  {
    STDOUTLN("This is GeneralSurfaceExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid working in " << dimw << " space!");
  }


  /**
   * @brief default destructor, frees memory
   */
  ~GeneralSurfaceExtractor();


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
    typename VertexInfoMap::const_iterator it = this->_vtxInfo.find(this->index<dim>(v));
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
    typename ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->index<0>(e));
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

}; // end of class GeneralSurfaceExtractor



template<typename GV, int dimG>
GeneralSurfaceExtractor<GV, dimG>::~GeneralSurfaceExtractor()
{
  // only the objects that have been allocated manually have to be
  // deallocated manually again
  for (typename VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}


template<typename GV, int dimG>
void GeneralSurfaceExtractor<GV, dimG>::clear()
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




template<typename GV, int dimG>
void GeneralSurfaceExtractor<GV, dimG>::update(const FaceDescriptor<GV>& descr)
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
      Dune::GeometryType gt = elit->type();

      // remember the indices of the faces that shall become
      // part of the surface
      std::set<int> boundary_faces;

      // iterate over all intersections of codim 1 and test if the
      // boundary intersections are to be added to the surface
      for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
      {
        // only look at boundary faces
        if (is->boundary() && descr.contains(elit, is->numberInSelf()))
          boundary_faces.insert(is->numberInSelf());
      }

      // if some face is part of the surface add it!
      if (boundary_faces.size() != 0)
      {
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
        eindex = this->indexSet().template index<0>(*elit);
        this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit, 0);

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (typename std::set<int>::const_iterator sit = boundary_faces.begin(); sit != boundary_faces.end(); ++sit)
        {
          // get the corner count of this face
          const int face_corners = Dune::ReferenceElements<ctype, dim>::general(gt).size(*sit, 1, dim);

          // register the additional face(s)
          this->_elmtInfo[eindex]->num += (face_corners - 2);

          // now we only have to care about the 3D case, i.e. a triangle face can be
          // inserted directly whereas a quadrilateral face has to be divided into two triangles
          switch (face_corners)
          {
          case 3 :
            std::cout << "adding TRI..." << std::endl;
            // we have a triangle here

            // add a new face to the temporary collection
            temp_faces.push_back(FaceInfo(simplex_index, eindex, *sit, 0));

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < simplex_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = orientedSubface<dim>(gt, *sit, i);

              // get the vertex pointer and the index from the index set
              VertexPtr vptr(elit->template entity<dim>(vertex_number));
              IndexType vindex = this->index<dim>(*vptr);

              // remember the vertex' number in parent element's vertices
              temp_faces.back().corners[i].num = vertex_number;

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
              if (vimit == this->_vtxInfo.end())
              {
                // insert into the map
                this->_vtxInfo[vindex] = new VertexInfo(vertex_index, vptr);
                // remember the vertex as a corner of the current face in temp_faces
                temp_faces.back().corners[i].idx = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only insert the index into the simplices array
                temp_faces.back().corners[i].idx = vimit->second->idx;
              }
            }

            // now increase the current face index
            simplex_index++;
            break;
          case 4 :
            std::cout << "adding QUAD..." << std::endl;
            // we have a quadrilateral here
            VertexPtr* vptrs[4];
            unsigned int vertex_indices[4];
            unsigned int vertex_numbers[4];

            // get the vertex pointers for the quadrilateral's corner vertices
            // and try for each of them whether it is already inserted or not
            for (int i = 0; i < cube_corners; ++i)
            {
              // get the number of the vertex in the parent element
              vertex_numbers[i] = orientedSubface<dim>(gt, *sit, i);

              // get the vertex pointer and the index from the index set
              vptrs[i] = new VertexPtr(elit->template entity<dim>(vertex_numbers[i]));
              IndexType vindex = this->index<dim>(*(*vptrs[i]));

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->_vtxInfo.find(vindex);
              if (vimit == this->_vtxInfo.end())
              {
                // insert into the map
                this->_vtxInfo[vindex] = new VertexInfo(vertex_index, *vptrs[i]);
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

            // now introduce the two triangles subdividing the quadrilateral
            // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
            // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1

            // add a new face to the temporary collection for the first tri
            temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, 1));
            temp_faces.back().corners[0].idx = vertex_indices[0];
            temp_faces.back().corners[1].idx = vertex_indices[1];
            temp_faces.back().corners[2].idx = vertex_indices[2];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[0];
            temp_faces.back().corners[1].num = vertex_numbers[1];
            temp_faces.back().corners[2].num = vertex_numbers[2];

            // add a new face to the temporary collection for the second tri
            temp_faces.push_back(FaceInfo(simplex_index++, eindex, *sit, 2));
            temp_faces.back().corners[0].idx = vertex_indices[3];
            temp_faces.back().corners[1].idx = vertex_indices[2];
            temp_faces.back().corners[2].idx = vertex_indices[1];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[3];
            temp_faces.back().corners[1].num = vertex_numbers[2];
            temp_faces.back().corners[2].num = vertex_numbers[1];
            break;
          default :
            DUNE_THROW(Dune::NotImplemented, "the extractor does only work for triangle and quadrilateral faces");
            break;
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
  typename VertexInfoMap::const_iterator it1 = this->_vtxInfo.begin();
  for (; it1 != this->_vtxInfo.end(); ++it1)
  {
    // get a pointer to the associated info object
    CoordinateInfo* current = &this->_coords[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
    current->coord = it1->second->p->geometry().corner(0);
  }


  //	const char prefix[] = "GeneralSurfaceExtractor: ";
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


template<typename GV, int dimG>
inline void GeneralSurfaceExtractor<GV, dimG>::globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
{
  // only interpolate barycentric in the given triangle (=> for flat quads this is exact, not so for non-flat quads!)
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, wcoords, dimw);
}


template<typename GV, int dimG>
inline void GeneralSurfaceExtractor<GV, dimG>::localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const
{
  if (this->_faces[index].category == 0)
  {
    // this is a triangle face
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->numberInSelf(index);
    Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->type();
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimw>(gt, num_in_self, i);
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, ecoords, dimw);
  }
  else
  {
    // this is a quadrilateral face
    Coords wcoords;
    this->localAndGlobalCoords(index, bcoords, ecoords, wcoords);
  }
}


template<typename GV, int dimG>
inline void GeneralSurfaceExtractor<GV, dimG>::localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const
{
  this->globalCoords(index, bcoords, wcoords);
  // for rectangles avoid using world coordinates
  if (this->_faces[index].category == 0)
    // triangle face
    this->localCoords(index, bcoords, ecoords);
  else
    // arbitrarily distorted quadrilateral face (i.e. non-rectangular)
    ecoords = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().local(wcoords);
}


template<typename GV, int dimG>
template<typename CoordContainer>
void GeneralSurfaceExtractor<GV, dimG>::globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  for (int i = 0; i < size; ++i)
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], wcoords[i], dimw);

}


template<typename GV, int dimG>
template<typename CoordContainer>
void GeneralSurfaceExtractor<GV, dimG>::localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, int size) const
{
  if (this->_faces[index].category == 0)
  {
    // this is a triangle face
    Dune::array<Coords, simplex_corners> corners;
    unsigned int num_in_self = this->numberInSelf(index);
    Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->type();
    for (int i = 0; i < simplex_corners; ++i)
      corners[i] = cornerLocalInRefElement<ctype, dimw>(gt, num_in_self, i);
    for (int i = 0; i < size; ++i)
      interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], ecoords[i], dimw);
  }
  else
  {
    // this is a quadrilateral face
    CoordContainer wcoords;
    this->localAndGlobalCoords(index, bcoords, ecoords, wcoords, size);
  }
}


template<typename GV, int dimG>
template<typename CoordContainer>
void GeneralSurfaceExtractor<GV, dimG>::localAndGlobalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, CoordContainer &wcoords, int size) const
{
  this->globalCoords(index, bcoords, wcoords, size);
  // for triangles
  if (this->_faces[index].category == 0)
    this->localCoords(index, bcoords, ecoords, size);
  else
  {
    ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].parent)->second->p;
    for (int i = 0; i < size; ++i)
      ecoords[i] = eptr->geometry().local(wcoords[i]);
  }
}



/*   S P E C I A L I Z A T I O N   F O R   2 D   G R I D S   */
template<typename GV>
class GeneralSurfaceExtractor<GV, 2> : public GeneralEdgeExtractor<GV>
{
private:

  typedef GeneralEdgeExtractor<GV>  Base;


public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    cube_corners = 1 << (Base::dim-1)
  };


  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  GeneralSurfaceExtractor(const GV& gv) : Base(gv)
  {
    STDOUTLN("This is GeneralSurfaceExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid working in " << Base::dimw << " space expecting faces of type " << Dune::GeometryType(Dune::GeometryType::cube, Base::dim) << "!");
  }
};


#endif // GENERALSURFACEEXTRACTOR_HH_
