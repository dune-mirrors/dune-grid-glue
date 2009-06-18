// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    generaledgeextractor.hh
 *  Version:     1.0
 *  Created on:  Jun 7, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue-howto
 *  Description: <short_description>
 *  subversion:  $Id$
 *
 */
/**
 * @file generaledgeextractor.hh
 * @brief extractor class that extracts (edge) boundary of 2d grids
 * with arbitrary elements
 */

#ifndef GENERALEDGEEXTRACTOR_HH_
#define GENERALEDGEEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <set>
#include <algorithm>
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/array.hh>
#include <dune/grid/common/geometry.hh>
#include "SurfaceDescriptor.hh"



/**
 * @class GeneralEdgeExtractor
 * @brief extractor implementation that extracts (edge) boundary of 2d grids
 * with arbitrary elements
 */
template<typename GV>
class GeneralEdgeExtractor
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  enum
  {
    dimw = 2
  };

  enum
  {
    dim = 2
  };

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = 2
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

    FaceInfo(unsigned int index_, IndexType parent_, unsigned int num_in_parent_)
      :       parent(parent_), index(index_), num_in_parent(num_in_parent_)
    {}

    /// @brief the index of the parent element (from index set)
    IndexType parent;

    /// @brief the index of this face (in internal storage scheme) // NEEDED??
    unsigned int index : 28;

    /// @brief the number of the face in the parent element
    unsigned int num_in_parent : 4;

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
    ElementInfo(unsigned int idx_, ElementPtr& p_, unsigned int num_ = 1) : idx(idx_), num(num_), p(p_)
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
  GeneralEdgeExtractor(const GV& gv) : _gv(gv)
  {
    STDOUTLN("This is GeneralEdgeExtractor on a <" << GV::dimension << "," << GV::dimensionworld << "> grid working in " << dimw << " space!");
  }


  /**
   * @brief default destructor, frees memory
   */
  ~GeneralEdgeExtractor();


  /*  F U N C T I O N A L I T Y  */

  /**
   * ASSUMPTION:
   * dim == dimw
   *
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a decider function
   * that returns a boolean value for each codim 0 entity. If an entity is "in" that
   * means that all of its boundary faces are considered in the surface.
   * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are simplices.
   * @param descr a descriptor that "selects" the elements whose faces are to add to the surface
   */
  void update(const ElementDescriptor<GV>& descr);


  /**
   * ASSUMPTION:
   * dim == dimw
   *
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a decider function
   * that returns a boolean value for each codim 0 entity's every face.
   * Note that the function is only called for faces that are part of the domain boundary.
   * Inner faces cannot be part of the surface.
   * It is further assumed that only one geometric shape (simplex!) exists on the boundary.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If we the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
   * segments have are simplices.
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


  /**
   * @brief this is a dummy
   *
   * Only required to ensure uniform interfaces with other extractors.
   * @param index the index of the face
   * @param coords local face coords
   * @return @c coords
   */
  Coords toFaceCoords(unsigned int index, const Coords &coords) const
  {
    return coords;
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

}; // end of class GeneralEdgeExtractor



template<typename GV>
GeneralEdgeExtractor<GV>::~GeneralEdgeExtractor()
{
  // only the objects that have been allocated manually have to be
  // deallocated manually again
  // free all the manually allocated memory
  for (typename VertexInfoMap::iterator it = this->_vtxInfo.begin(); it != this->_vtxInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
  for (typename ElementInfoMap::iterator it = this->_elmtInfo.begin(); it != this->_elmtInfo.end(); ++it)
    if (it->second != NULL)
      delete it->second;
}



template<typename GV>
void GeneralEdgeExtractor<GV>::clear()
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



template<typename GV>
void GeneralEdgeExtractor<GV>::update(const ElementDescriptor<GV>& descr)
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

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<FaceInfo> temp_faces;

    // iterate over all codim 0 elemets on the grid
    for (ElementIter elit = this->_gv.template begin<0>(); elit != this->_gv.template end<0>(); ++elit)
    {
      // there are no unwanted geometric shapes, every 2D element is valid!

      // only do sth. if this element is "interesting"
      // implicit cast is done automatically
      if (descr.contains(elit))
      {
        // watch next element
        bool added = false;

        Dune::GeometryType gt = elit->geometry().type();

        // iterate over all intersections of codim 1
        for (IsIter is = this->_gv.ibegin(*elit); is != this->_gv.iend(*elit); ++is)
        {
          // only look at boundary faces
          if (is->boundary())
          {
            // if the first boundary of this element => add it
            if (!added)
            {
              // add an entry to the element info map, the index will be set properly later
              eindex = this->indexSet().template index<0>(*elit);
              this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit);
              added = true;
            }
            else
            {
              // the element info exists already, register the additional face
              this->_elmtInfo[eindex]->num++;
            }

            // get faces index in parent element
            const int num_in_parent = is->numberInSelf();

            // add a new face to the temporary collection
            temp_faces.push_back(FaceInfo(simplex_index, eindex, num_in_parent));

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < simplex_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = orientedSubface<dim>(gt, num_in_parent, i);

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
          }
        }                         // end loop over intersections
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



  //	const char prefix[] = "GeneralEdgeExtractor: ";
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


template<typename GV>
void GeneralEdgeExtractor<GV>::update(const FaceDescriptor<GV>& descr)
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
      // there are no unwanted geometric shapes, every 2d element is valid!

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
        this->_elmtInfo[eindex] = new ElementInfo(simplex_index, elit, boundary_faces.size());

        Dune::GeometryType gt = elit->geometry().type();

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (typename std::set<int>::const_iterator sit = boundary_faces.begin(); sit != boundary_faces.end(); ++sit)
        {
          // add a new face to the temporary collection
          temp_faces.push_back(FaceInfo(simplex_index, eindex, *sit));

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


  //	const char prefix[] = "GeneralEdgeExtractor: ";
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


template<typename GV>
inline void GeneralEdgeExtractor<GV>::globalCoords(unsigned int index, const Coords &bcoords, Coords &wcoords) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, wcoords, dimw);
}


template<typename GV>
inline void GeneralEdgeExtractor<GV>::localCoords(unsigned int index, const Coords &bcoords, Coords &ecoords) const
{
  Dune::array<Coords, simplex_corners> corners;
  unsigned int num_in_self = this->numberInSelf(index);
  Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().type();
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = cornerLocalInRefElement<ctype, dimw>(gt, num_in_self, i);
  interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords, ecoords, dimw);
}


template<typename GV>
inline void GeneralEdgeExtractor<GV>::localAndGlobalCoords(unsigned int index, const Coords &bcoords, Coords &ecoords, Coords &wcoords) const
{
  this->localCoords(index, bcoords, ecoords);
  //	wcoords = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().global(ecoords);
  this->globalCoords(index, bcoords, wcoords);
}


template<typename GV>
template<typename CoordContainer>
inline void GeneralEdgeExtractor<GV>::globalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &wcoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = this->_coords[this->_faces[index].corners[i].idx].coord;
  for (int i = 0; i < size; ++i)
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], wcoords[i], dimw);
}


template<typename GV>
template<typename CoordContainer>
inline void GeneralEdgeExtractor<GV>::localCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  unsigned int num_in_self = this->numberInSelf(index);
  Dune::GeometryType gt = this->_elmtInfo.find(this->_faces[index].parent)->second->p->geometry().type();
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = cornerLocalInRefElement<ctype, dimw>(gt, num_in_self, i);
  for (int i = 0; i < size; ++i)
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], ecoords[i], dimw);
}


template<typename GV>
template<typename CoordContainer>
inline void GeneralEdgeExtractor<GV>::localAndGlobalCoords(unsigned int index, const CoordContainer &bcoords, CoordContainer &ecoords, CoordContainer &wcoords, int size) const
{
  Dune::array<Coords, simplex_corners> corners;
  ElementPtr eptr = this->_elmtInfo.find(this->_faces[index].parent)->second->p;
  Dune::GeometryType gt = eptr->geometry().type();
  unsigned int num_in_self = this->numberInSelf(index);
  for (int i = 0; i < simplex_corners; ++i)
    corners[i] = cornerLocalInRefElement<ctype, dimw>(gt, num_in_self, i);
  for (int i = 0; i < size; ++i)
  {
    interpolateBarycentric<dimw, ctype, Dune::FieldVector<ctype, dimw> >(corners, bcoords[i], ecoords[i], dimw);
    wcoords[i] = eptr->geometry().global(ecoords[i]);
  }
}



#endif // GENERALEDGEEXTRACTOR_HH_
