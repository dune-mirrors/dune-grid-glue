// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GeneralSurfaceExtractor.hh
 *  Version:     1.0
 *  Created on:  Feb 3, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: grid extractor implementation for hybrid surface grids (not implemented!)
 *  subversion:  $Id$
 *
 */
/**
 * @file GeneralSurfaceExtractor.hh
 * @brief grid extractor implementation for hybrid surface grids (not implemented!)
 */

#ifndef GENERALSURFACEEXTRACTOR_HH_
#define GENERALSURFACEEXTRACTOR_HH_

#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>

#include "SurfaceDescriptor.hh"


template<typename GV>
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

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimw>                                    Coords;

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



  /************************** MEMBER VARIABLES ************************/

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief default constructor, initializes member variables with null values
   * @param g the grid object
   * @þaram descriptor object of a helper class telling which entities are part of the surface
   * @param eps a tolerance value used in coordinate comparison
   */
  GeneralSurfaceExtractor(const GV& gv)
  {
    DUNE_THROW(NotImplemented, "GeneralSurfaceExtractor not yet implemented");
  }


  //	/**
  //	 * @brief default destructor, frees memory
  //	 */
  //	~GeneralSurfaceExtractor();


  /*  F U N C T I O N A L I T Y  */

  //	/**
  //	 * ASSUMPTION:
  //	 * dim == dimw
  //	 *
  //	 * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
  //	 * with the topology of the surface written to them. The description of the
  //	 * surface part that is to be extracted is given in form of a mapper or set object
  //	 * @c m specifying an index set with codimension 0 entities near or on the boundary.
  //	 * It is assumed that only one geometric shape exists on the boundary.
  //	 * The template parameter n then denotes the number of corners per boundary element
  //	 * (e.g. n==3 for triangles in a 3D grid of tetrahedra).
  //	 *
  //	 * Assumed that we are in 2D the coords array will have the structure
  //	 * x0 y0 x1 y1 ... x(n-1) y(n-1)
  //	 * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
  //	 * index 1 is associated with the position x1. If we the surface consists of triangles
  //	 * we have always groups of 3 indices describing one triangle.
  //	 *
  //	 * Hint: The exception Dune::MathError is thrown if not all "interesting" boundary
  //	 * segments have are simplices.
  //	 */
  //	void update();
  //
  //
  //	/**
  //	 * @brief delete everything build up so far and free the memory
  //	 */
  //	void clear();
  //
  //
  //	/*  S E T T E R S  */
  //
  //	/**
  //	 * @brief setter for the surface descriptor
  //	 * @param value the new reference
  //	 */
  //	void setDescriptor(const D& d_)
  //	{
  ////		this->_descriptor = d_;
  //	}
  //
  //
  //	/*  G E T T E R S  */
  //
  //	/**
  //	 * @brief getter for the coordinates array
  //	 * It is strongly recommended not to modify its contents.
  //	 * Deallocation is done in this class.
  //	 * @return the _coordinates array
  //	 */
  //	const vector<ctype>& coords() const
  //	{
  ////		return this->_coords;
  //		return vector<ctype>(0);
  //	}
  //
  //	/**
  //		 * @brief getter for the count of coordinates
  //		 * @return the count
  //		 */
  //	unsigned int nCoords() const
  //	{
  ////		return this->_vtxIndex.size();
  //		return 0;
  //	}
  //
  //	/**
  //	 * @brief getter for the indices array
  //	 * It is strongly recommended not to modify its contents.
  //	 * Deallocation is done in this class.
  //	 * @return the _indices array
  //	 */
  //	const vector<unsigned int>& indices() const
  //	{
  ////		return this->_indices;
  //		return vector<unsigned int>(0);
  //	}
  //
  //
  //	/**
  //	 * @brief getter for internally used index set (grid's index set)
  //	 * @return the index set
  //	 */
  //	const IndexSet& indexSet() const
  //	{
  //		return this->_gv.indexSet();
  //	}
  //
  //
  //	/**
  //	 * @brief getter for the index of an entity of codim cc
  //	 * @return the index specified by the grid's index set
  //	 */
  //	template<int cc>
  //	IndexType index(const typename GV::Traits::template Codim<cc>::Entity& e) const
  //	{
  //		return this->indexSet().template index<cc>(e);
  //	}
  //
  //
  //	/**
  //	 * @brief gets index of coordinate in _coords associated with given vertex
  //	 * @return the index if possible, -1 else
  //	 */
  //	int coordinateIndex(const VertexPtr& p) const
  //	{
  ////		typename VertexInfoMap::const_iterator it = this->_vtxInfo.find(this->index<dim>(*p));
  ////		if (it == this->_vtxInfo.end())
  ////			return -1;
  ////		else
  ////			return it->second->idx;
  //		return -1;
  //	}
  //
  //
  //	/**
  //	 * @brief gets index of coordinate in _coords associated with given vertex
  //	 * @return the index if possible, -1 else
  //	 */
  //	int firstFaceIndex(const ElementPtr& p) const
  //	{
  ////		typename ElementInfoMap::const_iterator it = this->_elmtInfo.find(this->index<0>(*p));
  ////		if (it == this->_elmtInfo.end())
  ////			return -1;
  ////		else
  ////			return it->second->idx;
  //		return -1;
  //	}
  //
  //
  //	/**
  //	 * @brief gets index of coordinate in _coords associated with given vertex
  //	 * @param p the element
  //	 * @param first will contain the first index if found, else -1
  //	 * @param count will contain the number of faces if found, else 0
  //	 * @return success
  //	 */
  //	bool faceIndices(const ElementPtr& p, int& first, int& count) const
  //	{
  ////		first = this->firstFaceIndex(p);
  ////		count = 0;
  ////		if (first < 0)
  ////			return false;
  ////		else
  ////			while (this->_faceInfo[first].parent == this->_faceInfo[++count + first].parent);
  ////		return true;
  //		return false;
  //	}
  //
  //
  //	/**
  //	 * @brief gets the parent element for a given face index
  //	 * @param index the index of the face
  //	 * @param p if successful the reference will be set
  //	 * @return success
  //	 */
  //	bool element(unsigned int index, ElementPtr& p) const
  //	{
  ////		if (index < this->_faceInfo.size())
  ////		{
  ////			p = (this->_elmtInfo.find(this->_faceInfo[index].parent))->second->p;
  ////			return true;
  ////		}
  ////		else
  //			return false;
  //	}
  //
  //
  //	/**
  //	 * @brief gets the vertex for a given coordinate index
  //	 * @param index the index of the coordinate
  //	 * @param p if successful the reference will be set
  //	 * @return success
  //	 */
  //	bool vertex(unsigned int index, VertexPtr& p) const
  //	{
  ////		if (index < this->_vtxIndex.size())
  ////		{
  ////			p = (this->_vtxInfo.find(this->_vtxIndex[index]))->second->p;
  ////			return true;
  ////		}
  ////		else
  //			return false;
  //	}
  //
  //
  //	/**
  //	 * @brief gets the indices of all faces with the given coordinate as corner
  //	 * @param index the index of the coordinate
  //	 * @return array with #elements in 1st element if given index was legal, NULL else
  //	 * DO NOT MODIFY THE ARRAY'S CONTENT!
  //	 */
  //	bool parentFaces(unsigned int index,  unsigned int const*& parents) const
  //	{
  ////		if (index >= this->_vtxIndex.size())
  ////			return false;
  ////		// index valid
  ////		parents = this->_vtxFaces[index];
  ////		return true;
  //		return false;
  //	}

}; // end of class GeneralSurfaceExtractor

//
//
//template<typename GV>
//GeneralSurfaceExtractor<GV>::~GeneralSurfaceExtractor()
//{
//}
//
//
//
//template<typename GV>
//void GeneralSurfaceExtractor<GV>::clear()
//{
//}
//
//
//
//template<typename GV>
//void GeneralSurfaceExtractor<GV>::update()
//{
//	// free everything there is in this object
//	this->clear();
//}


#endif // GENERALSURFACEEXTRACTOR_HH_
