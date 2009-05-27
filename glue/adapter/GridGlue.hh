// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlue.hh
 *  Version:     1.0
 *  Created on:  Feb 2, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Central component of the module implementing the coupling of two grids.
 *  subversion:  $Id$
 *
 */
/**
 * @file GridGlue.hh
 * @brief Central component of the module implementing the coupling of two grids.
 */


#ifndef GRIDGLUE_HH_
#define GRIDGLUE_HH_

#include "../../dune/common/fixedarray.hh"

#ifdef GRID_GLUE_USE_CONCEPTS
#include "../misc/conceptchecking.h"
#include "../merging/SurfaceMerge.hh"
#endif
#include "../surfaces/GridExtractor.hh"
#include "../merging/ContactMappingSurfaceMerge.hh"
#include "../surfaces/VtkSurfaceWriter.hh"
#include "../surfaces/GridExtractionTraits.hh"
#include "../surfaces/ExtractorSelector.hh"
#include "../surfaces/SurfaceDescriptor.hh"
#include "SimplexGeometry.hh"
#include "RemoteIntersection.hh"
#include "RemoteIntersectionIterators.hh"
#include "CoordinateTransformation.hh"


using namespace MeshClassification;
using namespace ExtractorClassification;

/**
 * @class GridGlue
 * @brief sequential adapter to couple two grids at specified close together boundaries
 *
 *
 * @param GET1 a first traits class to determine the type of surface extractor to use for the domain grid,
 * has to be a model of GridExtractionTraitsConcept
 * @param GET2 a second traits class to determine the type of surface extractor to use for the target grid,
 * has to be a model of GridExtractionTraitsConcept
 * @param SM the surface merging class, has to be a model of the SurfaceMergeConcept
 */
template<
    typename GET1,
    typename GET2,
    typename SM =
      ContactMappingSurfaceMerge<
        GET1::GridView::dimensionworld + static_cast<int>(GET1::GridView::dimensionworld == static_cast<int>(GET1::dimS)),
        typename GET1::GridView::Grid::ctype
        >
    >
class GridGlue
{
private:
  /*   C H E C K   C O N C E P T S   */

  typedef GET1 DomGridExtractionTraits;
#ifdef GRID_GLUE_USE_CONCEPTS
  CLASS_REQUIRE(DomGridExtractionTraits, GridExtractionTraitsConcept);
#endif

  typedef GET2 TarGridExtractionTraits;
#ifdef GRID_GLUE_USE_CONCEPTS
  CLASS_REQUIRE(TarGridExtractionTraits, GridExtractionTraitsConcept);
#endif


  /*   P R I V A T E   T Y P E S   */

  typedef GridGlue<GET1, GET2, SM> This;

  template<typename BSET1, typename BSET2, typename BSM, ExtractorType type1, ExtractorType type2>
  class BuilderImpl;

  class RemoteIntersectionImpl;

  class RemoteIntersectionIteratorImpl;

  class DomainIntersectionIteratorImpl;

  class TargetIntersectionIteratorImpl;

  typedef RemoteIntersectionInterface::RemoteIntersectionIterators<
      RemoteIntersectionImpl,
      RemoteIntersectionIteratorImpl,
      DomainIntersectionIteratorImpl,
      TargetIntersectionIteratorImpl
      > IntersectionIterators;


public:

  /*   P U B L I C   T Y P E S   A N D   C O N S T A N T S   */

  typedef typename DomGridExtractionTraits::GridView DomainGridView;

  typedef typename DomainGridView::Grid DomainGridType;

  typedef typename ExtractorSelector<DomGridExtractionTraits>::ExtractorType DomainExtractor;

  enum { domdim = DomainExtractor::dim };


  typedef typename TarGridExtractionTraits::GridView TargetGridView;

  typedef typename TargetGridView::Grid TargetGridType;

  typedef typename ExtractorSelector<TarGridExtractionTraits>::ExtractorType TargetExtractor;

  enum { tardim = TargetExtractor::dim };


  /// @brief export the world dimension
  /// (must be the same for both extractors!)
  enum { dimw = DomainExtractor::dimw };


  // TODO maybe use traits class to decide which has more precision (DomainGridType::ctype or TargetGridType::ctype) and then take this one
  typedef typename DomainGridType::ctype ctype;

  typedef FieldVector<ctype, dimw>                   Coords;

  typedef CoordinateTransformation<dimw, ctype>      Transformation;


  typedef typename DomainGridView::Traits::template Codim<0>::Entity DomainElement;

  typedef typename DomainGridView::Traits::template Codim<0>::EntityPointer DomainElementPtr;

  typedef typename DomainGridView::Traits::template Codim<DomainGridType::dimension>::Entity DomainVertex;

  typedef typename DomainGridView::Traits::template Codim<DomainGridType::dimension>::EntityPointer DomainVertexPtr;


  typedef typename TargetGridView::Traits::template Codim<0>::Entity TargetElement;

  typedef typename TargetGridView::Traits::template Codim<0>::EntityPointer TargetElementPtr;

  typedef typename TargetGridView::Traits::template Codim<TargetGridType::dimension>::Entity TargetVertex;

  typedef typename TargetGridView::Traits::template Codim<TargetGridType::dimension>::EntityPointer TargetVertexPtr;


#ifdef GRID_GLUE_USE_CONCEPTS
  typedef typename SurfaceMerge<SM>::ModelType Matcher;
#else
  typedef SM Matcher;
#endif

#ifdef GRID_GLUE_USE_CONCEPTS
  typedef BuilderImpl<GET1, GET2, SM, DomainExtractor::type, TargetExtractor::type>  Builder;
#else
  typedef BuilderImpl<GET1, GET2, SM,
      ExtractorClassifier<
          DomainExtractor::GridView::dimension,
          DomainExtractor::GridView::dimensionworld,
          DomainExtractor::dimw>::type,
      ExtractorClassifier<
          TargetExtractor::GridView::dimension,
          TargetExtractor::GridView::dimensionworld,
          TargetExtractor::dimw>::type>
  Builder;
#endif

  typedef RemoteIntersectionInterface::RemoteIntersection<RemoteIntersectionImpl>    RemoteIntersection;

  typedef typename IntersectionIterators::RemoteIntersectionIterator RemoteIntersectionIterator;

  typedef typename IntersectionIterators::RemoteIntersectionDomainIterator DomainIntersectionIterator;

  typedef typename IntersectionIterators::RemoteIntersectionTargetIterator TargetIntersectionIterator;


private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief the "domain" grid view
  const DomainGridView&        _domgv;

  /// @brief the "target" grid view
  const TargetGridView&        _targv;

  /// @brief the domain surface extractor
  DomainExtractor _domext;

  /// @brief the target surface extractor
  TargetExtractor _tarext;

  /// @brief the surface merging utility
#ifdef GRID_GLUE_USE_CONCEPTS
  SurfaceMerge<Matcher>        _sm;
#else
  Matcher                      &_sm;
#endif

  /// @brief the builder utility
  Builder _builder;


  /// @brief an invalid interesection object used as dummy and
  /// also as recognizable end object of iterations over intersections
  mutable RemoteIntersectionImpl NULL_INTERSECTION;

  /// @brief a vector with intersection elements
  mutable vector<RemoteIntersectionImpl>   _intersections;


protected:

  /**
   * @brief after building the merged grid the intersection can be updated
   * through this method (for internal use)
   */
  void updateIntersections()
  {
    // build the interesections array again
    this->_intersections.resize(this->_sm.nSimplices(), this->NULL_INTERSECTION);
    for (unsigned int i = 0; i < this->_sm.nSimplices(); ++i)
      this->_intersections[i] = RemoteIntersectionImpl(this, i);

    STDOUTLN("GridGlue::updateIntersections : The number of overlaps is " << this->_sm.nSimplices());
  }


public:

  /*   C O N S T R U C T O R S   A N D   D E S T R U C T O R S   */

  /**
   * @brief constructor
   *
   * Initializes components but does not "glue" the surfaces. The surfaces
   * are extracted from the grids here though.
   * @param gv1 the domain grid view
   * @param gv2 the target grid view
   * @param matcher The matcher object that is used to compute the merged grid. This class has
   * to be a model of the SurfaceMergeConcept.
   */
  GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Matcher &matcher);


  /*  S E T T E R S  */


  /*   G E T T E R S   */

  /**
   * @brief getter for the domain grid view
   * @return the object
   */
  const DomainGridView& domainGridView() const
  {
    return this->_domgv;
  }


  /**
   * @brief getter for the target grid view
   * @return the object
   */
  const TargetGridView& targetGridView() const
  {
    return this->_targv;
  }


  /**
   * @brief getter for the builder
   * @return the object
   */
  Builder& builder()
  {
    return this->_builder;
  }


  /**
   * @brief getter for the surface matcher utility
   *
   * This grants access to the surface matcher. This Matcher class has to be a model
   * of the SurfaceMergeConcept, but since different implementations of matchers
   * may provide different configuration possibilities this part has to remain
   * implementation specific. Through this getter one can retrieve the internal
   * matcher and configure it before its "build" member is called.
   * @return a (non-const) reference the object
   */
  Matcher& matcher()
  {
#ifdef GRID_GLUE_USE_CONCEPTS
    return this->_sm.getImplementation();
#else
    return this->_sm;
#endif
  }


  /*   F U N C T I O N A L I T Y   */

  //	/**
  //	 * @brief tells whether a vertex could be mapped
  //	 * @param v the vertex
  //	 * @return TRUE <=> if mapped successfully
  //	 */
  //	bool vertexMappedToTarget(const DomainVertex& v) const;


  //	/**
  //	 * @brief tells whether a vertex could be mapped
  //	 * @param v the vertex
  //	 * @return TRUE <=> if mapped successfully
  //	 */
  //	bool vertexMappedToDomain(const TargetVertex& v) const;


  /**
   * @brief tells whether a codim 0 entity's face(s) (or at least a part)
   * could be mapped
   *
   * For the given entity  could be mapped, if the latter is the case the number of the
   * face in the parent element is returned.
   * Note: Calling this function with only @c e given and checking the return value
   * (-1 or not) is the easiest way to determine whether any of the entity's faces were mapped.
   * @param e the element
   * @param index number of the first face that is checked, faces with lower index are ignored
   * @return -1 if there is no face mapped with number >=@c index, else next face's number
   */
  int domainEntityNextFace(const DomainElement& e, int index = 0) const;


  /**
   * @brief tells whether a codim 0 entity's face(s) (or at least a part)
   * could be mapped
   *
   * For the given entity  could be mapped, if the latter is the case the number of the
   * face in the parent element is returned.
   * Note: Calling this function with only @c e given and checking the return value
   * (-1 or not) is the easiest way to determine whether any of the entity's faces were mapped.
   * @param e the element
   * @param index number of the first face that is checked, faces with lower index are ignored
   * @return -1 if there is no face mapped with number >=@c index, else next face's number
   */
  int targetEntityNextFace(const TargetElement& e, int index = 0) const;


  //	/**
  //	 * @brief tells how many faces of an element could be (partially) mapped
  //	 * @param e the element
  //	 * @param faces An ALLOCATED array to which the indices of the faces in the Dune
  //	 * element are written starting at index 0. The user has to choose its size
  //	 * appropriatly (e.g. the number of faces for this element)
  //	 * @return the number of (partially) mapped faces for this element
  //	 */
  //	int facesMappedToTarget(const DomainElement& e, unsigned int* faces) const;
  //
  //
  //	/**
  //	 * @brief tells how many faces of an element could be (partially) mapped
  //	 * @param e the element
  //	 * @param faces An ALLOCATED array to which the indices of the faces in the Dune
  //	 * element are written starting at index 0. The user has to choose its size
  //	 * appropriatly (e.g. the number of faces for this element)
  //	 * @return the number of (partially) mapped faces for this element
  //	 */
  //	int facesMappedToDomain(const TargetElement& e, unsigned int* faces) const;


  /*   I N T E R S E C T I O N S   A N D   I N T E R S E C T I O N   I T E R A T O R S   */

  /**
   * @brief gets an iterator over all remote intersections in the merged grid between domain and target
   *
   * @return the iterator
   */
  RemoteIntersectionIterator iremotebegin() const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  DomainIntersectionIterator idomainbegin(const DomainElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @return the iterator
   */
  DomainIntersectionIterator idomainbegin(const DomainElement& e) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  TargetIntersectionIterator itargetbegin(const TargetElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @return the iterator
   */
  TargetIntersectionIterator itargetbegin(const TargetElement& e) const;


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  RemoteIntersectionIterator iremoteend() const
  {
    return RemoteIntersectionIterator(RemoteIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  DomainIntersectionIterator idomainend() const
  {
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }


  /**
   * @brief gets the (general) end-iterator for iterations over target codim 0 entities' faces
   *
   * @return the iterator
   */
  TargetIntersectionIterator itargetend() const
  {
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }



};


/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

template<typename GET1, typename GET2, typename SM>
GridGlue<GET1, GET2, SM>::GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Matcher &matcher)
  : _domgv(gv1), _targv(gv2),
    _domext(gv1), _tarext(gv2), _sm(matcher),
    _builder(*const_cast<GridGlue<GET1, GET2, SM>*>(this)),
    NULL_INTERSECTION(this, -1), _intersections(0, NULL_INTERSECTION)
{
  STDOUTLN("GridGlue: Constructor succeeded!");
}


//template<typename GET1, typename GET2, typename SM>
//inline bool GridGlue<GET1, GET2, SM>::vertexMappedToTarget(const DomainVertex& v) const
//{
//	int index = this->_domext.coordinateIndex(v);
//	if (index < 0)
//		return false;
//	return this->_sm.domainVertexMatched(index);
//}
//
//
//template<typename GET1, typename GET2, typename SM>
//inline bool GridGlue<GET1, GET2, SM>::vertexMappedToDomain(const TargetVertex& v) const
//{
//	int index = this->_tarext.coordinateIndex(v);
//	if (index < 0)
//		return false;
//	return this->_sm.targetVertexMatched(index);
//}


template<typename GET1, typename GET2, typename SM>
int GridGlue<GET1, GET2, SM>::domainEntityNextFace(const DomainElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->_domext.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (this->_domext.numberInSelf(first) < index || !this->_sm.domainSimplexMatched(first)))
    first++;
  if (first == count)
    return -1;             // no more faces
  else
    return this->_domext.numberInSelf(first);             // found, return the face's number
}


template<typename GET1, typename GET2, typename SM>
int GridGlue<GET1, GET2, SM>::targetEntityNextFace(const TargetElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->_tarext.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (this->_tarext.numberInSelf(first) < index || !this->_sm.targetSimplexMatched(first)))
    first++;
  if (first == count)
    return -1;             // no more faces
  else
    return this->_tarext.numberInSelf(first);             // found, return the face's number
}


//template<typename GET1, typename GET2, typename SM>
//int GridGlue<GET1, GET2, SM>::facesMappedToTarget(const DomainElement& e, unsigned int* faces) const
//{
//	// first check if the element forms a part of the extracted surface
//	int first, count;
//	bool in_surface = this->_domext.faceIndices(e, first, count);
//	if (!in_surface)
//		return 0;
//
//	// now we have found at least one face, but if there is more
//	// than one face in the map check whether the right one could be mapped
//	int result = 0;
//	for (count += first; first < count; ++first)
//	{
//		if (this->_sm.domainSimplexMatched(first))
//			faces[result++] = this->_domext.numberInSelf(first);
//	}
//
//	// no matter how many faces of the element were matched,
//	// now we have the right index
//	return result;
//}
//
//
//template<typename GET1, typename GET2, typename SM>
//int GridGlue<GET1, GET2, SM>::facesMappedToDomain(const TargetElement& e, unsigned int* faces) const
//{
//	// first check if the element forms a part of the extracted surface
//	int first, count;
//	bool in_surface = this->_tarext.faceIndices(e, first, count);
//	if (!in_surface)
//		return 0;
//
//	// now we have found at least one face, but if there is more
//	// than one face in the map check whether the right one could be mapped
//	int result = 0;
//	for (count += first; first < count; ++first)
//	{
//		if (this->_sm.targetSimplexMatched(first))
//			faces[result++] = this->_tarext.numberInSelf(first);
//	}
//
//	// no matter how many faces of the element were matched,
//	// now we have the right index
//	return result;
//}


template<typename GET1, typename GET2, typename SM>
typename GridGlue<GET1, GET2, SM>::RemoteIntersectionIterator GridGlue<GET1, GET2, SM>::iremotebegin() const
{
  return RemoteIntersectionIterator(RemoteIntersectionIteratorImpl(this->_intersections[0]));
}


template<typename GET1, typename GET2, typename SM>
typename GridGlue<GET1, GET2, SM>::DomainIntersectionIterator GridGlue<GET1, GET2, SM>::idomainbegin(const DomainElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_domext.faceIndices(e, first, count);
  if (!in_surface)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));

  count += first;
  while (first < count)
  {
    if (this->_domext.numberInSelf(first) == num && this->_sm.domainSimplexMatched(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      vector<unsigned int> global_results;
      vector<unsigned int> local_results;

      // get the remote intersections
      this->_sm.domainSimplexRefined(first, global_results);
      while (++first < count && this->_domext.numberInSelf(first) == num && this->_sm.domainSimplexRefined(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2, typename SM>
typename GridGlue<GET1, GET2, SM>::DomainIntersectionIterator GridGlue<GET1, GET2, SM>::idomainbegin(const DomainElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_domext.faceIndices(e, first, count);
  if (!in_surface)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));


  // now accumulate all remote intersections of the element's faces
  vector<unsigned int> global_results(0, 0);
  vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->_sm.domainSimplexRefined(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));
  else
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2, typename SM>
typename GridGlue<GET1, GET2, SM>::TargetIntersectionIterator GridGlue<GET1, GET2, SM>::itargetbegin(const TargetElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_tarext.faceIndices(e, first, count);
  if (!in_surface)
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));

  count += first;
  while (first < count)
  {
    if (this->_tarext.numberInSelf(first) == num && this->_sm.targetSimplexMatched(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      vector<unsigned int> global_results;
      vector<unsigned int> local_results;

      // get the remote intersections
      this->_sm.targetSimplexRefined(first, global_results);
      while (++first < count && this->_tarext.numberInSelf(first) == num && this->_sm.targetSimplexRefined(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2, typename SM>
typename GridGlue<GET1, GET2, SM>::TargetIntersectionIterator GridGlue<GET1, GET2, SM>::itargetbegin(const TargetElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_tarext.faceIndices(e, first, count);
  if (!in_surface)
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));


  // now accumulate all remote intersections of the element's faces
  vector<unsigned int> global_results(0, 0);
  vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->_sm.targetSimplexRefined(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));
  else
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
}



// include implementation of subclass BuilderImpl
#include "GridGlueBuilderImpl.hh"

// include implementation of subclass RemoteIntersectionImpl
#include "GridGlueRemoteIntersectionImpl.hh"

// include implementation of subclasses DomainIntersectionIteratorImpl and TargetIntersectionIteratorImpl
#include "GridGlueRemoteIntersectionIteratorImpl.hh"


#endif // GRIDGLUE_HH_
