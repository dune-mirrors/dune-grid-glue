// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SurfaceMergeConcept.hh
 *  Version:     1.0
 *  Created on:  Jan 25, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: concept class for the SurfaceMergeConcept
 *  subversion:  $Id$
 *
 */
/**
 * @file SurfaceMergeConcept.hh
 * @brief Can be used to check an implementation of the SurfaceMerge concept.
 */

#ifndef SURFACEMERGECONCEPT_HH_
#define SURFACEMERGECONCEPT_HH_

#include <dune/common/fvector.hh>
#include <vector>

/**
 * @class SurfaceMergeConcept
 *
 * @brief a concept class defining the syntactical requirements of valid wrapper class
 * used to encapsulate surface merging functionality used in the grid glue adapter
 */
template<class T>
struct SurfaceMergeConcept
{

  void constraints()
  {
    /* BUILDING THE MERGED GRID */

    /*
     * @brief builds the merged grid
     *
     * Note that the indices are used consequently throughout the whole class interface just like they are
     * introduced here.
     *
     * @param const std::vector<ctype>& the domain vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
     * @param const std::vector<unsigned int>& array with all domain simplices represented as corner indices into @c domain_coords;
     * the simplices are just written to this array one after another
     * @param const std::vector<ctype>& the target vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
     * @param const std::vector<unsigned int>& just like with the domain_simplices and domain_coords
     * @return TRUE <=> build successful and merged grid not empty
     */
    this->_b = this->_sm.build(this->_coords, this->_indices, this->_coords, this->_indices);


    /* QUESTIONING THE MERGED GRID */

    // @brief get the number of simplices in the merged grid
    // The indices are then in 0..nSimplices()-1
    this->_index = this->_sm.nSimplices();

    /*
     * @brief check if the domain resp. target simplex at _index could be matched in the merged grid
     * @param unsigned int the index of the domain/target simplex
     * @return TRUE <=> refined in merged grid
     */
    this->_b = this->_sm.domainSimplexMatched(this->_index);
    this->_b = this->_sm.targetSimplexMatched(this->_index);

    //          /*
    //           * @brief check if the domain resp. target vertex at _index could be matched in the merged grid
    //           * @param unsigned int the index of the domain/target vertex
    //           * @return TRUE <=> contained in merged grid
    //           */
    //          this->_b = this->_sm.domainVertexMatched(this->_index);
    //          this->_b = this->_sm.targetVertexMatched(this->_index);


    /* MAPPING ON INDEX BASIS */

    /*
     * @brief get index of domain/target parent simplex for given merged grid simplex
     * @param unsigned int index of the merged grid simplex
     * @return index of the domain/target parent simplex
     */
    this->_index = this->_sm.domainParent(this->_index);
    this->_index = this->_sm.targetParent(this->_index);

    /*
     * @brief get the merged grid simplices refining a given domain/target simplex
     * @param unsigned int index of domain/target simplex
     * @param std::vector<unsigned int>& will be resized first and then filled with the refining simplices
     * @return TRUE <=> given simplex could be matched and is part of the merged grid
     */
    this->_b = this->_sm.domainSimplexRefined(this->_index, this->_indices);
    this->_b = this->_sm.targetSimplexRefined(this->_index, this->_indices);


    /* GEOMETRICAL INFORMATION */

    /*
     * @brief get the domain/target parent's simplex local coordinates for a particular merged grid simplex corner
     * (parent's index can be obtained via "domainParent" resp. "targetParent")
     *
     * The order of corners is expected to be consistent with the order of corners of the parent (in domain
     * resp. target), which means that the orientation of the subface should be the same as the orientation of the face itself.
     * @param unsigned int the index of the merged grid simplex
     * @param unsigned int the index of the simplex' corner
     * @return barycentric coordinates in parent domain/target simplex
     */
    this->_c = this->_sm.domainParentLocal(this->_index , this->_index);
    this->_c = this->_sm.targetParentLocal(this->_index , this->_index);


    /*
     * @brief get the simplex local coordinates in target/domain view (oriented)
     * for a point in domain/target view local coordinates
     * @param idx the index of the merged grid simplex
     * @param local barycentric coordinates in oriented domain/target view of the simplex
     * @return barycentric coordinates in oriented target/domain view of the simplex
     */
    this->_c = this->_sm.targetLocals(this->_index, this->_c);
    this->_c = this->_sm.domainLocals(this->_index, this->_c);

    // implicit: the model class has to have a default destructor
  }

private:
  // These types and constants must be declared in the model's interface
  enum { dimw = T::dimw };

  typedef typename T::ctype ctype;

  // helpful types used here
  typedef typename T::Coords Coords;


  // dummy member variables used in interface definition in "constraints()"

  T _sm;

  unsigned int _index;

  bool _b;

  ctype _scalar;

  std::vector<unsigned int>        _indices;

  std::vector<ctype>               _coords;

  Coords _c;

  Dune::FieldVector<Coords, dimw>  _corners;
};




#endif // SURFACEMERGECONCEPT_HH_
