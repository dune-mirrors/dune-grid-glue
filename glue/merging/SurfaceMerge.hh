// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    SurfaceMerge.hh
 *  Version:     1.0
 *  Created on:  Apr 30, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: a wrapper class for models of the SurfaceMergeConcept class
 *  subversion:  $Id$
 *
 */
/**
 * @file SurfaceMerge.hh
 * @brief Class that wraps a model class of the SurfaceMerge concept
 * to ensure uniform access to a part of the functionality of different
 * implementations.
 */

#ifndef SURFACEMERGE_HH_
#define SURFACEMERGE_HH_

#include <vector>
#include "../misc/conceptchecking.hh"
#include "SurfaceMergeConcept.hh"


using namespace std;


/**
 * @class SurfaceMerge
 * @brief wraps a model class of the SurfaceMerge concept
 *
 * This wrapper applies the SurfaceMergeConcept to the given implemening
 * type thus guaranteeing that the class meets all requirements to a valid
 * implementation of the concept.
 * It is to be expected that the implementing class offers additional,
 * specific functionality that does not fit into the uniform interface.
 * That's why even non-const access to this implementation is provided.
 */
template<typename Implementation>
class SurfaceMerge
{
public:

  /// @brief the actual implementation behind this interface
  typedef Implementation ModelType;


private:

  /*   C H E C K   C O N C E P T S   */

  CLASS_REQUIRE(ModelType, SurfaceMergeConcept);


public:

  /*   E X P O R T E D   T Y P E S   A N D   C O N S T A N T S   */

  /// @brief the dimension via compile time constant
  enum { dimw = ModelType::dimw };

  /// @brief the numeric type used in this interface
  typedef typename ModelType::ctype ctype;

  /// @brief the coordinate type used in this interface
  typedef typename ModelType::Coords Coords;


private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief maximum distance between two matched points in the mapping
  ModelType              &realimpl;


public:

  /*   C O N S T R U C T O R S   A N D   D E S T R U C T O R S   */

  SurfaceMerge(ModelType &implementation) : realimpl(implementation) {}


  ~SurfaceMerge() {}


  /*   A C C E S S   T O   T H E   A C T U A L   I M P L E M E N T A T I O N   */

  ModelType& getImplementation() const
  {
    return this->realimpl;
  }



  /*   C O N C E P T   I M P L E M E N T I N G   I N T E R F A C E   */

  /**
   * @brief builds the merged grid
   *
   * Note that the indices are used consequently throughout the whole class interface just like they are
   * introduced here.
   *
   * @param domain_coords the domain vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param domain_simplices array with all domain simplices represented as corner indices into @c domain_coords;
   * the simplices are just written to this array one after another
   * @param target_coords the target vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param target_simplices just like with the domain_simplices and domain_coords
   * @param epsilon the estimate maximum deformation for the contact oracle
   * @param obsDirections If given this function is used to compute the normals for the vertices of domain and target.
   * Otherwise the default algorithm is used to compute them.
   * @return TRUE <=> build successful and merged grid not empty
   */
  bool build(
    const std::vector<ctype>& domain_coords,
    const std::vector<unsigned int>& domain_simplices,
    const std::vector<ctype>& target_coords,
    const std::vector<unsigned int>& target_simplices
    )
  {
    return this->realimpl.build(domain_coords, domain_simplices, target_coords, target_simplices);
  }


  /*   Q U E S T I O N I N G   T H E   M E R G E D   G R I D   */

  /// @brief get the number of simplices in the merged grid
  /// The indices are then in 0..nSimplices()-1
  unsigned int nSimplices() const
  {
    return this->realimpl.nSimplices();
  }


  /**
   * @brief check if given domain simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a domain simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "domainLocalToMerged" or "domainGlobalToMerged".
   * @param idx the index of the domain simplex
   * @return TRUE <=> refined in merged grid
   */
  bool domainSimplexMatched(unsigned int idx) const
  {
    return this->realimpl.domainSimplexMatched(idx);
  }


  /**
   * @brief check if given target simplex could be matched in the merged grid
   *
   * The result of this member even is positive if a target simplex only is
   * partially refined! That means the simplex is not necessarily completely
   * covered in the merged grid. Whether or not a particular point in the simplex
   * was mapped can be asked via "targetLocalToMerged" or "targetGlobalToMerged".
   * @param idx the index of the target simplex
   * @return TRUE <=> refined in merged grid
   */
  bool targetSimplexMatched(unsigned int idx) const
  {
    return this->realimpl.targetSimplexMatched(idx);
  }


  //	/**
  //	 * @brief check if given domain vertex could be matched in the merged grid
  //	 * @param idx the index of the domain vertex
  //	 * @return TRUE <=> contained in merged grid
  //	 */
  //	bool domainVertexMatched(unsigned int idx) const;
  //	{
  //		return this->realimpl.domainVertexMatched(idx);
  //	}
  //
  //
  //	/**
  //	 * @brief check if given target vertex could be matched in the merged grid
  //	 * @param idx the index of the target vertex
  //	 * @return TRUE <=> contained in merged grid
  //	 */
  //	bool targetVertexMatched(unsigned int idx) const
  //	{
  //		return this->realimpl.targetVertexMatched(idx);
  //	}


  /*   M A P P I N G   O N   I N D E X   B A S I S   */

  /**
   * @brief get index of domain parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the domain parent simplex
   */
  unsigned int domainParent(unsigned int idx) const
  {
    return this->realimpl.domainParent(idx);
  }


  /**
   * @brief get index of target parent simplex for given merged grid simplex
   * @param idx index of the merged grid simplex
   * @return index of the target parent simplex
   */
  unsigned int targetParent(unsigned int idx) const
  {
    return this->realimpl.targetParent(idx);
  }


  /**
   * @brief get the merged grid simplices refining a given domain simplex
   * @param idx index of domain simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool domainSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
  {
    return this->realimpl.domainSimplexRefined(idx, indices);
  }


  /**
   * @brief get the merged grid simplices refining a given target simplex
   * @param idx index of target simplex
   * @param indices will be resized first and then filled with the refining simplices
   * @return TRUE <=> given simplex could be matched and is part of the merged grid
   */
  bool targetSimplexRefined(unsigned int idx, std::vector<unsigned int>& indices) const
  {
    return this->realimpl.targetSimplexRefined(idx, indices);
  }


  /*   G E O M E T R I C A L   I N F O R M A T I O N   */

  /**
   * @brief get the domain parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "domainParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return barycentric coordinates in parent domain simplex
   */
  Coords domainParentLocal(unsigned int idx, unsigned int corner) const
  {
    return this->realimpl.domainParentLocal(idx, corner);
  }


  /**
   * @brief get the target parent's simplex local coordinates for a particular merged grid simplex corner
   * (parent's index can be obtained via "targetParent")
   * @param idx the index of the merged grid simplex
   * @param corner the index of the simplex' corner
   * @return barycentric coordinates in parent target simplex
   */
  Coords targetParentLocal(unsigned int idx, unsigned int corner) const
  {
    return this->realimpl.targetParentLocal(idx, corner);
  }


  /**
   * @brief get the simplex local coordinates in target view (oriented)
   * for a point in domain view local coordinates
   * @param idx the index of the merged grid simplex
   * @param local barycentric coordinates in oriented domain view of the simplex
   * @return barycentric coordinates in oriented target view of the simplex
   */
  Coords targetLocals(unsigned int idx, const Coords &local) const
  {
    return this->realimpl.targetLocals(idx, local);
  }


  /**
   * @brief get the simplex local coordinates in domain view (oriented)
   * for a point in target view local coordinates
   * @param idx the index of the merged grid simplex
   * @param local barycentric coordinates in oriented target view of the simplex
   * @return barycentric coordinates in oriented domain view of the simplex
   */
  Coords domainLocals(unsigned int idx, const Coords &local) const
  {
    return this->realimpl.domainLocals(idx, local);
  }

};

#endif // SURFACEMERGE_HH_
