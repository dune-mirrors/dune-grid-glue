// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef PSURFACE_EXTERN
#include "config.h"
#endif

#include <dune/grid-glue/merging/psurfacemerge.hh>
#include <psurface/ContactMapping.h>

#if HAVE_PSURFACE

template<int dim, int dimworld, typename T>
void PSurfaceMerge<dim, dimworld, T>::build(const std::vector<Dune::FieldVector<T,dimworld> >& domain_coords,
                                            const std::vector<unsigned int>& domain_elements,
                                            const std::vector<Dune::GeometryType>& domain_element_types,
                                            const std::vector<Dune::FieldVector<T,dimworld> >& target_coords,
                                            const std::vector<unsigned int>& target_elements,
                                            const std::vector<Dune::GeometryType>& target_element_types
                                            )
{
  PSURFACE_NAMESPACE ContactMapping<dim+1,ctype> cm;

  // //////////////////////////////////////////
  //   check input data for consistency
  // //////////////////////////////////////////
#ifndef NDEBUG
  // first the domain side
  unsigned int expectedCorners = 0;

  for (size_t i=0; i<domain_element_types.size(); i++) {

    // Check whether GeometryType has the correct dimension
    if (domain_element_types[i].dim() != dim)
      DUNE_THROW(Dune::GridError, "You cannot hand a " << domain_element_types[i]
                                                       << " to a " << dim << "-dimensional PSurfaceMerge!");

    expectedCorners += Dune::ReferenceElements<ctype,dim>::general(domain_element_types[i]).size(dim);

  }

  if (domain_elements.size() != expectedCorners)
    DUNE_THROW(Dune::GridError, domain_elements.size() << " element corners were handed over, "
               " but " << expectedCorners << " were expected!");

  // now the same for the target side
  expectedCorners = 0;

  for (size_t i=0; i<target_element_types.size(); i++) {

    // Check whether GeometryType has the correct dimension
    if (target_element_types[i].dim() != dim)
      DUNE_THROW(Dune::GridError, "You cannot hand a " << target_element_types[i]
                                                       << " to a " << dim << "-dimensional PSurfaceMerge!");

    expectedCorners += Dune::ReferenceElements<ctype,dim>::general(target_element_types[i]).size(dim);

  }

  if (target_elements.size() != expectedCorners)
    DUNE_THROW(Dune::GridError, target_elements.size() << " element corners were handed over, "
               " but " << expectedCorners << " were expected!");
#endif

  // //////////////////////////////////////////////////////////////////////////////////
  //   copy domain and target simplices to internal arrays
  //   split quadrilaterals to triangles because psurface can only handle triangles
  // //////////////////////////////////////////////////////////////////////////////////

  // First count the number of simplices we expect
  size_t numDomainSimplices = 0;
  size_t numTargetSimplices = 0;

  if (dim==1) {
    numDomainSimplices = domain_element_types.size();
    numTargetSimplices = target_element_types.size();
  } else {

    for (size_t i=0; i<domain_element_types.size(); i++)
      if (domain_element_types[i].isSimplex())
        numDomainSimplices++;
      else if (domain_element_types[i].isQuadrilateral())
        numDomainSimplices += 2;
      else
        DUNE_THROW(Dune::GridError, "Unknown element " << domain_element_types[i] << " found!");

    for (size_t i=0; i<target_element_types.size(); i++)
      if (target_element_types[i].isSimplex())
        numTargetSimplices++;
      else if (target_element_types[i].isQuadrilateral())
        numTargetSimplices += 2;
      else
        DUNE_THROW(Dune::GridError, "Unknown element " << target_element_types[i] << " found!");

  }

  domi_.resize(numDomainSimplices);
  tari_.resize(numTargetSimplices);

  // copy and split the domain surface
  Dune::BitSetVector<1> domainIsSecondTriangle(numDomainSimplices, false);
  std::vector<int> unsplitDomainElementNumbers(numDomainSimplices);

  size_t vertexCounter = 0;
  size_t simplexCounter = 0;

  for (size_t i = 0; i < domain_element_types.size(); ++i) {

    // dim is known at compile time, hence if dim==1 the following conditional is removed
    if (dim == 1 || domain_element_types[i].isSimplex()) {
      for (int j=0; j<dim+1; j++)
        domi_[simplexCounter][j] = domain_elements[vertexCounter++];

      unsplitDomainElementNumbers[simplexCounter] = i;

      simplexCounter++;

    } else {
      // quadrilateral: split it in two triangles
      domi_[simplexCounter][0] = domain_elements[vertexCounter];
      domi_[simplexCounter][1] = domain_elements[vertexCounter+1];
      domi_[simplexCounter][2] = domain_elements[vertexCounter+2];

      unsplitDomainElementNumbers[simplexCounter] = i;
      simplexCounter++;

      domi_[simplexCounter][0] = domain_elements[vertexCounter+3];
      domi_[simplexCounter][1] = domain_elements[vertexCounter+2];
      domi_[simplexCounter][2] = domain_elements[vertexCounter+1];

      domainIsSecondTriangle[simplexCounter] = true;
      unsplitDomainElementNumbers[simplexCounter] = i;

      simplexCounter++;
      vertexCounter += 4;
    }

  }

  // copy and split the target surface
  Dune::BitSetVector<1> targetIsSecondTriangle(numTargetSimplices, false);
  std::vector<int> unsplitTargetElementNumbers(numTargetSimplices);

  vertexCounter  = 0;
  simplexCounter = 0;

  for (size_t i = 0; i < target_element_types.size(); ++i) {

    // dim is known at compile time, hence if dim==1 the following conditional is removed
    if (dim==1 || target_element_types[i].isSimplex()) {
      for (size_t j=0; j<dim+1; j++)
        tari_[simplexCounter][j] = target_elements[vertexCounter++];

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;
    } else {
      // quadrilateral: split it in two triangles
      tari_[simplexCounter][0] = target_elements[vertexCounter];
      tari_[simplexCounter][1] = target_elements[vertexCounter+1];
      tari_[simplexCounter][2] = target_elements[vertexCounter+2];

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;

      tari_[simplexCounter][0] = target_elements[vertexCounter+3];
      tari_[simplexCounter][1] = target_elements[vertexCounter+2];
      tari_[simplexCounter][2] = target_elements[vertexCounter+1];

      targetIsSecondTriangle[simplexCounter] = true;

      unsplitTargetElementNumbers[simplexCounter] = i;
      simplexCounter++;
      vertexCounter += 4;
    }

  }

  // ////////////////////////////////////////////////////////////
  //   copy the coordinates to internal arrays of coordinates
  // ////////////////////////////////////////////////////////////
  domc_.resize(domain_coords.size());
  for (size_t i = 0; i < domc_.size(); ++i)
    for (size_t j = 0; j < dimworld; ++j)
      domc_[i][j] = domain_coords[i][j];

  tarc_.resize(target_coords.size());
  for (size_t i = 0; i < tarc_.size(); ++i)
    for (size_t j = 0; j < dimworld; ++j)
      tarc_[i][j] = target_coords[i][j];

  // psurface doesn't actually support the case dim==dimworld.  Therefore we
  // use a trick: we just embed everything in a space of dimension dimworld+1
  if (dim==dimworld) {
    for (size_t i = 0; i < domc_.size(); ++i)
      domc_[i][dim] = 0;

    for (size_t i = 0; i < domc_.size(); ++i)
      tarc_[i][dim] = 1;
  }

  std::cout << "PSurfaceMerge building merged grid..." << std::endl;

  if (dim==dimworld) {

    ConstantDirection<+1> positiveDirection;
    ConstantDirection<-1> negativeDirection;

    // compute the merged grid using the psurface library
    cm.build(domc_, domi_,tarc_, tari_,
             &positiveDirection, &negativeDirection);
  } else {

    // compute the merged grid using the psurface library
    cm.build(domc_, domi_,tarc_, tari_,
             domainDirections_, targetDirections_);

  }

  std::cout << "Finished building merged grid!" << std::endl;

  // get the representation from the contact mapping object
  std::vector<PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype> > overlaps;
  cm.getOverlaps(overlaps);

  // /////////////////////////////////////////////////////////////////////////////
  //  If overlaps refer to triangular elements that have been created
  //  by splitting a quadrilateral we have to reorder the element references
  //  and the local coordinates.
  // /////////////////////////////////////////////////////////////////////////////

  for (size_t i=0; i<overlaps.size(); i++) {

    //
    if (domainIsSecondTriangle[overlaps[i].tris[0]][0]) {

      assert(dim==2);

      // loop over the intersection corners
      for (int j=0; j<dim+1; j++) {

        // kinda umstaendlich: go from barycentric to reference coords, do the transformation
        // there, then go back to barycentric coordinates.  Streamlining this is left for later.
        Dune::FieldVector<double,dim+1> barycentric;
        barycentric[0] = overlaps[i].localCoords[0][j][0];
        barycentric[1] = overlaps[i].localCoords[0][j][1];
        barycentric[2] = 1 - barycentric[0] - barycentric[1];

        Dune::FieldVector<double,dim> ref = barycentricToReference(barycentric);

        // this is the actual transformation.  (0,0) becomes (1,1), (1,0) becomes (0,1),  (0,1) becomes (1,0)
        ref *= -1;
        ref +=  1;

        // back to barycentric
        barycentric = referenceToBarycentric(ref);
        overlaps[i].localCoords[0][j][0] = barycentric[0];
        overlaps[i].localCoords[0][j][1] = barycentric[1];

      }

    }

    if (targetIsSecondTriangle[overlaps[i].tris[1]][0]) {

      assert(dim==2);

      // loop over the intersection corners
      for (size_t j=0; j<dim+1; j++) {

        // kinda umstaendlich: go from barycentric to reference coords, do the transformation
        // there, then go back to barycentric coordinates.  Streamlining this is left for later.
        Dune::FieldVector<double,dim+1> barycentric;
        barycentric[0] = overlaps[i].localCoords[1][j][0];
        barycentric[1] = overlaps[i].localCoords[1][j][1];
        barycentric[2] = 1 - barycentric[0] - barycentric[1];

        Dune::FieldVector<double,dim> ref = barycentricToReference(barycentric);

        // this is the actual transformation.  (0,0) becomes (1,1), (1,0) becomes (0,1),  (0,1) becomes (1,0)
        ref *= -1;
        ref +=  1;

        // back to barycentric
        barycentric = referenceToBarycentric(ref);
        overlaps[i].localCoords[1][j][0] = barycentric[0];
        overlaps[i].localCoords[1][j][1] = barycentric[1];

      }


    }

    // The numbers in overlaps[].tri refer to the split elements.
    // Replace that with the actual numbers
    overlaps[i].tris[0] = unsplitDomainElementNumbers[overlaps[i].tris[0]];
    overlaps[i].tris[1] = unsplitTargetElementNumbers[overlaps[i].tris[1]];

  }

  // //////////////////////////////////////////////
  // initialize the merged grid overlap manager
  // //////////////////////////////////////////////

  this->olm_.setOverlaps(overlaps);

  valid = true;
}


template<int dim, int dimworld, typename T>
typename PSurfaceMerge<dim, dimworld, T>::LocalCoords PSurfaceMerge<dim, dimworld, T>::grid1ParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& ip = this->olm_.domain(idx);

  // read the local coordinates from the overlap's struct,
  LocalCoords result;

  if (dim==1) {
    result[0] = ip.localCoords[0][corner][0];
  } else {

    // The psurface local coordinates are strange: The coordinates x[0], x[1] are the first two
    // components of the three barycentric coordinates.  The third one is implied to be 1-x[0]-x[1].
    // Maybe I should change that some day.  OS
    assert(dim==2);
    Dune::FieldVector<double,dim+1> barycentric;
    barycentric[0] = ip.localCoords[0][corner][0];
    barycentric[1] = ip.localCoords[0][corner][1];
    barycentric[2] = 1 - barycentric[0] - barycentric[1];

    result = barycentricToReference(barycentric);
  }

  return result;
}


template<int dim, int dimworld, typename T>
typename PSurfaceMerge<dim, dimworld, T>::LocalCoords PSurfaceMerge<dim, dimworld, T>::grid2ParentLocal(unsigned int idx, unsigned int corner) const
{
  // get the simplex overlap
  const PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype>& ip = this->olm_.domain(idx);

  // read the local coordinates from the overlap's struct,
  LocalCoords result;
  if (dim==1) {
    result[0] = ip.localCoords[1][corner][0];
  } else {

    // The psurface local coordinates are strange: The coordinates x[0], x[1] are the first two
    // components of the three barycentric coordinates.  The third one is implied to be 1-x[0]-x[1].
    // Maybe I should change that some day.  OS
    assert(dim==2);
    Dune::FieldVector<double,dim+1> barycentric;
    barycentric[0] = ip.localCoords[1][corner][0];
    barycentric[1] = ip.localCoords[1][corner][1];
    barycentric[2] = 1 - barycentric[0] - barycentric[1];

    result = barycentricToReference(barycentric);
  }
  return result;
}


/* IMPLEMENTATION OF   O V E R L A P  M A N A G E R  SUBCLASS */

template<int dim, int dimworld, typename T>
void PSurfaceMerge<dim, dimworld, T>::OverlapManager::setOverlaps(const std::vector<PSURFACE_NAMESPACE IntersectionPrimitive<dim,ctype> >& unordered)
{
  this->domOrder.clear();
  this->tarOrder.clear();
  this->domOrder.resize(unordered.size());
  this->tarOrder.resize(unordered.size(), NULL);

  // initialize the pointer arrays
  for (unsigned int i = 0; i < unordered.size(); ++i)
  {
    this->domOrder[i] = unordered[i];
  }

  // now sort to ascending domain parent index order
  std::sort(this->domOrder.begin(), this->domOrder.end(), &domainParentSmaller);
  // fill the 2nd array with the pointers to the elements from the 1st array
  for (unsigned int i = 0; i < this->domOrder.size(); ++i)
    this->tarOrder[i] = &this->domOrder[i];
  // now sort only the pointers in O(N*logN) comparisons
  std::sort(this->tarOrder.begin(), this->tarOrder.end(), &targetParentSmaller);

  // get the base pointer of the array
  this->baseptr = &this->domOrder[0];

}


template<int dim, int dimworld, typename T>
unsigned int PSurfaceMerge<dim, dimworld, T>::OverlapManager::firstDomainParent(unsigned int parent) const
{
  // perform a binary search for a simplex overlap with the given parent
  unsigned int first = 0, last = this->domOrder.size(), p = last/2;
  bool continuing = true;
  // search until either found or until not found with the search interval's size decreased to 1
  while (((p = this->domain((first + last)/2).tris[0]) != parent) && (continuing = (last > first + 1)))
  {
    if (p > parent)
      last = (first + last) / 2;
    else
      first = (first + last) / 2;
  }

  // if the parent is not found return a value too high
  if (!continuing)
    return this->domOrder.size();

  // else go back if the found element is not the first one
  // (note that p is now the INDEX of the overlap)
  p = (first + last) / 2;
  while (p > 0 && static_cast<unsigned int>(this->domain(p-1).tris[0]) == parent)
    p--;
  return p;
}


template<int dim, int dimworld, typename T>
unsigned int PSurfaceMerge<dim, dimworld, T>::OverlapManager::firstTargetParent(unsigned int parent) const
{
  // perform a binary search for a simplex overlap with the given parent
  unsigned int first = 0, last = this->domOrder.size(), p = last/2;
  bool continuing = true;
  // search until either found or until not found with the search interval's size decreased to 1
  while (((p = this->target((first + last)/2).tris[1]) != parent) && (continuing = (last > first + 1)))
  {
    if (p > parent)
      last = (first + last) / 2;
    else
      first = (first + last) / 2;
  }

  // if the parent is not found return a value too high
  if (!continuing)
    return this->domOrder.size();

  // else go back if the found element is not the first one
  // (note that p is now the INDEX of the overlap)
  p = (first + last) / 2;
  while (p > 0 && static_cast<unsigned int>(this->target(p-1).tris[1]) == parent)
    p--;
  return p;
}

// Explicit instantiation
#ifndef PSURFACE_EXTERN
template class PSurfaceMerge<1,1,double>;
template class PSurfaceMerge<1,2,double>;
template class PSurfaceMerge<2,2,double>;
template class PSurfaceMerge<2,3,double>;
#endif

#endif // HAVE_PSURFACE
