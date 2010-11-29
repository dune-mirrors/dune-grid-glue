// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

template<typename P0, typename P1>
#if HAVE_MPI
GridGlue<P0, P1>::GridGlue(const Grid0Patch& gp1, const Grid1Patch& gp2, Merger* merger, MPI_Comm m)
#else
GridGlue<P0, P1>::GridGlue(const Grid0Patch & gp1, const Grid1Patch & gp2, Merger* merger)
#endif
  : domgv_(gp1.gridView()), targv_(gp2.gridView()),
    patch0_(gp1), patch1_(gp2), merger_(merger)
#if HAVE_MPI
    , mpicomm(m)
#endif
{
  std::cout << "GridGlue: Constructor succeeded!" << std::endl;
}

template<typename P0, typename P1>
void GridGlue<P0, P1>::build()
{
  // clear the contents from the current intersections array
  {
    std::vector<IntersectionData> dummy;
    intersections_.swap(dummy);
  }

  std::vector<Dune::FieldVector<ctype, dimworld> > domcoords;
  std::vector<unsigned int> domfaces;
  std::vector<Dune::GeometryType> domainElementTypes;
  std::vector<Dune::FieldVector<ctype,dimworld> > tarcoords;
  std::vector<unsigned int> tarfaces;
  std::vector<Dune::GeometryType> targetElementTypes;

  /*
   * extract global surface grids
   */

  // retrieve the coordinate and topology information from the extractors
  // and apply transformations if necessary
  extractGrid(patch0_, domcoords, domfaces, domainElementTypes);
  extractGrid(patch1_, tarcoords, tarfaces, targetElementTypes);

#ifdef WRITE_TO_VTK
  const int dimw = Parent::dimworld;
  const char prefix[] = "GridGlue::Builder::build() : ";
  char domainsurf[256];
  sprintf(domainsurf, "/tmp/vtk-domain-test-%i", mpi_rank);
  char targetsurf[256];
  sprintf(targetsurf, "/tmp/vtk-target-test-%i", mpi_rank);

  std::cout << prefix << "Writing domain surface to '" << domainsurf << ".vtk'...\n";
  VtkSurfaceWriter vtksw(domainsurf);
  vtksw.writeSurface(domcoords, domfaces, dimw, dimw);
  std::cout << prefix << "Done writing domain surface!\n";

  std::cout << prefix << "Writing target surface to '" << targetsurf << ".vtk'...\n";
  vtksw.setFilename(targetsurf);
  vtksw.writeSurface(tarcoords, tarfaces, dimw, dimw);
  std::cout << prefix << "Done writing target surface!\n";
#endif // WRITE_TO_VTK


  // start the actual build process
  merger_->build(domcoords, domfaces, domainElementTypes,
                 tarcoords, tarfaces, targetElementTypes);

  // the intersections need to be recomputed
  updateIntersections();
}


template<typename P0, typename P1>
template<typename Extractor>
void GridGlue<P0, P1>::extractGrid (const Extractor & extractor,
                                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                                    std::vector<unsigned int> & faces,
                                    std::vector<Dune::GeometryType>& geometryTypes) const
{
  std::vector<typename Extractor::Coords> tempcoords;
  std::vector<typename Extractor::VertexVector> tempfaces;

  extractor.getCoords(tempcoords);
  coords.clear();
  coords.reserve(tempcoords.size());

  for (unsigned int i = 0; i < tempcoords.size(); ++i)
  {
    assert(int(dimworld) == int(Extractor::dimworld));
    coords.push_back(Dune::FieldVector<ctype, dimworld>());
    for (size_t j = 0; j <dimworld; ++j)
      coords.back()[j] = tempcoords[i][j];
  }

  extractor.getFaces(tempfaces);
  faces.clear();

  for (unsigned int i = 0; i < tempfaces.size(); ++i) {
    for (unsigned int j = 0; j < tempfaces[i].size(); ++j)
      faces.push_back(tempfaces[i][j]);
  }

  // get the list of geometry types from the extractor
  extractor.getGeometryTypes(geometryTypes);

}


template<typename P0, typename P1>
int GridGlue<P0, P1>::domainEntityNextFace(const DomainElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!patch0_.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count &&    (patch0_.indexInInside(first) < index || !merger_->template simplexMatched<0>(first)))
    first++;
  if (first == count)
    return -1;     // no more faces
  else
    return patch0_.indexInInside(first);     // found, return the face's number
}


template<typename P0, typename P1>
int GridGlue<P0, P1>::targetEntityNextFace(const TargetElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!patch1_.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (patch1_.indexInInside(first) < index || !merger_->template simplexMatched<1>(first)))
    first++;
  if (first == count)
    return -1;     // no more faces
  else
    return patch1_.indexInInside(first);     // found, return the face's number
}


template<typename P0, typename P1>
typename GridGlue<P0, P1>::IntersectionIterator
GridGlue<P0, P1>::ibegin() const
{
  return IntersectionIterator(this, 0);
}


template<typename P0, typename P1>
typename GridGlue<P0, P1>::Grid0IntersectionIterator
GridGlue<P0, P1>::idomainbegin(const DomainElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = patch0_.faceIndices(e, first, count);
  if (!in_surface)
                      #warning Not Implemented
    return Grid0IntersectionIterator(this, index__sz);

  count += first;
  while (first < count)
  {
    if (patch0_.indexInInside(first) == num && merger_->template simplexMatched<0>(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      merger_->template simplexRefined<0>(first, global_results);
      while (++first < count && patch0_.indexInInside(first) == num && merger_->template simplexRefined<0>(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
                                            #warning Not Implemented
        assert(false);
      //                 return Grid0IntersectionIterator(intersections_[global_results[0]], global_results);

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
    #warning Not Implemented
  return Grid0IntersectionIterator(this, index__sz);
}


template<typename P0, typename P1>
typename GridGlue<P0, P1>::Grid0IntersectionIterator
GridGlue<P0, P1>::idomainbegin(const DomainElement& e) const
{
  // first check if the element has at least one extracted subEntity
  int first, count;
  bool hasExtractedSubEntity = patch0_.faceIndices(e, first, count);
  if (!hasExtractedSubEntity)
                                 #warning Not Implemented
    return Grid0IntersectionIterator(this, index__sz);


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplex refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (merger_->template simplexRefined<0>(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
                    #warning Not Implemented
    assert(false);
  //        return Grid0IntersectionIterator(this, intersections_[global_results[0]], global_results);
  else
          #warning Not Implemented
    return Grid0IntersectionIterator(this, index__sz);
}


template<typename P0, typename P1>
typename GridGlue<P0, P1>::Grid1IntersectionIterator
GridGlue<P0, P1>::itargetbegin(const TargetElement& e, int num) const
{
  // first check if the element has at least one extracted subEntity
  int first, count;
  bool hasExtractedSubEntity = patch1_.faceIndices(e, first, count);
  if (!hasExtractedSubEntity) return itargetend();

  count += first;
  while (first < count)
  {
    if (patch1_.indexInInside(first) == num && merger_->template simplexMatched<1>(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      merger_->template simplexRefined<1>(first, global_results);
      while (++first < count && patch1_.indexInInside(first) == num && merger_->template simplexRefined<1>(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
                                            #warning Not Implemented
        assert(false);
      // return Grid1IntersectionIterator(this, intersections_[global_results[0]], global_results);

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
    #warning Not Implemented
  return Grid1IntersectionIterator(this, index__sz);
}


template<typename P0, typename P1>
typename GridGlue<P0, P1>::Grid1IntersectionIterator
GridGlue<P0, P1>::itargetbegin(const TargetElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = patch1_.faceIndices(e, first, count);
  if (!in_surface) return itargetend();


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (merger_->template simplexRefined<1>(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
                    #warning Not Implemented
    assert(false);
  // return Grid1IntersectionIterator(this, intersections_[global_results[0]], global_results);
  else
          #warning Not Implemented
    return Grid1IntersectionIterator(this, index__sz);
}
