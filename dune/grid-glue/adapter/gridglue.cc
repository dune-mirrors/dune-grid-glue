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


// template<typename P0, typename P1>
// int GridGlue<P0, P1>::domainEntityNextFace(const DomainElement& e, int index) const
// {
//     int first, count;
//     // first check if the element forms a part of the extracted surface
//     if (!patch0_.faceIndices(e, first, count))
//         return -1;

//     // check all mapped faces and accept the first one with number >=index
//     count += first;
//     while (first < count &&    (patch0_.indexInInside(first) < index || !merger_->template simplexMatched<0>(first)))
//         first++;
//     if (first == count)
//         return -1; // no more faces
//     else
//         return patch0_.indexInInside(first); // found, return the face's number
// }


// template<typename P0, typename P1>
// int GridGlue<P0, P1>::targetEntityNextFace(const TargetElement& e, int index) const
// {
//     int first, count;
//     // first check if the element forms a part of the extracted surface
//     if (!patch1_.faceIndices(e, first, count))
//         return -1;

//     // check all mapped faces and accept the first one with number >=index
//     count += first;
//     while (first < count && (patch1_.indexInInside(first) < index || !merger_->template simplexMatched<1>(first)))
//         first++;
//     if (first == count)
//         return -1; // no more faces
//     else
//         return patch1_.indexInInside(first); // found, return the face's number
// }


template<typename P0, typename P1>
template<int P>
bool
GridGlue<P0, P1>::getIntersectionIndices(const typename GridGlueView<P0,P1,P>::GridElement& e, std::vector<unsigned int> & indices) const
{
  // first check if the element has at least one extracted subEntity
  int p_first, p_cnt;
  bool hasExtractedSubEntity = patch<P>().faceIndices(e, p_first, p_cnt);
  if (!hasExtractedSubEntity)
    return false;

  // now accumulate all remote intersections of the element's faces
  indices.clear();
  assert(indices.size() == 0);
  std::vector<unsigned int> tmp;

  // iterate over all simplices to check if there is more than one simplex refining the face
  for (int p = 0; p < p_cnt; p++)
  {
    if (merger_->template simplexRefined<P>(p_first+p, tmp))
    {
      for (unsigned int i = 0; i < tmp.size(); ++i)
      {
        indices.push_back(tmp[i]);
      }
    }
  }

#ifndef NDEBUG
  std::cout << "INTERSECTIONS: " << p_first << "/" << p_cnt << std::endl;
  std::cout << "INTERSECTIONS:";
  for (unsigned int j = 0; j < indices.size(); j++)
  {
    int is = intersections_[indices[j]].index_;
    int idx = merger_->template parent<P>(is);
    typedef typename GridGlueView<P0,P1,P>::Patch::GridView::template Codim<0>::EntityPointer EPtr;
    assert(idx >= p_first);
    assert(idx < p_first+p_cnt);
    EPtr ep(e);
    EPtr ex = patch<P>().element(idx);
    assert(ep == ex);
    std::cout << "\t" << is;
  }
  std::cout << std::endl;
#endif

  // add end iterator
  indices.push_back(index__sz);

  return (indices.size() > 0);
}


template<typename P0, typename P1>
template<int I>
typename GridGlueView<P0,P1,I>::CellIntersectionIterator
GridGlue<P0, P1>::ibegin(const typename GridGlueView<P0,P1,I>::GridElement& e) const
{
  typedef typename GridGlueView<P0,P1,I>::CellIntersectionIterator CellIntersectionIterator;
  // first check if the element has at least one intersection
  std::vector<unsigned int> indices;
  bool hasExtractedSubEntity = getIntersectionIndices<I>(e, indices);
  if (!hasExtractedSubEntity)
    return CellIntersectionIterator(this);

  return CellIntersectionIterator(this, indices);
}


template<typename P0, typename P1>
template<int I>
typename GridGlueView<P0,P1,I>::CellIntersectionIterator
GridGlue<P0, P1>::iend(const typename GridGlueView<P0,P1,I>::GridElement& e) const
{
  typedef typename GridGlueView<P0,P1,I>::CellIntersectionIterator CellIntersectionIterator;
  return CellIntersectionIterator(this);
}
