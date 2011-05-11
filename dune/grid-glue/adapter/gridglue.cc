// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

#include "intersection.hh"

template<typename P0, typename P1>
#if HAVE_MPI
GridGlue<P0, P1>::GridGlue(const Grid0Patch& gp1, const Grid1Patch& gp2, Merger* merger, MPI_Comm m)
#else
GridGlue<P0, P1>::GridGlue(const Grid0Patch & gp1, const Grid1Patch & gp2, Merger* merger)
#endif
  :
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

  std::vector<Dune::FieldVector<ctype, dimworld> > patch0coords;
  std::vector<unsigned int> patch0entities;
  std::vector<Dune::GeometryType> patch0types;
  std::vector<Dune::FieldVector<ctype,dimworld> > patch1coords;
  std::vector<unsigned int> patch1entities;
  std::vector<Dune::GeometryType> patch1types;

  /*
   * extract global surface patchs
   */

  // retrieve the coordinate and topology information from the extractors
  // and apply transformations if necessary
  extractGrid(patch0_, patch0coords, patch0entities, patch0types);
  extractGrid(patch1_, patch1coords, patch1entities, patch1types);

#ifdef WRITE_TO_VTK
  const int dimw = Parent::dimworld;
  const char prefix[] = "GridGlue::Builder::build() : ";
  char patch0surf[256];
  sprintf(patch0surf, "/tmp/vtk-patch0-test-%i", mpicomm.rank());
  char patch1surf[256];
  sprintf(patch1surf, "/tmp/vtk-patch1-test-%i", mpicomm.rank());

  std::cout << prefix << "Writing patch0 surface to '" << patch0surf << ".vtk'...\n";
  VtkSurfaceWriter vtksw(patch0surf);
  vtksw.writeSurface(patch0coords, patch0entities, dimw, dimw);
  std::cout << prefix << "Done writing patch0 surface!\n";

  std::cout << prefix << "Writing patch1 surface to '" << patch1surf << ".vtk'...\n";
  vtksw.setFilename(patch1surf);
  vtksw.writeSurface(patch1coords, patch1entities, dimw, dimw);
  std::cout << prefix << "Done writing patch1 surface!\n";
#endif // WRITE_TO_VTK

  // clear old intersection list
  intersections_.clear();

#if HAVE_MPI
  // setup parallel indexset
  domain_is.beginResize();
  target_is.beginResize();
#endif

  int myrank = 0;
  int commsize = 1;
#if HAVE_MPI
  MPI_Comm_rank(mpicomm, &myrank);
  MPI_Comm_size(mpicomm, &commsize);
#endif

  // merge local patches and add to intersection list
  mergePatches(patch0coords, patch0entities, patch0types, myrank,
               patch1coords, patch1entities, patch1types, myrank);

  if (commsize > 1)
  {
#warning implement ring comm of patches
    DUNE_THROW(Dune::Exception, "TODO: implement ring comm of patches");

    // communicate max patch size
    int patchCoords_maxSize;
    int patchEntites_maxSize;

    // allocate remote buffers (maxsize to avoid reallocation)
    std::vector<Dune::FieldVector<ctype, dimworld> > remotePatch0coords ( patchCoords_maxSize );
    std::vector<unsigned int> remotePatch0entities ( patchEntites_maxSize );
    std::vector<Dune::GeometryType> remotePatch0types ( patchEntites_maxSize );
    std::vector<Dune::FieldVector<ctype,dimworld> > remotePatch1coords ( patchCoords_maxSize );
    std::vector<unsigned int> remotePatch1entities ( patchEntites_maxSize );
    std::vector<Dune::GeometryType> remotePatch1types ( patchEntites_maxSize );

    // copy local patches to remote patch buffers
#warning todo

    // allocate tmp buffers (maxsize to avoid reallocation)
    std::vector<Dune::FieldVector<ctype, dimworld> > tmpPatchCoords ( patchCoords_maxSize );
    std::vector<unsigned int> tmpPatchEntities ( patchEntites_maxSize );
    std::vector<Dune::GeometryType> tmpPatchTypes ( patchEntites_maxSize );

    // communicate patches in the ring
    for (int i=0; i<commsize; i++)
    {
      int remoterank = (myrank - i) % commsize;
      int rightrank  = (myrank + 1) % commsize;
      int leftrank   = (myrank - 1) % commsize;

      // communicate actual patch sizes
      int patch0coordsSize;
      int patch0entitiesSize;
      int patch1coordsSize;
      int patch1entitiesSize;

#warning todo
      // int MPI_Send( void *buf, int count, MPI_Datatype datatype, int dest,
      //     int tag, MPI_Comm comm );

      /* patch0 */
      tmpPatchCoords.resize(patch0coordsSize);
      tmpPatchEntities.resize(patch0entitiesSize);
      tmpPatchTypes.resize(patch0entitiesSize);

      // send remote patch to right neighbor
#warning todo

      // receive remote patch from left neighbor
#warning todo

      // finish communication
#warning todo

      // swap communcation buffers
      remotePatch0coords.swap(tmpPatchCoords);
      remotePatch0entities.swap(tmpPatchEntities);
      remotePatch0types.swap(tmpPatchTypes);

      /* patch1 */
      tmpPatchCoords.resize(patch1coordsSize);
      tmpPatchEntities.resize(patch1entitiesSize);
      tmpPatchTypes.resize(patch1entitiesSize);

      // send remote patch to right neighbor
#warning todo

      // receive remote patch from left neighbor
#warning todo

      // finish communication
#warning todo

      // swap communcation buffers
      remotePatch1coords.swap(tmpPatchCoords);
      remotePatch1entities.swap(tmpPatchEntities);
      remotePatch1types.swap(tmpPatchTypes);

      /* merging */
      // merge local & remote patches
      mergePatches(patch0coords, patch0entities, patch0types, myrank,
                   remotePatch1coords, remotePatch1entities, remotePatch1types, remoterank);
      mergePatches(remotePatch0coords, remotePatch0entities, remotePatch0types, remoterank,
                   patch1coords, patch1entities, patch1types, myrank);
    }
  }

#if HAVE_MPI
  ////// finalize ParallelIndexSet & RemoteIndices
  domain_is.endResize();
  target_is.endResize();

  // setup remote index information
  remoteIndices.setIncludeSelf(false);
  remoteIndices.setIndexSets(domain_is, target_is, mpicomm) ;
  remoteIndices.rebuild<true>();
#endif

}

template<typename P0, typename P1>
void GridGlue<P0, P1>::mergePatches(
  const std::vector<Dune::FieldVector<ctype,dimworld> >& patch0coords,
  const std::vector<unsigned int>& patch0entities,
  const std::vector<Dune::GeometryType>& patch0types,
  const int patch0rank,
  const std::vector<Dune::FieldVector<ctype,dimworld> >& patch1coords,
  const std::vector<unsigned int>& patch1entities,
  const std::vector<Dune::GeometryType>& patch1types,
  const int patch1rank)
{

#warning howto handle overlap etc?

  int myrank = 0;
#if HAVE_MPI
  MPI_Comm_rank(mpicomm, &myrank);
#endif

  // which patches are local?
  const bool patch0local = (myrank == patch0rank);
  const bool patch1local = (myrank == patch1rank);

  // remember the number of previous remote intersections
  const unsigned int offset = intersections_.size();

  // start the actual build process
  merger_->build(patch0coords, patch0entities, patch0types,
                 patch1coords, patch1entities, patch1types);

  // append to intersections list
  intersections_.resize(merger_->nSimplices() + offset + 1);
  for (unsigned int i = 0; i < merger_->nSimplices(); ++i)
  {
    // currently we only support local merging!
    IntersectionData data(*this, offset+i, patch0local, patch1local);
    intersections_[offset+i] = data;
  }

  index__sz = intersections_.size() - 1;

  std::cout
  << "GridGlue::mergePatches : "
  << "The number of remote intersections is " << intersections_.size()-1 << std::endl;

#if HAVE_MPI
  for (unsigned int i = 0; i < merger_->nSimplices(); i++)
  {
    const IntersectionData & it = intersections_[i];
    GlobalId gid;
    gid.first.first = patch0rank;
    gid.first.second = patch1rank;
    gid.second = offset+i;
    if (it.grid0local_)
    {
      Dune::PartitionType ptype = patch0_.element(it.grid0index_)->partitionType();
      domain_is.add (gid, LocalIndex(offset+i, ptype) );
    }
    if (it.grid1local_)
    {
      Dune::PartitionType ptype = patch1_.element(it.grid1index_)->partitionType();
      target_is.add (gid, LocalIndex(offset+i, ptype) );
    }
  }
#endif

  // cleanup the merger
  merger_->clear();
}

template<typename P0, typename P1>
template<typename Extractor>
void GridGlue<P0, P1>::extractGrid (const Extractor & extractor,
                                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                                    std::vector<unsigned int> & entities,
                                    std::vector<Dune::GeometryType>& geometryTypes) const
{
  std::vector<typename Extractor::Coords> tempcoords;
  std::vector<typename Extractor::VertexVector> tempentities;

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

  extractor.getFaces(tempentities);
  entities.clear();

  for (unsigned int i = 0; i < tempentities.size(); ++i) {
    for (unsigned int j = 0; j < tempentities[i].size(); ++j)
      entities.push_back(tempentities[i][j]);
  }

  // get the list of geometry types from the extractor
  extractor.getGeometryTypes(geometryTypes);

}
