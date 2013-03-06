// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

#include "intersection.hh"
#include <vector>
#include <iterator>
#include <dune/common/mpitraits.hh>

#include "../common/multivector.hh"

/** \todo Implement MPI Status check with exception handling */
#define CheckMPIStatus(A,B) {}

#if HAVE_MPI
namespace {
  template<typename T>
  struct MPITypeInfo {};

  template<>
  struct MPITypeInfo< int >
  {
    static const unsigned int size = 1;
    static inline MPI_Datatype getType()
    {
      return MPI_INT;
    }
    static const int tag = 1234560;
  };

  template<typename K, int N>
  struct MPITypeInfo< Dune::FieldVector<K,N> >
  {
    static const unsigned int size = N;
    static inline MPI_Datatype getType()
    {
      return Dune::MPITraits<K>::getType();
    }
    static const int tag = 1234561;
  };

  template<>
  struct MPITypeInfo< unsigned int >
  {
    static const unsigned int size = 1;
    static inline MPI_Datatype getType()
    {
      return MPI_UNSIGNED;
    }
    static const int tag = 1234562;
  };

  template<>
  struct MPITypeInfo< Dune::GeometryType >
  {
    static const unsigned int size = 1;
    static inline MPI_Datatype getType()
    {
      return Dune::MPITraits< Dune::GeometryType >::getType();
    }
    static const int tag = 1234563;
  };

  /**
     Send std::vector<T> in the ring

   * data is sent to rankright
   * from rankleft tmp is received and swapped with data

   */
  template<typename T>
  void MPI_SendVectorInRing(
    std::vector<T> & data,
    std::vector<T> & tmp,
    int leftsize,
    int rightrank,
    int leftrank,
    MPI_Comm comm
    )
  {
    // mpi status stuff
    int result = 0;
    MPI_Status status;
    typedef MPITypeInfo<T> Info;
    // alloc buffer
    unsigned int tmpsize = tmp.size();
    tmp.resize(leftsize);
    // send data
    int rank; MPI_Comm_rank(comm, &rank);
    // std::cout << rank << " send " << data.size() << " to " << rightrank << std::endl;
    // std::cout << rank << " recv " << tmp.size() << " from " << leftrank << std::endl;
    if (leftsize > 0 && data.size() > 0)
    {
      // send & receive
      result =
        MPI_Sendrecv(
          &(data[0]), Info::size*data.size(), Info::getType(), rightrank, Info::tag,
          &(tmp[0]),  Info::size*tmp.size(),  Info::getType(), leftrank,  Info::tag,
          comm, &status);
    }
    if (leftsize == 0 && data.size() > 0)
    {
      // send
      result =
        MPI_Send(
          &(data[0]), Info::size*data.size(), Info::getType(), rightrank, Info::tag,
          comm);
    }
    if (leftsize > 0 && data.size() == 0)
    {
      // receive
      result =
        MPI_Recv(
          &(tmp[0]),  Info::size*tmp.size(),  Info::getType(), leftrank,  Info::tag,
          comm, &status);
    }
    // check result
    CheckMPIStatus(result, status);
    // swap buffers
    data.swap(tmp);
    // resize tmp buffer
    tmp.resize(tmpsize);
  }

  /** \todo Please doc me! */
  struct MaxPatchSizes
  {
    int coords, entities;
  };

  /** \todo Please doc me! */
  struct PatchSizes
  {
    int patch0coords, patch0entities, patch1coords, patch1entities;
  };
}
#endif // HAVE_MPI

template<typename P0, typename P1>
#if HAVE_MPI
GridGlue<P0, P1>::GridGlue(const Grid0Patch& gp1, const Grid1Patch& gp2, Merger* merger, MPI_Comm m)
#else
GridGlue<P0, P1>::GridGlue(const Grid0Patch & gp1, const Grid1Patch & gp2, Merger* merger)
#endif // HAVE_MPI
  :
    patch0_(gp1), patch1_(gp2), merger_(merger)
#if HAVE_MPI
    , mpicomm(m)
#endif // HAVE_MPI
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

  int myrank = 0;
#if HAVE_MPI
  int commsize = 1;
  MPI_Comm_rank(mpicomm, &myrank);
  MPI_Comm_size(mpicomm, &commsize);

  // setup parallel indexset
  domain_is.beginResize();
  target_is.beginResize();
#endif // HAVE_MPI

  // merge local patches and add to intersection list
  mergePatches(patch0coords, patch0entities, patch0types, myrank,
               patch1coords, patch1entities, patch1types, myrank);

#if HAVE_MPI

  // status variables of communication
  int mpi_result;
  MPI_Status mpi_status;

  if (commsize > 1)
  {
    // communicate max patch size
    MaxPatchSizes maxPatchSizes;
    {
      MaxPatchSizes myMaxPatchSizes;
      myMaxPatchSizes.coords = std::max(patch0coords.size(), patch1coords.size());
      myMaxPatchSizes.entities = std::max(patch0entities.size(), patch1entities.size());
      mpi_result = MPI_Allreduce(&myMaxPatchSizes, &maxPatchSizes,
                                 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      CheckMPIStatus(mpi_result, 0);
    }

    /**
       \todo Use vector<struct> for message buffer and MultiVector to copy these
     */
    // allocate remote buffers (maxsize to avoid reallocation)
    std::vector<Dune::FieldVector<ctype, dimworld> > remotePatch0coords ( maxPatchSizes.coords );
    std::vector<unsigned int> remotePatch0entities ( maxPatchSizes.entities );
    std::vector<Dune::GeometryType> remotePatch0types ( maxPatchSizes.entities );
    std::vector<Dune::FieldVector<ctype,dimworld> > remotePatch1coords ( maxPatchSizes.coords );
    std::vector<unsigned int> remotePatch1entities ( maxPatchSizes.entities );
    std::vector<Dune::GeometryType> remotePatch1types ( maxPatchSizes.entities );

    // copy local patches to remote patch buffers
    remotePatch0coords.clear();
    std::copy(patch0coords.begin(), patch0coords.end(), std::back_inserter(remotePatch0coords));
    remotePatch0entities.clear();
    std::copy(patch0entities.begin(), patch0entities.end(), std::back_inserter(remotePatch0entities));
    remotePatch0types.clear();
    std::copy(patch0types.begin(), patch0types.end(), std::back_inserter(remotePatch0types));

    remotePatch1coords.clear();
    std::copy(patch1coords.begin(), patch1coords.end(), std::back_inserter(remotePatch1coords));
    remotePatch1entities.clear();
    std::copy(patch1entities.begin(), patch1entities.end(), std::back_inserter(remotePatch1entities));
    remotePatch1types.clear();
    std::copy(patch1types.begin(), patch1types.end(), std::back_inserter(remotePatch1types));

    // allocate tmp buffers (maxsize to avoid reallocation)
    std::vector<Dune::FieldVector<ctype, dimworld> > tmpPatchCoords ( maxPatchSizes.coords );
    std::vector<unsigned int> tmpPatchEntities ( maxPatchSizes.entities );
    std::vector<Dune::GeometryType> tmpPatchTypes ( maxPatchSizes.entities );

    // communicate patches in the ring
    for (int i=0; i<commsize; i++)
    {
      int remoterank = (myrank - i) % commsize;
      int rightrank  = (myrank + 1) % commsize;
      int leftrank   = (myrank - 1) % commsize;

      // communicate actual patch sizes
      PatchSizes patchSizes;

      {
        patchSizes.patch0coords = remotePatch0coords.size();
        patchSizes.patch0entities = remotePatch0entities.size();
        patchSizes.patch1coords = remotePatch1coords.size();
        patchSizes.patch1entities = remotePatch1entities.size();
        // send to right neighbor, receive from left neighbor
        mpi_result =
          MPI_Sendrecv_replace(
            &patchSizes, 4, MPI_INT,
            rightrank, MPITypeInfo<int>::tag,
            leftrank,  MPITypeInfo<int>::tag,
            mpicomm, &mpi_status);
        CheckMPIStatus(mpi_result, mpi_status);
      }

      /* send remote patch to right neighbor and receive from left neighbor */

      // patch0coords
      MPI_SendVectorInRing(
        remotePatch0coords, tmpPatchCoords, patchSizes.patch0coords,
        rightrank, leftrank, mpicomm);

      // patch0entities
      MPI_SendVectorInRing(
        remotePatch0entities, tmpPatchEntities, patchSizes.patch0entities,
        rightrank, leftrank, mpicomm);

      // patch0types
      MPI_SendVectorInRing(
        remotePatch0types, tmpPatchTypes, patchSizes.patch0entities,
        rightrank, leftrank, mpicomm);

      // patch1coords
      MPI_SendVectorInRing(
        remotePatch1coords, tmpPatchCoords, patchSizes.patch1coords,
        rightrank, leftrank, mpicomm);

      // patch1entities
      MPI_SendVectorInRing(
        remotePatch1entities, tmpPatchEntities, patchSizes.patch1entities,
        rightrank, leftrank, mpicomm);

      // patch1types
      MPI_SendVectorInRing(
        remotePatch1types, tmpPatchTypes, patchSizes.patch1entities,
        rightrank, leftrank, mpicomm);

      /* merging */
      // merge local & remote patches
      // domain_is and target_is are updated automatically
      mergePatches(patch0coords, patch0entities, patch0types, myrank,
                   remotePatch1coords, remotePatch1entities, remotePatch1types, remoterank);
      mergePatches(remotePatch0coords, remotePatch0entities, remotePatch0types, remoterank,
                   patch1coords, patch1entities, patch1types, myrank);
    }
  }

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

  // howto handle overlap etc?

  int myrank = 0;
#if HAVE_MPI
  MPI_Comm_rank(mpicomm, &myrank);
#endif // HAVE_MPI

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
  // update remote index sets
#if DUNE_VERSION_NEWER_REV(DUNE_COMMON,2,2,1)
  assert(Dune::RESIZE == domain_is.state());
  assert(Dune::RESIZE == target_is.state());
#endif
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
#endif // HAVE_MPI

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
