// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

#include "intersection.hh"
#include <vector>
#include <iterator>
#include "gridglue.hh"

#include "../common/multivector.hh"
#include "../extractors/vtksurfacewriter.hh"

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

    MPI_Barrier(comm);
  }

  /** \brief struct to simplify communication of the patch data sizes */
  struct PatchSizes
  {
    PatchSizes() :
      patch0coords(0), patch0entities(0), patch0types(0),
      patch1coords(0), patch1entities(0), patch1types(0) {}

    //! initialize patch sizes
    PatchSizes(unsigned int c0, unsigned int e0, unsigned int t0,
               unsigned int c1, unsigned int e1, unsigned int t1) :
      patch0coords(c0), patch0entities(e0), patch0types(t0),
      patch1coords(c1), patch1entities(e1), patch1types(t1) {}

    //! initialize patch sizes using the data containers
    template<typename C, typename E, typename T>
    PatchSizes(const C & c0, const E &  e0, const T & t0,
               const C & c1, const E & e1, const T & t1) :
      patch0coords(c0.size()), patch0entities(e0.size()), patch0types(t0.size()),
      patch1coords(c1.size()), patch1entities(e1.size()), patch1types(t1.size()) {}

    unsigned int patch0coords, patch0entities, patch0types,
                 patch1coords, patch1entities, patch1types;

    unsigned int maxCoords() const { return std::max(patch0coords, patch1coords); }
    unsigned int maxEntities() const { return std::max(patch0entities, patch1entities); }
    unsigned int maxTypes() const { return std::max(patch0types, patch1types); }
  };
}
#endif // HAVE_MPI

template<typename P0, typename P1>
GridGlue<P0, P1>::GridGlue(const Grid0Patch& gp0, const Grid1Patch& gp1, Merger* merger) :
  patch0_(gp0), patch1_(gp1), merger_(merger)
{
#if HAVE_MPI
  // if we have only seq. meshes don't use parallel glueing
  if (gp0.gridView().comm().size() == 1
      && gp1.gridView().comm().size() == 1)
    mpicomm_ = MPI_COMM_SELF;
  else
    mpicomm_ = MPI_COMM_WORLD;
#endif // HAVE_MPI
  std::cout << "GridGlue: Constructor succeeded!" << std::endl;
}

template<typename P0, typename P1>
void GridGlue<P0, P1>::build()
{
  int myrank = 0;
#if HAVE_MPI
  int commsize = 1;
  MPI_Comm_rank(mpicomm_, &myrank);
  MPI_Comm_size(mpicomm_, &commsize);
#endif // HAVE_MPI

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

  std::cout << ">>>> rank " << myrank << " coords: "
            << patch0coords.size() << " and " << patch1coords.size() << std::endl;
  std::cout << ">>>> rank " << myrank << " entities: "
            << patch0entities.size() << " and " << patch1entities.size() << std::endl;
  std::cout << ">>>> rank " << myrank << " types: "
            << patch0types.size() << " and " << patch1types.size() << std::endl;

#ifdef WRITE_TO_VTK
  const char prefix[] = "GridGlue::Builder::build() : ";
  char patch0surf[256];
  sprintf(patch0surf, "/tmp/vtk-patch0-test-%i", myrank);
  char patch1surf[256];
  sprintf(patch1surf, "/tmp/vtk-patch1-test-%i", myrank);

  std::cout << prefix << "Writing patch0 surface to '" << patch0surf << ".vtk'...\n";
  VtkSurfaceWriter vtksw(patch0surf);
  vtksw.writeSurface(patch0coords, patch0entities, grid0dim, dimworld);
  std::cout << prefix << "Done writing patch0 surface!\n";

  std::cout << prefix << "Writing patch1 surface to '" << patch1surf << ".vtk'...\n";
  vtksw.setFilename(patch1surf);
  vtksw.writeSurface(patch1coords, patch1entities, grid1dim, dimworld);
  std::cout << prefix << "Done writing patch1 surface!\n";
#endif // WRITE_TO_VTK

#if HAVE_MPI
  if (commsize > 1)
  {
    // setup parallel indexset
    domain_is_.beginResize();
    target_is_.beginResize();
  }
#endif // HAVE_MPI

  // merge local patches and add to intersection list
  if (patch0entities.size() > 0 && patch1entities.size() > 0)
    mergePatches(patch0coords, patch0entities, patch0types, myrank,
                 patch1coords, patch1entities, patch1types, myrank);

#ifdef CALL_MERGER_TWICE
  if (patch0entities.size() > 0 && patch1entities.size() > 0)
    mergePatches(patch0coords, patch0entities, patch0types, myrank,
                 patch1coords, patch1entities, patch1types, myrank);
#endif

#if HAVE_MPI

  // status variables of communication
  int mpi_result;
  MPI_Status mpi_status;

  if (commsize > 1)
  {
    // get patch sizes
    PatchSizes patchSizes (patch0coords, patch0entities, patch0types,
                           patch1coords, patch1entities, patch1types);

    // communicate max patch size
    PatchSizes maxPatchSizes;
    mpi_result = MPI_Allreduce(&patchSizes, &maxPatchSizes,
                               6, MPI_UNSIGNED, MPI_MAX, MPI_COMM_WORLD);
    CheckMPIStatus(mpi_result, 0);

    /**
       \todo Use vector<struct> for message buffer and MultiVector to copy these
     */
    // allocate remote buffers (maxsize to avoid reallocation)
    std::vector<Dune::FieldVector<ctype, dimworld> > remotePatch0coords ( maxPatchSizes.patch0coords );
    std::vector<unsigned int> remotePatch0entities ( maxPatchSizes.patch0entities );
    std::vector<Dune::GeometryType> remotePatch0types ( maxPatchSizes.patch0types );
    std::vector<Dune::FieldVector<ctype,dimworld> > remotePatch1coords ( maxPatchSizes.patch1coords );
    std::vector<unsigned int> remotePatch1entities ( maxPatchSizes.patch1entities );
    std::vector<Dune::GeometryType> remotePatch1types ( maxPatchSizes.patch1types );

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
    std::vector<Dune::FieldVector<ctype, dimworld> > tmpPatchCoords ( maxPatchSizes.maxCoords() );
    std::vector<unsigned int> tmpPatchEntities ( maxPatchSizes.maxEntities() );
    std::vector<Dune::GeometryType> tmpPatchTypes ( maxPatchSizes.maxTypes() );

    // communicate patches in the ring
    for (int i=1; i<commsize; i++)
    {
      int remoterank = (myrank - i + commsize) % commsize;
      int rightrank  = (myrank + 1 + commsize) % commsize;
      int leftrank   = (myrank - 1 + commsize) % commsize;

      // communicate current patch sizes
      // patchsizes were initialized before
      {
        // send to right neighbor, receive from left neighbor
        mpi_result =
          MPI_Sendrecv_replace(
            &patchSizes, 6, MPI_UNSIGNED,
            rightrank, MPITypeInfo<unsigned int>::tag,
            leftrank,  MPITypeInfo<unsigned int>::tag,
            mpicomm_, &mpi_status);
        CheckMPIStatus(mpi_result, mpi_status);
      }

      /* send remote patch to right neighbor and receive from left neighbor */

      // patch0coords
      // std::cout << myrank << " patch0coords" << std::endl;
      MPI_SendVectorInRing(
        remotePatch0coords, tmpPatchCoords, patchSizes.patch0coords,
        rightrank, leftrank, mpicomm_);

      // patch0entities
      // std::cout << myrank << " patch0entities" << std::endl;
      MPI_SendVectorInRing(
        remotePatch0entities, tmpPatchEntities, patchSizes.patch0entities,
        rightrank, leftrank, mpicomm_);

      // patch0types
      // std::cout << myrank << " patch0types" << std::endl;
      MPI_SendVectorInRing(
        remotePatch0types, tmpPatchTypes, patchSizes.patch0types,
        rightrank, leftrank, mpicomm_);

      // patch1coords
      // std::cout << myrank << " patch1coords" << std::endl;
      MPI_SendVectorInRing(
        remotePatch1coords, tmpPatchCoords, patchSizes.patch1coords,
        rightrank, leftrank, mpicomm_);

      // patch1entities
      // std::cout << myrank << " patch1entities" << std::endl;
      MPI_SendVectorInRing(
        remotePatch1entities, tmpPatchEntities, patchSizes.patch1entities,
        rightrank, leftrank, mpicomm_);

      // patch1types
      // std::cout << myrank << " patch1types" << std::endl;
      MPI_SendVectorInRing(
        remotePatch1types, tmpPatchTypes, patchSizes.patch1types,
        rightrank, leftrank, mpicomm_);

      /* merging */
      // merge local & remote patches
      // domain_is_ and target__is are updated automatically
      if (remotePatch1entities.size() > 0 && patch0entities.size() > 0)
        mergePatches(patch0coords, patch0entities, patch0types, myrank,
                     remotePatch1coords, remotePatch1entities, remotePatch1types, remoterank);
      if (remotePatch0entities.size() > 0 && patch1entities.size() > 0)
        mergePatches(remotePatch0coords, remotePatch0entities, remotePatch0types, remoterank,
                     patch1coords, patch1entities, patch1types, myrank);

      std::cout << "Sync processes" << std::endl;
      MPI_Barrier(mpicomm_);
      std::cout << "...done" << std::endl;
    }
  }

  if (commsize > 1)
  {
    // finalize ParallelIndexSet & RemoteIndices
    domain_is_.endResize();
    target_is_.endResize();

    // setup remote index information
    remoteIndices_.setIncludeSelf(false);
    remoteIndices_.setIndexSets(domain_is_, target_is_, mpicomm_) ;
    remoteIndices_.rebuild<true>();
  }
#endif

}

template<typename T>
void printVector(const std::vector<T> & v, std::string name)
{
  std::cout << name << std::endl;
  for (size_t i=0; i<v.size(); i++)
  {
    std::cout << v[i] << "   ";
  }
  std::cout << std::endl;
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
  int commsize = 1;
#if HAVE_MPI
  MPI_Comm_rank(mpicomm_, &myrank);
  MPI_Comm_size(mpicomm_, &commsize);
#endif // HAVE_MPI

  // which patches are local?
  const bool patch0local = (myrank == patch0rank);
  const bool patch1local = (myrank == patch1rank);

  // remember the number of previous remote intersections
  const unsigned int offset = intersections_.size();

  std::cout << myrank
            << " GridGlue::mergePatches : rank " << patch0rank << " / " << patch1rank << std::endl;

  // start the actual build process
  merger_->build(patch0coords, patch0entities, patch0types,
                 patch1coords, patch1entities, patch1types);

  // append to intersections list
  intersections_.resize(merger_->nSimplices() + offset + 1);
  for (unsigned int i = 0; i < merger_->nSimplices(); ++i)
  {
    // currently we only support local merging!
    IntersectionData data(*this, i, offset, patch0local, patch1local);
    intersections_[offset+i] = data;
  }

  index__sz = intersections_.size() - 1;

  std::cout << myrank
            << " GridGlue::mergePatches : "
            << "The number of remote intersections is " << intersections_.size()-1 << std::endl;

  // printVector(patch0coords,"patch0coords");
  // printVector(patch0entities,"patch0entities");
  // printVector(patch0types,"patch0types");
  // printVector(patch1coords,"patch1coords");
  // printVector(patch1entities,"patch1entities");
  // printVector(patch1types,"patch1types");

#if HAVE_MPI
  if (commsize > 1)
  {
    // update remote index sets
#if DUNE_VERSION_NEWER_REV(DUNE_COMMON,2,2,1)
    assert(Dune::RESIZE == domain_is_.state());
    assert(Dune::RESIZE == target_is_.state());
#endif
    for (unsigned int i = 0; i < merger_->nSimplices(); i++)
    {
#warning only handle the newest intersections / merger info
      const IntersectionData & it = intersections_[i];
      GlobalId gid;
      gid.first.first = patch0rank;
      gid.first.second = patch1rank;
      gid.second = offset+i;
      if (it.grid0local_)
      {
        Dune::PartitionType ptype = patch0_.element(it.grid0index_)->partitionType();
        domain_is_.add (gid, LocalIndex(offset+i, ptype) );
      }
      if (it.grid1local_)
      {
        Dune::PartitionType ptype = patch1_.element(it.grid1index_)->partitionType();
        target_is_.add (gid, LocalIndex(offset+i, ptype) );
      }
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
