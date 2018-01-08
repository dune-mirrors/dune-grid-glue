// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

/** \todo Implement MPI Status check with exception handling */
#define CheckMPIStatus(A,B) {}

#if HAVE_MPI
#include <mpi.h>
#include <functional>

#include <dune/common/hybridutilities.hh>
#include <dune/common/std/utility.hh>
#include <dune/common/std/apply.hh>

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

  // template<typename K, int N>
  // struct MPITypeInfo< Dune::FieldVector<K,N> >
  // {
  //   static const unsigned int size = N;
  //   static inline MPI_Datatype getType()
  //   {
  //     return Dune::MPITraits<K>::getType();
  //   }
  //   static const int tag = 1234561;
  // };

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

  // template<>
  // struct MPITypeInfo< Dune::GeometryType >
  // {
  //   static const unsigned int size = 1;
  //   static inline MPI_Datatype getType()
  //   {
  //     return Dune::MPITraits< Dune::GeometryType >::getType();
  //   }
  //   static const int tag = 1234563;
  // };

  template<typename T>
  void MPI_SetVectorSize(
    std::vector<T> & data,
    MPI_Status & status)
  {
    typedef MPITypeInfo<T> Info;
    int sz;
    MPI_Get_count(&status, Info::getType(), &sz);
    data.resize(sz);
  }

  /**
     Send std::vector<T> in the ring

     \return a pair of MPI_Request, as this operation is asynchroneous
     * data is sent to rankright
     * from rankleft tmp is received and swapped with data

   */
  template<typename T>
  void MPI_SendVectorInRing(
    std::vector<T> & data,
    std::vector<T> & next,
    int tag,
    int rightrank,
    int leftrank,
    MPI_Comm comm,
    MPI_Request& r_send,
    MPI_Request& r_recv
    )
  {
    // mpi status stuff
    int result DUNE_UNUSED;
    result = 0;
    typedef MPITypeInfo<T> Info;
    // resize next buffer to maximum size
    next.resize(next.capacity());
    // send data (explicitly send data.size elements)
    result =
      MPI_Isend(
        &(data[0]), Info::size*data.size(), Info::getType(), rightrank, tag,
        comm, &r_send);
    // receive up to maximum size. The acutal size is stored in the status
    result =
      MPI_Irecv(
        &(next[0]),  Info::size*next.size(),  Info::getType(), leftrank,  tag,
        comm, &r_recv);
    // // check result
    // MPI_Status status;
    // CheckMPIStatus(result, status);
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
  };
}
#endif // HAVE_MPI

template<typename OP, std::size_t... Indices, typename... Args>
void MPI_AllApply_impl(MPI_Comm mpicomm,
  OP && op,
  std::index_sequence<Indices...> indices,
  const Args&... data)
{
  constexpr std::size_t N = sizeof...(Args);
  int myrank = 0;
  int commsize = 0;
#if HAVE_MPI
  MPI_Comm_rank(mpicomm, &myrank);
  MPI_Comm_size(mpicomm, &commsize);
  // status variables of communication
  int mpi_result;
#endif // HAVE_MPI

  if (commsize > 1)
  {
#ifdef DEBUG_GRIDGLUE_PARALLELMERGE
    std::cout << myrank << " Start Communication, size " << commsize << std::endl;
#endif

    // get data sizes
    std::array<unsigned int, N> size({ (unsigned int)data.size()... });

    // communicate max data size
    std::array<unsigned int, N> maxSize;
    mpi_result = MPI_Allreduce(&size, &maxSize,
      size.size(), MPI_UNSIGNED, MPI_MAX, mpicomm);
    CheckMPIStatus(mpi_result, 0);
#ifdef DEBUG_GRIDGLUE_PARALLELMERGE
    std::cout << myrank << " maxSize " << "done... " << std::endl;
#endif

    // allocate receiving buffers with maxsize to ensure sufficient buffer size for communication
    std::tuple<Args...> remotedata { maxSize[Indices]... };
    {
      int dummy[sizeof...(Indices)] =
        { (std::cout << myrank << ": size " << std::get<Indices>(remotedata).size() << std::endl, 0)... };
    }

    // copy local data to receiving buffer
    {
      [](std::initializer_list<auto>){}(
        { (std::get<Indices>(remotedata) = data)... });
    }

    // allocate second set of receiving buffers necessary for async communication
    std::tuple<Args...> nextdata { maxSize[Indices]... };

    // communicate data in the ring
    int rightrank  = (myrank + 1 + commsize) % commsize;
    int leftrank   = (myrank - 1 + commsize) % commsize;

    std::cout << myrank << ": size = " << commsize << std::endl;
    std::cout << myrank << ": left = " << leftrank
              << " right = " << rightrank << std::endl;

    // currently the remote data is our own data
    int remoterank = myrank;

    for (int i=1; i<commsize; i++)
    {
      // in this iteration we will receive data from nextrank
      int nextrank = (myrank - i + commsize) % commsize;

      std::cout << myrank << ": next = " << nextrank << std::endl;

      // send remote data to right neighbor and receive from left neighbor
      std::array<MPI_Request,2*N> requests;

      int tag = 12345678;
      Dune::Hybrid::forEach(indices,
        [&](auto i){
          MPI_SendVectorInRing(
            std::get<i>(remotedata),
            std::get<i>(nextdata),
            tag+i,
            rightrank, leftrank, mpicomm,
            requests[2*i],
            requests[2*i+1]);
        });

      // apply operator
      op(remoterank,std::get<Indices>(remotedata)...);

      // wait for communication to finalize
      std::array<MPI_Status,2*N> status;
      MPI_Waitall(2*N,&requests[0],&status[0]);

      // we finished receiving from nextrank and thus remoterank = nextrank
      remoterank = nextrank;

      // get current data sizes
      // and resize vectors
      Dune::Hybrid::forEach(indices,
        [&](auto i){
          MPI_SetVectorSize(std::get<i>(nextdata),status[2*i+1]);
        });

      // swap the communication buffers
      std::swap(remotedata,nextdata);
    }

    // last apply (or the only one in the case of sequential application)
    op(remoterank,std::get<Indices>(remotedata)...);
  }
  else // sequential
  {
    op(myrank,data...);
  }
}

template<typename OP, typename... Args>
void MPI_AllApply(MPI_Comm mpicomm,
  OP && op,
  const Args& ... data)
{
  MPI_AllApply_impl(
    mpicomm,
    std::forward<OP>(op),
    std::make_index_sequence<sizeof...(Args)>(),
    data...
    );
}
