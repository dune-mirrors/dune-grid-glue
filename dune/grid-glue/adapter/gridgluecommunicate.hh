// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUECOMMUNICATE_HH
#define DUNE_GRIDGLUECOMMUNICATE_HH

/**@file
   @author Christian Engwer
   @brief Describes the parallel communication interface class for Dune::GridGlue
 */

#include <dune/common/bartonnackmanifcheck.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/istl/communicator.hh>

namespace Dune
{

  /**
     \brief describes the features of a data handle for
     communication in parallel runs using the GridGlue::communicate methods.
     Here the Barton-Nackman trick is used to interprete data handle objects
     as its interface. Therefore usable data handle classes need to be
     derived from this class.

     \tparam DataHandleImp implementation of the users data handle
     \tparam DataTypeImp type of data that are going to be communicated which is exported as \c DataType (for example double)
     \ingroup GICollectiveCommunication
   */
  template <class DataHandleImp, class DataTypeImp>
  class GridGlueCommDataHandleIF
  {
  public:
    //! data type of data to communicate
    typedef DataTypeImp DataType;

  protected:
    // one should not create an explicit instance of this inteface object
    GridGlueCommDataHandleIF() {}

  public:

    /*! how many objects of type DataType have to be sent for a given intersection
       Note: Both sender and receiver side need to know this size.
     */
    template<class RISType>
    size_t size (RISType& i) const
    {
      CHECK_INTERFACE_IMPLEMENTATION((asImp().size(i)));
      return asImp().size(i);
    }

    /** @brief pack data from user to message buffer
        @param buff message buffer provided by the grid
        @param e entity for which date should be packed to buffer
        @param i RemoteIntersection for which date should be packed to buffer
     */
    template<class MessageBufferImp, class EntityType, class RISType>
    void gather (MessageBufferImp& buff, const EntityType& e, const RISType & i) const
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().gather(buffIF,e,i)));
    }

    /*! unpack data from message buffer to user
       n is the number of objects sent by the sender
       @param buff message buffer provided by the grid
       @param e entity for which date should be unpacked from buffer
       @param n number of data written to buffer for this entity before
     */
    template<class MessageBufferImp, class EntityType, class RISType>
    void scatter (MessageBufferImp& buff, const EntityType& e, const RISType & i, size_t n)
    {
      MessageBufferIF<MessageBufferImp> buffIF(buff);
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION((asImp().scatter(buffIF,e,i,n)));
    }

  private:
    //!  Barton-Nackman trick
    DataHandleImp& asImp () {
      return static_cast<DataHandleImp &> (*this);
    }
    //!  Barton-Nackman trick
    const DataHandleImp& asImp () const
    {
      return static_cast<const DataHandleImp &>(*this);
    }
  }; // end class CommDataHandleIF

  /**
   * Streaming MessageBuffer for the GridGlue communication
   */
  template<typename DT>
  class GridGlueMessageBuffer {
  public:
    typedef DT value_type;

    // Constructor
    GridGlueMessageBuffer (DT *p)
    {
      a=p;
      i=0;
      j=0;
    }

    // write data to message buffer, acts like a stream !
    template<class Y>
    void write (const Y& data)
    {
      dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
      a[i++] = data;
    }

    // read data from message buffer, acts like a stream !
    template<class Y>
    void read (Y& data) const
    {
      dune_static_assert(( is_same<DT,Y>::value ), "DataType missmatch");
      data = a[j++];
    }

    size_t counter() const { return i; }

    void clear()
    {
      i = 0;
      j = 0;
    }

    // we need access to these variables in an assertion
#ifdef NDEBUG
  private:
#endif
    DT *a;
    size_t i;
    mutable size_t j;
  };

  //// CommPolicy
  template <typename GG, class DataHandleImp, class DataTypeImp>
  struct GridGlueCommInfo
  {
    typedef DataTypeImp value_type;
    typedef GG GridGlue;
    typedef DataTypeImp DataType;

    GridGlueCommInfo() : buffer(100), mbuffer(&buffer[0])
    {}

    // tunnel information to the policy and the operators
    const GridGlue * gridglue;
    GridGlueCommDataHandleIF<DataHandleImp, DataTypeImp> * data;

    // state variables
    std::vector<DataType> buffer;
    mutable GridGlueMessageBuffer<DataType> mbuffer;
    size_t currentsize;
  };

#if HAVE_MPI
  template<typename GG, class DataHandleImp, class DataTypeImp>
  struct CommPolicy< GridGlueCommInfo<GG, DataHandleImp, DataTypeImp> >
  {
    /**
     * @brief The type of the GridGlueCommInfo
     */
    typedef GridGlueCommInfo<GG, DataHandleImp, DataTypeImp> Type;

    /**
     * @brief The datatype that should be communicated.
     */
    typedef DataTypeImp IndexedType;

    /**
     * @brief Each intersection can communicate a different number of objects.
     */
    // typedef SizeOne IndexedTypeFlag;
    typedef VariableSize IndexedTypeFlag;

    /**
     * @brief Get the number of objects at an intersection.
     */
    static size_t getSize(const Type& commInfo, size_t i)
    {
      // get RemoteIntersection
      typedef typename Type::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->getIntersection(i));

      // ask data handle for size
      return commInfo.data->size(ris);
    }
  };
#endif

  // forward gather scatter to use defined class
  template<int dir>
  class GridGlueCommunicationOperator
  {
  public:
    template<class CommInfo>
    static const typename CommInfo::DataType& gather(const CommInfo& commInfo, size_t i, size_t j = 0)
    {
      // get RemoteIntersection
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->getIntersection(i));

      // fill buffer if we have a new intersection
      if (j == 0)
      {
        commInfo.mbuffer.clear();
        if (dir == Dune::ForwardCommunication)
        {
          // read from domain
          assert(ris.hasDomain());
          commInfo.data->gather(commInfo.mbuffer, ris.entityDomain(), ris);
        }
        else     // (dir == Dune::BackwardCommunication)
        {
          // read from target
          assert(ris.hasTarget());
          commInfo.data->gather(commInfo.mbuffer, ris.entityTarget(), ris);
        }
      }

      // return the j'th value in the buffer
      assert(j < commInfo.mbuffer.i);
      return commInfo.buffer[j];
    }

    template<class CommInfo>
    static void scatter(CommInfo& commInfo, const typename CommInfo::DataType& v, std::size_t i, std::size_t j = 0)
    {
      // extract GridGlue objects...
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->getIntersection(i));

      // get size if we have a new intersection
      if (j == 0)
      {
        commInfo.mbuffer.clear();
        commInfo.currentsize = commInfo.data->size(ris);
      }

      // write entry to buffer
      commInfo.buffer[j] = v;

      // write back the buffer if we are at the end of this intersection
      if (j == commInfo.currentsize-1)
      {
        if (dir == Dune::ForwardCommunication)
        {
          // write to target
          assert(ris.hasTarget());
          commInfo.data->scatter(commInfo.mbuffer, ris.entityTarget(), ris, commInfo.currentsize);
        }
        else     // (dir == Dune::BackwardCommunication)
        {
          // write to domain
          assert(ris.hasDomain());
          commInfo.data->scatter(commInfo.mbuffer, ris.entityDomain(), ris, commInfo.currentsize);
        }
        assert(commInfo.mbuffer.j <= commInfo.currentsize);
      }
    }
  };

  typedef GridGlueCommunicationOperator<Dune::ForwardCommunication> GridGlueForwardOperator;
  typedef GridGlueCommunicationOperator<Dune::BackwardCommunication> GridGlueBackwardOperator;

} // end namespace Dune
#endif
