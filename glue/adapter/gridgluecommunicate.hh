// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUECOMMUNICATE_HH
#define DUNE_GRIDGLUECOMMUNICATE_HH

/**@file
   @author Christian Engwer
   @brief Describes the parallel communication interface class for Dune::GridGlue
 */

#include <dune/grid/common/datahandleif.hh>
#include <dune/common/bartonnackmanifcheck.hh>

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

  private:
    DT *a;
    size_t i;
    mutable int j;
  };

  //// CommPolicy
  template <typename GG, class DataHandleImp, class DataTypeImp>
  struct GridGlueCommInfo
  {
    typedef DataTypeImp value_type;
    typedef GG GridGlue;
    typedef DataTypeImp DataType;

    GridGlueCommInfo() : buffer(100), mbuffer(&buffer[0]) {}

    // tunnel information to the policy and the operators
    const GridGlue * gridglue;
    GridGlueCommDataHandleIF<DataHandleImp, DataTypeImp> * data;

    // state variables
    std::vector<DataType> buffer;
    mutable GridGlueMessageBuffer<DataType> mbuffer;
    size_t currentsize;
  };

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
    typedef SizeOne IndexedTypeFlag;
    //typedef VariableSize IndexedTypeFlag;

    /**
     * @brief Get the number of objects at an intersection.
     */
    static int getSize(const Type& commInfo, size_t i)
    {
      assert(false);
      // get RemoteIntersection
      typedef typename Type::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->_intersections[i]);

      assert(commInfo.data->size(ris) == 1);
      // ask data handle for size
      return commInfo.data->size(ris);
    }
  };

  // forward gather scatter to use defined class
  class GridGlueForwardOperator
  {
  public:
    template<class CommInfo>
    static const typename CommInfo::value_type& gather(const CommInfo& commInfo, std::size_t i)
    {
      typedef typename CommInfo::value_type DataType;
      typedef Dune::GridGlueMessageBuffer<DataType> MessageBuffer;
      static DataType buffer[10];
      GridGlueMessageBuffer<DataType> mbuffer(buffer);

      // extract GridGlue objects...
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->_intersections[i]);

      // read from domain
      assert(ris.hasDomain());
      commInfo.data->gather(mbuffer, ris.entityDomain(), ris);

      // return _the_ value
      return buffer[0];
    }

    template<class CommInfo, class DataType>
    static void scatter(CommInfo& commInfo, const DataType& v, std::size_t i)
    {
      typedef typename CommInfo::value_type DataType;
      typedef Dune::GridGlueMessageBuffer<DataType> MessageBuffer;
      DataType buffer[10];
      GridGlueMessageBuffer<DataType> mbuffer(buffer);

      // extract GridGlue objects...
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->_intersections[i]);

      // fill buffer
      mbuffer.write(v);

      // write to target
      assert(ris.hasTarget());
      commInfo.data->scatter(mbuffer, ris.entityTarget(), ris, 1);
    }
  };

  class GridGlueBackwardOperator
  {
  public:
    template<class CommInfo>
    static const typename CommInfo::value_type& gather(const CommInfo& commInfo, std::size_t i)
    {
      typedef typename CommInfo::value_type DataType;
      typedef Dune::GridGlueMessageBuffer<DataType> MessageBuffer;
      static DataType buffer[10];
      GridGlueMessageBuffer<DataType> mbuffer(buffer);

      // extract GridGlue objects...
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->_intersections[i]);

      // read from target
      assert(ris.hasTarget());
      commInfo.data->gather(mbuffer, ris.entityTarget(), ris);

      // return _the_ value
      return buffer[0];
    }

    template<class CommInfo, class DataType>
    static void scatter(CommInfo& commInfo, const DataType& v, std::size_t i)
    {
      typedef typename CommInfo::value_type DataType;
      typedef Dune::GridGlueMessageBuffer<DataType> MessageBuffer;
      DataType buffer[10];
      GridGlueMessageBuffer<DataType> mbuffer(buffer);

      // extract GridGlue objects...
      typedef typename CommInfo::GridGlue::RemoteIntersection RemoteIntersection;
      RemoteIntersection ris(commInfo.gridglue->_intersections[i]);

      // fill buffer
      mbuffer.write(v);

      // write to domain
      assert(ris.hasDomain());
      commInfo.data->scatter(mbuffer, ris.entityDomain(), ris, 1);
    }
  };

} // end namespace Dune
#endif
