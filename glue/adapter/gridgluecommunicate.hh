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
  class GridGlueCommDataHandleIF : public CommDataHandleIF<DataHandleImp, DataTypeImp>
  {
  protected:
    // one should not create an explicit instance of this inteface object
    GridGlueCommDataHandleIF() {}

  public:

    /*! how many objects of type DataType have to be sent for a given intersection
       Note: Only the sender side needs to know this size.
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

} // end namespace Dune
#endif
