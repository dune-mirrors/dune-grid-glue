// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
 * @file
 * @brief Central component of the module implementing the coupling of two grids.
 * @author Gerrit Buse, Christian Engwer
 */


#ifndef DUNE_GRIDGLUE_GRIDGLUE_HH
#define DUNE_GRIDGLUE_GRIDGLUE_HH

#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>
#include <dune/common/shared_ptr.hh>

#include "adapter/gridgluecommunicate.hh"
#include <dune/grid-glue/merging/merger.hh>

#if DUNE_VERSION_NEWER_REV(DUNE_COMMON,2,3,0)
#include <dune/common/parallel/mpitraits.hh>
#include <dune/common/parallel/mpicollectivecommunication.hh>
#else
#include <dune/common/mpitraits.hh>
#include <dune/common/mpicollectivecommunication.hh>
#endif
#include <dune/common/parallel/indexset.hh>
#include <dune/common/parallel/plocalindex.hh>
#include <dune/common/parallel/remoteindices.hh>
#include <dune/common/parallel/communicator.hh>
#include <dune/common/parallel/interface.hh>

namespace Dune {
namespace GridGlue {

/** Document the relation between the old grid names and the new numbering */
enum GridOrdering {
  Domain = 0,
  Target = 1
};

// forward declarations
template<typename P0, typename P1>
class GridGlue;

template<typename P0, typename P1>
class IntersectionData;

template<typename P0, typename P1, int inside, int outside>
class Intersection;

template<typename P0, typename P1, int inside, int outside>
class IntersectionIterator;

template<typename P0, typename P1>
class IntersectionIndexSet;

template<typename P0, typename P1, int P>
struct GridGlueView;

template<typename P0, typename P1>
struct GridGlueView<P0,P1,0>
{
  typedef P0 Patch;
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> IntersectionIterator;
  typedef typename Patch::GridView::template Codim<0>::Entity GridElement;
  static const P0& patch(const GridGlue<P0,P1>& g)
  {
    return *g.patch0_;
  }
};

template<typename P0, typename P1>
struct GridGlueView<P0,P1,1>
{
  typedef P1 Patch;
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,1,0> IntersectionIterator;
  typedef typename Patch::GridView::template Codim<0>::Entity GridElement;
  static const P1& patch(const GridGlue<P0,P1>& g)
  {
    return *g.patch1_;
  }
};

/**
 * @class GridGlue
 * @brief sequential adapter to couple two grids at specified close together boundaries
 *
 * @tparam P0 patch (extractor) to use for grid 0
 * @tparam P1 patch (extractor) to use for grid 1
 *
 * \todo adapt member names according to style guide
 */
template<typename P0, typename P1>
class GridGlue
{
private:

  /*   F R I E N D S   */

  friend class IntersectionData<P0,P1>;
  friend class Intersection<P0,P1,0,1>;
  friend class Intersection<P0,P1,1,0>;
  friend class IntersectionIterator<P0,P1,0,1>;
  friend class IntersectionIterator<P0,P1,1,0>;
  friend class IntersectionIndexSet<P0,P1>;
  friend struct GridGlueView<P0,P1,0>;
  friend struct GridGlueView<P0,P1,1>;

  /*   P R I V A T E   T Y P E S   */

  /** \brief GlobalId type of an intersection (used for communication) */
  typedef ::Dune::GridGlue::GlobalId GlobalId;

  /** \brief LocalIndex type of an intersection (used for communication) */
  typedef Dune::ParallelLocalIndex <Dune::PartitionType> LocalIndex;

  /** \brief ParallelIndexSet type (used for communication communication) */
  typedef Dune::ParallelIndexSet <GlobalId, LocalIndex> PIndexSet;

public:

  /*   P U B L I C   T Y P E S   A N D   C O N S T A N T S   */

  /** \brief GridView of grid 0 (aka domain grid) */
  typedef typename P0::GridView Grid0View;

  /** \brief Grid 0 type */
  typedef typename Grid0View::Grid Grid0;

  /** \brief Grid 0 type
      \deprecated please use typede Grid0
   */
  typedef Grid0 DomainGridType DUNE_DEPRECATED;

  /** \brief Coupling patch of grid 0 */
  typedef P0 Grid0Patch;

  /** \brief Dimension of the grid 0 extractor */
  enum {
    /** \brief Dimension of the grid 0 extractor */
    grid0dim = Grid0Patch::dim,
    domdim = Grid0Patch::dim,
    /** \brief World dimension of the grid 0 extractor */
    grid0dimworld = Grid0Patch::dimworld,
    domdimworld = Grid0Patch::dimworld
  };

  /** \brief GridView of grid 1 (aka target grid) */
  typedef typename P1::GridView Grid1View;

  /** \brief Grid 1 type */
  typedef typename Grid1View::Grid Grid1;

  /** \brief Grid 1 type
      \deprecated please use typede Grid0
   */
  typedef Grid1 TargetGridType DUNE_DEPRECATED;

  /** \brief Coupling patch of grid 1 */
  typedef P1 Grid1Patch;

  /** \todo */
  typedef unsigned int IndexType;

  /** \brief Dimension of the grid 1 extractor */
  enum {
    /** \brief Dimension of the grid 1 extractor */
    tardim = Grid1Patch::dim,
    grid1dim = Grid1Patch::dim,
    /** \brief World dimension of the grid 1 extractor */
    tardimworld = Grid1Patch::dimworld,
    grid1dimworld = Grid1Patch::dimworld
  };


  /** \brief export the world dimension */
  enum {
    /** \brief export the world dimension :
        maximum of the two extractor world dimensions */
    dimworld = ((int)Grid0Patch::dimworld > (int)Grid1Patch::dimworld) ? (int)Grid0Patch::dimworld : (int)Grid1Patch::dimworld

  };

  /** \brief The type used for coordinates
      \todo maybe use traits class to decide which has more precision (Grid0View::ctype or Grid1View::ctype) and then take this one
   */
  typedef typename Grid0View::ctype ctype;

  /** \brief The type used for coordinate vectors */
  typedef Dune::FieldVector<ctype, dimworld>                   Coords;

  /** \brief The type of the Grid0 elements */
  typedef typename Grid0View::Traits::template Codim<0>::Entity Grid0Element;

  /** \brief The type of the Grid0 elements
      \deprecated please use Grid0Element
   */
  typedef typename Grid0View::Traits::template Codim<0>::Entity DomainElement DUNE_DEPRECATED;

  /** \brief Pointer type to Grid0 elements */
  typedef typename Grid0View::Traits::template Codim<0>::EntityPointer Grid0ElementPtr;

  /** \brief Pointer type to Grid0 elements
      \deprecated please use Grid0ElementPtr
   */
  typedef typename Grid0View::Traits::template Codim<0>::EntityPointer DomainElementPtr DUNE_DEPRECATED;

  /** \brief The type of the Grid0 vertices */
  typedef typename Grid0View::Traits::template Codim<Grid0::dimension>::Entity Grid0Vertex;

  /** \brief The type of the Grid0 vertices
      \deprecated please use Grid0ElementPtr
   */
  typedef typename Grid0View::Traits::template Codim<Grid0::dimension>::Entity DomainVertex DUNE_DEPRECATED;

  /** \brief Pointer type to Grid0 vertices */
  typedef typename Grid0View::Traits::template Codim<Grid0::dimension>::EntityPointer Grid0VertexPtr;

  /** \brief Pointer type to Grid0 vertices
      \deprecated please use Grid0VertexPtr
   */
  typedef typename Grid0View::Traits::template Codim<Grid0::dimension>::EntityPointer DomainVertexPtr DUNE_DEPRECATED;

  /** \brief The type of the Grid1 elements */
  typedef typename Grid1View::Traits::template Codim<0>::Entity Grid1Element;

  /** \brief The type of the Grid1 elements
      \deprecated please use Grid1Element
   */
  typedef typename Grid1View::Traits::template Codim<0>::Entity TargetElement DUNE_DEPRECATED;

  /** \brief Pointer type to Grid1 elements */
  typedef typename Grid1View::Traits::template Codim<0>::EntityPointer Grid1ElementPtr;

  /** \brief Pointer type to Grid1 elements
      \deprecated please use Grid1ElementPtr
   */
  typedef typename Grid1View::Traits::template Codim<0>::EntityPointer TargetElementPtr DUNE_DEPRECATED;

  /** \brief The type of the Grid1 vertices */
  typedef typename Grid1View::Traits::template Codim<Grid1::dimension>::Entity Grid1Vertex;

  /** \brief The type of the Grid1 vertices
      \deprecated please use Grid1Vertex
   */
  typedef typename Grid1View::Traits::template Codim<Grid1::dimension>::Entity TargetVertex DUNE_DEPRECATED;

  /** \brief Pointer type to Grid1 vertices */
  typedef typename Grid1View::Traits::template Codim<Grid1::dimension>::EntityPointer Grid1VertexPtr;

  /** \brief Pointer type to Grid1 vertices
      \deprecated please use Grid1VertexPtr
   */
  typedef typename Grid1View::Traits::template Codim<Grid1::dimension>::EntityPointer TargetVertexPtr DUNE_DEPRECATED;

  /** \brief Instance of a Merger */
  typedef ::Merger<ctype,
      Grid0::dimension - Grid0Patch::codim,
      Grid1::dimension - Grid1Patch::codim,
      dimworld> Merger;

  /** \brief Type of remote intersection objects */
  typedef Dune::GridGlue::Intersection<P0,P1,0,1> Intersection;

  /** \brief Type of remote intersection indexSet */
  typedef Dune::GridGlue::IntersectionIndexSet<P0,P1> IndexSet;

  /** \brief Type of the iterator that iterates over remove intersections */

  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> Grid0IntersectionIterator;
  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,1,0> Grid1IntersectionIterator;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief the patch0 description
  const std::shared_ptr<const Grid0Patch> patch0_;

  /// @brief the patch1 description
  const std::shared_ptr<const Grid1Patch> patch1_;

  /// @brief the surface merging utility
  const std::shared_ptr<Merger> merger_;

  /// @brief number of intersections
  IndexType index__sz;

#if HAVE_MPI
  /// @brief MPI_Comm which this GridGlue is working on
  MPI_Comm mpicomm_;

  /// @brief parallel indexSet for the intersections with a local grid0 entity
  PIndexSet patch0_is_;

  /// @brief parallel indexSet for the intersections with a local grid1 entity
  PIndexSet patch1_is_;

  /// @brief keeps information about which process has which intersection
  Dune::RemoteIndices<PIndexSet> remoteIndices_;
#endif // HAVE_MPI

  /// \todo
  typedef Dune::GridGlue::IntersectionData<P0,P1> IntersectionData;

  /// @brief a vector with intersection elements
  mutable std::vector<IntersectionData>   intersections_;

protected:

  /**
   * @brief after building the merged grid the intersection can be updated
   * through this method (for internal use)
   *
   * @param patch0coords the patch0 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param patch0entities array with all patch0 entities represented as corner indices into @c patch0coords;
   * the entities are just written to this array one after another
   * @param patch0types array with all patch0 entities types
   * @param patch0rank  rank of the process where patch0 was extracted
   *
   * @param patch1coords the patch2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param patch1entities just like with the patch0entities and patch0corners
   * @param patch1types array with all patch1 entities types
   * @param patch1rank  rank of the process where patch1 was extracted
   *
   */
  void mergePatches(const std::vector<Dune::FieldVector<ctype,dimworld> >& patch0coords,
                    const std::vector<unsigned int>& patch0entities,
                    const std::vector<Dune::GeometryType>& patch0types,
                    const int patch0rank,
                    const std::vector<Dune::FieldVector<ctype,dimworld> >& patch1coords,
                    const std::vector<unsigned int>& patch1entities,
                    const std::vector<Dune::GeometryType>& patch1types,
                    const int patch1rank);


  template<typename Extractor>
  void extractGrid (const Extractor & extractor,
                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                    std::vector<unsigned int> & faces,
                    std::vector<Dune::GeometryType>& geometryTypes) const;

public:

  /*   C O N S T R U C T O R S   A N D   D E S T R U C T O R S   */

  /**
   * @brief constructor
   *
   * Initializes components but does not "glue" the surfaces. The surfaces
   * are extracted from the grids here though.
   * @param gp0 the grid0 patch
   * @param gp1 the grid1 patch
   * @param merger The merger object that is used to compute the merged grid. This class has
   * to be a model of the SurfaceMergeConcept.
   */
  GridGlue(const Grid0Patch& gp0, const Grid1Patch& gp1, Merger* merger);
  GridGlue(const std::shared_ptr<const Grid0Patch> gp0, const std::shared_ptr<const Grid1Patch> gp1, const std::shared_ptr<Merger> merger);

  /*   G E T T E R S   */

  /** \todo Please doc me! */
  template<int P>
  const typename GridGlueView<P0,P1,P>::Patch & patch() const
  {
    return GridGlueView<P0,P1,P>::patch(*this);
  }

  /**
   * @brief getter for the GridView of patch P
   * @return the object
   */
  template<int P>
  const typename GridGlueView<P0,P1,P>::Patch::GridView & gridView() const
  {
    return GridGlueView<P0,P1,P>::patch(*this).gridView();
  }


  /*   F U N C T I O N A L I T Y   */

  void build();

  /*   I N T E R S E C T I O N S   A N D   I N T E R S E C T I O N   I T E R A T O R S   */

  /**
   * @brief gets an iterator over all remote intersections in the merged grid between grid0 and grid1
   * @tparam I select inside grid I=0 or I=1
   *
   * @return the iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::IntersectionIterator ibegin() const
  {
    return typename GridGlueView<P0,P1,I>::IntersectionIterator(this, 0);
  }


  /**
   * @brief gets the (general) end-iterator for grid glue iterations
   * @tparam I select inside grid I=0 or I=1
   *
   * @return the iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::IntersectionIterator iend() const
  {
    return typename GridGlueView<P0,P1,I>::IntersectionIterator(this, index__sz);
  }


  /*! \brief Communicate information on the MergedGrid of a GridGlue

     Template parameter is a model of Dune::GridGlue::CommDataHandle

     \param data GridGlueDataHandle
     \param iftype Interface for which the Communication should take place
     \param dir Communication direction (Forward means grid0 to grid1, Backward is the reverse)

     \todo fix mixed communication: seq->par use commSeq, par->seq use commPar
     \todo add directed version communicate<FROM,TO, DH,DT>(data,iftype,dir)
   */
  template<class DataHandleImp, class DataTypeImp>
  void communicate (Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> & data,
                    Dune::InterfaceType iftype, Dune::CommunicationDirection dir) const
  {
    typedef Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> DataHandle;
    typedef typename DataHandle::DataType DataType;

#if HAVE_MPI

    if (mpicomm_ != MPI_COMM_SELF)
    {
      /*
       * P A R A L L E L   V E R S I O N
       */
      // setup communication interfaces
      Dune::dinfo << "GridGlue: parallel communication" << std::endl;
      typedef Dune::EnumItem <Dune::PartitionType, Dune::InteriorEntity> InteriorFlags;
      typedef Dune::EnumItem <Dune::PartitionType, Dune::OverlapEntity>  OverlapFlags;
      typedef Dune::EnumRange <Dune::PartitionType, Dune::InteriorEntity, Dune::GhostEntity>  AllFlags;
      Dune::Interface interface;
      assert(remoteIndices_.isSynced());
      switch (iftype)
      {
      case Dune::InteriorBorder_InteriorBorder_Interface :
        interface.build (remoteIndices_, InteriorFlags(), InteriorFlags() );
        break;
      case Dune::InteriorBorder_All_Interface :
        if (dir == Dune::ForwardCommunication)
          interface.build (remoteIndices_, InteriorFlags(), AllFlags() );
        else
          interface.build (remoteIndices_, AllFlags(), InteriorFlags() );
        break;
      case Dune::Overlap_OverlapFront_Interface :
        interface.build (remoteIndices_, OverlapFlags(), OverlapFlags() );
        break;
      case Dune::Overlap_All_Interface :
        if (dir == Dune::ForwardCommunication)
          interface.build (remoteIndices_, OverlapFlags(), AllFlags() );
        else
          interface.build (remoteIndices_, AllFlags(), OverlapFlags() );
        break;
      case Dune::All_All_Interface :
        interface.build (remoteIndices_, AllFlags(), AllFlags() );
        break;
      default :
        DUNE_THROW(Dune::NotImplemented, "GridGlue::communicate for interface " << iftype << " not implemented");
      }

      // setup communication info (class needed to tunnel all info to the operator)
      typedef Dune::GridGlue::CommInfo<GridGlue,DataHandleImp,DataTypeImp> CommInfo;
      CommInfo commInfo;
      commInfo.dir = dir;
      commInfo.gridglue = this;
      commInfo.data = &data;

      // create communicator
      Dune::BufferedCommunicator bComm ;
      bComm.template build< CommInfo >(commInfo, commInfo, interface);

      // do communication
      // choose communication direction.
      if (dir == Dune::ForwardCommunication)
        bComm.forward< Dune::GridGlue::ForwardOperator >(commInfo, commInfo);
      else
        bComm.backward< Dune::GridGlue::BackwardOperator >(commInfo, commInfo);
    }
    else
#endif // HAVE_MPI
    {
      /*
       * S E Q U E N T I A L   V E R S I O N
       */
      Dune::dinfo << "GridGlue: sequential fallback communication" << std::endl;

      // get comm buffer size
      int ssz = size() * 10;       // times data per intersection
      int rsz = size() * 10;

      // allocate send/receive buffer
      DataType* sendbuffer = new DataType[ssz];
      DataType* receivebuffer = new DataType[rsz];

      // iterators
      Grid0IntersectionIterator rit = ibegin<0>();
      Grid0IntersectionIterator ritend = iend<0>();

      // gather
      Dune::GridGlue::StreamingMessageBuffer<DataType> gatherbuffer(sendbuffer);
      for (; rit != ritend; ++rit)
      {
        /*
           we need to have to variants depending on the communication direction.
         */
        if (dir == Dune::ForwardCommunication)
        {
          /*
             dir : Forward (grid0 -> grid1)
           */
          if (rit->self())
          {
            data.gather(gatherbuffer, rit->inside(), *rit);
          }
        }
        else         // (dir == Dune::BackwardCommunication)
        {
          /*
             dir : Backward (grid1 -> grid0)
           */
          if (rit->neighbor())
          {
            data.gather(gatherbuffer, rit->outside(), rit->flip());
          }
        }
      }

      assert(ssz == rsz);
      for (int i=0; i<ssz; i++)
        receivebuffer[i] = sendbuffer[i];

      // scatter
      Dune::GridGlue::StreamingMessageBuffer<DataType> scatterbuffer(receivebuffer);
      for (rit = ibegin<0>(); rit != ritend; ++rit)
      {
        /*
           we need to have to variants depending on the communication direction.
         */
        if (dir == Dune::ForwardCommunication)
        {
          /*
             dir : Forward (grid0 -> grid1)
           */
          if (rit->neighbor())
            data.scatter(scatterbuffer, rit->outside(), rit->flip(),
                         data.size(*rit));
        }
        else         // (dir == Dune::BackwardCommunication)
        {
          /*
             dir : Backward (grid1 -> grid0)
           */
          if (rit->self())
            data.scatter(scatterbuffer, rit->inside(), *rit,
                         data.size(*rit));
        }
      }

      // cleanup pointers
      delete[] sendbuffer;
      delete[] receivebuffer;
    }
  }

  /*
   * @brief return an IndexSet mapping from Intersection to IndexType
   */
  IndexSet indexSet() const
  {
    return IndexSet(this);
  }

#if QUICKHACK_INDEX
  // indexset size
  size_t indexSet_size() const
  {
    return index__sz;
  }

#endif

  Intersection getIntersection(int i) const
  {
    return Intersection(this, & intersections_[i]);
  }

  size_t size() const
  {
    return index__sz;
  }

};

} // end namespace GridGlue
} // end namespace Dune

#include "adapter/gridglue.cc"

#include "adapter/intersection.hh"
#include "adapter/intersectioniterator.hh"
#include "adapter/intersectionindexset.hh"

#include "adapter/rangegenerators.hh"

#endif // DUNE_GRIDGLUE_GRIDGLUE_HH
