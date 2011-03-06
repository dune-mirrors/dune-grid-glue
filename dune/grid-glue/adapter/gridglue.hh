// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlue.hh
 *  Version:     1.0
 *  Created on:  Feb 2, 2009
 *  Author:      Gerrit Buse, Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Central component of the module implementing the coupling of two grids.
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Central component of the module implementing the coupling of two grids.
 */


#ifndef GRIDGLUE_HH_
#define GRIDGLUE_HH_

#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/iteratorfacades.hh>

#define QUICKHACK_INDEX 1

#include "gridgluecommunicate.hh"
#include <dune/grid-glue/merging/merger.hh>

#define DUNE_COMMON_VERSION_NUMBER (DUNE_COMMON_VERSION_MAJOR * 10 + DUNE_COMMON_VERSION_MINOR)
#if DUNE_COMMON_VERSION_NUMBER > 20
  #include <dune/common/parallel/indexset.hh>
  #include <dune/common/parallel/plocalindex.hh>
  #include <dune/common/parallel/remoteindices.hh>
  #include <dune/common/parallel/communicator.hh>
  #include <dune/common/parallel/interface.hh>
#else
  #include <dune/istl/indexset.hh>
  #include <dune/istl/plocalindex.hh>
  #include <dune/istl/remoteindices.hh>
  #include <dune/istl/communicator.hh>
  #include <dune/istl/interface.hh>
#endif

/** Document the relation between the old grid names and the new numbering */
enum GridOrdering {
  Domain = 0,
  Target = 1
};

// forward declaration
template<typename P0, typename P1>
class GridGlue;

namespace Dune {
  namespace GridGlue {

    template<typename P0, typename P1>
    class IntersectionData;

    template<typename P0, typename P1, int inside, int outside>
    class Intersection;

    template<typename P0, typename P1, int inside, int outside>
    class IntersectionIterator;

    template<typename P0, typename P1, int inside, int outside>
    class CellIntersectionIterator;

  }
}

template<typename P0, typename P1, int P>
struct GridGlueView;

template<typename P0, typename P1>
struct GridGlueView<P0,P1,0>
{
  typedef P0 Patch;
  typedef Dune::GridGlue::CellIntersectionIterator<P0,P1,0,1> CellIntersectionIterator;
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> IntersectionIterator;
  typedef typename Patch::GridView::template Codim<0>::Entity GridElement;
  static const P0& patch(const GridGlue<P0,P1>& g)
  {
    return g.patch0_;
  }
};

template<typename P0, typename P1>
struct GridGlueView<P0,P1,1>
{
  typedef P1 Patch;
  typedef Dune::GridGlue::CellIntersectionIterator<P0,P1,1,0> CellIntersectionIterator;
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,1,0> IntersectionIterator;
  typedef typename Patch::GridView::template Codim<0>::Entity GridElement;
  static const P1& patch(const GridGlue<P0,P1>& g)
  {
    return g.patch1_;
  }
};

/**
 * @class GridGlue
 * @brief sequential adapter to couple two grids at specified close together boundaries
 *
 *
 * @tparam GET1 a first traits class to determine the type of surface extractor to use for the domain grid,
 * has to be a model of GridExtractionTraitsConcept
 * @tparam GET2 a second traits class to determine the type of surface extractor to use for the target grid,
 * has to be a model of GridExtractionTraitsConcept
 * @tparam SM the surface merging class, has to be a model of the SurfaceMergeConcept
 */
template<typename P0, typename P1>
class GridGlue
{
private:

  /*   P R I V A T E   T Y P E S   */

  /** \brief GlobalId type of an intersection (used for communication) */
  typedef unsigned int GlobalId;

  /** \brief LocalIndex type of an intersection (used for communication) */
  typedef Dune::ParallelLocalIndex<Dune::PartitionType> LocalIndex;

  /** \brief ParallelIndexSet type (used for communication communication) */
  typedef Dune::ParallelIndexSet <GlobalId, LocalIndex > PIndexSet;

public:

  /*   P U B L I C   T Y P E S   A N D   C O N S T A N T S   */

  /** \brief Grid view of the domain grid */
  typedef typename P0::GridView Grid0View;

  /** \brief Domain grid type */
  typedef typename Grid0View::Grid DomainGridType;

  /** \brief Coupling patch of grid 0 */
  typedef P0 Grid0Patch;

  /** \brief Dimension of the domain extractor */
  enum {
    /** \brief Dimension of the domain extractor */
    domdim = Grid0Patch::dim,
    /** \brief WOrld dimension of the domain extractor */
    domdimworld = Grid0Patch::dimworld
  };

  /** \brief Grid view of the target grid */
  typedef typename P1::GridView Grid1View;

  /** \brief Target grid type */
  typedef typename Grid1View::Grid TargetGridType;

  /** \brief Coupling patch of grid 1 */
  typedef P1 Grid1Patch;

  /** \todo */
  typedef unsigned int IndexType;

  /** \brief Dimension of the target extractor */
  enum {
    /** \brief Dimension of the target extractor */
    tardim = Grid1Patch::dim,
    /** \brief World dimension of the target extractor */
    tardimworld = Grid1Patch::dimworld
  };


  /** \brief export the world dimension */
  enum {
    /** \brief export the world dimension :
        maximum of the two extractor world dimensions */
    dimworld = ((int)Grid0Patch::dimworld > (int)Grid1Patch::dimworld) ? (int)Grid0Patch::dimworld : (int)Grid1Patch::dimworld

  };

  /** \brief The type used for coordinates
      \todo maybe use traits class to decide which has more precision (DomainGridType::ctype or TargetGridType::ctype) and then take this one
   */
  typedef typename DomainGridType::ctype ctype;

  /** \brief The type used for coordinate vectors */
  typedef Dune::FieldVector<ctype, dimworld>                   Coords;

  /** \brief The type of the domain grid elements */
  typedef typename Grid0View::Traits::template Codim<0>::Entity DomainElement;

  /** \brief Pointer type to domain grid elements */
  typedef typename Grid0View::Traits::template Codim<0>::EntityPointer DomainElementPtr;

  /** \brief The type of the domain grid vertices */
  typedef typename Grid0View::Traits::template Codim<DomainGridType::dimension>::Entity DomainVertex;

  /** \brief Pointer type to domain grid vertices */
  typedef typename Grid0View::Traits::template Codim<DomainGridType::dimension>::EntityPointer DomainVertexPtr;

  /** \brief The type of the target grid elements */
  typedef typename Grid1View::Traits::template Codim<0>::Entity TargetElement;

  /** \brief Pointer type to target grid elements */
  typedef typename Grid1View::Traits::template Codim<0>::EntityPointer TargetElementPtr;

  /** \brief The type of the target grid vertices */
  typedef typename Grid1View::Traits::template Codim<TargetGridType::dimension>::Entity TargetVertex;

  /** \brief Pointer type to target grid vertices */
  typedef typename Grid1View::Traits::template Codim<TargetGridType::dimension>::EntityPointer TargetVertexPtr;

  /** \todo Please doc me! */
  typedef ::Merger<typename DomainGridType::ctype,
      DomainGridType::dimension - Grid0Patch::codim,
      TargetGridType::dimension - Grid1Patch::codim,
      dimworld>                         Merger;

  /** \brief Type of remote intersection objects */
  typedef Dune::GridGlue::Intersection<P0,P1,0,1> Intersection;

  friend class Dune::GridGlue::IntersectionData<P0,P1>;
  friend class Dune::GridGlue::Intersection<P0,P1,0,1>;
  friend class Dune::GridGlue::Intersection<P0,P1,1,0>;
  friend class Dune::GridGlue::IntersectionIterator<P0,P1,0,1>;
  friend class Dune::GridGlue::IntersectionIterator<P0,P1,1,0>;
  friend class Dune::GridGlue::CellIntersectionIterator<P0,P1,0,1>;
  friend class Dune::GridGlue::CellIntersectionIterator<P0,P1,1,0>;
  friend class GridGlueView<P0,P1,0>;
  friend class GridGlueView<P0,P1,1>;

  /** \brief Type of the iterator that iterates over remove intersections */

  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> Grid0IntersectionIterator;
  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,1,0> Grid1IntersectionIterator;

  /** \todo Please doc me! */
  typedef Dune::GridGlue::CellIntersectionIterator<P0,P1,0,1> Grid0CellIntersectionIterator;
  /** \todo Please doc me! */
  typedef Dune::GridGlue::CellIntersectionIterator<P0,P1,1,0> Grid1CellIntersectionIterator;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief the domain surface extractor
  const Grid0Patch&       patch0_;

  /// @brief the target surface extractor
  const Grid1Patch&       patch1_;

  /// @brief the surface merging utility
  Merger*                 merger_;

  /// @brief number of intersections
  IndexType index__sz;

#if HAVE_MPI
  /// @brief MPI_Comm which this GridGlue is working on
  MPI_Comm mpicomm;

  /// @brief parallel indexSet for the intersections with a local domain entity
  PIndexSet domain_is;

  /// @brief parallel indexSet for the intersections with a local target entity
  PIndexSet target_is;

  /// @brief keeps information about which process has which intersection
  Dune::RemoteIndices<PIndexSet> remoteIndices;
#endif

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
   * @param patch0local boolean to indicate if this patch is local (i.e. the entity is accessible) or not
   *
   * @param patch1coords the patch2 vertices' coordinates ordered like e.g. in 3D x_0 y_0 z_0 x_1 y_1 ... y_(n-1) z_(n-1)
   * @param patch1entities just like with the patch0entities and patch0corners
   * @param patch1types array with all patch1 entities types
   * @param patch1local boolean to indicate if this patch is local (i.e. the entity is accessible) or not
   *
   */
  void mergePatches(const std::vector<Dune::FieldVector<ctype,dimworld> >& patch0coords,
                    const std::vector<unsigned int>& patch0entities,
                    const std::vector<Dune::GeometryType>& patch0types,
                    const bool patch0local,
                    const std::vector<Dune::FieldVector<ctype,dimworld> >& patch1coords,
                    const std::vector<unsigned int>& patch1entities,
                    const std::vector<Dune::GeometryType>& patch1types,
                    const bool patch1local);


  template<typename Extractor>
  void extractGrid (const Extractor & extractor,
                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                    std::vector<unsigned int> & faces,
                    std::vector<Dune::GeometryType>& geometryTypes) const;

  template<int P>
  bool getIntersectionIndices(const typename GridGlueView<P0,P1,P>::GridElement& e, std::vector<unsigned int> & indices) const;

public:

  /*   C O N S T R U C T O R S   A N D   D E S T R U C T O R S   */

  /**
   * @brief constructor
   *
   * Initializes components but does not "glue" the surfaces. The surfaces
   * are extracted from the grids here though.
   * @param gv1 the domain grid view
   * @param gv2 the target grid view
   * @param matcher The matcher object that is used to compute the merged grid. This class has
   * to be a model of the SurfaceMergeConcept.
   */
#if HAVE_MPI
  GridGlue(const Grid0Patch& gp1, const Grid1Patch& gp2, Merger* merger, MPI_Comm mpicomm = MPI_COMM_WORLD);
#else
  GridGlue(const Grid0Patch& gp1, const Grid1Patch& gp2, Merger* merger);
#endif
  /*   G E T T E R S   */

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
   * @brief gets an iterator over all remote intersections in the merged grid between domain and target
   *
   * @return the iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::IntersectionIterator ibegin() const
  {
    return typename GridGlueView<P0,P1,I>::IntersectionIterator(this, 0);
  }


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::IntersectionIterator iend() const
  {
    return typename GridGlueView<P0,P1,I>::IntersectionIterator(this, index__sz);
  }


  /**
   * @brief gets an iterator over the intersections of a given codim 0 entity in grid I
   *
   * @tparam I number of the inside grid
   * @param e codim 0 entity in grid I
   * @return the iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::CellIntersectionIterator
  ibegin(const typename GridGlueView<P0,P1,I>::GridElement& e) const;

  /**
   * @brief gets the end iterator for the intersections of a given codim 0 entity in grid I
   *
   * @tparam I number of the inside grid
   * @param e codim 0 entity in grid I
   * @return the end-iterator
   */
  template<int I>
  typename GridGlueView<P0,P1,I>::CellIntersectionIterator
  iend(const typename GridGlueView<P0,P1,I>::GridElement& e) const;


  /*! \brief Communicate information on the MergedGrid of a GridGlue

     Template parameter is a model of Dune::GridGlue::CommDataHandle

     \param data GridGlueDataHandle
     \param iftype Interface for which the Communication should take place
     \param dir Communication direction (Forward means Domain to Target, Backward is the reverse)
   */

  template<class DataHandleImp, class DataTypeImp>
  void communicate (Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> & data,
                    Dune::InterfaceType iftype, Dune::CommunicationDirection dir) const
  {
    typedef Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> DataHandle;
    typedef typename DataHandle::DataType DataType;

#if HAVE_MPI

    /*
     * P A R A L L E L   V E R S I O N
     */

    // setup communication interfaces
    typedef Dune::EnumItem <Dune::PartitionType, Dune::InteriorEntity> InteriorFlags;
    typedef Dune::EnumItem <Dune::PartitionType, Dune::OverlapEntity>  OverlapFlags;
    typedef Dune::EnumRange <Dune::PartitionType, Dune::InteriorEntity, Dune::GhostEntity>  AllFlags;
    Dune::Interface interface;
    switch (iftype)
    {
    case Dune::InteriorBorder_InteriorBorder_Interface :
      interface.build (remoteIndices, InteriorFlags(), InteriorFlags() );
      break;
    case Dune::InteriorBorder_All_Interface :
      if (dir == Dune::ForwardCommunication)
        interface.build (remoteIndices, InteriorFlags(), AllFlags() );
      else
        interface.build (remoteIndices, AllFlags(), InteriorFlags() );
      break;
    case Dune::Overlap_OverlapFront_Interface :
      interface.build (remoteIndices, OverlapFlags(), OverlapFlags() );
      break;
    case Dune::Overlap_All_Interface :
      if (dir == Dune::ForwardCommunication)
        interface.build (remoteIndices, OverlapFlags(), AllFlags() );
      else
        interface.build (remoteIndices, AllFlags(), OverlapFlags() );
      break;
    case Dune::All_All_Interface :
      interface.build (remoteIndices, AllFlags(), AllFlags() );
      break;
    default :
      DUNE_THROW(Dune::NotImplemented, "GridGlue::communicate for interface " << iftype << " not implemented");
    }

    // setup communication info (class needed to tunnel all info to the operator)
    typedef Dune::GridGlue::CommInfo<GridGlue,DataHandleImp,DataTypeImp> CommInfo;
    CommInfo commInfo;
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

#else
    /*
     * S E Q U E N T I A L   V E R S I O N
     */

    // get comm buffer size
    int ssz = indexSet_size() * 10;     // times data per intersection
    int rsz = indexSet_size() * 10;

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
           dir : Forward (domain -> target)
         */
        if (rit->hasDomain())
        {
          data.gather(gatherbuffer, rit->entityDomain(), *rit);
        }
      }
      else       // (dir == Dune::BackwardCommunication)
      {
        /*
           dir : Backward (target -> domain)
         */
        if (rit->hasTarget())
        {
          data.gather(gatherbuffer, rit->entityTarget(), *rit);
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
           dir : Forward (domain -> target)
         */
        if (rit->hasTarget())
          data.scatter(scatterbuffer, rit->entityTarget(), *rit,
                       data.size(*rit));
      }
      else       // (dir == Dune::BackwardCommunication)
      {
        /*
           dir : Backward (target -> domain)
         */
        if (rit->hasDomain())
          data.scatter(scatterbuffer, rit->entityDomain(), *rit,
                       data.size(*rit));
      }
    }

    // cleanup pointers
    delete[] sendbuffer;
    delete[] receivebuffer;
#endif
  }

#if QUICKHACK_INDEX
  /*
   * @brief return an IndexSet mapping from RemoteIntersection to IndexType
   */

  // indexset size
  size_t indexSet_size() const
  {
    return index__sz;
  }

#endif

  Intersection getIntersection(int i) const
  {
    return Intersection(*this, i);
  }

  size_t size() const
  {
    return index__sz;
  }

};

#include "gridglue.cc"

#include "intersection.hh"
#include "intersectioniterator.hh"

#endif // GRIDGLUE_HH_
