// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    GridGlue.hh
 *  Version:     1.0
 *  Created on:  Feb 2, 2009
 *  Author:      Gerrit Buse
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

  }
}

template<typename P0, typename P1, int P>
struct GridGlueView;

template<typename P0, typename P1>
struct GridGlueView<P0,P1,0>
{
  typedef P0 Patch;
  static const P0& patch(const GridGlue<P0,P1>& g)
  {
    return g.patch0_;
  }
};

template<typename P0, typename P1>
struct GridGlueView<P0,P1,1>
{
  typedef P1 Patch;
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
  friend class GridGlueView<P0,P1,0>;
  friend class GridGlueView<P0,P1,1>;

  /** \todo */
  typedef Dune::GridGlue::Intersection<P0,P1,0,1> Grid0Intersection;
  typedef Dune::GridGlue::Intersection<P0,P1,1,0> Grid1Intersection;

  /** \brief Type of the iterator that iterates over remove intersections */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> IntersectionIterator;

  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,0,1> Grid0IntersectionIterator;
  /** \todo Please doc me! */
  typedef Dune::GridGlue::IntersectionIterator<P0,P1,1,0> Grid1IntersectionIterator;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief the "domain" grid view
  const Grid0View&        domgv_;

  /// @brief the "target" grid view
  const Grid1View&        targv_;

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
   */
  void updateIntersections()
  {
    // store the number of remote intersection for later use
    index__sz = merger_->nSimplices();

    // build the intersections array again
    intersections_.resize(merger_->nSimplices() + 1);
    for (unsigned int i = 0; i < merger_->nSimplices(); ++i)
    {
      IntersectionData data(*this, i);
      intersections_[i] = data;
    }

    std::cout << "GridGlue::updateIntersections : "
    "The number of remote intersections is " << index__sz << std::endl;

    ////// create ParallelIndexSet & RemoteIndices
#if HAVE_MPI && 0
    // setup parallel indexset
    domain_is.beginResize();
    target_is.beginResize();
    IntersectionIterator rit = ibegin();
    IntersectionIterator ritend = iend();
    for (; rit != ritend; ++rit)
    {
            #error update this!
      if (rit->hasGrid0())
      {
        domain_is.add (rit->globalIndex(),
                       LocalIndex(rit->index(), rit->entityGrid0()->partitionType()) ) ;
      }
      if (rit->hasGrid1())
      {
        target_is.add (rit->globalIndex(),
                       LocalIndex(rit->index(), rit->entityGrid1()->partitionType()) ) ;
      }
    }
    domain_is.endResize();
    target_is.endResize();

    // setup remote index information
    remoteIndices.setIndexSets(domain_is, target_is, mpicomm) ;
    remoteIndices.rebuild<true>();
#endif
  }

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
   * @brief getter for the domain grid view
   * @return the object
   */
    #warning
  const Grid0View& domainGridView() const
  {
    return domgv_;
  }


  /**
   * @brief getter for the target grid view
   * @return the object
   */
  const Grid1View& targetGridView() const
  {
    return targv_;
  }


  /**
   * @brief getter for the surface matcher utility
   *
   * This grants access to the surface merger. This Merger class has to be a model
   * of the SurfaceMergeConcept, but since different implementations of matchers
   * may provide different configuration possibilities this part has to remain
   * implementation specific. Through this getter one can retrieve the internal
   * matcher and configure it before its "build" member is called.
   * @return a (non-const) reference the object
   */
  const Merger* merger() const
  {
    return merger_;
  }

  /*   F U N C T I O N A L I T Y   */

  void build();

  /**
   * @brief tells whether a codim 0 entity's face(s) (or at least a part)
   * could be mapped
   *
   * For the given entity  could be mapped, if the latter is the case the number of the
   * face in the parent element is returned.
   * Note: Calling this function with only @c e given and checking the return value
   * (-1 or not) is the easiest way to determine whether any of the entity's faces were mapped.
   * @param e the element
   * @param index number of the first face that is checked, faces with lower index are ignored
   * @return -1 if there is no face mapped with number >=@c index, else next face's number
   */
  int domainEntityNextFace(const DomainElement& e, int index = 0) const;


  /**
   * @brief tells whether a codim 0 entity's face(s) (or at least a part)
   * could be mapped
   *
   * For the given entity  could be mapped, if the latter is the case the number of the
   * face in the parent element is returned.
   * Note: Calling this function with only @c e given and checking the return value
   * (-1 or not) is the easiest way to determine whether any of the entity's faces were mapped.
   * @param e the element
   * @param index number of the first face that is checked, faces with lower index are ignored
   * @return -1 if there is no face mapped with number >=@c index, else next face's number
   */
  int targetEntityNextFace(const TargetElement& e, int index = 0) const;

  /*   I N T E R S E C T I O N S   A N D   I N T E R S E C T I O N   I T E R A T O R S   */

  /**
   * @brief gets an iterator over all remote intersections in the merged grid between domain and target
   *
   * @return the iterator
   */
  IntersectionIterator ibegin() const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  Grid0IntersectionIterator idomainbegin(const DomainElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @return the iterator
   */
  Grid0IntersectionIterator idomainbegin(const DomainElement& e) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  Grid1IntersectionIterator itargetbegin(const TargetElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @return the iterator
   */
  Grid1IntersectionIterator itargetbegin(const TargetElement & e) const;


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  IntersectionIterator iend() const
  {
    return IntersectionIterator(this, index__sz);
  }


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  Grid0IntersectionIterator idomainend() const
  {
        #warning
    return Grid0IntersectionIterator(this, index__sz);
  }


  /**
   * @brief gets the (general) end-iterator for iterations over target codim 0 entities' faces
   *
   * @return the iterator
   */
  Grid1IntersectionIterator itargetend() const
  {
        #warning
    return Grid1IntersectionIterator(this, index__sz);
  }

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
    RemoteIntersectionIterator rit = iremotebegin();
    RemoteIntersectionIterator ritend = iremoteend();

    // gather
    Dune::GridGlue::StreamingMessageBuffer<DataType> gatherbuffer(sendbuffer);
    for (rit = iremotebegin(); rit != ritend; ++rit)
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
    for (rit = iremotebegin(); rit != ritend; ++rit)
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

};

#include "gridglue.cc"

#include "intersection.hh"
#include "intersectioniterator.hh"

#endif // GRIDGLUE_HH_
