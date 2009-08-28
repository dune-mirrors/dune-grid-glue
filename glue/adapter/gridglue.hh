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
 * @file GridGlue.hh
 * @brief Central component of the module implementing the coupling of two grids.
 */


#ifndef GRIDGLUE_HH_
#define GRIDGLUE_HH_

#include <dune/common/array.hh>

#define QUICKHACK_INDEX 1

#ifdef GRID_GLUE_USE_CONCEPTS
#include "../misc/conceptchecking.hh"
#endif
#include "../extractors/gridextractor.hh"
#include "../merging/psurfacemerge.hh"
#include "../extractors/vtksurfacewriter.hh"
#include "../extractors/gridextractiontraits.hh"
#include "../extractors/extractorselector.hh"
#include "../extractors/surfacedescriptor.hh"
#include "simplexgeometry.hh"
#include "remoteintersection.hh"
#include "remoteintersectioniterators.hh"
#include "coordinatetransformation.hh"
#include "gridgluecommunicate.hh"

#include <dune/istl/indexset.hh>
#include <dune/istl/plocalindex.hh>
#include <dune/istl/remoteindices.hh>
#include <dune/istl/communicator.hh>

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
template<typename GET1, typename GET2>
class GridGlue
{
private:

  /*   C H E C K   C O N C E P T S   */

  typedef GET1 DomGridExtractionTraits;
#ifdef GRID_GLUE_USE_CONCEPTS
  CLASS_REQUIRE(DomGridExtractionTraits, GridExtractionTraitsConcept);
#endif

  typedef GET2 TarGridExtractionTraits;
#ifdef GRID_GLUE_USE_CONCEPTS
  CLASS_REQUIRE(TarGridExtractionTraits, GridExtractionTraitsConcept);
#endif


  /*   P R I V A T E   T Y P E S   */

  typedef GridGlue<GET1, GET2> This;

  /** \brief GlobalId type of an intersection (used for communication) */
  typedef unsigned int GlobalId;

  /** \brief LocalIndex type of an intersection (used for communication) */
  typedef Dune::ParallelLocalIndex<Dune::PartitionType> LocalIndex;

  /** \brief ParallelIndexSet type (used for communication communication) */
  typedef Dune::ParallelIndexSet <GlobalId, LocalIndex > PIndexSet;

public:
  /** \todo Please doc me! */
  template<typename BSET1, typename BSET2, int codim1, int codim2>
  class BuilderImpl;

  class RemoteIntersectionImpl;

  class RemoteIntersectionIteratorImpl;

  class DomainIntersectionIteratorImpl;

  class TargetIntersectionIteratorImpl;

public:

  /*   P U B L I C   T Y P E S   A N D   C O N S T A N T S   */

  /** \brief Grid view of the domain grid */
  typedef typename DomGridExtractionTraits::GridView DomainGridView;

  /** \brief Domain grid type */
  typedef typename DomainGridView::Grid DomainGridType;

  /** \brief Extractor used for the domain grid */
  typedef typename ExtractorSelector<DomGridExtractionTraits>::ExtractorType DomainExtractor;

  /** \brief Dimension of the domain extractor */
  enum {
    /** \brief Dimension of the domain extractor */
    domdim = DomainExtractor::dim
  };

  /** \brief Grid view of the target grid */
  typedef typename TarGridExtractionTraits::GridView TargetGridView;

  /** \brief Target grid type */
  typedef typename TargetGridView::Grid TargetGridType;

  /** \brief Extractor used for the target grid */
  typedef typename ExtractorSelector<TarGridExtractionTraits>::ExtractorType TargetExtractor;

  /** \brief Dimension of the target extractor */
  enum {
    /** \brief Dimension of the target extractor */
    tardim = TargetExtractor::dim
  };


  /** \brief export the world dimension */
  enum {
    /** \brief export the world dimension :
        maximum of the two extractor world dimensions */
    dimworld = ((int)DomainExtractor::dimworld > (int)TargetExtractor::dimworld) ? (int)DomainExtractor::dimworld : (int)TargetExtractor::dimworld

  };

  /** \brief The type used for coordinates
      \todo maybe use traits class to decide which has more precision (DomainGridType::ctype or TargetGridType::ctype) and then take this one
   */
  typedef typename DomainGridType::ctype ctype;

  /** \brief The type used for coordinate vectors */
  typedef Dune::FieldVector<ctype, dimworld>                   Coords;

  /** \brief The type of transformation used for the domain grid*/
  typedef CoordinateTransformation<DomainExtractor::dimworld, dimworld, ctype>      DomainTransformation;

  /** \brief The type of transformation used for the domain grid*/
  typedef CoordinateTransformation<TargetExtractor::dimworld, dimworld, ctype>      TargetTransformation;

  /** \brief The type of the domain grid elements */
  typedef typename DomainGridView::Traits::template Codim<0>::Entity DomainElement;

  /** \brief Pointer type to domain grid elements */
  typedef typename DomainGridView::Traits::template Codim<0>::EntityPointer DomainElementPtr;

  /** \brief The type of the domain grid vertices */
  typedef typename DomainGridView::Traits::template Codim<DomainGridType::dimension>::Entity DomainVertex;

  /** \brief Pointer type to domain grid vertices */
  typedef typename DomainGridView::Traits::template Codim<DomainGridType::dimension>::EntityPointer DomainVertexPtr;

  /** \brief The type of the target grid elements */
  typedef typename TargetGridView::Traits::template Codim<0>::Entity TargetElement;

  /** \brief Pointer type to target grid elements */
  typedef typename TargetGridView::Traits::template Codim<0>::EntityPointer TargetElementPtr;

  /** \brief The type of the target grid vertices */
  typedef typename TargetGridView::Traits::template Codim<TargetGridType::dimension>::Entity TargetVertex;

  /** \brief Pointer type to target grid vertices */
  typedef typename TargetGridView::Traits::template Codim<TargetGridType::dimension>::EntityPointer TargetVertexPtr;

  /** \todo Please doc me! */
  typedef ::Merger<typename DomainGridType::ctype,
      DomainGridType::dimension - DomainExtractor::codim,
      TargetGridType::dimension - TargetExtractor::codim,
      dimworld>                         Merger;

  /** \todo Please doc me! */
  typedef BuilderImpl<GET1, GET2, DomainExtractor::codim, TargetExtractor::codim>                    Builder;

  /** \brief Type of remote intersection objects */
  typedef RemoteIntersectionInterface::RemoteIntersection<RemoteIntersectionImpl>    RemoteIntersection;

  /** \brief Type of the iterator that iterates over remove intersections */
  typedef RemoteIntersectionInterface::RemoteIntersectionIterator<RemoteIntersectionImpl, RemoteIntersectionIteratorImpl>
  RemoteIntersectionIterator;

  /** \todo Please doc me! */
  typedef RemoteIntersectionInterface::RemoteIntersectionIterator<RemoteIntersectionImpl, DomainIntersectionIteratorImpl>
  DomainIntersectionIterator;

  /** \todo Please doc me! */
  typedef RemoteIntersectionInterface::RemoteIntersectionIterator<RemoteIntersectionImpl, TargetIntersectionIteratorImpl>
  TargetIntersectionIterator;

private:

  /*   M E M B E R   V A R I A B L E S   */

  /// @brief the "domain" grid view
  const DomainGridView&        _domgv;

  /// @brief the "target" grid view
  const TargetGridView&        _targv;

  /// @brief the domain surface extractor
  DomainExtractor _domext;

  /// @brief the target surface extractor
  TargetExtractor _tarext;

  /// @brief the surface merging utility
  Merger* _merg;

  /// @brief the builder utility
  Builder _builder;

  /// @brief number of intersections
  unsigned int _index_sz;

  /// @brief number of domain intersections
  unsigned int _dindex_sz;

  /// @brief number of target intersections
  unsigned int _tindex_sz;

  /// @brief an invalid intersection object used as dummy and
  /// also as recognizable end object of iterations over intersections
  mutable RemoteIntersectionImpl NULL_INTERSECTION;

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

#warning HACK
public:
  /// @brief a vector with intersection elements
  mutable std::vector<RemoteIntersectionImpl>   _intersections;

protected:

  /**
   * @brief after building the merged grid the intersection can be updated
   * through this method (for internal use)
   */
  void updateIntersections()
  {
    // build the intersections array again
    this->_intersections.resize(this->_merg->nSimplices(), this->NULL_INTERSECTION);
    _index_sz = 0;
    _dindex_sz = 0;
    _tindex_sz = 0;
    for (unsigned int i = 0; i < this->_merg->nSimplices(); ++i)
    {
      RemoteIntersectionImpl ri(this, i);
#if HAVE_MPI
      if ((ri.hasTarget() || ri.hasDomain()))
      {
        if (ri.hasDomain())
          ri.domainIndex() = _dindex_sz++;
        if (ri.hasTarget())
          ri.targetIndex() = _tindex_sz++;
        ri.index() = _index_sz;
        this->_intersections[_index_sz++] = ri;
      }
#endif
    }

    std::cout << "GridGlue::updateIntersections : The number of overlaps is " << _index_sz
              << " with " << _dindex_sz << " domain and "
              << " with " << _tindex_sz << " target entities" << std::endl;

    ////// create ParallelIndexSet & RemoteIndices
#if HAVE_MPI
    // setup parallel indexset
    domain_is.beginResize();
    target_is.beginResize();
    RemoteIntersectionIterator rit = iremotebegin();
    RemoteIntersectionIterator ritend = iremoteend();
    for (; rit != ritend; ++rit)
    {
      if (rit->hasDomain())
      {
        domain_is.add (rit->globalIndex(),
                       LocalIndex(rit->index(), rit->entityDomain()->partitionType()) ) ;
      }
      if (rit->hasTarget())
      {
        target_is.add (rit->globalIndex(),
                       LocalIndex(rit->index(), rit->entityTarget()->partitionType()) ) ;
      }
    }
    domain_is.endResize();
    target_is.endResize();

    // setup remote index information
    remoteIndices.setIndexSets(domain_is, target_is, mpicomm) ;
    remoteIndices.rebuild<true>();
#endif
  }


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
  GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Merger* merger, MPI_Comm mpicomm = MPI_COMM_WORLD);
#else
  GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Merger* merger);
#endif
  /*  S E T T E R S  */


  /*   G E T T E R S   */

  /**
   * @brief getter for the domain grid view
   * @return the object
   */
  const DomainGridView& domainGridView() const
  {
    return this->_domgv;
  }


  /**
   * @brief getter for the target grid view
   * @return the object
   */
  const TargetGridView& targetGridView() const
  {
    return this->_targv;
  }


  /**
   * @brief getter for the builder
   * @return the object
   */
  Builder& builder()
  {
    return this->_builder;
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
  const Merger* merger()
  {
    return this->_merg;
  }

  /*   F U N C T I O N A L I T Y   */

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
  RemoteIntersectionIterator iremotebegin() const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  DomainIntersectionIterator idomainbegin(const DomainElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the domain grid
   *
   * @param e codim 0 entity in the domain grid
   * @return the iterator
   */
  DomainIntersectionIterator idomainbegin(const DomainElement& e) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 1 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @param num the index of the face (codim 1 entity) in @c e entity, ignored if only one face in the surface
   * @return the iterator
   */
  TargetIntersectionIterator itargetbegin(const TargetElement& e, int num) const;


  /**
   * @brief gets an iterator over the remote intersections of a given codim 0 entity in the target grid
   *
   * @param e codim 0 entity in the target grid
   * @return the iterator
   */
  TargetIntersectionIterator itargetbegin(const TargetElement& e) const;


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  RemoteIntersectionIterator iremoteend() const
  {
    return RemoteIntersectionIterator(RemoteIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  DomainIntersectionIterator idomainend() const
  {
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }


  /**
   * @brief gets the (general) end-iterator for iterations over target codim 0 entities' faces
   *
   * @return the iterator
   */
  TargetIntersectionIterator itargetend() const
  {
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
  }

  /*! \brief Communicate information on the MergedGrid of a GridGlue

     Template parameter is a model of Dune::GridGlueCommDataHandleIF

     \param data GridGlueDataHandle
     \param iftype Interface for which the Communication should take place
     \param dir Communication direction (Forward means Domain to Target, Backward is the reverse)
   */

  template<class DataHandleImp, class DataTypeImp>
  void communicate (Dune::GridGlueCommDataHandleIF<DataHandleImp,DataTypeImp> & data,
                    Dune::InterfaceType iftype, Dune::CommunicationDirection dir) const
  {
    typedef Dune::GridGlueCommDataHandleIF<DataHandleImp,DataTypeImp> DataHandle;
    typedef typename DataHandle::DataType DataType;

#if HAVE_MPI

    /*
     * P A R A L L E L   V E R S I O N
     */

    // setup communication interfaces
    typedef Dune::EnumItem <Dune::PartitionType, Dune::InteriorEntity> InteriorFlags;
    typedef Dune::EnumItem <Dune::PartitionType, Dune::OverlapEntity>  OverlapFlags;
    typedef Dune::EnumRange <Dune::PartitionType, Dune::InteriorEntity, Dune::GhostEntity>  AllFlags;
    Dune::Interface < PIndexSet > interface;
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
    typedef Dune::GridGlueCommInfo<GridGlue,DataHandleImp,DataTypeImp> CommInfo;
    CommInfo commInfo;
    commInfo.gridglue = this;
    commInfo.data = &data;

    // create communicator
    Dune::BufferedCommunicator<PIndexSet> bComm ;
    bComm.template build< CommInfo >(commInfo, commInfo, interface);

    // do communication
    // choose communication direction.
    if (dir == Dune::ForwardCommunication)
      bComm.forward< Dune::GridGlueForwardOperator >(commInfo, commInfo);
    else
      bComm.backward< Dune::GridGlueBackwardOperator >(commInfo, commInfo);

#else
    /*
     * S E Q U E N T I A L   V E R S I O N
     */

    // get comm buffer size
    int ssz = this->indexSet_size() * 10;     // times data per intersection
    int rsz = this->indexSet_size() * 10;

    // allocate send/receive buffer
    DataType* sendbuffer = new DataType[ssz];
    DataType* receivebuffer = new DataType[rsz];

    // iterators
    RemoteIntersectionIterator rit = iremotebegin();
    RemoteIntersectionIterator ritend = iremoteend();

    // gather
    Dune::GridGlueMessageBuffer<DataType> gatherbuffer(sendbuffer);
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
    Dune::GridGlueMessageBuffer<DataType> scatterbuffer(receivebuffer);
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
    return _index_sz;
  }

  size_t domainIndexSet_size() const
  {
    return _dindex_sz;
  }

  size_t targetIndexSet_size() const
  {
    return _tindex_sz;
  }
#endif
};

/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

template<typename GET1, typename GET2>
#if HAVE_MPI
GridGlue<GET1, GET2>::GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Merger* merger, MPI_Comm m)
#else
GridGlue<GET1, GET2>::GridGlue(const DomainGridView & gv1, const TargetGridView & gv2, Merger* merger)
#endif
  : _domgv(gv1), _targv(gv2),
    _domext(gv1), _tarext(gv2), _merg(merger),
    _builder(*const_cast<GridGlue<GET1, GET2>*>(this)),
    NULL_INTERSECTION(this),
#if HAVE_MPI
    mpicomm(m),
#endif
    _intersections(0, NULL_INTERSECTION)
{
  std::cout << "GridGlue: Constructor succeeded!" << std::endl;
}


template<typename GET1, typename GET2>
int GridGlue<GET1, GET2>::domainEntityNextFace(const DomainElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->_domext.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (this->_domext.indexInInside(first) < index || !this->_merg->domainSimplexMatched(first)))
    first++;
  if (first == count)
    return -1;             // no more faces
  else
    return this->_domext.indexInInside(first);             // found, return the face's number
}


template<typename GET1, typename GET2>
int GridGlue<GET1, GET2>::targetEntityNextFace(const TargetElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->_tarext.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (this->_tarext.indexInInside(first) < index || !this->_merg->targetSimplexMatched(first)))
    first++;
  if (first == count)
    return -1;             // no more faces
  else
    return this->_tarext.indexInInside(first);             // found, return the face's number
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::RemoteIntersectionIterator GridGlue<GET1, GET2>::iremotebegin() const
{
  return RemoteIntersectionIterator(RemoteIntersectionIteratorImpl(this->_intersections[0]));
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::DomainIntersectionIterator GridGlue<GET1, GET2>::idomainbegin(const DomainElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_domext.faceIndices(e, first, count);
  if (!in_surface)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));

  count += first;
  while (first < count)
  {
    if (this->_domext.indexInInside(first) == num && this->_merg->domainSimplexMatched(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      this->_merg->domainSimplexRefined(first, global_results);
      while (++first < count && this->_domext.indexInInside(first) == num && this->_merg->domainSimplexRefined(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::DomainIntersectionIterator GridGlue<GET1, GET2>::idomainbegin(const DomainElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_domext.faceIndices(e, first, count);
  if (!in_surface)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->_merg->domainSimplexRefined(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));
  else
    return DomainIntersectionIterator(DomainIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::TargetIntersectionIterator GridGlue<GET1, GET2>::itargetbegin(const TargetElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_tarext.faceIndices(e, first, count);
  if (!in_surface) return itargetend();

  count += first;
  while (first < count)
  {
    if (this->_tarext.indexInInside(first) == num && this->_merg->targetSimplexMatched(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      this->_merg->targetSimplexRefined(first, global_results);
      while (++first < count && this->_tarext.indexInInside(first) == num && this->_merg->targetSimplexRefined(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::TargetIntersectionIterator GridGlue<GET1, GET2>::itargetbegin(const TargetElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->_tarext.faceIndices(e, first, count);
  if (!in_surface) return itargetend();


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->_merg->targetSimplexRefined(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->_intersections[global_results[0]], global_results));
  else
    return TargetIntersectionIterator(TargetIntersectionIteratorImpl(this->NULL_INTERSECTION));
}

// include implementation of subclass BuilderImpl
#include "gridgluebuilderimpl.hh"

// include implementation of subclass RemoteIntersectionImpl
#include "gridglueremoteintersectionimpl.hh"

// include implementation of subclasses DomainIntersectionIteratorImpl and TargetIntersectionIteratorImpl
#include "gridglueremoteintersectioniteratorimpl.hh"


#endif // GRIDGLUE_HH_
