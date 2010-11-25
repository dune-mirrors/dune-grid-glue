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

#include "../extractors/gridextractiontraits.hh"
#include "../extractors/extractorselector.hh"
#include "../extractors/extractorpredicate.hh"
#include "coordinatetransformation.hh"
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

  typedef GET1 DomGridExtractionTraits;
  typedef GET2 TarGridExtractionTraits;

  /*   P R I V A T E   T Y P E S   */

  typedef GridGlue<GET1, GET2> This;

  /** \brief GlobalId type of an intersection (used for communication) */
  typedef unsigned int GlobalId;

  /** \brief LocalIndex type of an intersection (used for communication) */
  typedef Dune::ParallelLocalIndex<Dune::PartitionType> LocalIndex;

  /** \brief ParallelIndexSet type (used for communication communication) */
  typedef Dune::ParallelIndexSet <GlobalId, LocalIndex > PIndexSet;

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
    domdim = DomainExtractor::dim,
    /** \brief WOrld dimension of the domain extractor */
    domdimworld = DomainExtractor::dimworld
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
    tardim = TargetExtractor::dim,
    /** \brief World dimension of the target extractor */
    tardimworld = TargetExtractor::dimworld
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

  /** \brief Type of remote intersection objects */
  class RemoteIntersection;

  /** \brief Type of the iterator that iterates over remove intersections */
  class RemoteIntersectionIterator;

  /** \todo Please doc me! */
  class DomainIntersectionIterator;

  /** \todo Please doc me! */
  class TargetIntersectionIterator;

private:

  typedef ExtractorPredicate<DomainGridView, DomainExtractor::codim> DomainDescriptor;

  typedef ExtractorPredicate<TargetGridView, TargetExtractor::codim> TargetDescriptor;

  /*   M E M B E R   V A R I A B L E S   */

  const DomainDescriptor*                   domelmntdescr_;

  const TargetDescriptor*                   tarelmntdescr_;

  const DomainTransformation*    domtrafo_;

  const TargetTransformation*    tartrafo_;

  /// @brief the "domain" grid view
  const DomainGridView&        domgv_;

  /// @brief the "target" grid view
  const TargetGridView&        targv_;

  /// @brief the domain surface extractor
  DomainExtractor domext_;

  /// @brief the target surface extractor
  TargetExtractor tarext_;

  /// @brief the surface merging utility
  Merger* merger_;

  /// @brief number of intersections
  unsigned int index__sz;

  /// @brief number of domain intersections
  unsigned int _dindex_sz;

  /// @brief number of target intersections
  unsigned int _tindex_sz;

  /// @brief an invalid intersection object used as dummy and
  /// also as recognizable end object of iterations over intersections
  mutable RemoteIntersection NULL_INTERSECTION;

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

  /// @brief a vector with intersection elements
  mutable std::vector<RemoteIntersection>   intersections_;

protected:

  /**
   * @brief after building the merged grid the intersection can be updated
   * through this method (for internal use)
   */
  void updateIntersections()
  {
    // build the intersections array again
    this->intersections_.resize(this->merger_->nSimplices(), this->NULL_INTERSECTION);
    _dindex_sz = 0;
    _tindex_sz = 0;
    for (unsigned int i = 0; i < this->merger_->nSimplices(); ++i)
    {

      RemoteIntersection ri(this, i);
      if (ri.hasDomain())
        ri.domainIndex() = _dindex_sz++;
      if (ri.hasTarget())
        ri.targetIndex() = _tindex_sz++;
      ri.index() = i;
      this->intersections_[i] = ri;

    }

    std::cout << "GridGlue::updateIntersections : The number of remote intersections is " << merger_->nSimplices()
              << " with " << _dindex_sz << " domain entities"
              << " and " << _tindex_sz << " target entities" << std::endl;

    // store the number of remote intersection for later use
    index__sz = merger_->nSimplices();

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

  template<typename Extractor>
  void extractGrid (const Extractor & extractor,
                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                    std::vector<unsigned int> & faces,
                    std::vector<Dune::GeometryType>& geometryTypes,
                    const CoordinateTransformation<Extractor::dimworld, dimworld, ctype>* trafo) const;

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

  DomainExtractor & domainExtractor() { return domext_; }

  TargetExtractor & targetExtractor() { return tarext_; }

  /*   G E T T E R S   */

  /**
   * @brief getter for the domain grid view
   * @return the object
   */
  const DomainGridView& domainGridView() const
  {
    return this->domgv_;
  }


  /**
   * @brief getter for the target grid view
   * @return the object
   */
  const TargetGridView& targetGridView() const
  {
    return this->targv_;
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
    return this->merger_;
  }

  void setDomainDescriptor(const DomainDescriptor& descr)
  {
    domelmntdescr_ = &descr;
  }


  void setTargetDescriptor(const TargetDescriptor& descr)
  {
    tarelmntdescr_ = &descr;
  }

  void setDomainTransformation(const DomainTransformation* trafo)
  {
    domtrafo_ = trafo;
  }


  void setTargetTransformation(const TargetTransformation* trafo)
  {
    tartrafo_ = trafo;
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
    return RemoteIntersectionIterator(this->NULL_INTERSECTION);
  }


  /**
   * @brief gets the (general) end-iterator for iterations over domain codim 0 entities' faces
   *
   * @return the iterator
   */
  DomainIntersectionIterator idomainend() const
  {
    return DomainIntersectionIterator(this->NULL_INTERSECTION);
  }


  /**
   * @brief gets the (general) end-iterator for iterations over target codim 0 entities' faces
   *
   * @return the iterator
   */
  TargetIntersectionIterator itargetend() const
  {
    return TargetIntersectionIterator(this->NULL_INTERSECTION);
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
    int ssz = this->indexSet_size() * 10;     // times data per intersection
    int rsz = this->indexSet_size() * 10;

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

  size_t domainIndexSet_size() const
  {
    return _dindex_sz;
  }

  size_t targetIndexSet_size() const
  {
    return _tindex_sz;
  }
#endif

  RemoteIntersection getIntersection(int i) const
  {
    return RemoteIntersection(intersections_[i]);
  }

};

/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

template<typename GET1, typename GET2>
#if HAVE_MPI
GridGlue<GET1, GET2>::GridGlue(const DomainGridView& gv1, const TargetGridView& gv2, Merger* merger, MPI_Comm m)
#else
GridGlue<GET1, GET2>::GridGlue(const DomainGridView & gv1, const TargetGridView & gv2, Merger* merger)
#endif
  : domelmntdescr_(NULL), tarelmntdescr_(NULL),
    domtrafo_(NULL), tartrafo_(NULL),
    domgv_(gv1), targv_(gv2),
    domext_(gv1), tarext_(gv2), merger_(merger),
    NULL_INTERSECTION(this),
#if HAVE_MPI
    mpicomm(m),
#endif
    intersections_(0, NULL_INTERSECTION)
{
  std::cout << "GridGlue: Constructor succeeded!" << std::endl;
}

template<typename GET1, typename GET2>
void GridGlue<GET1, GET2>::build()
{
  // setup the domain surface extractor
  if (domelmntdescr_ != NULL)
    domext_.update(*domelmntdescr_);
  else
    DUNE_THROW(Dune::Exception, "GridGlue::Builder : no domain surface descriptor set");

  // setup the target surface extractor
  if (tarelmntdescr_ != NULL)
    tarext_.update(*tarelmntdescr_);
  else
    DUNE_THROW(Dune::Exception, "GridGlue::Builder : no target surface descriptor set");

  // clear the contents from the current intersections array
  {
    std::vector<RemoteIntersection> dummy(0, NULL_INTERSECTION);
    intersections_.swap(dummy);
  }

  std::vector<Dune::FieldVector<ctype, dimworld> > domcoords;
  std::vector<unsigned int> domfaces;
  std::vector<Dune::GeometryType> domainElementTypes;
  std::vector<Dune::FieldVector<ctype,dimworld> > tarcoords;
  std::vector<unsigned int> tarfaces;
  std::vector<Dune::GeometryType> targetElementTypes;

  /*
   * extract global surface grids
   */

  // retrieve the coordinate and topology information from the extractors
  // and apply transformations if necessary
  extractGrid(domext_, domcoords, domfaces, domainElementTypes, domtrafo_);
  extractGrid(tarext_, tarcoords, tarfaces, targetElementTypes, tartrafo_);

#ifdef WRITE_TO_VTK
  const int dimw = Parent::dimworld;
  const char prefix[] = "GridGlue::Builder::build() : ";
  char domainsurf[256];
  sprintf(domainsurf, "/tmp/vtk-domain-test-%i", mpi_rank);
  char targetsurf[256];
  sprintf(targetsurf, "/tmp/vtk-target-test-%i", mpi_rank);

  std::cout << prefix << "Writing domain surface to '" << domainsurf << ".vtk'...\n";
  VtkSurfaceWriter vtksw(domainsurf);
  vtksw.writeSurface(domcoords, domfaces, dimw, dimw);
  std::cout << prefix << "Done writing domain surface!\n";

  std::cout << prefix << "Writing target surface to '" << targetsurf << ".vtk'...\n";
  vtksw.setFilename(targetsurf);
  vtksw.writeSurface(tarcoords, tarfaces, dimw, dimw);
  std::cout << prefix << "Done writing target surface!\n";
#endif // WRITE_TO_VTK


  // start the actual build process
  merger_->build(domcoords, domfaces, domainElementTypes,
                 tarcoords, tarfaces, targetElementTypes);

  // the intersections need to be recomputed
  updateIntersections();
}


template<typename GET1, typename GET2>
template<typename Extractor>
void GridGlue<GET1, GET2>::extractGrid (const Extractor & extractor,
                                        std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                                        std::vector<unsigned int> & faces,
                                        std::vector<Dune::GeometryType>& geometryTypes,
                                        const CoordinateTransformation<Extractor::dimworld, dimworld, ctype>* trafo) const
{
  std::vector<typename Extractor::Coords> tempcoords;
  std::vector<typename Extractor::VertexVector> tempfaces;

  extractor.getCoords(tempcoords);
  coords.clear();
  coords.reserve(tempcoords.size());

  if (trafo != NULL)
  {
    std::cout << "GridGlue::Builder : apply trafo\n";
    for (size_t i = 0; i < tempcoords.size(); ++i)
    {
      Coords temp = (*trafo)(tempcoords[i]);
      coords.push_back(temp);
    }
  }
  else
  {
    for (unsigned int i = 0; i < tempcoords.size(); ++i)
    {
      assert(int(dimworld) == int(Extractor::dimworld));
      coords.push_back(Dune::FieldVector<ctype, dimworld>());
      for (size_t j = 0; j <dimworld; ++j)
        coords.back()[j] = tempcoords[i][j];
    }
  }

  extractor.getFaces(tempfaces);
  faces.clear();

  for (unsigned int i = 0; i < tempfaces.size(); ++i) {
    for (unsigned int j = 0; j < tempfaces[i].size(); ++j)
      faces.push_back(tempfaces[i][j]);
  }

  // get the list of geometry types from the extractor
  extractor.getGeometryTypes(geometryTypes);

}


template<typename GET1, typename GET2>
int GridGlue<GET1, GET2>::domainEntityNextFace(const DomainElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->domext_.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count &&    (this->domext_.indexInInside(first) < index || !this->merger_->template simplexMatched<0>(first)))
    first++;
  if (first == count)
    return -1;     // no more faces
  else
    return this->domext_.indexInInside(first);     // found, return the face's number
}


template<typename GET1, typename GET2>
int GridGlue<GET1, GET2>::targetEntityNextFace(const TargetElement& e, int index) const
{
  int first, count;
  // first check if the element forms a part of the extracted surface
  if (!this->tarext_.faceIndices(e, first, count))
    return -1;

  // check all mapped faces and accept the first one with number >=index
  count += first;
  while (first < count && (this->tarext_.indexInInside(first) < index || !this->merger_->template simplexMatched<1>(first)))
    first++;
  if (first == count)
    return -1;     // no more faces
  else
    return this->tarext_.indexInInside(first);     // found, return the face's number
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::RemoteIntersectionIterator GridGlue<GET1, GET2>::iremotebegin() const
{
  return (this->intersections_.size() > 0)
         ? RemoteIntersectionIterator(intersections_[0])
         : RemoteIntersectionIterator(NULL_INTERSECTION);
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::DomainIntersectionIterator GridGlue<GET1, GET2>::idomainbegin(const DomainElement& e, int num) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->domext_.faceIndices(e, first, count);
  if (!in_surface)
    return DomainIntersectionIterator(this->NULL_INTERSECTION);

  count += first;
  while (first < count)
  {
    if (this->domext_.indexInInside(first) == num && this->merger_->template simplexMatched<0>(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      this->merger_->template simplexRefined<0>(first, global_results);
      while (++first < count && this->domext_.indexInInside(first) == num && this->merger_->template simplexRefined<0>(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return DomainIntersectionIterator(this->intersections_[global_results[0]], global_results);

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return DomainIntersectionIterator(this->NULL_INTERSECTION);
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::DomainIntersectionIterator GridGlue<GET1, GET2>::idomainbegin(const DomainElement& e) const
{
  // first check if the element has at least one extracted subEntity
  int first, count;
  bool hasExtractedSubEntity = this->domext_.faceIndices(e, first, count);
  if (!hasExtractedSubEntity)
    return DomainIntersectionIterator(this->NULL_INTERSECTION);


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplex refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->merger_->template simplexRefined<0>(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return DomainIntersectionIterator(this->intersections_[global_results[0]], global_results);
  else
    return DomainIntersectionIterator(this->NULL_INTERSECTION);
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::TargetIntersectionIterator GridGlue<GET1, GET2>::itargetbegin(const TargetElement& e, int num) const
{
  // first check if the element has at least one extracted subEntity
  int first, count;
  bool hasExtractedSubEntity = this->tarext_.faceIndices(e, first, count);
  if (!hasExtractedSubEntity) return itargetend();

  count += first;
  while (first < count)
  {
    if (this->tarext_.indexInInside(first) == num && this->merger_->template simplexMatched<1>(first))
    {
      // perfect candidate found! done searching bec. of consecutive order of extracted simplices!
      std::vector<unsigned int> global_results;
      std::vector<unsigned int> local_results;

      // get the remote intersections
      this->merger_->template simplexRefined<1>(first, global_results);
      while (++first < count && this->tarext_.indexInInside(first) == num && this->merger_->template simplexRefined<1>(first, local_results))
      {
        for (unsigned int i = 0; i < local_results.size(); ++i)
          global_results.push_back(local_results[i]);
      }

      // if sth. has been found, return the iterator
      if (global_results.size() > 0)
        return TargetIntersectionIterator(this->intersections_[global_results[0]], global_results);

      // else leave the loop
      break;
    }
    first++;
  }

  // nothing has been found
  return TargetIntersectionIterator(this->NULL_INTERSECTION);
}


template<typename GET1, typename GET2>
typename GridGlue<GET1, GET2>::TargetIntersectionIterator GridGlue<GET1, GET2>::itargetbegin(const TargetElement& e) const
{
  // first check if the element forms a part of the extracted surface
  int first, count;
  bool in_surface = this->tarext_.faceIndices(e, first, count);
  if (!in_surface) return itargetend();


  // now accumulate all remote intersections of the element's faces
  std::vector<unsigned int> global_results(0, 0);
  std::vector<unsigned int> local_results;

  // iterate over all simplices to check if there is more than one simplix refining the face
  bool found_sth = false;
  count += first;
  while (first < count)
  {
    if (this->merger_->template simplexRefined<1>(first, local_results))
    {
      if (local_results.size() > 0)
        found_sth = true;
      for (unsigned int i = 0; i < local_results.size(); ++i)
        global_results.push_back(local_results[i]);
    }
    first++;
  }

  if (found_sth)
    return TargetIntersectionIterator(this->intersections_[global_results[0]], global_results);
  else
    return TargetIntersectionIterator(this->NULL_INTERSECTION);
}

// include implementation of subclass RemoteIntersection
#include "gridglueremoteintersectionimpl.hh"

// include implementation of subclasses DomainIntersectionIterator and TargetIntersectionIterator
#include "gridglueremoteintersectioniteratorimpl.hh"


#endif // GRIDGLUE_HH_
