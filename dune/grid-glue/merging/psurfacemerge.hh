// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/**
   @file
   @author Gerrit Buse, Christian Engwer, Oliver Sander
   @brief implementation of the SurfaceMerge concept based on libpsurface

   This implementation of the SurfaceMerge concept can be used to
   compute 1d, 2d and 3d intersections. It uses psurface routines to
   compute the merged grid and provides access to it via the interface
   specified by the concept.
 */

#ifndef DUNE_GRIDGLUE_MERGING_PSURFACEMERGE_HH
#define DUNE_GRIDGLUE_MERGING_PSURFACEMERGE_HH

#warning PSurfaceMerge is deprecated. Please use ContactMerge or OverlappingMerge instead.

#include <memory>

#include <dune/common/shared_ptr.hh>

#include <dune/grid-glue/merging/contactmerge.hh>
#include <dune/grid-glue/merging/overlappingmerge.hh>

#if HAVE_PSURFACE
#include <psurface/DirectionFunction.h>
#else
// forward declaration of PSurface classes
template <int dim, typename ctype> class DirectionFunction;
// switch off the macro that contains (in certain versions) the psurface namespace prefix
#define PSURFACE_NAMESPACE
#endif

namespace Dune {
namespace GridGlue {

#if HAVE_PSURFACE
namespace Implementation {

template<typename Vector>
class PSurfaceDirectionFunctionAdapter
  : public Dune::VirtualFunction<Vector, Vector>
{
  using ctype = typename Vector::field_type;
  const static unsigned int dim = Vector::dimension;
  using DF = typename PSURFACE_NAMESPACE DirectionFunction<dim, ctype>;
  using ADF = typename PSURFACE_NAMESPACE AnalyticDirectionFunction<dim, ctype>;

  std::shared_ptr<const ADF> m_direction;

public:
  PSurfaceDirectionFunctionAdapter(const std::shared_ptr<const DF> direction)
    : m_direction(std::dynamic_pointer_cast<const ADF>(direction))
    {
      if (!m_direction)
        DUNE_THROW(Dune::Exception, "Only psurface's AnalyticDirectionFunction is supported.");
    }

  void evaluate(const Vector& x, Vector& y) const override
    {
      PSURFACE_NAMESPACE StaticVector<ctype, dim> x_, y_;
      for (std::size_t i = 0; i < dim; ++i)
        x_[i] = x[i];
      y_ = (*m_direction)(x_);
      for (std::size_t i = 0; i < dim; ++i)
        y[i] = y_[i];
    }
};

} /* namespace Implementation */
#endif

/** \brief Standard implementation of the SurfaceMerge concept using the psurface library.

   \tparam dim Grid dimension of the coupling grids.  Must be the same for both sides
   \tparam dimworld  Dimension of the world coordinates.  Must be equal to dim or to dim+1
   \tparam T Type used for coordinates
 */
template<int dim, int dimworld, typename T = double>
class PSurfaceMerge
  : public ContactMerge<dimworld, T>
{
  using Base = ContactMerge<dimworld, T>;

  static_assert(dim+1 == dimworld, "The PSurfaceMerger only supports dim==dimworld and dim+1==dimworld");

public:
  /// @brief the numeric type used in this interface
  typedef T ctype;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dimworld>  WorldCoords;

  /// @brief the coordinate type used in this interface
  typedef Dune::FieldVector<T, dim>  LocalCoords;

private:
#if HAVE_PSURFACE
  /** \brief Vector field on the domain surface which prescribes the direction
      in which the domain surface is projected onto the target surface
   */
  std::shared_ptr<const Dune::VirtualFunction<WorldCoords, WorldCoords> > psurfaceDomainDirections_;

  /** \brief Vector field on the target surface which prescribes a 'forward'
      direction.

      PSurface uses the normals of the target side to increase projection
      robustness.  If these cannot be computed from the surface directly
      (e.g. because it is not properly oriented), they can be given
      explicitly through the targetDirections field.
   */
  std::shared_ptr<const Dune::VirtualFunction<WorldCoords, WorldCoords> > psurfaceTargetDirections_;
#endif

public:

  PSurfaceMerge(const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype>* domainDirections,
                const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype>* targetDirections)
    : PSurfaceMerge(Dune::stackobject_to_shared_ptr(*domainDirections), Dune::stackobject_to_shared_ptr(*targetDirections))
    { /* Nothing. */ }

  PSurfaceMerge(std::shared_ptr<const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype> > domainDirections = nullptr,
                std::shared_ptr<const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype> > targetDirections = nullptr)
    {
      this->minNormalAngle(0.0);
      this->enableFallback(true);
      setSurfaceDirections(domainDirections, targetDirections);
    }

  using Base::setSurfaceDirections;

  /**
   * @brief Set surface direction functions
   *
   * The matching of the geometries offers the possibility to specify a function for
   * the exact evaluation of domain surface normals. If no such function is specified
   * (default) normals are interpolated.
   * @param domainDirections the new function for the outer normal of grid0 (domain) (or NULL to unset the function)
   * @param targetDirections the new function for the outer normal of grid1 (domain) (or NULL to unset the function)
   */
  void setSurfaceDirections(const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype>* domainDirections,
                            const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype>* targetDirections)
    {
      setSurfaceDirections(Dune::stackobject_to_shared_ptr(*domainDirections), Dune::stackobject_to_shared_ptr(*targetDirections));
    }

  void setSurfaceDirections(std::shared_ptr<const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype> > domainDirections,
                            std::shared_ptr<const PSURFACE_NAMESPACE DirectionFunction<dimworld,ctype> > targetDirections)
    {
#if HAVE_PSURFACE
      using Adapter = Implementation::PSurfaceDirectionFunctionAdapter<WorldCoords>;
      psurfaceDomainDirections_ = domainDirections ? std::make_shared<Adapter>(domainDirections) : nullptr;
      psurfaceTargetDirections_ = targetDirections ? std::make_shared<Adapter>(domainDirections) : nullptr;

      this->setSurfaceDirections(psurfaceDomainDirections_.get(), psurfaceTargetDirections_.get());
#endif
    }
};

template<int dimworld, typename T>
class PSurfaceMerge<dimworld, dimworld, T>
  : public OverlappingMerge<dimworld, dimworld, dimworld, T>
{
  /* Nothing. */
};

} /* namespace GridGlue */
} /* namespace Dune */

#endif // DUNE_GRIDGLUE_MERGING_PSURFACEMERGE_HH
