#ifndef DUNE_GRIDGLUE_ADAPTER_RANGEGENERATORS_HH
#define DUNE_GRIDGLUE_ADAPTER_RANGEGENERATORS_HH

#include <dune/common/version.hh>

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)
#  include <dune/common/iteratorrange.hh>
#endif

namespace Dune {
namespace GridGlue {

template<bool reverse>
struct Reverse
  : std::integral_constant<bool, reverse>
{
  typedef Reverse type;

  constexpr
  Reverse<!reverse> operator!() const
    { return {}; }
};

namespace {
const Reverse<true> reversed;
} /* namespace */

#if DUNE_VERSION_NEWER(DUNE_COMMON, 2, 4)

template<typename P0, typename P1, bool reverse = false>
IteratorRange<typename GridGlueView<P0, P1, reverse ? 1 : 0>::IntersectionIterator>
intersections(const GridGlue<P0, P1>& glue, const Reverse<reverse>& = {})
{
  const static int side = reverse ? 1 : 0;
  return {glue.template ibegin<side>(), glue.template iend<side>()};
}

#endif

} /* namespace GridGlue */
} /* namespace Dune */

#endif
