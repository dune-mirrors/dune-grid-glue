Changes in dune-grid-glue v2.4.0
================================

Incompatible changes in dune-grid-glue v2.4.0
---------------------------------------------

* All interfaces have been moved into the `Dune::GridGlue` namespace.

* The `contains()` method of the `Codim0Extractor` and `Codim1Extractor`
  classes now take a codim-0-`Entity` instead of a codim-0-`EntityPointer`.
  New code can also use a `std::function` instead of the extractor classes.

* Methods that returned an `EntityPointer` now return an `Entity` instead
  when `dune-grid-glue` is built against version 2.4 or later of the DUNE
  core modules.

Major changes in dune-grid-glue v2.4.0
--------------------------------------

* This is the first release that supports version 2.4 of the DUNE core modules.
  The older 2.2 and 2.3 releases are also still supported.

* `dune-grid-glue` now requires C++11 support.

* Both `ContactMerge` and `OverlappingMerge` had large changes and were
  partially rewritten.  The newer version should both be faster and hopefully
  have less bugs as well.

* Support for the range-based for statement has been added.  It is
  possible to iterate over all intersections of a `GridGlue` object by
  code like
  ```
  GridGlue<...> glue;
  for (const auto& in : intersections(glue)) ...;
  ```
  See the documentation for further details.

* `PSurfaceMerge` has been replaced by a wrapper around `ContactMerge` and
  `OverlappingMerge`.  It is now also deprecated and will be removed in the
  next release.
