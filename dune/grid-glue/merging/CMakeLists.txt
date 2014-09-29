#install headers
install(FILES
  conformingmerge.hh
  contactmerge.hh
  merger.hh
  overlappingmerge.cc
  overlappingmerge.hh
  psurfacemerge.cc
  psurfacemerge.hh
  standardmerge.cc
  standardmerge.hh
  DESTINATION include/dune/grid_glue/merging)

dune_add_library(dunegridglue
  psurfacemerge.cc
  standardmerge.cc
  ADD_LIBS ${DUNE_LIBS})
add_dune_psurface_flags(dunegridglue)
