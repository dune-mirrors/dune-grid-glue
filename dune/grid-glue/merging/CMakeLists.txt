#install headers
install(FILES
  conformingmerge.hh
  contactmerge.cc
  contactmerge.hh
  merger.hh
  overlappingmerge.cc
  overlappingmerge.hh
  psurfacemerge.hh
  standardmerge.cc
  standardmerge.hh
  DESTINATION include/dune/grid-glue/merging)

dune_add_library(dunegridglue
  standardmerge.cc
  ADD_LIBS ${DUNE_LIBS})
add_dune_psurface_flags(dunegridglue)
