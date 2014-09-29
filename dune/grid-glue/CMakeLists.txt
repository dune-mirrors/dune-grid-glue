add_subdirectory(adapter)
add_subdirectory(merging)
add_subdirectory(common)
add_subdirectory(extractors)
add_subdirectory(test EXCLUDE_FROM_ALL)
install(FILES
  gridglue.hh
  DESTINATION include/dune/grid-glue)
