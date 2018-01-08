if(UG_PARALLEL STREQUAL "yes")
  set(HAVE_UG_PARALLEL 1)
else()
  set(HAVE_UG_PARALLEL 0)
endif()

dune_add_test(SOURCES callmergertwicetest.cc)
dune_add_test(SOURCES ringcommtest.cc
              MPI_RANKS 1 2 3 4 5 6 7 8 9 10
              TIMEOUT 300)
dune_add_test(SOURCES computecyclicordertest.cc)
dune_add_test(SOURCES disconnectedtest.cc
              CMAKE_GUARD UG_FOUND "NOT HAVE_UG_PARALLEL")
dune_add_test(SOURCES mixeddimcouplingtest.cc)
dune_add_test(SOURCES mixeddimoverlappingtest.cc)
dune_add_test(SOURCES mixeddimscalingtest.cc)
dune_add_test(SOURCES nonoverlappingcouplingtest.cc
              MPI_RANKS 1 2 4
              TIMEOUT 300
              COMPILE_DEFINITIONS "CALL_MERGER_TWICE")
dune_add_test(SOURCES overlappingcouplingtest.cc
              COMPILE_DEFINITIONS "HAVE_UG_PARALLEL=${HAVE_UG_PARALLEL}"
              COMPILE_FLAGS "-frounding-math")
dune_add_test(SOURCES projectiontest.cc)
