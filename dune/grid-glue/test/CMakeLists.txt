if(MPI_C_FOUND)
  set(MPI_PROGRAMS nonoverlappingcouplingtest_mpi)
endif(MPI_C_FOUND)

# tests where program to build and program to run are equal
set(TESTS
 callmergertwicetest
 computecyclicordertest
 mixeddimcouplingtest
 mixeddimoverlappingtest
 mixeddimscalingtest
 multivectortest
 nonoverlappingcouplingtest
 overlappingcouplingtest
 projectiontest
 ${MPI_PROGRAMS})

add_executable(callmergertwicetest callmergertwicetest.cc)
add_executable(computecyclicordertest computecyclicordertest.cc)
add_executable(mixeddimcouplingtest mixeddimcouplingtest.cc)
add_executable(mixeddimoverlappingtest mixeddimoverlappingtest.cc)
add_executable(mixeddimscalingtest mixeddimscalingtest.cc)
add_executable(multivectortest multivectortest.cc)
add_executable(nonoverlappingcouplingtest nonoverlappingcouplingtest.cc)
set_property(TARGET nonoverlappingcouplingtest
  APPEND PROPERTY COMPILE_DEFINITIONS "CALL_MERGER_TWICE")
add_executable(overlappingcouplingtest overlappingcouplingtest.cc)
set_property(TARGET overlappingcouplingtest
  APPEND PROPERTY COMPILE_FLAGS "-frounding-math")
add_executable(projectiontest projectiontest.cc)

if(MPI_C_FOUND)
  add_executable(nonoverlappingcouplingtest_mpi nonoverlappingcouplingtest.cc)
  set_property(TARGET nonoverlappingcouplingtest_mpi
    APPEND PROPERTY COMPILE_F "WRITE_TO_VTK")
  add_dune_mpi_flags(nonoverlappingcouplingtest_mpi)
endif(MPI_C_FOUND)

if(UG_FOUND)
  add_executable(disconnectedtest disconnectedtest.cc)
  add_dune_ug_flags(disconnectedtest)
  list(APPEND TESTS disconnectedtest)
endif()

foreach(_exe ${TESTS})
  add_dune_psurface_flags(${_exe})
  target_link_libraries(${_exe} "dunegridglue" ${DUNE_LIBS})
  add_test(${_exe} ${_exe})
endforeach(_exe ${TESTS})

set(SOURCES
  callmergertwicetest.cc
  communicationtest.hh
  couplingtest.hh
  mixeddimcouplingtest.cc
  mixeddimoverlappingtest.cc
  mixeddimscalingtest.cc
  multivectortest.cc
  nonoverlappingcouplingtest.cc
  overlappingcouplingtest.cc)
