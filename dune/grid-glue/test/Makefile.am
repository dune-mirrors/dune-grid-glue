# -*- tab-width: 4; indent-tabs-mode: nil -*-

TESTPROGS = \
    callmergertwicetest        \
    mixeddimcouplingtest       \
    mixeddimoverlappingtest    \
    multivectortest            \
    nonoverlappingcouplingtest \
    overlappingcouplingtest

if MPI
TESTPROGS += nonoverlappingcouplingtest_mpi
endif
if UG
TESTPROGS += disconnectedtest
endif

# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)

# programs to run when "make check" is used
TESTS = $(TESTPROGS)

AM_CPPFLAGS += $(PSURFACE_CPPFLAGS)
AM_LDFLAGS  += $(PSURFACE_LDFLAGS) $(AMIRAMESH_LDFLAGS)

# If UG is available we can test pyramids and prisms
AM_CPPFLAGS += $(UG_CPPFLAGS)
AM_LDFLAGS  += $(UG_LDFLAGS)
LDADD       += $(UG_LIBS)

# define the programs
callmergertwicetest_SOURCES = callmergertwicetest.cc
nonoverlappingcouplingtest_SOURCES = nonoverlappingcouplingtest.cc
nonoverlappingcouplingtest_CPPFLAGS = $(AM_CPPFLAGS) # -DCALL_MERGER_TWICE
if UG
disconnectedtest_SOURCES = disconnectedtest.cc
endif
if MPI
nonoverlappingcouplingtest_mpi_SOURCES = nonoverlappingcouplingtest.cc
nonoverlappingcouplingtest_mpi_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS) -DWRITE_TO_VTK
nonoverlappingcouplingtest_mpi_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS)
nonoverlappingcouplingtest_mpi_LDADD = $(LDADD) $(DUNEMPILIBS)
endif
mixeddimcouplingtest_SOURCES = mixeddimcouplingtest.cc
mixeddimoverlappingtest_SOURCES = mixeddimoverlappingtest.cc
multivectortest_SOURCES = multivectortest.cc
overlappingcouplingtest_SOURCES = overlappingcouplingtest.cc
overlappingcouplingtest_CPPFLAGS = $(AM_CPPFLAGS) -frounding-math

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules
