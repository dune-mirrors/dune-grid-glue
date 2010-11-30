# -*- tab-width: 4; indent-tabs-mode: nil -*-

TESTPROGS = nonoverlappingcouplingtest overlappingcouplingtest mixeddimcouplingtest

if MPI
TESTPROGS += nonoverlappingcouplingtest_mpi
endif

if UG
TESTPROGS += orientedsubfacetest
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
nonoverlappingcouplingtest_SOURCES = nonoverlappingcouplingtest.cc
if MPI
nonoverlappingcouplingtest_mpi_SOURCES = nonoverlappingcouplingtest.cc
nonoverlappingcouplingtest_mpi_CPPFLAGS = $(AM_CPPFLAGS) $(DUNEMPICPPFLAGS)
nonoverlappingcouplingtest_mpi_LDFLAGS = $(AM_LDFLAGS) $(DUNEMPILDFLAGS)
endif
mixeddimcouplingtest_SOURCES = mixeddimcouplingtest.cc
overlappingcouplingtest_SOURCES = overlappingcouplingtest.cc
orientedsubfacetest_SOURCES = orientedsubfacetest.cc

include $(top_srcdir)/am/global-rules
