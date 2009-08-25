# -*- tab-width: 4; indent-tabs-mode: nil -*-

TESTPROGS = couplingtest overlappingcouplingtest mixeddimcouplingtest orientedsubfacetest


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)


AM_CPPFLAGS += $(PSURFACE_CPPFLAGS)
AM_LDFLAGS += $(LOCAL_LIBS) $(PSURFACE_LDFLAGS) $(PSURFACE_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS)

# If UG is available we can test pyramids and prisms
AM_CPPFLAGS += $(UG_CPPFLAGS)
AM_LDFLAGS  += $(UG_LDFLAGS) $(UG_LIBS)

# define the programs
couplingtest_SOURCES = couplingtest.cc
mixeddimcouplingtest_SOURCES = mixeddimcouplingtest.cc
overlappingcouplingtest_SOURCES = overlappingcouplingtest.cc
orientedsubfacetest_SOURCES = orientedsubfacetest.cc

include $(top_srcdir)/am/global-rules
