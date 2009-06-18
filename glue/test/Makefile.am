# -*- tab-width: 4; indent-tabs-mode: nil -*-

TESTPROGS = couplingtest


# programs just to build when "make check" is used
check_PROGRAMS = $(TESTPROGS)


AM_CPPFLAGS += $(PSURFACE_CPPFLAGS)
AM_LDFLAGS += $(LOCAL_LIBS) $(PSURFACE_LDFLAGS) $(PSURFACE_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS)

# define the programs
couplingtest_SOURCES = couplingtest.cc

include $(top_srcdir)/am/global-rules
