
# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

SUBDIRS = m4 glue doc

# possible options
#LDADD = $(UG_LDFLAGS) $(AMIRAMESH_LDFLAGS) $(UG_LIBS) $(AMIRAMESH_LIBS)
#AM_CPPFLAGS = $(UG_CPPFLAGS) $(AMIRAMESH_CPPFLAGS)


dune_grid_glue_CXXFLAGS = $(MPI_CPPFLAGS) $(UG_CPPFLAGS) $(PSURFACE_CPPFLAGS) $(AMIRAMESH_CPPFLAGS) $(ALBERTA_CPPFLAGS)
dune_grid_glue_LDADD = $(MPI_LDFLAGS) $(ALBERTA_LDFLAGS) $(ALBERTA_LIBS) $(PSURFACE_LDFLAGS) $(PSURFACE_LIBS) $(AMIRAMESH_LDFLAGS) $(AMIRAMESH_LIBS) $(UG_LDFLAGS) $(UG_LIBS) $(MPI_LDFLAGS) $(DUNE_LDFLAGS) $(DUNE_LIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"


include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules

