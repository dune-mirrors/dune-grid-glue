
# we need the module file to be able to build via dunecontrol
EXTRA_DIST= \
  CHANGELOG.md \
  CHANGES \
  CMakeLists.txt \
  COPYING \
  README.md \
  TODO \
  dune-grid-glue.pc.in \
  dune.module

SUBDIRS = m4 dune doc lib

# don't follow the full GNU-standard
# we need automake 1.9 or newer
AUTOMAKE_OPTIONS = foreign 1.9
# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) --with-dune-geometry=$(DUNE_GEOMETRY_ROOT) --with-dune-grid=$(DUNE_GRID_ROOT) --with-dune-istl=$(DUNE_ISTL_ROOT) --with-dune-localfunctions=$(DUNE_LOCALFUNCTIONS_ROOT) --with-dune-pdelab=$(DUNE_PDELAB_ROOT)  CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
include $(top_srcdir)/am/top-rules

configdir = $(datadir)/dune-grid-glue
dist_config_DATA = config.h.cmake

include $(top_srcdir)/am/cmake-pkg-config
