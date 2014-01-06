# $Id: $

lib_LTLIBRARIES = libdunegridglue.la

libdunegridglue_la_SOURCES = \
  ../dune/grid-glue/merging/psurfacemerge.cc \
  ../dune/grid-glue/merging/standardmerge.cc
libdunegridglue_la_CPPFLAGS = $(AM_CPPFLAGS) $(PSURFACE_CPPFLAGS)
libdunegridglue_la_LDFLAGS  = $(AM_LDFLAGS) $(PSURFACE_LDFLAGS) $(AMIRAMESH_LDFLAGS)
libdunegridglue_la_LIBADD = $(DUNE_LIBS) $(PSURFACE_LIBS) $(AMIRAMESH_LIBS)

#  ../dune/grid-glue/merging/conformingmerge.cc

EXTRA_DIST = CMakeLists.txt

include $(top_srcdir)/am/global-rules
