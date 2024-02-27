# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception

#install headers
install(FILES
  computeintersection.hh
  computeintersection.cc
  conformingmerge.hh
  contactmerge.cc
  contactmerge.hh
  intersectionlist.hh
  merger.hh
  overlappingmerge.cc
  overlappingmerge.hh
  simplexintersection.cc
  standardmerge.hh
  DESTINATION include/dune/grid-glue/merging)

if(dune-common_VERSION VERSION_GREATER_EQUAL 2.10.0)
    dune_add_library(dunegridglue standardmerge.cc EXPORT_NAME GridGlue)
    target_link_libraries(dunegridglue PUBLIC ${DUNE_LIBS})
else()
    dune_add_library(dunegridglue
            standardmerge.cc
            ADD_LIBS ${DUNE_LIBS})
endif()
