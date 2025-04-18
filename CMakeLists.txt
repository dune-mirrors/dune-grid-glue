# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception

cmake_minimum_required(VERSION 3.13)
project(dune-grid-glue CXX)

if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

#include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

dune_register_package_flags(LIBRARIES dunegridglue)

add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("examples")

# finalize the dune project, e.g. generating config.h etc.
if(dune-common_VERSION VERSION_GREATER_EQUAL 2.10.0)
    finalize_dune_project()
else()
    finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
endif()
