# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception

add_subdirectory(adapter)
add_subdirectory(merging)
add_subdirectory(common)
add_subdirectory(extractors)
add_subdirectory(test)
install(FILES
  gridglue.hh
  DESTINATION include/dune/grid-glue)
