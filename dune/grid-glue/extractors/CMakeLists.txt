# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception

#install headers
install(FILES
  codim0extractor.hh
  codim1extractor.hh
  extractor.hh
  vtksurfacewriter.hh
  DESTINATION include/dune/grid-glue/extractors)
