# SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
# SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception

#install headers
install(FILES
  areawriter.hh
  areawriter_impl.hh
  crossproduct.hh
  projection.hh
  projection_impl.hh
  projectionwriter.hh
  projectionwriter_impl.hh
  ringcomm.hh
  DESTINATION include/dune/grid-glue/common)
