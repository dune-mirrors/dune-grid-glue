// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception
#ifndef DUNE_GRIDGLUE_COMMON_PROJECTIONWRITER_HH
#define DUNE_GRIDGLUE_COMMON_PROJECTIONWRITER_HH

#include <iostream>
#include <string>

#include <dune/grid-glue/common/projection.hh>

namespace Dune {
namespace GridGlue {

/**
 * \brief write projection in VTK format
 *
 * This file writes the projection information in VTK format to the given output
 * stream. It is intended to be used for debugging so the details of the output
 * might change at any time.
 *
 * \note This currently only works in 3D!
 *
 * \param projection projection result
 * \param corners corners of the projected triangles
 * \param normals normals of the projected triangles
 * \param out output stream
 */
template<typename Coordinate, typename Corners, typename Normals>
void write(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals,
           std::ostream& out);

/**
 * \brief write projection in VTK format
 *
 * This function works the same as \ref write(const Projection<Coordinate>&, const Corners&, const Normals&, std::ostream&)
 * but creates a new file with the given name.
 */
template<typename Coordinate, typename Corners, typename Normals>
void write(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals,
           const std::string& filename);
/**
 * \brief Print information about the projection to std::cout stream
 *
 * This method is mainly intended to be used for debugging.
 *
 * \param projection projection result
 * \param corners corners of the projected triangles
 * \param normals normals of the projected triangles
 */
template<typename Coordinate, typename Corners, typename Normals>
void print(const Projection<Coordinate>& projection,
           const Corners& corners,
           const Normals& normals);

} /* namespace GridGlue */
} /* namespace Dune */

#include "projectionwriter_impl.hh"

#endif
