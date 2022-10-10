// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception
#ifndef DUNE_GRIDGLUE_COMMON_CROSSPRODUCT_HH
#define DUNE_GRIDGLUE_COMMON_CROSSPRODUCT_HH 1

namespace Dune {
namespace GridGlue {

/**
 * \brief compute cross product
 *
 * \return <code>a × b</code>
 */
template <class T, int dim>
static Dune::FieldVector<T,dim> crossProduct(const Dune::FieldVector<T,dim>& a,
                                             const Dune::FieldVector<T,dim>& b)
{
  if (dim!=3)
    DUNE_THROW(Dune::NotImplemented, "crossProduct does not work for dimension " << dim);

  Dune::FieldVector<T,dim> c;
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  return c;
}

} /* namespace GridGlue */
} /* namespace Dune */

#endif
