<!--
SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-GPL-2.0-only-with-dune-grid-glue-exception
-->

# dune-grid-glue

The `dune-grid-glue` module provides infrastructure for the coupling of two unrelated Dune grids.
The coupling may be overlapping or nonoverlapping, conforming or nonconforming.
The two grids are not requested to be of the same type, and they may even be of different dimensions.  Here are a few possible scenarios:

<img src="/doc/gfx/coupling_scenarios.png" class="center-block" style="width:80%; max-width:600px">

Couplings are described as sets of remote intersections.
Conceptually, these remote intersections are very close to what the regular intersections
in the Dune grid interface are, with the difference that the `inside` and `outside` entities
are taken from different grids.
You can iterate over the global set of remote intersections or over the ones
of a given element.  This allows to assemble the terms needed, for example,
for mortar methods and other domain decomposition methods.

<img src="/doc/gfx/nonoverlapping_coupling.png" style="width:80%">

The _dune-grid-glue_module follows the usual Dune philosophy:

* It defines abstract interfaces to general grid coupling mechanisms, allowing to implement most existing domain decomposition algorithms.
* It allows and encourages the use of existing coupling implementations as legacy backends.
* The code should be efficient, using generic programming where appropriate.

The actual computation of the remote intersections is handled by exchangeable backends. Currently, three backends are available:

* `OverlappingMerge`: For overlapping couplings in 1d, 2d, and 3d
* `ContactMerge`: For nonoverlapping couplings, including contact problems, where there is a positive distance between the two contact boundaries.
* `ConformingMerge`: A fast implementation for conforming nonoverlapping couplings.

All three backends are based on the optimal-time advancing front algorithm by Martin Gander and Caroline Japhet.

* [M. Gander, C. Japhet, An Algorithm for Non-Matching Grid Projections with Linear Complexity, In 'Domain Decomposition Methods in Science and Engineering XVIII', Springer, 2009, pp. 185-192](https://dx.doi.org/10.1007/978-3-642-02677-5_19)

## Releases

Releases of _dune-grid-glue_ usually appear together with releases of the
Dune core and staging modules.  Release tarballs can be download from the
main [Dune project size](https://www.dune-project.org).  Recent releases are

<table>
<tr>
  <th>version</th>
  <th>source</th>
  <th>signature</th>
</tr>
<tr>
  <td>2.9.0</td>
  <td><a href="https://dune-project.org/download/dune-grid-glue/dune-grid-glue-2.9.0.tar.gz" download>dune-grid-glue-2.9.0.tar.gz</a></td>
  <td><a href="https://dune-project.org/download/dune-grid-glue/dune-grid-glue-2.9.0.tar.gz.asc" download>dune-grid-glue-2.9.0.tar.gz.asc</a></td>
</tr>
<tr>
  <td>2.8.0</td>
  <td><a href="https://dune-project.org/download/dune-grid-glue/dune-grid-glue-2.8.0.tar.gz" download>dune-grid-glue-2.8.0.tar.gz</a></td>
  <td><a href="https://dune-project.org/download/dune-grid-glue/dune-grid-glue-2.8.0.tar.gz.asc" download>dune-grid-glue-2.8.0.tar.gz.asc</a></td>
</tr>
</table>

You may also check your favorite Linux distribution.  Some of them have _dune-grid-glue_
available as binary packages.

## Installation

Please see the [general instructions for building DUNE modules](https://www.dune-project.org/doc/installation-notes.html) for detailed instructions on how to build the module.

## Documentation

There is also a class-documentation created by [doxygen](http://www.stack.nl/~dimitri/doxygen/):

* [API documentation for _dune-grid-glue_ master](https://dune-project.org//doxygen/dune-grid-glue/master)

## Maintainers
_dune-grid-glue_ has been written by [Christian Engwer](https://www.uni-muenster.de/AMM/engwer/team/engwer.shtml) and [Oliver Sander](https://tu-dresden.de/mn/math/numerik/sander), based on work by Gerrit Buse (Universität Stuttgart & TU München). Lots of help has come from Ansgar Burchardt, Katja Hanowski, and Jonathan Youett.

## Publications

The concepts of _dune-grid-glue_ are presented in the following publication:

* [P. Bastian, G. Buse, O. Sander: Infrastructure for the Coupling of Dune Grids, In 'Proceedings of ENUMATH 2009', Springer, 2010, pp. 107-114](https://dx.doi.org/10.1007/978-3-642-11795-4_10)
* [C. Engwer, S. Müthing: Concepts for flexible parallel multi-domain simulations, In 'Domain Decomposition Methods in Science and Engineering XXII', Springer](https://dx.doi.org/10.1007/978-3-319-18827-0_17)

## License

The `dune-grid-glue` module is licensed under the GNU Lesser General Public License, version 3 or later, or the GNU General Public License, version 2, with a special runtime exception.

Please see the COPYING file for details.
