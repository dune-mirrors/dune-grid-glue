// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    orientedsubface.hh
 *  Version:     1.0
 *  Created on:  Jan 28, 2009
 *  Author:      Oliver Sander
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: Renumber vertices of a face so they are correctly oriented
 *  subversion:  $Id$
 *
 */
/**
 * @file
 * @brief Renumber vertices of a face so they are correctly oriented
 */

#ifndef ORIENTED_SUBFACE_HH
#define ORIENTED_SUBFACE_HH

#include <dune/common/geometrytype.hh>
#include <dune/grid/common/genericreferenceelements.hh>

template <int dim>
int orientedSubface(const Dune::GeometryType& type, int face, int vertex)
{
  const Dune::GenericReferenceElement<double,dim>& refElement =
    Dune::GenericReferenceElements<double, dim>::general(type);

  // edges
  if (dim == 1)
    return vertex;

  // Triangle
  if (type.isTriangle() && face==1)
    return refElement.subEntity(face,1, (vertex+1)%2, dim);

  // Quadrilateral
  if (type.isQuadrilateral() && (face==0 || face==3))
    return refElement.subEntity(face,1, (vertex+1)%2, dim);

  // Tetrahedron, face 0 and 2:  swap local vertices 1 and 2
  if (type.isTetrahedron() && (face==0 || face==2) && vertex>=1)
    return refElement.subEntity(face,1, (vertex%2)+1, dim);

  // Pyramid, quad face (0):  swap 1 and 2
  if (type.isPyramid() && face==0 && vertex>=1 && vertex<=2)
    return refElement.subEntity(face,1, (vertex%2)+1, dim);

  // Pyramid, faces 1 and 4: swap 0 and 1
  if (type.isPyramid() && (face==1 || face==4) && vertex<2)
    return refElement.subEntity(face,1, (vertex+1)%2, dim);

  // Prism, face 3: swap vertices 1 and 2
  if (type.isPrism() && face==3 && vertex>=1)
    return refElement.subEntity(face,1, (vertex%2)+1, dim);

  // Prism: face 1: swap vertices 0 and 1; and 2 and 3
  if (type.isPrism() && face==1)
    return (vertex<2)
           ? refElement.subEntity(face,1, (vertex+1)%2, dim)
           : refElement.subEntity(face,1, ((vertex-1)%2)+2, dim);

  // Hexahedron:
  // face 0: 0 2 4 6 --> 2 0 6 4
  // face 3: 2 3 6 7 --> 3 2 7 6
  if (type.isHexahedron() && (face==0 || face==3))
    return refElement.subEntity(face,1, (vertex<2) ? (vertex+1)%2 : ((vertex-1)%2)+2, dim);


  // face 4: 0 1 2 3 --> 0 2 1 3
  if (type.isHexahedron() && face==4) {
    const int renumber[4] = {0,2,1,3};
    return refElement.subEntity(face,1, renumber[vertex], dim);
  }

  return refElement.subEntity(face,1,vertex,dim);

}

#endif // ORIENTED_SUBFACE_HH
