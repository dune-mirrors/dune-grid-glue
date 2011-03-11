// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef REMOTE_SIMPLICIAL_INTERSECTION_HH
#define REMOTE_SIMPLICIAL_INTERSECTION_HH

#include <dune/common/fvector.hh>
#include <dune/common/array.hh>

template <class T, int grid1Dim, int grid2Dim, int dimworld>
struct RemoteSimplicialIntersection
{
  // Local coordinates in the grid1 entity
  Dune::array<Dune::FieldVector<T,grid1Dim>, dimworld+1> grid1Local_;

  // Local coordinates in the grid1 entity
  Dune::array<Dune::FieldVector<T,grid2Dim>, dimworld+1> grid2Local_;

  //
  int grid1Entity_;

  int grid2Entity_;

};

#endif
