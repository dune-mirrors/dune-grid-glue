// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    codim1extractor.hh
 *  Version:     1.0
 *  Created on:  Jun 23, 2009
 *  Author:      Oliver Sander, Christian Engwer
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: class for grid extractors extracting surface grids
 *
 */
/**
 * @file
 * @brief Grid extractor class for codim 1 subgrids
 */

#ifndef DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH
#define DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH

#include "extractor.hh"
#include "extractorpredicate.hh"

#include <deque>

#include <dune/grid-glue/common/projectionhelper.hh>

namespace Dune {

  namespace GridGlue {

template<typename GV>
class Codim1Extractor : public Extractor<GV,1>
{
public:

  /*  E X P O R T E D  T Y P E S   A N D   C O N S T A N T S  */

  using Extractor<GV,1>::dimworld;
  using Extractor<GV,1>::dim;
  using Extractor<GV,1>::codim;
  using Extractor<GV,1>::cube_corners;
  typedef typename Extractor<GV,1>::IndexType IndexType;

  /// @brief compile time number of corners of surface simplices
  enum
  {
    simplex_corners = dim
  };

  typedef GV GridView;

  typedef typename GV::Grid::ctype ctype;
  typedef Dune::FieldVector<ctype, dimworld>                       Coords;

  typedef typename GV::Traits::template Codim<dim>::EntityPointer VertexPtr;
  typedef typename GV::Traits::template Codim<dim>::Entity Vertex;
  typedef typename GV::Traits::template Codim<0>::EntityPointer ElementPtr;
  typedef typename GV::Traits::template Codim<0>::Entity Element;

  static const Dune::PartitionIteratorType PType = Dune::Interior_Partition;
  typedef typename GV::Traits::template Codim<0>::template Partition<PType>::Iterator ElementIter;

  typedef typename GV::IntersectionIterator IsIter;

  // import typedefs from base class
  typedef typename Extractor<GV,1>::SubEntityInfo SubEntityInfo;
  typedef typename Extractor<GV,1>::ElementInfo ElementInfo;
  typedef typename Extractor<GV,1>::VertexInfo VertexInfo;
  typedef typename Extractor<GV,1>::CoordinateInfo CoordinateInfo;
  typedef typename Extractor<GV,1>::VertexInfoMap VertexInfoMap;

public:

  /*  C O N S T R U C T O R S   A N D   D E S T R U C T O R S  */

  /**
   * @brief Constructor
   * @param gv the grid view object to work with
   * @param descr a predicate to mark entities for extraction (unary functor returning bool)
   */
  Codim1Extractor(const GV& gv, const ExtractorPredicate<GV,1>& descr)
    :  Extractor<GV,1>(gv)
  {
    std::cout << "This is Codim1Extractor on a <" << dim
              << "," << dimworld << "> grid!"
              << std::endl;
    update(descr);
  }

private:

  /**
   * Extracts a codimension 1 surface from the grid @c g and builds up two arrays
   * with the topology of the surface written to them. The description of the
   * surface part that is to be extracted is given in form of a predicate class.
   *
   * Assumed that we are in 2D the coords array will have the structure
   * x0 y0 x1 y1 ... x(n-1) y(n-1)
   * Values in the @c _indices array then refer to the indices of the coordinates, e.g.
   * index 1 is associated with the position x1. If the surface consists of triangles
   * we have always groups of 3 indices describing one triangle.
   *
   * @param descr a predicate class that "selects" the faces to add to the surface
   */
  void update(const ExtractorPredicate<GV,1>& descr);

};


template<typename GV>
void Codim1Extractor<GV>::update(const ExtractorPredicate<GV,1>& descr)
{
  // free everything there is in this object
  this->clear();

  // In this first pass iterate over all entities of codim 0.
  // For each codim 1 intersection check if it is part of the boundary and if so,
  // get its corner vertices, find resp. store them together with their associated index,
  // and remember the indices of the boundary faces' corners.
  {
    // several counter for consecutive indexing are needed
    int simplex_index = 0;
    int vertex_index = 0;
    IndexType eindex = 0;     // supress warning

    // needed later for insertion into a std::set which only
    // works with const references

    // a temporary container where newly acquired face
    // information can be stored at first
    std::deque<SubEntityInfo> temp_faces;

    // iterate over interior codim 0 elements on the grid
    for (ElementIter elit = this->gv_.template begin<0, PType>();
         elit != this->gv_.template end<0, PType>(); ++elit)
    {
      const auto& elmt = *elit;
      Dune::GeometryType gt = elmt.type();

      // if some face is part of the surface add it!
      if (elit->hasBoundaryIntersections())
      {
        // add an entry to the element info map, the index will be set properly later,
        // whereas the number of faces is already known
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
        eindex = this->cellMapper_.index(elmt);
#else
        eindex = this->cellMapper_.map(elmt);
#endif
        this->elmtInfo_[eindex] = new ElementInfo(simplex_index, elmt, 0);

        // now add the faces in ascending order of their indices
        // (we are only talking about 1-4 faces here, so O(n^2) is ok!)
        for (IsIter is = this->gv_.ibegin(elmt); is != this->gv_.iend(elmt); ++is)
        {
          // Stop only at selected boundary faces
          if (!is->boundary() or !descr.contains(elmt, is->indexInInside()))
            continue;

#if DUNE_VERSION_NEWER(DUNE_GEOMETRY,2,3)
          const Dune::ReferenceElement<ctype, dim>& refElement = Dune::ReferenceElements<ctype, dim>::general(gt);
#else
          const Dune::GenericReferenceElement<ctype, dim>& refElement = Dune::GenericReferenceElements<ctype, dim>::general(gt);
#endif
          // get the corner count of this face
          const int face_corners = refElement.size(is->indexInInside(), 1, dim);

          // now we only have to care about the 3D case, i.e. a triangle face can be
          // inserted directly whereas a quadrilateral face has to be divided into two triangles
          switch (face_corners)
          {
          case 2 :
          case 3:
          {
            // we have a simplex here

            // register the additional face(s)
            this->elmtInfo_[eindex]->faces++;

            // add a new face to the temporary collection
            temp_faces.push_back(SubEntityInfo(eindex, is->indexInInside(),
                                               Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)));

            std::vector<FieldVector<ctype,dimworld> > cornerCoords(face_corners);

            // try for each of the faces vertices whether it is already inserted or not
            for (int i = 0; i < face_corners; ++i)
            {
              // get the number of the vertex in the parent element
              int vertex_number = refElement.subEntity(is->indexInInside(), 1, i, dim);

              // get the vertex pointer and the index from the index set
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
              const Vertex vertex = elit->template subEntity<dim>(vertex_number);
#else
              const Vertex& vertex = *elit->template subEntity<dim>(vertex_number);
#endif
              cornerCoords[i] = vertex.geometry().corner(0);

              IndexType vindex = this->gv_.indexSet().template index<dim>(vertex);

              // remember the vertex' number in parent element's vertices
              temp_faces.back().corners[i].num = vertex_number;

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->vtxInfo_.find(vindex);
              if (vimit == this->vtxInfo_.end())
              {
                // insert into the map
                this->vtxInfo_[vindex] = new VertexInfo(vertex_index, vertex);
                // remember the vertex as a corner of the current face in temp_faces
                temp_faces.back().corners[i].idx = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only insert the index into the simplices array
                temp_faces.back().corners[i].idx = vimit->second->idx;
              }
            }

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want them to be oriented such that all boundary edges
            // point in the counterclockwise direction.  Therefore, we check the orientation of the
            // new face and possibly switch the two vertices.
            FieldVector<ctype,dimworld> realNormal = is->centerUnitOuterNormal();

            // Compute segment normal
            FieldVector<ctype,dimworld> reconstructedNormal;
            if (dim==2)  // boundary face is a line segment
            {
              reconstructedNormal[0] =  cornerCoords[1][1] - cornerCoords[0][1];
              reconstructedNormal[1] =  cornerCoords[0][0] - cornerCoords[1][0];
            } else {    // boundary face is a triangle
              FieldVector<ctype,dimworld> segment1 = cornerCoords[1] - cornerCoords[0];
              FieldVector<ctype,dimworld> segment2 = cornerCoords[2] - cornerCoords[0];
              reconstructedNormal = Projection::crossProduct(segment1, segment2);
            }
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);

            // now increase the current face index
            simplex_index++;
            break;
          }
          case 4 :
          {
            assert(dim == 3);
            // we have a quadrilateral here
            unsigned int vertex_indices[4];
            unsigned int vertex_numbers[4];

            // register the additional face(s) (2 simplices)
            this->elmtInfo_[eindex]->faces += 2;

            Dune::array<FieldVector<ctype,dimworld>, 4> cornerCoords;

            // get the vertex pointers for the quadrilateral's corner vertices
            // and try for each of them whether it is already inserted or not
            for (int i = 0; i < cube_corners; ++i)
            {
              // get the number of the vertex in the parent element
              vertex_numbers[i] = refElement.subEntity(is->indexInInside(), 1, i, dim);

              // get the vertex pointer and the index from the index set
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
              const Vertex vertex = elit->template subEntity<dim>(vertex_numbers[i]);
#else
              const Vertex &vertex = *elit->template subEntity<dim>(vertex_numbers[i]);
#endif
              cornerCoords[i] = vertex.geometry().corner(0);

              IndexType vindex = this->gv_.indexSet().template index<dim>(vertex);

              // if the vertex is not yet inserted in the vertex info map
              // it is a new one -> it will be inserted now!
              typename VertexInfoMap::iterator vimit = this->vtxInfo_.find(vindex);
              if (vimit == this->vtxInfo_.end())
              {
                // insert into the map
                this->vtxInfo_[vindex] = new VertexInfo(vertex_index, vertex);
                // remember this vertex' index
                vertex_indices[i] = vertex_index;
                // increase the current index
                vertex_index++;
              }
              else
              {
                // only remember the vertex' index
                vertex_indices[i] = vimit->second->idx;
              }
            }

            // now introduce the two triangles subdividing the quadrilateral
            // ATTENTION: the order of vertices given by "orientedSubface" corresponds to the order
            // of a Dune quadrilateral, i.e. the triangles are given by 0 1 2 and 3 2 1

            // add a new face to the temporary collection for the first tri
            temp_faces.push_back(SubEntityInfo(eindex, is->indexInInside(),
                                               Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)));
            temp_faces.back().corners[0].idx = vertex_indices[0];
            temp_faces.back().corners[1].idx = vertex_indices[1];
            temp_faces.back().corners[2].idx = vertex_indices[2];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[0];
            temp_faces.back().corners[1].num = vertex_numbers[1];
            temp_faces.back().corners[2].num = vertex_numbers[2];

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want the triangle vertices on counterclockwise order,
            // when viewed from the outside of the grid. Therefore, we check the orientation of the
            // new face and possibly switch two vertices.
            FieldVector<ctype,dimworld> realNormal = is->centerUnitOuterNormal();

            // Compute segment normal
            FieldVector<ctype,dimworld> reconstructedNormal = Projection::crossProduct(cornerCoords[1] - cornerCoords[0],
                                                                                       cornerCoords[2] - cornerCoords[0]);
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);


            // add a new face to the temporary collection for the second tri
            temp_faces.push_back(SubEntityInfo(eindex, is->indexInInside(),
                                               Dune::GeometryType(Dune::GeometryType::simplex,dim-codim)));
            temp_faces.back().corners[0].idx = vertex_indices[3];
            temp_faces.back().corners[1].idx = vertex_indices[2];
            temp_faces.back().corners[2].idx = vertex_indices[1];
            // remember the vertices' numbers in parent element's vertices
            temp_faces.back().corners[0].num = vertex_numbers[3];
            temp_faces.back().corners[1].num = vertex_numbers[2];
            temp_faces.back().corners[2].num = vertex_numbers[1];

            // Now we have the correct vertices in the last entries of temp_faces, but they may
            // have the wrong orientation.  We want the triangle vertices on counterclockwise order,
            // when viewed from the outside of the grid. Therefore, we check the orientation of the
            // new face and possibly switch two vertices.
            // Compute segment normal
            reconstructedNormal = Projection::crossProduct(cornerCoords[2] - cornerCoords[3],
                                                           cornerCoords[1] - cornerCoords[3]);
            reconstructedNormal /= reconstructedNormal.two_norm();

            if (realNormal * reconstructedNormal < 0.0)
              std::swap(temp_faces.back().corners[0], temp_faces.back().corners[1]);

            simplex_index+=2;
            break;
          }
          default :
            DUNE_THROW(Dune::NotImplemented, "the extractor does only work for triangle and quadrilateral faces (" << face_corners << " corners)");
            break;
          }
        }         // end loop over found surface parts
      }
    }     // end loop over elements

    std::cout << "added " << simplex_index << " subfaces\n";

    // allocate the array for the face specific information...
    this->subEntities_.resize(simplex_index);
    // ...and fill in the data from the temporary containers
    copy(temp_faces.begin(), temp_faces.end(), this->subEntities_.begin());
  }


  // now first write the array with the coordinates...
  this->coords_.resize(this->vtxInfo_.size());
  typename VertexInfoMap::const_iterator it1 = this->vtxInfo_.begin();
  for (; it1 != this->vtxInfo_.end(); ++it1)
  {
    // get a pointer to the associated info object
    CoordinateInfo* current = &this->coords_[it1->second->idx];
    // store this coordinates index // NEEDED?
    current->index = it1->second->idx;
    // store the vertex' index for the index2vertex mapping
    current->vtxindex = it1->first;
    // store the vertex' coordinates under the associated index
    // in the coordinates array
#if DUNE_VERSION_NEWER(DUNE_GRID, 2, 4)
    const auto vtx = this->grid().entity(it1->second->p);
#else
    const auto& vtx = *this->grid().entityPointer(it1->second->p);
#endif
    current->coord = vtx.geometry().corner(0);
  }

}

}  // namespace GridGlue

}  // namespace Dune

#endif // DUNE_GRIDGLUE_EXTRACTORS_CODIM1EXTRACTOR_HH
