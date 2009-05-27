// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*
 *  Filename:    VtkSurfaceWriter.hh
 *  Version:     1.0
 *  Created on:  Jan 16, 2009
 *  Author:      Gerrit Buse
 *  ---------------------------------
 *  Project:     dune-grid-glue
 *  Description: helper class for graphical output of grids in generic representation
 *  subversion:  $Id$
 *
 */
/**
 * @file VtkSurfaceWriter.hh
 * @brief helper class for graphical output of grids in generic representation
 */

#ifndef VTKSURFACEWRITER_HH_
#define VTKSURFACEWRITER_HH_

#include <fstream>
#include <iomanip>
#include <vector>

#include "../adapter/GridGlueVtkWriter.hh"

using namespace std;


class VtkSurfaceWriter
{
private:




public:


  VtkSurfaceWriter(const char* filename) : _filename(filename)
  {}

  ~VtkSurfaceWriter()
  {}

  void setFilename(const char* name)
  {
    if (strlen(name) > 0)
      this->_filename = name;
  }


  template<typename K>
  void writeSurface(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, int dim)
  {
    ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->_filename);
    fos.open(buffer);
    fos << setprecision(8) << setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    fos.close();
  }


  template<typename K, typename T>
  void writeSurfaceElementData(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, const std::vector<T>& data, const char* dataname, int dim)
  {
    ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->_filename);
    fos.open(buffer);
    fos << setprecision(8) << setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    this->writeCellData(data, dataname, dim, fos);
    fos.close();
  }


  template<typename K, typename T>
  void writeSurfaceVertexData(const std::vector<K>& coords, const std::vector<unsigned int>& indices, int corners, const std::vector<T>& data, const char* dataname, int dim)
  {
    ofstream fos;
    char buffer[64];
    sprintf(buffer, "%s.vtk", this->_filename);
    fos.open(buffer);
    fos << std::setprecision(8) << setw(1);
    // write preamble
    fos << "# vtk DataFile Version 2.0\nFilename: " << buffer << "\nASCII" << endl;
    this->writePoints(coords, dim, fos);
    const int polycount = indices.size()/corners;
    int corner_count[polycount];
    for (int i = 0; i < polycount; ++i)
      corner_count[i] = corners;
    this->writePolygons(indices, corner_count, polycount, dim, fos);
    this->writePointData(data, dataname, dim, fos);
    fos.close();
  }

protected:

  template<typename K>
  void writePoints(const std::vector<K>& coords, int dim, ofstream& fos)
  {
    int coord_count = coords.size() / dim;
    fos << "DATASET POLYDATA\nPOINTS " << coord_count*(dim == 2 ? 2 : 1) << " " << TypeNames[Nametraits<K>::nameidx] << endl;
    const K* current = &coords[0];
    for (int i = 0; i < coord_count; ++i)
    {
      fos << *current;
      if (dim == 2)
        fos << " " << *(current+1) << " 0 \n" << *current << " " << *(current+1) << " 0.01" << std::endl;
      else                   // dim == 3
        fos << " " << *(current+1) << " "  << *(current+2) << std::endl;
      // move pointer
      current += dim;
    }
  }

  void writePolygons(const std::vector<unsigned int>& indices, const int* corners, int ncorners, int dim, ofstream& fos)
  {
    if (dim == 2)
    {
      fos << "POLYGONS " << indices.size()/2 << " " << 5*(indices.size() / 2) << endl;
      for (unsigned int i = 0; i < indices.size(); i += 2)
        fos << "4 " << 2*indices[i] << " " << 2*indices[i+1] << " " << 2*indices[i+1]+1 << " "<< 2*indices[i]+1 << std::endl;

      // arbitrary shapes - ignored here!
      //			int sum = ncorners;
      //			for (int i = 0; i < ncorners; ++i)
      //				sum += (corners[i] > 2 ? corners[i] : 3);
      //
      //			fos << "POLYGONS " << ncorners << " " << sum << endl;
      //			int index = 0;
      //			for (int i = 0; i < ncorners; ++i)
      //			{
      //				// write the first index twice if it is an egde
      //				// => triangle instead of edge - paraview can display it then
      //				if (corners[i] > 2)
      //					fos << corners[i];
      //				else
      //					fos << "3 " << indices[index];
      //
      //				for (int j = 0; j < corners[i]; ++j)
      //					fos << " " << indices[index++];
      //				fos << endl;
      //			}
    }
    else
    {
      int sum = ncorners;
      for (int i = 0; i < ncorners; ++i)
        sum += corners[i];
      fos << "POLYGONS " << ncorners << " " << sum << endl;
      int index = 0;
      for (int i = 0; i < ncorners; ++i)
      {
        fos << corners[i];
        for (int j = 0; j < corners[i]; ++j)
          fos << " " << indices[index++];
        fos << endl;
      }
    }
  }

  template<typename T>
  void writeCellData(const std::vector<T>& data, const char* dataname, int dim, ofstream& fos)
  {
    fos << "CELL_DATA " << data.size()*(dim == 2 ? 2 : 1) << endl;
    fos << "SCALARS " << dataname << " " << TypeNames[Nametraits<T>::nameidx] << " 1" << endl;
    fos << "LOOKUP_TABLE default" << endl;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      fos << data[i] << endl;
      if (dim == 2)
        fos << data[i] << endl;
    }
  }

  template<typename T>
  void writePointData(const std::vector<T>& data, const char* dataname, int dim, ofstream& fos)
  {
    fos << "POINT_DATA " << data.size()*(dim == 2 ? 2 : 1) << endl;
    fos << "SCALARS " << dataname << " " << TypeNames[Nametraits<T>::nameidx] << " 1" << endl;
    fos << "LOOKUP_TABLE default" << endl;
    for (unsigned int i = 0; i < data.size(); ++i)
    {
      fos << data[i] << endl;
      if (dim == 2)
        fos << data[i] << endl;
    }
  }


private:
  const char*  _filename;
};


#endif // VTKSURFACEWRITER_HH_
