// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __MIR_MESH_H__
#define __MIR_MESH_H__


#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "MIRMeshTypes.hpp"

// C/C++ includes
#include <cmath>          // for definition of M_PI, exp()
#include <vector>
#include <string>
#include <fstream>

namespace numerics = axom::numerics;
namespace slam = axom::slam;


//--------------------------------------------------------------------------------
namespace axom
{
namespace mir
{
  struct EdgeClipInfo
  {
    int vertexOne;
    int vertexTwo;
    int colorOne;
    int colorTwo;
    float t;        // The percent distance from vertexOne to VertexTwo where the clip occurs.
  };

  //--------------------------------------------------------------------------------

  class MIRMesh
  {
    /****************************************************************
     *                       MESH FUNCTIONS
     ****************************************************************/
    public:
      MIRMesh();
      MIRMesh(MIRMesh* _mesh);   // copy constructor
      ~MIRMesh();

      void        InitializeMesh(std::vector<PosType> _evInds, std::vector<PosType> _evBegins, std::vector<PosType> _veInds, std::vector<PosType> _veBegins, VertSet _verts, ElemSet _elems, int _numMaterials);

      void        constructMeshRelations();                                                 /// Constructs the mesh boundary and coboundary relations   
      void        constructMeshVolumeFractionMaps(std::vector<axom::float64*> materialVolumeFractionsData); /// Constructs the volume fraction maps on the vertices
      void        constructMeshVolumeFractionsVertex(std::vector<std::vector<axom::float64> > vertexVF);  /// Constructs the vertex volume fraction map
      void        constructVertexPositionMap(Point2* data);                                /// Constucts the positions map on the vertices

      void        printElementScalarMap(ScalarMap& elements, std::string prefix);           /// Prints out the map values for each element
      void        printVertexScalarMap(ScalarMap& vertices, std::string prefix);            /// Prints out the map values for each vertex

      void        print();                                                                  // Print out a variety of useful information about the mesh

      void        readMeshFromFile();                                                       /// Reads in and constructs a mesh from a file
      void        writeMeshToFile(std::string filename);                                    /// Writes out the current mesh to a file
      void        attachVertexMap();                                                        /// Adds a map of data to the mesh

    /****************************************************************
     *                        VARIABLES
     ****************************************************************/
    public:
      // support data for mesh connectivity
      std::vector<PosType> evInds;
      std::vector<PosType> evBegins;
      std::vector<PosType> veInds;
      std::vector<PosType> veBegins;

    public:
      // Mesh Set Definitions
      VertSet verts;  // the set of vertices in the mesh
      ElemSet elems;  // the set of elements in the mesh

    public:
      // Mesh Relation Definitions
      ElemToVertRelation bdry;      // Boundary relation from elements to vertices
      VertToElemRelation cobdry;   // Coboundary relation from vertices to elements

    public:
      // Mesh Map Definitions
      PointMap vertexPositions;   // vertex position
      std::vector<ScalarMap> materialVolumeFractionsElement;  // the volume fractions of each material for each element
      std::vector<ScalarMap> materialVolumeFractionsVertex;   // the volume fractions of each material for each vertex

    public:
      int numMaterials;
  };

  //--------------------------------------------------------------------------------
}
}
#endif