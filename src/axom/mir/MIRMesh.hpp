// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef __MIR_MESH_H__
#define __MIR_MESH_H__


#include "axom/core.hpp"  // for axom macros
#include "axom/slam.hpp"  // unified header for slam classes and functions

#include "MIRMeshTypes.hpp"
#include "CellData.hpp"

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

  #define NULL_MAT -1

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
      void        constructMeshVolumeFractionsMaps(std::vector<std::vector<axom::float64> > elementVF);  /// Constructs the element and vertex volume fraction maps
      void        constructMeshVolumeFractionsVertex(std::vector<std::vector<axom::float64> > vertexVF);  /// Constructs the vertex volume fraction map
      void        constructVertexPositionMap(Point2* data);                                /// Constucts the positions map on the vertices
      void        constructElementParentMap(int* cellParents);                                /// Constructs the map of cells to their original cell parent
      void        constructElementDominantMaterialMap(std::vector<int> dominantMaterials);                      /// Constructs the map of elements to their dominant materials

      void        print();                                                                  // Print out a variety of useful information about the mesh

      void        readMeshFromFile(std::string filename);                                   /// Reads in and constructs a mesh from a file
      void        writeMeshToFile(std::string filename);                                    /// Writes out the current mesh to a file
      void        attachVertexMap();                                                        /// Adds a map of data to the mesh

      std::vector<std::vector<axom::float64> >        computeOriginalElementVolumeFractions();
    
    private:
      axom::float64 computeTriangleArea(Point2 p0, Point2 p1, Point2 p2);
      axom::float64 computeQuadArea(Point2 p0, Point2 p1, Point2 p2, Point2 p3);

    /****************************************************************
     *                        VARIABLES
     ****************************************************************/
    public:
      // Mesh Set Definitions
      VertSet verts;  // the set of vertices in the mesh
      ElemSet elems;  // the set of elements in the mesh

    public:
      // Mesh Relation Definitions
      ElemToVertRelation bdry;      // Boundary relation from elements to vertices
      VertToElemRelation cobdry;    // Coboundary relation from vertices to elements

    public:
      // Mesh Map Definitions
      PointMap vertexPositions;                               // vertex position for each vertex
      std::vector<ScalarMap> materialVolumeFractionsElement;  // the volume fractions of each material for each element
      std::vector<ScalarMap> materialVolumeFractionsVertex;   // the volume fractions of each material for each vertex

    public:
      int numMaterials;

    public:
      IntMap elementParentIDs;              // the ID of the parent element from the original mesh
      IntMap elementDominantMaterials;      // the dominant material of the cell (a processed mesh should have only one material per cell)

      public:
      CellData data;                        // Contains mesh connectivity data, volume fraction data, vertex position data
  };

  //--------------------------------------------------------------------------------
}
}
#endif