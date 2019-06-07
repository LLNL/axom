// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MIRMesh.hpp"

namespace axom
{
namespace mir
{
  MIRMesh::MIRMesh()
  {

  }

  MIRMesh::MIRMesh(MIRMesh* _mesh)  // copy constructor
  {
      evInds = _mesh->evInds;
      evBegins = _mesh->evBegins;
      veInds = _mesh->veInds;
      veBegins = _mesh->veBegins;
      verts = _mesh->verts;
      elems = _mesh->elems;
      bdry = _mesh->bdry;
      cobdry = _mesh->cobdry;
      vertexPositions = _mesh->vertexPositions;
      materialVolumeFractionsElement = _mesh->materialVolumeFractionsElement;
      materialVolumeFractionsVertex = _mesh->materialVolumeFractionsVertex;
      numMaterials = _mesh->numMaterials;
  }

  MIRMesh::~MIRMesh()
  {

  }

  /// Initializes a mesh with the given topology.
  void MIRMesh::InitializeMesh(std::vector<PosType> _evInds, std::vector<PosType> _evBegins, std::vector<PosType> _veInds, std::vector<PosType> _veBegins, VertSet _verts, ElemSet _elems, int _numMaterials)
  {
    evInds = _evInds;
    evBegins = _evBegins;
    veInds = _veInds;
    veBegins = _veBegins;
    verts = _verts;
    elems = _elems;
    numMaterials = _numMaterials;

    // Check validity of the sets
    SLIC_ASSERT_MSG( verts.isValid(), "Vertex set is not valid.");
    SLIC_ASSERT_MSG( elems.isValid(), "Element set is not valid.");
  }

//--------------------------------------------------------------------------------

  /// Constructs the mesh boundary and coboundary relations
  void MIRMesh::constructMeshRelations()
  {
    // construct boundary relation from elements to vertices using variable cardinality
    {
      using RelationBuilder = ElemToVertRelation::RelationBuilder;
      bdry = RelationBuilder()
              .fromSet( &elems )
              .toSet( &verts )
              .begins( RelationBuilder::BeginsSetBuilder()
                        .size( elems.size() ) 
                        .data( evBegins.data() )  )
              .indices ( RelationBuilder::IndicesSetBuilder() 
                          .size( evInds.size() ) 
                          .data( evInds.data() ) );
    }


    {
      // _quadmesh_example_construct_cobdry_relation_start
      // construct coboundary relation from vertices to elements
      using RelationBuilder = VertToElemRelation::RelationBuilder;
      cobdry = RelationBuilder()
               .fromSet( &verts )
               .toSet( &elems )
               .begins( RelationBuilder::BeginsSetBuilder()
                        .size( verts.size() )
                        .data( veBegins.data() ) )
               .indices( RelationBuilder::IndicesSetBuilder()
                         .size( veInds.size() )
                         .data( veInds.data() ) );
      // _quadmesh_example_construct_cobdry_relation_end
    }

    SLIC_ASSERT_MSG( bdry.isValid(), "Boundary relation is not valid.");
    SLIC_ASSERT_MSG( cobdry.isValid(), "Coboundary relation is not valid.");

    SLIC_INFO("Elem-Vert relation has size " << bdry.totalSize());
    SLIC_INFO("Vert-Elem relation has size " << cobdry.totalSize());
  }

//--------------------------------------------------------------------------------

  /// Constructs the volume fraction maps on the elements and vertices given element volume fraction data.
  void MIRMesh::constructMeshVolumeFractionMaps(std::vector<axom::float64*> materialVolumeFractionsData)
  {
    // Initialize the maps for all of the materials with the input volume fraction data for each material
    for (int matID = 0; matID < materialVolumeFractionsData.size(); ++matID)
    {
      // Initialize the map for the current material
      materialVolumeFractionsElement.push_back(ScalarMap( &elems ));

      // Copy the data for the current material
      for (int eID = 0; eID < elems.size(); ++eID)
      {
        materialVolumeFractionsElement[matID][eID] = materialVolumeFractionsData[matID][eID];
      }

      SLIC_ASSERT_MSG( materialVolumeFractionsElement[matID].isValid(), "Element volume fraction map is not valid.");
    }

    // Initialize the maps for all of the vertex volume fractions
    for (int matID = 0; matID < materialVolumeFractionsData.size(); ++matID)
    {
      // Initialize the new map for the volume fractions
      materialVolumeFractionsVertex.push_back(ScalarMap( &verts ) );

      // Calculate the average volume fraction value for the current vertex for the current material
      for (int vID = 0; vID < verts.size(); ++vID)
      {
        // Compute the per vertex volume fractions for the green material
        axom::float64 sum = 0;
        auto vertexElements = cobdry[vID];

        for (int i = 0; i < vertexElements.size(); ++i)
        {
          auto eID = vertexElements[i];
          sum += materialVolumeFractionsElement[matID][eID];
        }

        materialVolumeFractionsVertex[matID][vID] = sum / vertexElements.size();
      }

      SLIC_ASSERT_MSG( materialVolumeFractionsVertex[matID].isValid(), "Vertex volume fraction map is not valid.");
    }
  }

//--------------------------------------------------------------------------------  

/// Construct the vertex volume fraction data given vertex volume fraction data
void MIRMesh::constructMeshVolumeFractionsVertex(std::vector<std::vector<axom::float64> > vertexVF)
{
  // Initialize the maps for all of the materials with the input volume fraction data for each vertex
  for (int matID = 0; matID < numMaterials; ++matID)
  {
    // Initialize the map for the current material
    materialVolumeFractionsVertex.push_back(ScalarMap( &verts ));

    // Copy the data for the current material
    for (int vID = 0; vID < verts.size(); ++vID)
    {
      materialVolumeFractionsVertex[matID][vID] = vertexVF[matID][vID];
    }

    SLIC_ASSERT_MSG( materialVolumeFractionsVertex[matID].isValid(), "Vertex volume fraction map is not valid.");
  } 
}

//--------------------------------------------------------------------------------  

  /// Constucts the positions map on the vertices
  void MIRMesh::constructVertexPositionMap(Point2* data)
  {
    // construct the position map on the vertices
    vertexPositions = PointMap( &verts );

    for (int vID = 0; vID < verts.size(); ++vID)
      vertexPositions[vID] = data[vID];

    SLIC_ASSERT_MSG( vertexPositions.isValid(), "Position map is not valid.");
  }

//--------------------------------------------------------------------------------

  /// Prints out the map values for each element
  void MIRMesh::printElementScalarMap(ScalarMap& elements, std::string prefix)
  {
    std::cout << prefix;
    for (int eID = 0; eID < elems.size(); ++eID)
    {
      printf("Element %d: %f\n", eID, elements[eID]);
    }
  }

//--------------------------------------------------------------------------------

  /// Prints out the map values for each vertex
  void MIRMesh::printVertexScalarMap(ScalarMap& vertices, std::string prefix)
  {
    std::cout << prefix;
    for (int vID = 0; vID < verts.size(); ++vID)
    {
      printf("Vertex %d: %f\n", vID, vertices[vID]);
    }
  }

//--------------------------------------------------------------------------------

  void MIRMesh::print()
  {
    printf("------------------------Printing Mesh Information:------------------------\n");
    printf("number of vertices: %d\n", verts.size());
    printf("number of elements: %d\n", elems.size());
    printf("number of materials: %d\n", numMaterials);

    printf("evInds: { ");
    for (int i = 0; i < evInds.size(); i++)
    {
      printf("%d ", evInds[i]);
    }
    printf("}\n");

    printf("evBegins: { ");
    for (int i = 0; i < evBegins.size(); i++)
    {
      printf("%d ", evBegins[i]);
    }
    printf("}\n");

    printf("veInds: { ");
    for (int i = 0; i < veInds.size(); i++)
    {
      printf("%d ", veInds[i]);
    }
    printf("}\n");

    printf("veBegins: { ");
    for (int i = 0; i < veBegins.size(); i++)
    {
      printf("%d ", veBegins[i]);
    }
    printf("}\n");

    printf("vertexPositions: { ");
    for (int i = 0; i < vertexPositions.size(); ++i)
    {
      printf("{%.2f, %.2f} ", vertexPositions[i].m_x, vertexPositions[i].m_y);
    }

    printf("}\n");
    printf("--------------------------------------------------------------------------\n");
  }

//--------------------------------------------------------------------------------

  /// Reads in a constructs a mesh from the given file
  void MIRMesh::readMeshFromFile()
  {
    printf("Mesh writing functionality not implemented yet.");
    
  }

//--------------------------------------------------------------------------------

  /// Writes out the mesh to a file
  void MIRMesh::writeMeshToFile(std::string filename)
  {
    std::ofstream meshfile;
    meshfile.open(filename);
    std::ostream_iterator<PosType> out_it(meshfile, " ");

    // write header
    meshfile << "# vtk DataFile Version 3.0\n"
             << "vtk output\n"
             << "ASCII\n"
             << "DATASET UNSTRUCTURED_GRID\n\n"
             << "POINTS " << verts.size() << " double\n";

    // write positions
    for (int vID = 0; vID < verts.size(); ++vID)
    {      
      meshfile << vertexPositions[vID].m_x << " " << vertexPositions[vID].m_y << " 0\n"; // must always set all 3 coords; set Z=0 for 2D
    }

    // // write elem-to-vert boundary relation
    // meshfile << "\nCELLS " << elems.size() << " " << 5 * elems.size();  // TODO: This will not always be 5
    // for(auto e: elems)
    // {
    //   meshfile<<"\n4 ";                                                 // TODO: This will not always be 4
    //   std::copy ( bdry.begin(e), bdry.end(e), out_it );
    // }

    // // write element types ( 9 == VTK_QUAD, 5 == VTK_TRIANGLE )
    // meshfile << "\n\nCELL_TYPES " << elems.size() << "\n";
    // for(int i=0 ; i< elems.size() ; ++i)
    // {
    //   meshfile << "9 ";                                                 // TODO: This will not always be 9, but will be 5 for any triangle elements
    // }    

    // // write element ids
    // meshfile << "\n\nCELL_DATA " << elems.size()
    //          << "\nSCALARS cellIds int 1"
    //          << "\nLOOKUP_TABLE default \n";
    // for(int i=0 ; i< elems.size() ; ++i)
    // {
    //   meshfile << elems[i] <<" ";
    // }

    // // write vertex ids
    // meshfile << "\n\nPOINT_DATA " << verts.size()
    //          << "\nSCALARS vertIds int 1"
    //          << "\nLOOKUP_TABLE default \n";
    // for(int i=0 ; i< verts.size() ; ++i)
    // {
    //   meshfile << verts[i] <<" ";
    // }
    meshfile <<"\n";
  }

//--------------------------------------------------------------------------------

}
}