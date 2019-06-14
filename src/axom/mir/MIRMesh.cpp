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
      data.evInds = _mesh->data.evInds;
      data.evBegins = _mesh->data.evBegins;
      data.veInds = _mesh->data.veInds;
      data.veBegins = _mesh->data.veBegins;
      verts = _mesh->verts;
      elems = _mesh->elems;
      bdry = _mesh->bdry;
      cobdry = _mesh->cobdry;
      vertexPositions = _mesh->vertexPositions;
      materialVolumeFractionsElement = _mesh->materialVolumeFractionsElement;
      materialVolumeFractionsVertex = _mesh->materialVolumeFractionsVertex;
      elementParentIDs = _mesh->elementParentIDs;
      elementDominantMaterials = _mesh->elementDominantMaterials;
      numMaterials = _mesh->numMaterials;
  }

  MIRMesh::~MIRMesh()
  {

  }

  /// Initializes a mesh with the given topology.
  void MIRMesh::InitializeMesh(std::vector<PosType> _evInds, std::vector<PosType> _evBegins, std::vector<PosType> _veInds, std::vector<PosType> _veBegins, VertSet _verts, ElemSet _elems, int _numMaterials)
  {
    data.evInds = _evInds;
    data.evBegins = _evBegins;
    data.veInds = _veInds;
    data.veBegins = _veBegins;
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
                        .data( data.evBegins.data() )  )
              .indices ( RelationBuilder::IndicesSetBuilder() 
                          .size( data.evInds.size() ) 
                          .data( data.evInds.data() ) );
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
                        .data( data.veBegins.data() ) )
               .indices( RelationBuilder::IndicesSetBuilder()
                         .size( data.veInds.size() )
                         .data( data.veInds.data() ) );
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

  /// Constructs the elementParentsID map of the ID of the parent element in the original mesh for each generated element
  void MIRMesh::constructElementParentMap(int* elementParents)
  {
    // Initialize the map for the elements' parent IDs
    elementParentIDs = IntMap( &elems );

    // Copy the data for the elements
    for (int eID = 0; eID < elems.size(); ++eID)
      elementParentIDs[eID] = elementParents[eID];

    SLIC_ASSERT_MSG( elementParentIDs.isValid(), "Element parent map is not valid.");
  }

//--------------------------------------------------------------------------------  

  /// Constructs the elementDominantMaterials map of each element's single most dominant material
  void MIRMesh::constructElementDominantMaterialMap(std::vector<int> dominantMaterials)
  {
    // Initialize the map for the elements' dominant colors
    elementDominantMaterials = IntMap( &elems );

    // Copy the dat for the elements
    for (int eID = 0; eID < elems.size(); ++eID)
      elementDominantMaterials[eID] = dominantMaterials[eID];

    SLIC_ASSERT_MSG( elementDominantMaterials.isValid(), "Element dominant materials map is not valid.");
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
    printf("\n------------------------Printing Mesh Information:------------------------\n");
    printf("number of vertices: %d\n", verts.size());
    printf("number of elements: %d\n", elems.size());
    printf("number of materials: %d\n", numMaterials);

    printf("evInds: { ");
    for (int i = 0; i < data.evInds.size(); i++)
    {
      printf("%d ", data.evInds[i]);
    }
    printf("}\n");

    printf("evBegins: { ");
    for (int i = 0; i < data.evBegins.size(); i++)
    {
      printf("%d ", data.evBegins[i]);
    }
    printf("}\n");

    printf("veInds: { ");
    for (int i = 0; i < data.veInds.size(); i++)
    {
      printf("%d ", data.veInds[i]);
    }
    printf("}\n");

    printf("veBegins: { ");
    for (int i = 0; i < data.veBegins.size(); i++)
    {
      printf("%d ", data.veBegins[i]);
    }
    printf("}\n");

    printf("vertexPositions: { ");
    for (int i = 0; i < verts.size(); ++i)      // TODO: This was previously vertexPositions.size() and was working...
    {
      printf("{%.2f, %.2f} ", vertexPositions[i].m_x, vertexPositions[i].m_y);
    }
    printf("}\n");

    printf("elementParentIDs: { ");
    for (int i = 0; i < elems.size(); ++i)
    {
      printf("%d ", elementParentIDs[i]);
    }
    printf("}\n");

    printf("elementDominantMaterials: { ");
    for (int i = 0; i < elems.size(); ++i)
    {
      printf("%d ", elementDominantMaterials[i]);
    }
    printf("}\n");

    printf("vertexVolumeFractions: { \n");
    for (int i = 0; i < materialVolumeFractionsVertex.size(); ++i)
    {
      printf("  { ");
      for (int j = 0; j < verts.size(); ++j)
      {
        printf("%.3f, ", materialVolumeFractionsVertex[i][j]);
      }
      printf("}\n");
    }
    printf("}\n");
    printf("--------------------------------------------------------------------------\n");
  }

//--------------------------------------------------------------------------------

  /// Reads in and constructs a mesh from the given file
  /// Note: Must currently be an ASCII, UNSTRUCTURED_GRID .vtk file
  void MIRMesh::readMeshFromFile(std::string filename)
  {
    printf("Mesh reading functionality not implemented yet.");
    
    // Read in header

    // Read in POINTS

    // Read in CELLS

    // Read in CELL_TYPES

    // Read in CELL_DATA (element volume fractions)
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
             << "DATASET UNSTRUCTURED_GRID\n"
             << "POINTS " << verts.size() << " double\n";

    // write positions
    for (int vID = 0; vID < verts.size(); ++vID)
    {      
      meshfile << vertexPositions[vID].m_x << " " << vertexPositions[vID].m_y << " 0\n"; // must always set all 3 coords; set Z=0 for 2D
    }

    meshfile << "\nCELLS " << elems.size() << " " << data.evInds.size() + elems.size();
    for (int i = 0; i < elems.size(); ++i)
    {
      int nVerts = data.evBegins[i+1] - data.evBegins[i];
      meshfile << "\n" << nVerts;
      for (int j = 0; j < nVerts; ++j)
      {
        int startIndex = data.evBegins[i];
        meshfile << " " << data.evInds[startIndex + j];
      }
    } 

    meshfile << "\n\nCELL_TYPES " << elems.size() << "\n";
    for (int i = 0; i < elems.size(); ++i)
    {
      int nVerts = data.evBegins[i + 1] - data.evBegins[i];
      if (nVerts == 3)
        meshfile << "5\n";
      else if (nVerts == 4)
        meshfile << "9\n";
    }

    // write element materials
    meshfile << "\n\nCELL_DATA " << elems.size()
             << "\nSCALARS cellIds int 1"
             << "\nLOOKUP_TABLE default \n";
    for(int i=0 ; i< elems.size() ; ++i)
    {
      meshfile << elementDominantMaterials[i] << " ";
    }

    meshfile <<"\n";
  }

//--------------------------------------------------------------------------------

/// Computes the volume fractions of the elements of the original mesh,
void MIRMesh::computeOriginalElementVolumeFractions()
{
  std::map<int, axom::float64> totalAreaOriginalElements; // the total area of the original elements
  std::map<int, axom::float64> newElementAreas;           // the area of each of the generated child elements
  std::map<int, int> numChildren;                         // the number of child elements generated from one of the original elements

  // Compute the total area of each element of the original mesh and also the area of each new element
  for (int eID = 0; eID < elems.size(); ++eID)
  {
    int numVertices = data.evBegins[eID + 1] - data.evBegins[eID];
    if (numVertices == 3)
    {
      Point2 trianglePoints[3];
      for (int i = 0; i < 3; ++i)
        trianglePoints[i] = vertexPositions[data.evInds[ data.evBegins[eID] + i] ];
      newElementAreas[eID] = computeTriangleArea(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
    }
    if (numVertices == 4) 
    {
      Point2 trianglePoints[4];
      for (int i = 0; i < 4; ++i)
        trianglePoints[i] = vertexPositions[data.evInds[ data.evBegins[eID] + i] ];
      newElementAreas[eID] = computeQuadArea(trianglePoints[0], trianglePoints[1], trianglePoints[2], trianglePoints[3]);
    }
   
    totalAreaOriginalElements[ elementParentIDs[eID] ] += newElementAreas[eID];
    numChildren[ elementParentIDs[eID] ] += .4;
  }
  
  // Intialize the element volume fraction vectors
  std::vector<std::vector<axom::float64> > elementVolumeFractions;    // indexed as: elementVolumeFractions[material][originalElementID] = volumeFraction
  elementVolumeFractions.resize( numMaterials );
  for (int i = 0; i < elementVolumeFractions.size(); ++i)
    elementVolumeFractions[i].resize( totalAreaOriginalElements.size(), 0.0 );

  // Compute the volume fractions for each of the original mesh elements
  for (auto itr = newElementAreas.begin(); itr != newElementAreas.end(); itr++)
  {
    int parentElementID = elementParentIDs[itr->first];
    int materialID = elementDominantMaterials[itr->first];
    elementVolumeFractions[materialID][parentElementID] += (itr->second / totalAreaOriginalElements[parentElementID]);
  }

  // Print out the results // TODO: Return the values and use them.
  printf("elementVolumeFractions: {\n");
  for (int matID = 0; matID < elementVolumeFractions.size(); ++matID)
  {
    printf("Material %d: {", matID);
    for (int eID = 0; eID < elementVolumeFractions[matID].size(); ++eID)
    {
      printf(" %f,", elementVolumeFractions[matID][eID]);
    }
    printf("}\n");
  }
  printf("}\n");
}

//--------------------------------------------------------------------------------

/// Computes the area of the triangle defined by the given three vertex positions. Uses Heron's formula.
axom::float64 MIRMesh::computeTriangleArea(Point2 p0, Point2 p1, Point2 p2)
{
  axom::float64 a = sqrt( ((p1.m_x - p0.m_x) * (p1.m_x - p0.m_x)) + ((p1.m_y - p0.m_y) * (p1.m_y - p0.m_y)) );// the distance from p0 to p1
  axom::float64 b = sqrt( ((p2.m_x - p1.m_x) * (p2.m_x - p1.m_x)) + ((p2.m_y - p1.m_y) * (p2.m_y - p1.m_y)) );// the distance from p1 to p2
  axom::float64 c = sqrt( ((p0.m_x - p2.m_x) * (p0.m_x - p2.m_x)) + ((p0.m_y - p2.m_y) * (p0.m_y - p2.m_y)) );// the distance from p2 to p0

  axom::float64 s = (a + b + c) / 2;  // the semi-perimeter of the triangle

  return sqrt(s * (s - a) * (s - b) * (s - c));
}

//--------------------------------------------------------------------------------

/// Computes the area of the quad defined by the given four vertex positions.
/// Note: It is assumed the points are given in consecutive order.
axom::float64 MIRMesh::computeQuadArea(Point2 p0, Point2 p1, Point2 p2, Point2 p3)
{
  return computeTriangleArea(p0, p1, p2) + computeTriangleArea(p2, p3, p0);
}

//--------------------------------------------------------------------------------

}
}