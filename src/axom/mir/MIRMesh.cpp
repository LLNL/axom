// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "MIRMesh.hpp"

namespace axom
{
namespace mir
{

//--------------------------------------------------------------------------------

MIRMesh::MIRMesh()
{

}

//--------------------------------------------------------------------------------

MIRMesh::MIRMesh(MIRMesh* _mesh)
{
    m_meshTopology.m_evInds = _mesh->m_meshTopology.m_evInds;
    m_meshTopology.m_evBegins = _mesh->m_meshTopology.m_evBegins;
    m_meshTopology.m_veInds = _mesh->m_meshTopology.m_veInds;
    m_meshTopology.m_veBegins = _mesh->m_meshTopology.m_veBegins;
    m_verts = _mesh->m_verts;
    m_elems = _mesh->m_elems;
    m_bdry = _mesh->m_bdry;
    m_cobdry = _mesh->m_cobdry;
    m_vertexPositions = _mesh->m_vertexPositions;
    m_materialVolumeFractionsElement = _mesh->m_materialVolumeFractionsElement;
    m_materialVolumeFractionsVertex = _mesh->m_materialVolumeFractionsVertex;
    m_elementParentIDs = _mesh->m_elementParentIDs;
    m_elementDominantMaterials = _mesh->m_elementDominantMaterials;
    m_numMaterials = _mesh->m_numMaterials;
}

//--------------------------------------------------------------------------------

MIRMesh::~MIRMesh()
{

}

//--------------------------------------------------------------------------------

void MIRMesh::initializeMesh(VertSet _verts, ElemSet _elems, int _numMaterials, CellTopologyData _topology, CellMapData _mapData, std::vector<std::vector<axom::float64> > _elementVF)
{
  // Initialize the vertex and element sets
  m_verts = _verts;
  m_elems = _elems;

  m_numMaterials = _numMaterials;

  // Intialize the mesh topology
  m_meshTopology.m_evInds = _topology.m_evInds;
  m_meshTopology.m_evBegins = _topology.m_evBegins;
  m_meshTopology.m_veInds = _topology.m_veInds;
  m_meshTopology.m_veBegins = _topology.m_veBegins;

  // Initialize the mesh relations
  constructMeshRelations();

  // Initialize the mesh's data maps
  constructVertexPositionMap(_mapData.m_vertexPositions.data());
  constructElementParentMap(_mapData.m_elementParents.data());
  constructElementDominantMaterialMap(_mapData.m_elementDominantMaterials);

  // Initialize the element and vertex volume fraction maps
  if (_elementVF.size() > 0)
    constructMeshVolumeFractionsMaps(_elementVF);

  // Check validity of the sets
  SLIC_ASSERT_MSG( m_verts.isValid(), "Vertex set is not valid.");
  SLIC_ASSERT_MSG( m_elems.isValid(), "Element set is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::constructMeshRelations()
{
  // construct boundary relation from elements to vertices using variable cardinality
  {
    using RelationBuilder = ElemToVertRelation::RelationBuilder;
    m_bdry = RelationBuilder()
            .fromSet( &m_elems )
            .toSet( &m_verts )
            .begins( RelationBuilder::BeginsSetBuilder()
                      .size( m_elems.size() ) 
                      .data( m_meshTopology.m_evBegins.data() )  )
            .indices ( RelationBuilder::IndicesSetBuilder() 
                        .size( m_meshTopology.m_evInds.size() ) 
                        .data( m_meshTopology.m_evInds.data() ) );
  }


  {
    // _quadmesh_example_construct_cobdry_relation_start
    // construct coboundary relation from vertices to elements
    using RelationBuilder = VertToElemRelation::RelationBuilder;
    m_cobdry = RelationBuilder()
              .fromSet( &m_verts )
              .toSet( &m_elems )
              .begins( RelationBuilder::BeginsSetBuilder()
                      .size( m_verts.size() )
                      .data( m_meshTopology.m_veBegins.data() ) )
              .indices( RelationBuilder::IndicesSetBuilder()
                        .size( m_meshTopology.m_veInds.size() )
                        .data( m_meshTopology.m_veInds.data() ) );
    // _quadmesh_example_construct_cobdry_relation_end
  }

  SLIC_ASSERT_MSG( m_bdry.isValid(), "Boundary relation is not valid.");
  SLIC_ASSERT_MSG( m_cobdry.isValid(), "Coboundary relation is not valid.");

  SLIC_INFO("Elem-Vert relation has size " << m_bdry.totalSize());
  SLIC_INFO("Vert-Elem relation has size " << m_cobdry.totalSize());
}

//--------------------------------------------------------------------------------  

void MIRMesh::constructMeshVolumeFractionsMaps(std::vector<std::vector<axom::float64> > elementVF)
{
  // Clear the old maps
  m_materialVolumeFractionsElement.clear();
  m_materialVolumeFractionsVertex.clear();

  // Initialize the maps for all of the materials with the input volume fraction data for each material
  for (int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the map for the current material
    m_materialVolumeFractionsElement.push_back(ScalarMap( &m_elems ));

    // Copy the data for the current material
    for (int eID = 0; eID < m_elems.size(); ++eID)
    {
      m_materialVolumeFractionsElement[matID][eID] = elementVF[matID][eID];
    }

    SLIC_ASSERT_MSG( m_materialVolumeFractionsElement[matID].isValid(), "Element volume fraction map is not valid.");
  }

  // Initialize the maps for all of the vertex volume fractions
  for (int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the new map for the volume fractions
    m_materialVolumeFractionsVertex.push_back(ScalarMap( &m_verts ) );

    // Calculate the average volume fraction value for the current vertex for the current material
    for (int vID = 0; vID < m_verts.size(); ++vID)
    {
      // Compute the per vertex volume fractions for the green material
      axom::float64 sum = 0;
      auto vertexElements = m_cobdry[vID];

      for (int i = 0; i < vertexElements.size(); ++i)
      {
        auto eID = vertexElements[i];
        sum += m_materialVolumeFractionsElement[matID][eID];
      }

      m_materialVolumeFractionsVertex[matID][vID] = sum / vertexElements.size();
    }

    SLIC_ASSERT_MSG( m_materialVolumeFractionsVertex[matID].isValid(), "Vertex volume fraction map is not valid.");
  }
}

//--------------------------------------------------------------------------------  

void MIRMesh::constructMeshVolumeFractionsVertex(std::vector<std::vector<axom::float64> > vertexVF)
{
  // Initialize the maps for all of the materials with the input volume fraction data for each vertex
  for (int matID = 0; matID < m_numMaterials; ++matID)
  {
    // Initialize the map for the current material
    m_materialVolumeFractionsVertex.push_back(ScalarMap( &m_verts ));

    // Copy the data for the current material
    for (int vID = 0; vID < m_verts.size(); ++vID)
    {
      m_materialVolumeFractionsVertex[matID][vID] = vertexVF[matID][vID];
    }

    SLIC_ASSERT_MSG( m_materialVolumeFractionsVertex[matID].isValid(), "Vertex volume fraction map is not valid.");
  } 
}

//--------------------------------------------------------------------------------  

void MIRMesh::constructVertexPositionMap(Point2* data)
{
  // construct the position map on the vertices
  m_vertexPositions = PointMap( &m_verts );

  for (int vID = 0; vID < m_verts.size(); ++vID)
    m_vertexPositions[vID] = data[vID];

  SLIC_ASSERT_MSG( m_vertexPositions.isValid(), "Position map is not valid.");
}

//--------------------------------------------------------------------------------  

void MIRMesh::constructElementParentMap(int* elementParents)
{
  // Initialize the map for the elements' parent IDs
  m_elementParentIDs = IntMap( &m_elems );

  // Copy the data for the elements
  for (int eID = 0; eID < m_elems.size(); ++eID)
    m_elementParentIDs[eID] = elementParents[eID];

  SLIC_ASSERT_MSG( m_elementParentIDs.isValid(), "Element parent map is not valid.");
}

//--------------------------------------------------------------------------------  

void MIRMesh::constructElementDominantMaterialMap(std::vector<int> dominantMaterials)
{
  // Initialize the map for the elements' dominant colors
  m_elementDominantMaterials = IntMap( &m_elems );

  // Copy the dat for the elements
  for (int eID = 0; eID < m_elems.size(); ++eID)
    m_elementDominantMaterials[eID] = dominantMaterials[eID];

  SLIC_ASSERT_MSG( m_elementDominantMaterials.isValid(), "Element dominant materials map is not valid.");
}

//--------------------------------------------------------------------------------

void MIRMesh::print()
{
  printf("\n------------------------Printing Mesh Information:------------------------\n");
  printf("number of vertices: %d\n", m_verts.size());
  printf("number of elements: %d\n", m_elems.size());
  printf("number of materials: %d\n", m_numMaterials);

  printf("evInds: { ");
  for (unsigned long i = 0; i < m_meshTopology.m_evInds.size(); i++)
  {
    printf("%d ", m_meshTopology.m_evInds[i]);
  }
  printf("}\n");

  printf("evBegins: { ");
  for (unsigned long i = 0; i < m_meshTopology.m_evBegins.size(); i++)
  {
    printf("%d ", m_meshTopology.m_evBegins[i]);
  }
  printf("}\n");

  printf("veInds: { ");
  for (unsigned long i = 0; i < m_meshTopology.m_veInds.size(); i++)
  {
    printf("%d ", m_meshTopology.m_veInds[i]);
  }
  printf("}\n");

  printf("veBegins: { ");
  for (unsigned long i = 0; i < m_meshTopology.m_veBegins.size(); i++)
  {
    printf("%d ", m_meshTopology.m_veBegins[i]);
  }
  printf("}\n");

  printf("vertexPositions: { ");
  for (int i = 0; i < m_verts.size(); ++i)
  {
    printf("{%.2f, %.2f} ", m_vertexPositions[i].m_x, m_vertexPositions[i].m_y);
  }
  printf("}\n");

  printf("elementParentIDs: { ");
  for (int i = 0; i < m_elems.size(); ++i)
  {
    printf("%d ", m_elementParentIDs[i]);
  }
  printf("}\n");

  printf("elementDominantMaterials: { ");
  for (int i = 0; i < m_elems.size(); ++i)
  {
    printf("%d ", m_elementDominantMaterials[i]);
  }
  printf("}\n");

  printf("vertexVolumeFractions: { \n");
  for (unsigned long i = 0; i < m_materialVolumeFractionsVertex.size(); ++i)
  {
    printf("  { ");
    for (int j = 0; j < m_verts.size(); ++j)
    {
      printf("%.3f, ", m_materialVolumeFractionsVertex[i][j]);
    }
    printf("}\n");
  }
  printf("}\n");
  printf("--------------------------------------------------------------------------\n");
}

//--------------------------------------------------------------------------------

void MIRMesh::readMeshFromFile(std::string filename)
{
  printf("Mesh reading functionality not implemented yet. Can't read file: %s", filename.c_str());
  
  // Read in header

  // Read in POINTS

  // Read in CELLS

  // Read in CELL_TYPES

  // Read in CELL_DATA (element volume fractions)
}

//--------------------------------------------------------------------------------

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
            << "POINTS " << m_verts.size() << " double\n";

  // write positions
  for (int vID = 0; vID < m_verts.size(); ++vID)
  {      
    meshfile << m_vertexPositions[vID].m_x << " " << m_vertexPositions[vID].m_y << " 0\n"; // must always set all 3 coords; set Z=0 for 2D
  }

  meshfile << "\nCELLS " << m_elems.size() << " " << m_meshTopology.m_evInds.size() + m_elems.size();
  for (int i = 0; i < m_elems.size(); ++i)
  {
    int nVerts = m_meshTopology.m_evBegins[i+1] - m_meshTopology.m_evBegins[i];
    meshfile << "\n" << nVerts;
    for (int j = 0; j < nVerts; ++j)
    {
      int startIndex = m_meshTopology.m_evBegins[i];
      meshfile << " " << m_meshTopology.m_evInds[startIndex + j];
    }
  } 

  meshfile << "\n\nCELL_TYPES " << m_elems.size() << "\n";
  for (int i = 0; i < m_elems.size(); ++i)
  {
    int nVerts = m_meshTopology.m_evBegins[i + 1] - m_meshTopology.m_evBegins[i];
    if (nVerts == 3)
      meshfile << "5\n";
    else if (nVerts == 4)
      meshfile << "9\n";
  }

  // write element materials
  meshfile << "\n\nCELL_DATA " << m_elems.size()
            << "\nSCALARS cellIds int 1"
            << "\nLOOKUP_TABLE default \n";
  for(int i=0 ; i< m_elems.size() ; ++i)
  {
    meshfile << m_elementDominantMaterials[i] << " ";
  }

  meshfile <<"\n";
}

//--------------------------------------------------------------------------------

std::vector<std::vector<axom::float64> > MIRMesh::computeOriginalElementVolumeFractions()
{
  std::map<int, axom::float64> totalAreaOriginalElements; // the total area of the original elements
  std::map<int, axom::float64> newElementAreas;           // the area of each of the generated child elements
  std::map<int, int> numChildren;                         // the number of child elements generated from one of the original elements

  // Compute the total area of each element of the original mesh and also the area of each new element
  for (int eID = 0; eID < m_elems.size(); ++eID)
  {
    int numVertices = m_meshTopology.m_evBegins[eID + 1] - m_meshTopology.m_evBegins[eID];
    if (numVertices == 3)
    {
      Point2 trianglePoints[3];
      for (int i = 0; i < 3; ++i)
        trianglePoints[i] = m_vertexPositions[m_meshTopology.m_evInds[ m_meshTopology.m_evBegins[eID] + i] ];
      newElementAreas[eID] = computeTriangleArea(trianglePoints[0], trianglePoints[1], trianglePoints[2]);
    }
    if (numVertices == 4) 
    {
      Point2 trianglePoints[4];
      for (int i = 0; i < 4; ++i)
        trianglePoints[i] = m_vertexPositions[m_meshTopology.m_evInds[ m_meshTopology.m_evBegins[eID] + i] ];
      newElementAreas[eID] = computeQuadArea(trianglePoints[0], trianglePoints[1], trianglePoints[2], trianglePoints[3]);
    }
   
    totalAreaOriginalElements[ m_elementParentIDs[eID] ] += newElementAreas[eID];
    numChildren[ m_elementParentIDs[eID] ] += .4;
  }
  
  // Intialize the element volume fraction vectors
  std::vector<std::vector<axom::float64> > elementVolumeFractions;    // indexed as: elementVolumeFractions[material][originalElementID] = volumeFraction
  elementVolumeFractions.resize( m_numMaterials );
  for (unsigned long i = 0; i < elementVolumeFractions.size(); ++i)
    elementVolumeFractions[i].resize( totalAreaOriginalElements.size(), 0.0 );

  // Compute the volume fractions for each of the original mesh elements
  for (auto itr = newElementAreas.begin(); itr != newElementAreas.end(); itr++)
  {
    int parentElementID = m_elementParentIDs[itr->first];
    int materialID = m_elementDominantMaterials[itr->first];
    elementVolumeFractions[materialID][parentElementID] += (itr->second / totalAreaOriginalElements[parentElementID]);
  }

  return elementVolumeFractions;
}

//--------------------------------------------------------------------------------

axom::float64 MIRMesh::computeTriangleArea(Point2 p0, Point2 p1, Point2 p2)
{
  axom::float64 a = sqrt( ((p1.m_x - p0.m_x) * (p1.m_x - p0.m_x)) + ((p1.m_y - p0.m_y) * (p1.m_y - p0.m_y)) );// the distance from p0 to p1
  axom::float64 b = sqrt( ((p2.m_x - p1.m_x) * (p2.m_x - p1.m_x)) + ((p2.m_y - p1.m_y) * (p2.m_y - p1.m_y)) );// the distance from p1 to p2
  axom::float64 c = sqrt( ((p0.m_x - p2.m_x) * (p0.m_x - p2.m_x)) + ((p0.m_y - p2.m_y) * (p0.m_y - p2.m_y)) );// the distance from p2 to p0

  axom::float64 s = (a + b + c) / 2;  // the semi-perimeter of the triangle

  return sqrt(s * (s - a) * (s - b) * (s - c));
}

//--------------------------------------------------------------------------------

axom::float64 MIRMesh::computeQuadArea(Point2 p0, Point2 p1, Point2 p2, Point2 p3)
{
  return computeTriangleArea(p0, p1, p2) + computeTriangleArea(p2, p3, p0);
}

//--------------------------------------------------------------------------------

}
}