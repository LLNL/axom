/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*!
 *******************************************************************************
 * \file
 * \brief Basic demo for the in/out octree.
 *        WIP towards point containment acceleration structure over surface.
 *******************************************************************************
 */

// ATK Toolkit includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "common/FileUtilities.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/Field.hpp"
#include "quest/FieldData.hpp"
#include "quest/FieldVariable.hpp"
#include "quest/Mesh.hpp"
#include "quest/Point.hpp"
#include "quest/STLReader.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/Triangle.hpp"
#include "quest/UniformMesh.hpp"
#include "quest/UnstructuredMesh.hpp"
#include "quest/Point.hpp"
#include "quest/SpatialOctree.hpp"
#include "quest/InOutOctree.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "slam/Utilities.hpp"
#include "slam/RangeSet.hpp"
#include "slam/Map.hpp"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

using namespace asctoolkit;

typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;

typedef quest::InOutOctree<3> Octree3D;

typedef Octree3D::GeometricBoundingBox GeometricBoundingBox;
typedef Octree3D::SpacePt SpacePt;
typedef Octree3D::SpaceVector SpaceVector;
typedef Octree3D::GridPt GridPt;
typedef Octree3D::BlockIndex BlockIndex;


//------------------------------------------------------------------------------
/**
 * \brief Computes the bounding box of the surface mesh
 */
GeometricBoundingBox compute_bounds( meshtk::Mesh* mesh)
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   GeometricBoundingBox meshBB;
   SpacePt pt;

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
   {
       mesh->getMeshNode( i, pt.data() );
       meshBB.addPoint( pt );
   }

   SLIC_ASSERT( meshBB.isValid() );

   return meshBB;
}


//------------------------------------------------------------------------------
void write_vtk( meshtk::Mesh* mesh, const std::string& fileName )
{
  SLIC_ASSERT( mesh != ATK_NULLPTR );

  std::ofstream ofs;
  ofs.open( fileName.c_str() );

  // STEP 0: Write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " Unstructured Mesh (";
  ofs << mesh->getBlockId() << ", " << mesh->getPartitionId() << ")\n";

  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // STEP 1: Write mesh nodes
  const int num_nodes = mesh->getMeshNumberOfNodes();
  ofs << "POINTS " << num_nodes << " double\n";
  for ( int node=0; node < num_nodes; ++node ) {

    ofs << mesh->getMeshNodeCoordinate( node, 0 ) << " ";
    ofs << mesh->getMeshNodeCoordinate( node, 1 ) << " ";
    if ( mesh->getDimension()==3 ) {

        ofs << mesh->getMeshNodeCoordinate( node, 2 );

    } else {

        ofs << "0.0";
    }

    ofs << std::endl;

  } // END for all nodes

  // STEP 2: Write mesh cell connectivity
  // TODO: note currently this does not work with mixed cell types
  const int ncells   = mesh->getMeshNumberOfCells();
  const int maxnodes = mesh->getMeshNumberOfCellNodes( 0 );
  ofs << "CELLS " << ncells << " " << ncells*(maxnodes+1) << std::endl;

  std::vector< int > cell;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {

    const int nnodes = mesh->getMeshNumberOfCellNodes( cellIdx );
    cell.resize( nnodes );
    mesh->getMeshCell( cellIdx, &cell[0] );

    ofs << nnodes << " ";
    for ( int i=0; i < nnodes; ++i ) {
      ofs << cell[ i ] << " ";
    } // END for all nodes
    ofs << std::endl;

  } // END for all cells

  // STEP 3: Write cell types
  ofs << "CELL_TYPES " << ncells << std::endl;
  for ( int cellIdx=0; cellIdx < ncells; ++cellIdx ) {
    int ctype    = mesh->getMeshCellType( cellIdx );
    int vtk_type = meshtk::cell::vtk_types[ ctype ];
    ofs << vtk_type << std::endl;
  } // END for all cells

  // STEP 4: Write Cell Data
  ofs << "CELL_DATA " << ncells << std::endl;
  meshtk::FieldData* CD = mesh->getCellFieldData();
  for ( int f=0; f < CD->getNumberOfFields(); ++f ) {

      meshtk::Field* field = CD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == meshtk::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          SLIC_ASSERT( dataPtr != ATK_NULLPTR );

          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < ncells; ++i ) {
             ofs << dataPtr[ i ] << std::endl;
          }
      }



  }

  // STEP 5: Write Point Data
  const int nnodes = mesh->getMeshNumberOfNodes();
  ofs << "POINT_DATA " << nnodes << std::endl;
  meshtk::FieldData* PD = mesh->getNodeFieldData();
  for ( int f=0; f < PD->getNumberOfFields(); ++f ) {

      meshtk::Field* field = PD->getField( f );

      ofs << "SCALARS " << field->getName() << " ";
      if ( field->getType() == meshtk::DOUBLE_FIELD_TYPE ) {

          double* dataPtr = field->getDoublePtr();
          ofs << "double\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      } else {

          int* dataPtr = field->getIntPtr();
          ofs << "int\n";
          ofs << "LOOKUP_TABLE default\n";
          for (int i=0; i < nnodes; ++i) {
              ofs << dataPtr[ i ] << std::endl;
          }

      }
  }

  ofs.close();
}



/**
 * \brief Computes some statistics about the surface mesh.
 *
 * Specifically, computes histograms (and ranges) of the edge lengths and triangle areas
 * on a lg scale and logs the results
 */
void print_surface_stats( meshtk::Mesh* mesh)
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   SpacePt pt;

   typedef quest::BoundingBox<double,1> MinMaxRange;
   typedef MinMaxRange::PointType LengthType;

   MinMaxRange meshEdgeLenRange;
   MinMaxRange meshTriAreaRange;
   const int nCells = mesh->getMeshNumberOfCells();
   typedef std::set<int> TriIdxSet;
   TriIdxSet badTriangles;

   // simple binning based on the exponent
   typedef std::map<int,int> LogHistogram;
   LogHistogram edgeLenHist;        // Create histogram of edge lengths (log scale)
   LogHistogram areaHist;           // Create histogram of triangle areas (log scale)

   typedef std::map<int,MinMaxRange> LogRangeMap;
   LogRangeMap edgeLenRangeMap;     // Tracks range of edge lengths at each scale
   LogRangeMap areaRangeMap;        // Tracks range of triangle areas at each scale

   typedef quest::Point<int,3> TriVertIndices;

   // Traverse mesh triangles and bin the edge lengths and areas
   for ( int i=0; i < nCells; ++i )
   {
      // Get the indices and positions of the triangle's three vertices
      TriVertIndices vertIndices;
      mesh->getMeshCell( i, vertIndices.data() );

      SpacePt vertPos[3];
      for(int j=0; j<3; ++j)
          mesh->getMeshNode( vertIndices[j], vertPos[j].data() );

      // Compute edge stats -- note edges are double counted
      for(int j=0; j<3; ++j)
      {
          double len = SpaceVector(vertPos[j],vertPos[(j+1)%3]).norm();
          if(len == 0)
          {
              badTriangles.insert(i);
          }
          else
          {
              LengthType edgeLen(len);
              meshEdgeLenRange.addPoint( edgeLen );
              int expBase2;
              std::frexp (len, &expBase2);
              edgeLenHist[expBase2]++;
              edgeLenRangeMap[expBase2].addPoint( edgeLen );
          }
      }

      // Compute triangle area stats
      double area = SpaceVector::cross_product(
                  SpaceVector(vertPos[0],vertPos[1])
                  , SpaceVector(vertPos[0],vertPos[2])
                  ).norm();
      if(area == 0.)
          badTriangles.insert(i);
      else
      {
          LengthType triArea(area);
          meshTriAreaRange.addPoint ( triArea );
          int expBase2;
          std::frexp (area, &expBase2);
          areaHist[expBase2]++;
          areaRangeMap[expBase2].addPoint( triArea);
      }
   }


   // Log the results
   const int nVerts = mesh->getMeshNumberOfNodes();
   SLIC_INFO("Mesh has " << nVerts << " vertices "
             <<"and " << nCells << " triangles.");

   SLIC_INFO("Edge length range is: "  << meshEdgeLenRange);
   SLIC_INFO("Triangle area range is: "  << meshTriAreaRange);

   std::stringstream edgeHistStr;
   edgeHistStr<<"\tEdge length histogram (lg-arithmic): ";
   for(LogHistogram::const_iterator it = edgeLenHist.begin()
           ; it != edgeLenHist.end()
           ; ++it)
   {
       edgeHistStr << "\n\t exp: " << it->first
                   <<"\t count: " << (it->second / 2)
                   <<"\tRange: " << edgeLenRangeMap[it->first];
   }
   SLIC_INFO(edgeHistStr.str());

   std::stringstream triHistStr;
   triHistStr<<"\tTriangle areas histogram (lg-arithmic): ";
   for(LogHistogram::const_iterator it =areaHist.begin()
           ; it != areaHist.end()
           ; ++it)
   {
       triHistStr<<"\n\t exp: " << it->first
                 <<"\t count: " << it->second
                 << "\tRange: " << areaRangeMap[it->first];
   }
   SLIC_INFO(triHistStr.str());

   if(! badTriangles.empty() )
   {
       std::stringstream badTriStr;
       badTriStr<<"The following triangle(s) have zero area/edge lengths:";
       for(TriIdxSet::const_iterator it = badTriangles.begin()
               ; it != badTriangles.end()
               ; ++it)
       {
           badTriStr<< "\n\tTriangle " << *it;
           TriVertIndices vertIndices;
           mesh->getMeshCell( *it, vertIndices.data() );

           SpacePt vertPos;
           for(int j=0; j<3; ++j)
           {
               mesh->getMeshNode( vertIndices[j], vertPos.data() );
               badTriStr<<"\n\t\t vId: " << vertIndices[j] <<" @ position: " << vertPos;
           }
       }
       SLIC_INFO(badTriStr.str());
   }
}

/**
 * Use octree index over mesh vertices to convert the 'triangle soup'
 * from the stl file into an indexed triangle mesh representation.
 * In particular, all vertices that are nearly coincident will be merged,
 * and degenerate triangles (where the three vertices do not have unique indices)
 * will be removed.
 */
void update_trimesh_verts(meshtk::Mesh*& mesh, Octree3D& octree)
{
    typedef int VertexIndex;
    static const VertexIndex NO_VERTEX = quest::InOutLeafData::NO_VERTEX;


    typedef asctoolkit::slam::PositionSet MeshVertsSet;
    typedef asctoolkit::slam::Map<VertexIndex> IndexMap;

    // Create a map from old vertex IDs to new vertex ids
    MeshVertsSet oldVerts( mesh->getMeshNumberOfNodes() );
    IndexMap vertexIndexMap( &oldVerts, NO_VERTEX);

    // Generate unique indices for mesh vertices
    int uniqueVertexCounter = 0;
    for(int i=0; i< oldVerts.size(); ++i)
    {
        // Get the coordinates of the vertex
        SpacePt pt;
        mesh->getMeshNode(i, pt.data());

        // Find the block and its indexed vertex in the octree
        BlockIndex leafBlock = octree.findLeafBlock(pt);
        SLIC_ASSERT( octree[leafBlock].hasVertex() );
        VertexIndex vInd = octree[ leafBlock ].vertexIndex();

        // If the indexed vertex doesn't have a new id, give it one
        if(vertexIndexMap[vInd] == NO_VERTEX)
        {
            vertexIndexMap[vInd] = uniqueVertexCounter++;
        }

        // If this is not the indexed vertex, grab that vertex's new IDX
        if(vInd != i)
        {
            vertexIndexMap[i] = vertexIndexMap[vInd];
        }
    }

    // Find coordinates for the new mesh vertices
    // Cache in temporary SLAM Map since these may be out of order w.r.t. new vertex indices
    typedef asctoolkit::slam::Map<SpacePt> MeshCoordsField;
    MeshVertsSet newVerts( uniqueVertexCounter );
    MeshCoordsField newVertCoords( &newVerts);
    for(int i=0; i< oldVerts.size(); ++i)
    {
        VertexIndex vInd = vertexIndexMap[i];
        mesh->getMeshNode(i, newVertCoords[vInd].data() );
    }

    // Create a new triangles mesh.
    // * add new vertices from octree to the mesh
    // * update vertex references in octree to new vertex indices
    TriangleMesh* newMesh = new TriangleMesh(3);
    for(int i = 0; i< newVerts.size(); ++i)
    {
        const SpacePt& pos = newVertCoords[i];
        newMesh->insertNode(pos[0], pos[1], pos[2]);

        BlockIndex leafBlock = octree.findLeafBlock(pos);
        SLIC_ASSERT( octree.isLeaf(leafBlock) && octree[leafBlock].hasVertex() );
        octree[ leafBlock ].setVertex(i);
    }

    // Add triangles from old mesh to new mesh using updated vertex ids
    typedef quest::Point<int,3> TriVertIndices;
    int numOldMeshTris =  mesh->getMeshNumberOfCells();
    for(int i=0; i< numOldMeshTris ; ++i)
    {
       TriVertIndices vertIndices;
       mesh->getMeshCell( i, vertIndices.data() );

       // Remap the vertex IDs
       for(int j=0; j< 3; ++j)
           vertIndices[j] = vertexIndexMap[ vertIndices[j] ];

       // Skip degenerate triangles -- need 3 unique vertex IDS
       if(    (vertIndices[0] != vertIndices[1])
           && (vertIndices[1] != vertIndices[2])
           && (vertIndices[2] != vertIndices[0]) )
       {
           newMesh->insertCell( vertIndices.data(), meshtk::LINEAR_TRIANGLE, 3);
       }
    }

    // Delete old mesh, redirect pointer to newly created mesh
    delete mesh;
    mesh = newMesh;
}


/**
 * \brief Finds the octree leaf containing the given query point, and optionally refines the leaf
 */
void refineAndPrint(Octree3D& octree, const SpacePt& queryPt, bool shouldRefine = true)
{
    BlockIndex leafBlock = octree.findLeafBlock(queryPt);

    if(shouldRefine)
    {
        octree.refineLeaf( leafBlock );
        leafBlock = octree.findLeafBlock(queryPt);
    }

    GeometricBoundingBox blockBB = octree.blockBoundingBox( leafBlock);
    bool containsPt = blockBB.contains(queryPt);

    SLIC_INFO("\t{gridPt: " << leafBlock.pt()
            <<"; lev: " << leafBlock.level()
            <<"} "
            <<" with bounds " << blockBB
            << (containsPt? " contains " : "does not contain ") << "query point.");
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  slic::UnitTestLogger logger;  // create & initialize logger

  bool hasInputArgs = argc > 1;

  // STEP 1: Get file from user or use default
  std::string stlFile;
  if(hasInputArgs)
  {
      stlFile = std::string( argv[1] );
  }
  else
  {
      const std::string defaultFileName = "plane_simp.stl";
      const std::string defaultDir = "src/components/quest/data/";

      stlFile = asctoolkit::utilities::filesystem::joinPath(defaultDir, defaultFileName);
  }

  stlFile = asctoolkit::slam::util::findFileInAncestorDirs(stlFile);
  SLIC_ASSERT( asctoolkit::utilities::filesystem::pathExists( stlFile));

  // STEP 2: read mesh file
  SLIC_INFO("Reading file: " << stlFile << "...");

  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( stlFile );
  reader->read();
  SLIC_INFO("done.");


  // STEP 3: create surface mesh
  meshtk::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );
  // dump mesh info
  SLIC_INFO("Mesh has "
          << surface_mesh->getMeshNumberOfNodes() << " nodes and "
          << surface_mesh->getMeshNumberOfCells() << " cells.");

  // STEP 4: Delete the reader
  delete reader;
  reader = ATK_NULLPTR;


  // STEP 5: Compute the bounding box and log some stats about the surface
  GeometricBoundingBox meshBB = compute_bounds( surface_mesh);
  SLIC_INFO( "Mesh bounding box: " << meshBB );
  print_surface_stats(surface_mesh);

  // STEP 6: Create octree over mesh's bounding box and query a point in space
  Octree3D octree(meshBB, surface_mesh);
  octree.generateIndex();

  // STEP 7: Update triangle mesh vertex indices based on octree partition
  SLIC_INFO("Updating mesh based on octree vertices.");
  update_trimesh_verts(surface_mesh, octree);

  print_surface_stats(surface_mesh);
  write_vtk(surface_mesh, "meldedTriMesh.vtk");

  // STEP 8: Insert triangles into the octree



  //

  asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Warning);

  // Other -- find leaf block of a given query point at various levels of resolution
  double alpha = 2./3.;
  SpacePt queryPt = SpacePt::lerp(meshBB.getMin(), meshBB.getMax(), alpha);

  SLIC_INFO("Finding associated grid point for query point: " << queryPt );
  for(int lev = 0; lev < octree.maxLeafLevel(); ++lev)
  {
      GridPt gridPt = octree.findGridCellAtLevel(queryPt, lev);
      SLIC_INFO("  @level " << lev
              <<":\n\t" <<  gridPt
              <<"\n\t[max gridPt: " << octree.maxGridCellAtLevel(lev)
              <<"; spacing" << octree.spacingAtLevel(lev)
              <<";\n\t bounding box " << octree.blockBoundingBox(gridPt, lev)
              <<"]");
  }


  SLIC_INFO("Recursively refining around query point: " << queryPt);
  refineAndPrint(octree, queryPt, false);
//  for(int i=0; i< octree.maxInternalLevel(); ++i)
//      refineAndPrint(octree, queryPt);

  return 0;
}
