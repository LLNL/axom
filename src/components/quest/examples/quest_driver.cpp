/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Dec 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

// ATK Toolkit includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "common/FileUtilities.hpp"
#include "common/Timer.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/BVHTree.hpp"
#include "quest/Field.hpp"
#include "quest/FieldData.hpp"
#include "quest/FieldVariable.hpp"
#include "quest/Mesh.hpp"
#include "quest/Point.hpp"
#include "quest/STLReader.hpp"
#include "quest/Triangle.hpp"
#include "quest/UniformMesh.hpp"
#include "quest/UnstructuredMesh.hpp"
#include "quest/SignedDistance.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

#include "slam/Utilities.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>

using namespace asctoolkit;

typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;

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

//------------------------------------------------------------------------------
quest::BoundingBox< double,3 > compute_bounds( meshtk::Mesh* mesh)
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   using namespace quest;

   quest::BoundingBox< double,3 > meshBB;
   quest::Point< double,3 > pt;

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
   {
       mesh->getMeshNode( i, pt.data() );
       meshBB.addPoint( pt );
   } // END for all nodes

   SLIC_ASSERT( meshBB.isValid() );

   return meshBB;
}


//------------------------------------------------------------------------------
void distance_field( meshtk::Mesh* surface_mesh, meshtk::UniformMesh* umesh )
{
  SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
  SLIC_ASSERT( umesh != ATK_NULLPTR );

  quest::SignedDistance< 3 > signedDistance( surface_mesh, 25, 10 );

#ifdef ATK_DEBUG
  // write the bucket tree to a file
  const quest::BVHTree< int, 3>* btree = signedDistance.getBVHTree();
  SLIC_ASSERT( btree != ATK_NULLPTR );

  btree->writeLegacyVtkFile( "bucket-tree.vtk" );

  // mark bucket IDs on surface mesh
  const int ncells = surface_mesh->getMeshNumberOfCells();
  meshtk::FieldData* CD = surface_mesh->getCellFieldData();
  CD->addField( new meshtk::FieldVariable<int>( "BucketID", ncells ) );
  int* bidx = CD->getField( "BucketID" )->getIntPtr();
  SLIC_ASSERT( bidx != ATK_NULLPTR );

  const int numObjects = btree->getNumberOfObjects();
  for ( int i=0; i < numObjects; ++i ) {

     const int idx = btree->getObjectBucketIndex( i );
     bidx[ i ] = idx;

  } // END for all objects

  write_vtk( surface_mesh, "partitioned_surface_mesh.vtk" );
#endif

  const int nnodes = umesh->getNumberOfNodes();
  meshtk::FieldData* PD = umesh->getNodeFieldData();
  SLIC_ASSERT( PD != ATK_NULLPTR );

  PD->addField( new meshtk::FieldVariable< double >("phi",nnodes) );
  double* phi = PD->getField( "phi" )->getDoublePtr();
  SLIC_ASSERT( phi != ATK_NULLPTR );

  for ( int inode=0; inode < nnodes; ++inode ) {

      quest::Point< double,3 > pt;
      umesh->getMeshNode( inode, pt.data() );

      phi[ inode ] = signedDistance.computeDistance( pt );

  } // END for all nodes

}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingMsgLevel( asctoolkit::slic::message::Debug );

  // Create a more verbose message for this application (only level and message)
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  slic::GenericOutputStream* defaultStream =
          new slic::GenericOutputStream(&std::cout);
  slic::GenericOutputStream* compactStream =
          new slic::GenericOutputStream(&std::cout, slicFormatStr);
  slic::addStreamToMsgLevel(defaultStream, asctoolkit::slic::message::Fatal) ;
  slic::addStreamToMsgLevel(defaultStream, asctoolkit::slic::message::Error);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Warning);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Info);
  slic::addStreamToMsgLevel(compactStream, asctoolkit::slic::message::Debug);

  bool hasInputArgs = argc > 1;

  // STEP 1: get file from user or use default
  std::string stlFile;
  if( hasInputArgs ) {

    stlFile = std::string( argv[1] );

  }
  else {

    const std::string defaultFileName = "plane_simp.stl";
    const std::string defaultDir = "src/components/quest/data/";

    stlFile =
       asctoolkit::utilities::filesystem::joinPath(defaultDir, defaultFileName);

  }
  stlFile = asctoolkit::slam::util::findFileInAncestorDirs(stlFile);
  SLIC_ASSERT( asctoolkit::utilities::filesystem::pathExists(stlFile));

  // STEP 2: read file
  SLIC_INFO( "Reading file: " << stlFile << "...");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( stlFile );
  reader->read();
  SLIC_INFO("done");

  // STEP 3: get surface mesh
  meshtk::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );
  SLIC_INFO("Mesh has "
          << surface_mesh->getMeshNumberOfNodes() << " nodes and "
          << surface_mesh->getMeshNumberOfCells() << " cells.");

  // STEP 4: Delete the reader
  delete reader;
  reader = ATK_NULLPTR;

  // STEP 5: write vtk file
  write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 6: compute bounds
  quest::BoundingBox< double,3 > meshBB = compute_bounds( surface_mesh);
  SLIC_INFO("Mesh bounding box: " << meshBB );

  // STEP 7: get dimensions from user
  int nx, ny, nz;
  if ( hasInputArgs ) {

    std::cout << "Enter Nx Ny Nz:";
    std::cin >> nx >> ny >> nz;

  } else {

    nx = ny = nz = 32;

  }

  quest::BoundingBox< double,3 > queryBounds = meshBB;
  // queryBounds.scale(2.);

  double h[3];
  const quest::Vector< double,3 >& bbDiff = queryBounds.range();
  h[0] = bbDiff[0] / nx;
  h[1] = bbDiff[1] / ny;
  h[2] = bbDiff[2] / nz;
  SLIC_INFO("grid cell size:(" << h[0] << "," << h[1] << "," << h[2] << ")\n" );

  int node_ext[6];
  node_ext[0] = 0;
  node_ext[1] = nx;
  node_ext[2] = 0;
  node_ext[3] = ny;
  node_ext[4] = 0;
  node_ext[5] = nz;

  // STEP 8: Construct uniform mesh
  meshtk::UniformMesh* umesh =
          new meshtk::UniformMesh(3, queryBounds.getMin().data(), h, node_ext);

  // STEP 9: Compute the distance field on the uniform mesh
  SLIC_INFO( "computing distance field..." );
  utilities::Timer timer;
  timer.start();
  distance_field( surface_mesh, umesh );
  timer.stop();
  SLIC_INFO( "Time to compute distance field: " << timer.elapsed() );

  // STEP 10: write the uniform mesh
  write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 11: clean up
  delete surface_mesh;
  surface_mesh = ATK_NULLPTR;

  delete umesh;
  umesh = ATK_NULLPTR;

  // STEP 12: Finalize SLIC environment
  slic::finalize();
  return 0;
}
