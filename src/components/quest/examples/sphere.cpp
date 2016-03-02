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
 * \file sphere.cpp
 *
 * \date Dec 16, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
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
#include "quest/HyperSphere.hpp"
#include "quest/Mesh.hpp"
#include "quest/Orientation.hpp"
#include "quest/Point.hpp"
#include "quest/Vector.hpp"
#include "quest/STLReader.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/Triangle.hpp"
#include "quest/UniformMesh.hpp"
#include "quest/UnstructuredMesh.hpp"
#include "quest/fuzzy_compare.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "slam/Utilities.hpp"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

using namespace asctoolkit;

static const int DIM = 3;
typedef meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE > TriangleMesh;


typedef quest::Point<int, DIM> GridPoint;
typedef quest::Point<double, DIM> SpacePoint;;
typedef quest::Vector<double, DIM> SpaceVector;

typedef quest::BoundingBox<double, DIM> SpatialBoundingBox;;

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
SpatialBoundingBox compute_bounds( meshtk::Mesh* mesh)
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   SpatialBoundingBox meshBB;
   SpacePoint pt;

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
   {
       mesh->getMeshNode( i, pt.data() );
       meshBB.addPoint( pt );
   } // END for all nodes

   SLIC_ASSERT( meshBB.isValid() );

   return meshBB;
}


//------------------------------------------------------------------------------
quest::BoundingBox<double,3> getCellBoundingBox( int cellIdx,meshtk::Mesh* surface_mesh )
{
   // Sanity checks
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( cellIdx >= 0 && cellIdx < surface_mesh->getMeshNumberOfCells());

   using namespace quest;

   int const DIM = 3;

   Point<int,DIM> cell;
   surface_mesh->getMeshCell( cellIdx, cell.data() );

   BoundingBox<double,DIM> bb;
   Point<double,DIM> pt;

   for ( int i=0; i < 3; ++i ) {
      surface_mesh->getMeshNode( cell[i], pt.data() );
      bb.addPoint( pt );

   } // END for all cell nodes

   return bb;
}

//------------------------------------------------------------------------------
void n2( meshtk::Mesh* surface_mesh, meshtk::UniformMesh* umesh )
{
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( umesh != ATK_NULLPTR );

   // STEP 1: Setup node-centered signed distance field on uniform mesh
   const int nnodes = umesh->getNumberOfNodes();
   meshtk::FieldData* PD = umesh->getNodeFieldData();
   SLIC_ASSERT( PD != ATK_NULLPTR );

   PD->addField( new meshtk::FieldVariable<double>("phi",nnodes) );
   double* phi = PD->getField( "phi" )->getDoublePtr();
   SLIC_ASSERT( phi != ATK_NULLPTR );


   // STEP 2: loop over uniform mesh nodes and compute distance field
   SLIC_INFO("Calculating distance field...");

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < nnodes; ++i ) {

      // get target node
      SpacePoint pnt;
      umesh->getNode( i, pnt.data() );
      SpacePoint Q = pnt;

      double unsignedMinDistSQ = std::numeric_limits< double >::max();
      int sign = 0;

      const int ncells = surface_mesh->getMeshNumberOfCells();
      for (int j=0; j < ncells; ++j ) {

          // find minimum distance from query point to triangle
          GridPoint closest_cell;
          surface_mesh->getMeshCell( j, closest_cell.data() );

          SpacePoint a,b,c;
          surface_mesh->getMeshNode( closest_cell[0], a.data() );
          surface_mesh->getMeshNode( closest_cell[1], b.data() );
          surface_mesh->getMeshNode( closest_cell[2], c.data() );
          quest::Triangle< double,3 > T( a,b,c);
          const double sqDist = quest::squared_distance( Q, T );
          if ( sqDist < unsignedMinDistSQ)  {

              bool negSide = quest::orientation( Q, T ) == quest::ON_NEGATIVE_SIDE;
              sign = negSide? -1: 1;
              unsignedMinDistSQ = sqDist;

          } // END if

      } // END for all cells on the surface mesh

      phi[i] = sign * std::sqrt( unsignedMinDistSQ );

   } // END for all nodes on the uniform mesh
   SLIC_INFO("done." );
}

//------------------------------------------------------------------------------
void expected_phi(meshtk::UniformMesh* umesh)
{
   SLIC_ASSERT( umesh != ATK_NULLPTR );

   // STEP 0: Construct sphere centered at (0.0,0.0,0.0) with radius 5.0
   quest::HyperSphere< double, 3 > sphere( 5.0 );

   // STEP 1: Add node field to stored exact distance field.
   const int nnodes = umesh->getNumberOfNodes();
   meshtk::FieldData* PD = umesh->getNodeFieldData();
   SLIC_ASSERT( PD != ATK_NULLPTR );

   PD->addField( new meshtk::FieldVariable<double>("expected_phi",nnodes) );
   double* phi = PD->getField( "expected_phi" )->getDoublePtr();
   SLIC_ASSERT( phi != ATK_NULLPTR );

   // STEP 2: loop over uniform mesh nodes and compute distance field
   SLIC_INFO( "Calculating analytic distance field...");
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < nnodes; ++i ) {

       double pnt[3];
       umesh->getNode( i, pnt );

       phi[ i ]  = sphere.getSignedDistance( pnt );
   }

   SLIC_INFO("done.");
}

//------------------------------------------------------------------------------
void l2norm( meshtk::UniformMesh* umesh )
{
   SLIC_ASSERT( umesh != ATK_NULLPTR );

   const int nnodes = umesh->getMeshNumberOfNodes();

   // STEP 0: Get computed & expected fields
   meshtk::FieldData* PD = umesh->getNodeFieldData();
   double* phi_computed = PD->getField( "phi" )->getDoublePtr();
   double* phi_expected = PD->getField( "expected_phi" )->getDoublePtr();

   // STEP 1: Add field to store error
   PD->addField( new meshtk::FieldVariable< double >( "error",nnodes ) );
   double* error = PD->getField( "error" )->getDoublePtr();
   SLIC_ASSERT( error != ATK_NULLPTR );

   // STEP 2: loop and calculate the error
   for ( int i=0; i < nnodes; ++i ) {

       // relative error
       const double delta = phi_computed[ i ]-phi_expected[ i ];
       error[ i ] = delta / phi_expected[ i ];

   } // END for all nodes

}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::UnitTestLogger logger;  // create & initialize logger

  bool hasInputArgs = argc > 1;

  // STEP 1: get file from user or use default
  std::string stlFile;
  if(hasInputArgs)
  {
      stlFile = std::string( argv[1] );
  }
  else
  {
      const std::string defaultFileName = "sphere.stl";
      const std::string defaultDir = "src/components/quest/data/";

      stlFile = asctoolkit::utilities::filesystem::joinPath(defaultDir, defaultFileName);
  }

  stlFile = asctoolkit::slam::util::findFileInAncestorDirs(stlFile);
  SLIC_ASSERT( asctoolkit::utilities::filesystem::pathExists( stlFile));

  // STEP 2: read file
  SLIC_INFO("Reading file: " << stlFile << "...");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( stlFile );
  reader->read();
  SLIC_INFO("done");

  // STEP 3: get surface mesh
  meshtk::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );

  // STEP 4: Delete the reader
  delete reader;
  reader = ATK_NULLPTR;

  // STEP 5: write vtk file
  write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 6: compute bounds
  SpatialBoundingBox meshBounds = compute_bounds( surface_mesh);

  SLIC_INFO("Mesh bounding  box: " << meshBounds << "\n");

  // Add inflation factor
  std::cout << "Inflate by N: \n";
  double f = 2.;
  if(hasInputArgs)
  {
      std::cin >> f;
  }
  meshBounds.expand(f);
  SLIC_INFO("Bounding  box after inflating by " << f <<": " << meshBounds );

  // STEP 7: get dimensions from user
  std::cout << "Enter Nx Ny Nz: \n";
  GridPoint gridRes(16);
  if(hasInputArgs)
  {
      std::cin >> gridRes[0] >> gridRes[1] >> gridRes[2];
  }
  SLIC_INFO("Distance field grid is: " << gridRes);


  SpaceVector h( meshBounds.getMin(), meshBounds.getMax());
  for(int i=0; i<DIM; ++i)
      h[i] /= gridRes[i];
  SLIC_DEBUG("h: " << h );

  int ext[6];
  ext[0] = 0;
  ext[1] = gridRes[0];
  ext[2] = 0;
  ext[3] = gridRes[1];
  ext[4] = 0;
  ext[5] = gridRes[2];

  // STEP 8: Construct uniform mesh
  meshtk::UniformMesh* umesh = new meshtk::UniformMesh(3,meshBounds.getMin().data(),h.data(),ext);

  // STEP 9: Flag boundary cells on uniform mesh
  //  flag_boundary( surface_mesh, umesh );
  n2(surface_mesh,umesh);
  expected_phi( umesh );
  l2norm( umesh );
  write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 10: clean up
  delete surface_mesh;
  surface_mesh = ATK_NULLPTR;

  delete umesh;
  umesh = ATK_NULLPTR;

  return 0;
}
