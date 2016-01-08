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
#include "quest/BoundingBox.hpp"
#include "quest/Field.hpp"
#include "quest/FieldData.hpp"
#include "quest/FieldVariable.hpp"
#include "quest/HyperSphere.hpp"
#include "quest/Mesh.hpp"
#include "quest/Orientation.hpp"
#include "quest/Point.hpp"
#include "quest/STLReader.hpp"
#include "quest/SquaredDistance.hpp"
#include "quest/Triangle.hpp"
#include "quest/UniformMesh.hpp"
#include "quest/UnstructuredMesh.hpp"
#include "quest/fuzzy_compare.hpp"
#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"

// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>

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
void compute_bounds( meshtk::Mesh* mesh, double min[3], double max[3] )
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   std::fill( min, min+3, std::numeric_limits< double >::max() );
   std::fill( max, max+3, std::numeric_limits< double >::min() );

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i ) {

       double pnt[3];
       mesh->getMeshNode( i, pnt );

       for (int dim=0; dim < 3; ++dim ) {

           if ( pnt[dim] < min[dim] ) {

               min[dim] = pnt[dim];

           } else if ( pnt[dim] > max[dim] ) {

               max[dim] = pnt[dim];

           }

       } // END for all dims


   } // END for all nodes


}


//------------------------------------------------------------------------------
quest::BoundingBox<double,3> getCellBoundingBox( int cellIdx,meshtk::Mesh* surface_mesh )
{
   // Sanity checks
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( cellIdx >= 0 && cellIdx < surface_mesh->getMeshNumberOfCells());

   int cell[3];
   surface_mesh->getMeshCell( cellIdx, cell );
   quest::BoundingBox<double,3> bb;

   for ( int i=0; i < 3; ++i ) {

      double pnt[3];
      surface_mesh->getMeshNode( cell[i], pnt );

      bb.addPoint( quest::Point<double,3>(pnt) );

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
   std::cout << "Calculating distance field...";
   std::cout.flush();
   typedef quest::Point3D Point3D;
   typedef quest::Point<int, 3> GridPt;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < nnodes; ++i ) {

      // get target node
      Point3D pnt;
      umesh->getNode( i, pnt.data() );
      Point3D Q = pnt;

      phi[ i ] = std::numeric_limits< double >::max();

      const int ncells = surface_mesh->getMeshNumberOfCells();
      for (int j=0; j < ncells; ++j ) {

          // compute exact distance to the cell
          GridPt closest_cell;
          surface_mesh->getMeshCell( j, closest_cell.data() );

          Point3D a,b,c;
          surface_mesh->getMeshNode( closest_cell[0], a.data() );
          surface_mesh->getMeshNode( closest_cell[1], b.data() );
          surface_mesh->getMeshNode( closest_cell[2], c.data() );
          quest::Triangle< double,3 > T( a,b,c);
          const double dist = std::sqrt( quest::squared_distance( Q, T ) );
          if ( dist < std::abs( phi[i] ) )  {

              double sign = 1.0;
              if ( quest::orientation( Q, T ) == quest::ON_NEGATIVE_SIDE ) {
                  sign = -1.0f;
              }

              phi[ i ] = sign*dist;

          } // END if

      } // END for all cells on the surface mesh

   } // END for all nodes on the uniform mesh
   std::cout << "[DONE]" << std::endl;

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
   std::cout << "Calculating analytic distance field...";
   std::cout.flush();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < nnodes; ++i ) {

       double pnt[3];
       umesh->getNode( i, pnt );

       phi[ i ]  = sphere.getSignedDistance( pnt );
   }

   std::cout << "[DONE]" << std::endl;
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
int main( int ATK_NOT_USED(argc), char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingMsgLevel( asctoolkit::slic::message::Debug );
  slic::addStreamToAllMsgLevels( new slic::GenericOutputStream(&std::cout) );

  // STEP 1: get file from user or use default
  std::string stlFile = std::string( argv[1] ) ;

  // STEP 2: read file
  std::cout << "Reading file: " << stlFile << "...";
  std::cout.flush();
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( stlFile );
  reader->read();
  std::cout << "[DONE]\n";
  std::cout.flush();

  // STEP 3: get surface mesh
  meshtk::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader-> getMesh( static_cast<TriangleMesh*>( surface_mesh ) );

  // STEP 4: Delete the reader
  delete reader;
  reader = ATK_NULLPTR;

  // STEP 5: write vtk file
  write_vtk( surface_mesh, "surface_mesh.vtk" );

  // STEP 6: compute bounds
  double min[3];
  double max[3];
  compute_bounds( surface_mesh, min, max );
  std::cout << "min: " << min[0] << ", " << min[1] << ", " << min[2] << "\n";
  std::cout << "max: " << max[0] << ", " << max[1] << ", " << max[2] << "\n";

  double f;
  std::cout << "Inflate by N: ";
  std::cin >> f;
  for (int i=0; i < 3; ++i ) {

      min[ i ] -= f;
      max[ i ] += f;
  }
  std::cout << "min: " << min[0] << ", " << min[1] << ", " << min[2] << "\n";
  std::cout << "max: " << max[0] << ", " << max[1] << ", " << max[2] << "\n";

  // STEP 7: get dimensions from user
  int nx, ny, nz;
  std::cout << "Enter Nx Ny Nz: ";
  std::cin >> nx >> ny >> nz;
  double h[3];
  h[0] = (max[0]-min[0]) / nx;
  h[1] = (max[1]-min[1]) / ny;
  h[2] = (max[2]-min[2]) / nz;
  std::cout << "h: " << h[0] << " " << h[1] << " " << h[2] << std::endl;
  std::cout.flush();

  int ext[6];
  ext[0] = 0;
  ext[1] = nx;
  ext[2] = 0;
  ext[3] = ny;
  ext[4] = 0;
  ext[5] = nz;

  // STEP 8: Construct uniform mesh
  meshtk::UniformMesh* umesh = new meshtk::UniformMesh(3,min,h,ext);

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

  // STEP 11: Finalize SLIC environment
  slic::finalize();
  return 0;
}
