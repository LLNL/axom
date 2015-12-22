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
#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"
//#include "slam/FileUtilities.hpp"

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
void compute_bounds( meshtk::Mesh* mesh, double minPt[3], double maxPt[3] )
{
   SLIC_ASSERT( mesh != ATK_NULLPTR );

   using namespace quest;

   int const DIM = 3;
   BoundingBox<double, DIM> meshBB;
   Point3D pt;

   for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
   {
       mesh->getMeshNode( i, pt.data() );
       meshBB.addPoint( pt );
   } // END for all nodes

   SLIC_ASSERT( meshBB.isValid() );

   // copy bounds out to function parameters
   meshBB.getMin().to_array(minPt);
   meshBB.getMax().to_array(maxPt);

}

//------------------------------------------------------------------------------
void flag_boundary( meshtk::Mesh* surface_mesh, meshtk::UniformMesh* umesh )
{
   /* Sanity checks */
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( umesh != ATK_NULLPTR );

   // STEP 0: Create field variable on uniform mesh
   const int umesh_ncells = umesh->getMeshNumberOfCells();
   meshtk::FieldData* CD = umesh->getCellFieldData();
   CD->addField( new meshtk::FieldVariable< int >("tag", umesh_ncells) );

   int* bndry = CD->getField( "tag" )->getIntPtr();
   SLIC_ASSERT( bndry != ATK_NULLPTR );
   for ( int i=0; i < umesh_ncells; ++i ) {
       bndry[ i ] = 0;
   }

   // STEP 1: get uniform mesh origin & spacing
   typedef quest::Point3D Point3D;
   Point3D origin;
   umesh->getOrigin( origin.data() );

   Point3D h;
   umesh->getSpacing( h.data() );

   // STEP 2: tag boundary elements
   typedef quest::BoundingBox<int, 3> GridBounds;
   typedef GridBounds::PointType GridPt;
   GridBounds globalGridBounds;
   const int ncells = surface_mesh->getMeshNumberOfCells();
   for ( int i=0; i < ncells; ++i )
   {
      GridPt cellids;
      surface_mesh->getMeshCell( i, cellids.data() );

      GridBounds ijkBounds;

      // find sub-extent that encapsulates the surface element
      for ( int cn=0; cn < 3; ++cn ) {

          const int nodeIdx = cellids[ cn ];

          Point3D pnt;
          surface_mesh->getMeshNode( nodeIdx, pnt.data() );

          const double dx = pnt[0]-origin[0];
          const double dy = pnt[1]-origin[1];
          const double dz = pnt[2]-origin[2];

          const int i = std::floor( dx/h[0] ) ;
          const int j = std::floor( dy/h[1] ) ;
          const int k = std::floor( dz/h[2] ) ;

          ijkBounds.addPoint( GridPt::make_point(i,j,k) );
          globalGridBounds.addPoint( GridPt::make_point(i,j,k) );

      } // END for all cell nodes;

      // flag all elements within the sub-extend as boundary elements
      const GridPt& ijkmin = ijkBounds.getMin();
      const GridPt& ijkmax = ijkBounds.getMax();

      // std::cout <<"\nGrid bounding box: " << ijkBounds << std::endl;

      for ( int ii=ijkmin[0]; ii <= ijkmax[0]; ++ii ) {
         for ( int jj=ijkmin[1]; jj <= ijkmax[1]; ++jj ) {
            for ( int kk=ijkmin[2]; kk <= ijkmax[2]; ++kk ) {

                const int idx_wrong = umesh->getCellLinearIndex( ii, jj, kk );
                const int RES = 32;
                const int idx = ii + RES * jj + RES*RES* kk;
                ++bndry[ idx_wrong ];


            } // END for all kk
         } // END for all jj
      } // END for all ii

   } // END for all surface elements


   std::cout<<"Overall grid bounds: " << globalGridBounds << std::endl;
}

/**
 * \brief Try to find a valid input file.
 * \return The file
 */
std::string getStlFileName(const std::string & inputFileName)
{
    std::string fileName = inputFileName;
    const int MAX_ATTEMPTS= 3;

    if( inputFileName == "")
    {
        // Try to use the default file name.
        // HACK: If it doesn't work, try to access from a parent directory
        //       Fix if/when we have better file utilities
        //       and replace with path to the root of the data repo when available
        const std::string defaultFileName = "plane.stl";
        const std::string defaultDir = "src/components/quest/data/";
        fileName = defaultDir + defaultFileName;

        std::ifstream meshFile(fileName.c_str());
        for(int attempts = 0; !meshFile && attempts < MAX_ATTEMPTS; ++attempts)
        {
          fileName = "../" + fileName;
          meshFile.open( fileName.c_str());
        }

        SLIC_ERROR_IF( !meshFile
            , "fstream error -- problem opening file: '"  << defaultDir + defaultFileName
                                                          << "' (also tried several"
                                                          << " ancestors up to '../" << fileName << "')."
                                                          << "\nThe current working directory is: '"
                                                          /*<< asctoolkit::slam::util::getCWD() << "'"*/);
        meshFile.close();
    }

    return fileName;
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

   // STEP 0: compute bounding boxes for each surface element O(n)
   const int ncells = surface_mesh->getMeshNumberOfCells();
   std::vector< quest::BoundingBox<double,3> > cell_boxes( ncells );

   std::cout << "Computing bounding boxes....";
   std::cout.flush();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < ncells; ++i ) {

       cell_boxes[ i ] = getCellBoundingBox( i, surface_mesh );
   }
   std::cout << "[DONE]" << std::endl;

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
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
   for ( int i=0; i < nnodes; ++i ) {

      // get target node
      double pnt[3];
      umesh->getNode( i, pnt );
      quest::Point< double, 3 > Q =
              quest::Point< double,3 >::make_point(pnt[0],pnt[1],pnt[2]);

      int cellIdx              = -1;
      double cell_box_distance = std::numeric_limits< double >::max();

      // bounding-box filter
      const int ncells = surface_mesh->getMeshNumberOfCells();
      for (int j=0; j < ncells; ++j ) {

          const double dist = quest::squared_distance( Q, cell_boxes[j] );
          if ( std::abs(dist) < cell_box_distance ) {
              cellIdx = j;
              cell_box_distance = std::abs(dist);
          }
      } // END for all cells on the surface mesh

      // compute exact distance to the cell
      SLIC_ASSERT( cellIdx >= 0 && cellIdx < ncells );
      int closest_cell[3];
      surface_mesh->getMeshCell( cellIdx, closest_cell );

      double a[3];
      double b[3];
      double c[3];
      surface_mesh->getMeshNode( closest_cell[0], a );
      surface_mesh->getMeshNode( closest_cell[1], b );
      surface_mesh->getMeshNode( closest_cell[2], c );
      quest::Point< double, 3 > A =
              quest::Point< double,3 >::make_point( a[0], a[1], a[2] );
      quest::Point< double,3 > B =
              quest::Point< double,3 >::make_point( b[0], b[1], b[2] );
      quest::Point< double,3 > C =
              quest::Point< double,3 >::make_point( c[0], c[1], c[2] );

      quest::Triangle< double,3 > T( A,B,C );

      phi[ i ] = quest::squared_distance( Q, T );
   } // END for all nodes on the uniform mesh
   std::cout << "[DONE]" << std::endl;

}
//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingLevel( asctoolkit::slic::message::Debug );
  slic::addStreamToAllLevels( new slic::GenericOutputStream(&std::cout) );

  bool hasInputArgs = argc > 1;

  // STEP 1: get file from user or use default
  std::string stlFile = getStlFileName(hasInputArgs ? std::string( argv[1] ) : "") ;

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
  // dump mesh info
  std::cout<<"Mesh has "
          << surface_mesh->getMeshNumberOfNodes() << " nodes and "
          << surface_mesh->getMeshNumberOfCells() << " cells."
          << std::endl;

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

  // STEP 7: get dimensions from user
  int nx, ny, nz;
  if(hasInputArgs > 1)
  {
      std::cout << "Enter Nx Ny Nz:";
      std::cin >> nx >> ny >> nz;
  }
  else
  {
      nx = ny = nz = 32;
  }

  double h[3];
  h[0] = (max[0]-min[0]) / (nx);
  h[1] = (max[1]-min[1]) / (ny);
  h[2] = (max[2]-min[2]) / (nz);
  std::cout << "h: " << h[0] << " " << h[1] << " " << h[2] << std::endl;
  std::cout.flush();

  int node_ext[6];
  node_ext[0] = 0;
  node_ext[1] = nx;
  node_ext[2] = 0;
  node_ext[3] = ny;
  node_ext[4] = 0;
  node_ext[5] = nz;

  // STEP 8: Construct uniform mesh
  meshtk::UniformMesh* umesh = new meshtk::UniformMesh(3,min,h,node_ext);

  // STEP 9: Flag boundary cells on uniform mesh
  flag_boundary( surface_mesh, umesh );
  write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 10: clean up
  delete surface_mesh;
  surface_mesh = ATK_NULLPTR;

  // STEP 11: Finalize SLIC environment
  slic::finalize();
  return 0;
}
