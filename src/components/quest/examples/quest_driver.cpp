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
#include "quest/BoundingBox.hpp"
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
   double origin[3];
   umesh->getOrigin( origin );

   double h[3];
   umesh->getSpacing( h );

   // STEP 2: tag boundary elements
   const int ncells = surface_mesh->getMeshNumberOfCells();
   for ( int i=0; i < ncells; ++i ) {

      int cellids[3];
      surface_mesh->getMeshCell( i, cellids );

      int ijkmin[3];
      std::fill( ijkmin, ijkmin+3, std::numeric_limits<int>::max() );

      int ijkmax[3];
      std::fill( ijkmax, ijkmax+3, std::numeric_limits<int>::min() );

      // find sub-extent that encapsulates the surface element
      for ( int cn=0; cn < 3; ++cn ) {

          const int nodeIdx = cellids[ cn ];

          double pnt[3];
          surface_mesh->getMeshNode( nodeIdx, pnt );

          const double dx = pnt[0]-origin[0];
          const double dy = pnt[1]-origin[1];
          const double dz = pnt[2]-origin[2];

          const int i = std::floor( dx/h[0] );
          const int j = std::floor( dy/h[1] );
          const int k = std::floor( dz/h[2] );

          // i dimension
          if ( i <= ijkmin[0] ) {
             ijkmin[0] = i;
          }
          if ( i >= ijkmax[0] ) {
             ijkmax[0] = i;
          }

          // j dimension
          if ( j <= ijkmin[1] ) {
              ijkmin[1] = j;
          }
          if ( j >= ijkmax[1] ) {
              ijkmax[1] = j;
          }

          // k dimension
          if ( k <= ijkmin[2] ) {
              ijkmin[2] = k;
          }
          if ( k >= ijkmax[2] ) {
              ijkmax[2] = k;
          }

      } // END for all cell nodes;

      // flag all elements within the sub-extend as boundary elements
      for ( int ii=ijkmin[0]; ii <= ijkmax[0]; ++ii ) {
         for ( int jj=ijkmin[1]; jj <= ijkmax[1]; ++jj ) {
            for ( int kk=ijkmin[2]; kk <= ijkmax[2]; ++kk ) {

                const int idx = umesh->getLinearIndex( ii, jj, kk );
                bndry[ idx ] = 1;
            } // END for all kk
         } // END for all jj
      } // END for all ii


   } // END for all surface elements

}

//------------------------------------------------------------------------------
quest::BoundingBox getCellBoundingBox( int cellIdx,meshtk::Mesh* surface_mesh )
{
   // Sanity checks
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( cellIdx >= 0 && cellIdx < surface_mesh->getMeshNumberOfCells());

   int cell[3];
   surface_mesh->getMeshCell( cellIdx, cell );

   double min[3];
   double max[3];
   std::fill( min, min+3, std::numeric_limits<double>::max() );
   std::fill( max, max+3, std::numeric_limits<double>::min() );

   for ( int i=0; i < 3; ++i ) {

      double pnt[3];
      surface_mesh->getMeshNode( cell[i], pnt );

      for ( int j=0; j < 3; ++j ) {

          if ( pnt[ j ] < min[ j ] ) {

              min[ j ] = pnt[ j ];

          }

          if ( pnt[ j ] > max[ j ] ) {

              max[ j ] = pnt[ j ];
          }

      } // END for each dimension

   } // END for all cell nodes

   quest::BoundingBox bb;
   bb.setMin( min );
   bb.setMax( max );
   return bb;
}

//------------------------------------------------------------------------------
void n2( meshtk::Mesh* surface_mesh, meshtk::UniformMesh* umesh )
{
   SLIC_ASSERT( surface_mesh != ATK_NULLPTR );
   SLIC_ASSERT( umesh != ATK_NULLPTR );

   // STEP 0: compute bounding boxes for each surface element O(n)
   const int ncells = surface_mesh->getMeshNumberOfCells();
   std::vector< quest::BoundingBox > cell_boxes( ncells );

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
int main( int ATK_NOT_USED(argc), char** argv )
{
  // STEP 0: Initialize SLIC Environment
  slic::initialize();
  slic::setLoggingLevel( asctoolkit::slic::message::Debug );
  slic::addStreamToAllLevels( new slic::GenericOutputStream(&std::cout) );

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

  double f;
  std::cout << "Inflate by N:";
  std::cin >> f;
  for (int i=0; i < 3; ++i ) {

      min[ i ] -= f;
      max[ i ] += f;
  }
  std::cout << "min: " << min[0] << ", " << min[1] << ", " << min[2] << "\n";
  std::cout << "max: " << max[0] << ", " << max[1] << ", " << max[2] << "\n";

  // STEP 7: get dimensions from user
  int nx, ny, nz;
  std::cout << "Enter Nx Ny Nz:";
  std::cin >> nx >> ny >> nz;
  double h[3];
  h[0] = (max[0]-min[0]) / nx;
  h[1] = (max[1]-min[1]) / ny;
  h[2] = (max[2]-min[2]) / nz;
  std::cout << "h: " << h[0] << " " << h[1] << " " << h[2] << std::endl;
  std::cout.flush();

  int ext[6];
  ext[0] = 0;
  ext[1] = nx-1;
  ext[2] = 0;
  ext[3] = ny-1;
  ext[4] = 0;
  ext[5] = nz-1;

  // STEP 8: Construct uniform mesh
  meshtk::UniformMesh* umesh = new meshtk::UniformMesh(3,min,h,ext);

  // STEP 9: Flag boundary cells on uniform mesh
//  flag_boundary( surface_mesh, umesh );
  n2(surface_mesh,umesh);
  write_vtk( umesh, "uniform_mesh.vtk" );

  // STEP 10: clean up
  delete surface_mesh;
  surface_mesh = ATK_NULLPTR;

  // STEP 11: Finalize SLIC environment
  slic::finalize();
  return 0;
}
