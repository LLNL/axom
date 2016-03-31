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
#include "common/Timer.hpp"

#include "quest/BoundingBox.hpp"
#include "quest/Field.hpp"
#include "quest/FieldData.hpp"
#include "quest/FieldVariable.hpp"
#include "quest/Mesh.hpp"
#include "quest/Orientation.hpp"
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


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>  // for setprecision()

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


void testIntersectionOnRegularGrid()
{
    static int const DIM = 3;
    typedef quest::Point< double,DIM >   PointType;
    typedef quest::Triangle< double,DIM > TriangleType;
    typedef quest::BoundingBox< double,DIM > BoundingBoxType;

    double xArr[3] = { 1., 0., 0.};
    double yArr[3] = { 0., 1., 0.};
    double zArr[3] = { 0., 0., 1.};

    PointType ptX(xArr);
    PointType ptY(yArr);
    PointType ptZ(zArr);

    TriangleType unitTri( ptX, ptY, ptZ );

    typedef meshtk::UnstructuredMesh<meshtk::MIXED> DebugMesh;
    DebugMesh* debugMesh = new DebugMesh(3);

    // Add triangle to mesh
    debugMesh->insertNode( ptX[0], ptX[1], ptX[2]);
    debugMesh->insertNode( ptY[0], ptY[1], ptY[2]);
    debugMesh->insertNode( ptZ[0], ptZ[1], ptZ[2]);

    int tArr[3] = {0,1,2};
    debugMesh->insertCell(tArr, meshtk::LINEAR_TRIANGLE, 3);

    PointType bbMin(-0.1);
    PointType bbMax(1.1);
    BoundingBoxType bbox( bbMin,bbMax );

    typedef quest::SpatialOctree<DIM, quest::BlockData> SpaceOctree;
    SpaceOctree oct( bbox);


    int lev =3;
    for(int i=0; i< 1<<lev; ++i)
    {
        for(int j=0; j< 1<<lev; ++j)
        {
            for(int k=0; k< 1<<lev; ++k)
            {
                SpaceOctree::BlockIndex block( quest::Point<int,3>::make_point(i,j,k), lev );
                SpaceOctree::GeometricBoundingBox blockBB = oct.blockBoundingBox(block);

                if( quest::intersect( unitTri, blockBB))
                {
                    // Add to debug mesh
                    int vStart = debugMesh->getMeshNumberOfNodes();

                    debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
                    debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
                    debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMin()[2]);
                    debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMin()[2]);

                    debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
                    debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
                    debugMesh->insertNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMax()[2]);
                    debugMesh->insertNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMax()[2]);

                    int data[8];
                    for(int i=0; i< 8; ++i)
                        data[i] = vStart + i;

                    debugMesh->insertCell( data, meshtk::LINEAR_HEX, 8);
                }
            }
        }
    }


    debugMesh->toVtkFile("gridIntersections.vtk");

    delete debugMesh;
}


void testContainmentOnRegularGrid(const Octree3D& inOutOctree
           , const GeometricBoundingBox& meshBounds
           , int gridRes)
{
    SpaceVector h( meshBounds.getMin(), meshBounds.getMax());
    for(int i=0; i<3; ++i)
        h[i] /= gridRes;

    int ext[6];
    ext[0] = ext[2] = ext[4] = 0;
    ext[1] = ext[3] = ext[5] = gridRes;

    meshtk::UniformMesh* umesh =
            new meshtk::UniformMesh(3,meshBounds.getMin().data(),h.data(),ext);


    const int nnodes = umesh->getNumberOfNodes();
    meshtk::FieldData* PD = umesh->getNodeFieldData();
    SLIC_ASSERT( PD != ATK_NULLPTR );

    PD->addField( new meshtk::FieldVariable< int >("containment",nnodes) );
    int* containment = PD->getField( "containment" )->getIntPtr();
    SLIC_ASSERT( containment != ATK_NULLPTR );


    asctoolkit::utilities::Timer timer(true);
    for ( int inode=0; inode < nnodes; ++inode )
    {
        quest::Point< double,3 > pt;
        umesh->getMeshNode( inode, pt.data() );

        containment[ inode ] = inOutOctree.within(pt) ? 1 : 0;
    }
    timer.stop();
    SLIC_INFO("\tQuerying "<< gridRes << "^3 containment field took " << timer.elapsed() << " seconds"
              << " (@ " << nnodes / timer.elapsed() << " queries per second)"
        );

  #ifdef DUMP_VTK_MESH
    std::stringstream sstr;
    sstr << "gridContainment_" << gridRes << ".vtk";
    write_vtk( umesh, sstr.str());
  #endif

    delete umesh;
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
   SLIC_DEBUG(edgeHistStr.str());

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
   SLIC_DEBUG(triHistStr.str());

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
       SLIC_DEBUG(badTriStr.str());
   }
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
  meshBB.scale(1.01);

  asctoolkit::slic::setLoggingMsgLevel( asctoolkit::slic::message::Info);


  testIntersectionOnRegularGrid();

  // STEP 6: Create octree over mesh's bounding box and query a point in space
  Octree3D octree(meshBB, surface_mesh);
  octree.generateIndex();


  print_surface_stats(surface_mesh);
  write_vtk(surface_mesh, "meldedTriMesh.vtk");

  SLIC_INFO("-- About to query the octree");
  // STEP 7: Insert triangles into the octree
  for(int i=1; i< 9 ; ++i)
      testContainmentOnRegularGrid( octree, meshBB, 1<<i);


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
