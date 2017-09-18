/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/*!
 * \file containment_driver.cpp
 * \brief Basic demo of point containment acceleration structure over surfaces.
 */

// axom includes
#include "axom/config.hpp"
#include "axom/Macros.hpp"
#include "axom/Types.hpp"
#include "axom_utils/FileUtilities.hpp"
#include "axom_utils/Timer.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"

#include "primal/orientation.hpp"
#include "primal/squared_distance.hpp"

#include "quest/STLReader.hpp"
#include "quest/SpatialOctree.hpp"
#include "quest/InOutOctree.hpp"

#include "mint/Field.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"
#include "mint/Mesh.hpp"
#include "mint/UniformMesh.hpp"
#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"

#include "slic/slic.hpp"
#include "slic/UnitTestLogger.hpp"

#include "slam/Utilities.hpp"

#include "fmt/format.h"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <fstream>
#include <iomanip>  // for setprecision()

using namespace axom;

typedef axom::mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;

typedef axom::quest::InOutOctree<3> Octree3D;

typedef axom::primal::Point<int,3> TriVertIndices;
typedef axom::primal::Triangle<double, 3> SpaceTriangle;

typedef Octree3D::GeometricBoundingBox GeometricBoundingBox;
typedef Octree3D::SpacePt SpacePt;
typedef Octree3D::SpaceVector SpaceVector;
typedef Octree3D::GridPt GridPt;
typedef Octree3D::BlockIndex BlockIndex;

#ifdef AXOM_DEBUG
const int MAX_CONTAINMENT_QUERY_LEVEL = 7;
#else
const int MAX_CONTAINMENT_QUERY_LEVEL = 9;
#endif


//------------------------------------------------------------------------------
/**
 * \brief Computes the bounding box of the surface mesh
 */
GeometricBoundingBox compute_bounds( axom::mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  GeometricBoundingBox meshBB;
  SpacePt pt;

  for ( int i=0 ; i < mesh->getMeshNumberOfNodes() ; ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  }

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}


void testIntersectionOnRegularGrid()
{
  static int const DIM = 3;
  typedef axom::primal::Point< double,DIM >   PointType;
  typedef axom::primal::Triangle< double,DIM > TriangleType;
  typedef axom::primal::BoundingBox< double,DIM > BoundingBoxType;

  double xArr[3] = { 1., 0., 0.};
  double yArr[3] = { 0., 1., 0.};
  double zArr[3] = { 0., 0., 1.};

  PointType ptX(xArr);
  PointType ptY(yArr);
  PointType ptZ(zArr);

  TriangleType unitTri( ptX, ptY, ptZ );

  typedef axom::mint::UnstructuredMesh< MINT_MIXED_CELL > DebugMesh;
  DebugMesh* debugMesh = new DebugMesh(3);

  // Add triangle to mesh
  debugMesh->addNode( ptX[0], ptX[1], ptX[2]);
  debugMesh->addNode( ptY[0], ptY[1], ptY[2]);
  debugMesh->addNode( ptZ[0], ptZ[1], ptZ[2]);

  int tArr[3] = {0,1,2};
  debugMesh->addCell(tArr, MINT_TRIANGLE, 3);

  PointType bbMin(-0.1);
  PointType bbMax(1.1);
  BoundingBoxType bbox( bbMin,bbMax );

  typedef axom::quest::SpatialOctree<DIM, quest::BlockData> SpaceOctree;
  SpaceOctree oct( bbox);


  int lev =3;
  for(int i=0 ; i< 1<<lev ; ++i)
  {
    for(int j=0 ; j< 1<<lev ; ++j)
    {
      for(int k=0 ; k< 1<<lev ; ++k)
      {
        SpaceOctree::BlockIndex block(
        axom::primal::Point<int,3>::make_point(i,j,k), lev );
        SpaceOctree::GeometricBoundingBox blockBB = oct.blockBoundingBox(block);

        if( axom::primal::intersect( unitTri, blockBB))
        {
          for(int k=0; k< 1<<lev; ++k)
          {
            SpaceOctree::BlockIndex block( axom::primal::Point<int,3>::make_point(i,j,k), lev );
            SpaceOctree::GeometricBoundingBox blockBB = oct.blockBoundingBox(block);

            if( axom::primal::intersect( unitTri, blockBB))
            {
              // Add to debug mesh
              int vStart = debugMesh->getMeshNumberOfNodes();

              debugMesh->addNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
              debugMesh->addNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMin()[2]);
              debugMesh->addNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMin()[2]);
              debugMesh->addNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMin()[2]);

              debugMesh->addNode( blockBB.getMin()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
              debugMesh->addNode( blockBB.getMax()[0], blockBB.getMin()[1], blockBB.getMax()[2]);
              debugMesh->addNode( blockBB.getMax()[0], blockBB.getMax()[1], blockBB.getMax()[2]);
              debugMesh->addNode( blockBB.getMin()[0], blockBB.getMax()[1], blockBB.getMax()[2]);

              int data[8];
              for(int i=0; i< 8; ++i)
              {
                data[i] = vStart + i;
              }

              debugMesh->addCell( data, MINT_HEX, 8);
            }
          }
        }
      }
    }
  }


  axom::mint::write_vtk(debugMesh, "gridIntersections.vtk");

  delete debugMesh;
}


void testContainmentOnRegularGrid(
  const Octree3D& inOutOctree,
  const GeometricBoundingBox& queryBounds,
  int gridRes)
{
  SpaceVector h( queryBounds.getMin(), queryBounds.getMax());
  for(int i=0 ; i<3 ; ++i)
    h[i] /= gridRes;

  int ext[6];
  ext[0] = ext[2] = ext[4] = 0;
  ext[1] = ext[3] = ext[5] = gridRes;

  axom::mint::UniformMesh* umesh =
    new axom::mint::UniformMesh(3,queryBounds.getMin().data(),h.data(),ext);


  const int nnodes = umesh->getNumberOfNodes();
  axom::mint::FieldData* PD = umesh->getNodeFieldData();
  SLIC_ASSERT( PD != AXOM_NULLPTR );

  PD->addField( new axom::mint::FieldVariable< int >("containment",nnodes) );
  int* containment = PD->getField( "containment" )->getIntPtr();
  SLIC_ASSERT( containment != AXOM_NULLPTR );


  axom::utilities::Timer timer(true);
  for ( int inode=0 ; inode < nnodes ; ++inode )
  {
    axom::primal::Point< double,3 > pt;
    umesh->getMeshNode( inode, pt.data() );

    containment[ inode ] = inOutOctree.within(pt) ? 1 : 0;
  }
  timer.stop();
  SLIC_INFO(
    fmt::format("\tQuerying {}^3 containment field "
                "took {} seconds (@ {} queries per second)",
                gridRes, timer.elapsed(), nnodes / timer.elapsed()));

  #ifdef DUMP_VTK_MESH
  std::stringstream sstr;
  sstr << "gridContainment_" << gridRes << ".vtk";
  axom::mint::write_vtk( umesh, sstr.str());
  #endif

  delete umesh;
}


/**
 * \brief Extracts the vertex indices of cell cellIndex from the mesh
 */
TriVertIndices getTriangleVertIndices(axom::mint::Mesh* mesh, int cellIndex)
{
  SLIC_ASSERT(mesh != AXOM_NULLPTR);
  SLIC_ASSERT(cellIndex >= 0 && cellIndex < mesh->getMeshNumberOfCells());

  TriVertIndices tvInd;
  mesh->getMeshCell( cellIndex, tvInd.data() );
  return tvInd;
}

/**
 * \brief Extracts the positions of a traingle's vertices from the mesh
 * \return The triangle vertex positions in a SpaceTriangle instance
 */
SpaceTriangle getMeshTriangle(axom::mint::Mesh* mesh,
                              const TriVertIndices& vertIndices )
{
  SLIC_ASSERT(mesh != AXOM_NULLPTR);

  SpaceTriangle tri;
  for(int i=0 ; i< 3 ; ++i)
    mesh->getMeshNode( vertIndices[i], tri[i].data() );

  return tri;
}

/**
 * \brief Computes some statistics about the surface mesh.
 *
 * Specifically, computes histograms (and ranges) of the edge lengths and
 * triangle areas son a logarithmic scale and logs the results
 */
void print_surface_stats( axom::mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpacePt pt;

  typedef axom::primal::BoundingBox<double,1> MinMaxRange;
  typedef MinMaxRange::PointType LengthType;

  MinMaxRange meshEdgeLenRange;
  MinMaxRange meshTriAreaRange;
  const int nCells = mesh->getMeshNumberOfCells();
  typedef std::set<int> TriIdxSet;
  TriIdxSet badTriangles;

  // simple binning based on the exponent
  typedef std::map<int,int> LogHistogram;
  LogHistogram edgeLenHist;         // Create histogram of edge lengths (log
                                    // scale)
  LogHistogram areaHist;            // Create histogram of triangle areas (log
                                    // scale)

  typedef std::map<int,MinMaxRange> LogRangeMap;
  LogRangeMap edgeLenRangeMap;      // Tracks range of edge lengths at each
                                    // scale
  LogRangeMap areaRangeMap;         // Tracks range of triangle areas at each
                                    // scale

  typedef axom::primal::Point<int,3> TriVertIndices;
  int expBase2;

  // Traverse mesh triangles and bin the edge lengths and areas
  for ( int i=0 ; i < nCells ; ++i )
  {
    // Get the indices and positions of the triangle's three vertices
    TriVertIndices vertIndices = getTriangleVertIndices(mesh, i);
    SpaceTriangle tri = getMeshTriangle(mesh, vertIndices);

    // Compute edge stats -- note edges are double counted
    for(int j=0 ; j<3 ; ++j)
    {
      double len = SpaceVector(tri[j],tri[(j+1)%3]).norm();
      if(axom::utilities::isNearlyEqual(len,0.))
      {
        badTriangles.insert(i);
      }
      else
      {
        LengthType edgeLen(len);
        meshEdgeLenRange.addPoint( edgeLen );
        std::frexp (len, &expBase2);
        edgeLenHist[expBase2]++;
        edgeLenRangeMap[expBase2].addPoint( edgeLen );
      }
    }

    // Compute triangle area stats
    double area = tri.area();
    if( axom::utilities::isNearlyEqual(area, 0.))
    {
      badTriangles.insert(i);
    }
    else
    {
      LengthType triArea(area);
      meshTriAreaRange.addPoint ( triArea );
      std::frexp (area, &expBase2);
      areaHist[expBase2]++;
      areaRangeMap[expBase2].addPoint( triArea);
    }
  }


  // Log the results
  const int nVerts = mesh->getMeshNumberOfNodes();
  SLIC_INFO(fmt::format("Mesh has {} vertices  and {} triangles.", nVerts,
                        nCells));

  SLIC_INFO("Edge length range: " << meshEdgeLenRange);
  SLIC_INFO("Triangle area range is: " << meshTriAreaRange);

  fmt::MemoryWriter edgeHistStr;
  edgeHistStr<<"Edge length histogram (lg-arithmic): ";
  for(LogHistogram::const_iterator it = edgeLenHist.begin()
      ; it != edgeLenHist.end()
      ; ++it)
  {
    edgeHistStr.write("\n\texp: {}\tcount: {}\tRange: {}",
                      it->first, it->second / 2, edgeLenRangeMap[it->first]);
  }
  SLIC_DEBUG(edgeHistStr.str());

  fmt::MemoryWriter triHistStr;
  triHistStr<<"Triangle areas histogram (lg-arithmic): ";
  for(LogHistogram::const_iterator it =areaHist.begin()
      ; it != areaHist.end()
      ; ++it)
  {
    triHistStr.write("\n\texp: {}\tcount: {}\tRange: {}",
                     it->first, it->second, areaRangeMap[it->first]);
  }
  SLIC_DEBUG(triHistStr.str());

  if(!badTriangles.empty() )
  {
    fmt::MemoryWriter badTriStr;
    badTriStr<<"The following triangle(s) have zero area/edge lengths:";
    for(TriIdxSet::const_iterator it = badTriangles.begin()
        ; it != badTriangles.end()
        ; ++it)
    {
      badTriStr<< "\n\tTriangle " << *it;
      TriVertIndices vertIndices;
      mesh->getMeshCell( *it, vertIndices.data() );

      SpacePt vertPos;
      for(int j=0 ; j<3 ; ++j)
      {
        mesh->getMeshNode( vertIndices[j], vertPos.data() );
        badTriStr.write("\n\t\t vId: {} @ position: {}",
                        vertIndices[j], vertPos);
      }
    }
    SLIC_DEBUG(badTriStr.str());
  }
}

/**
 * \brief Finds the octree leaf containing the given query point, and optionally
 *  refines the leaf
 */
void refineAndPrint(Octree3D& octree, const SpacePt& queryPt,
                    bool shouldRefine = true)
{
  BlockIndex leafBlock = octree.findLeafBlock(queryPt);

  if(shouldRefine)
  {
    octree.refineLeaf( leafBlock );
    leafBlock = octree.findLeafBlock(queryPt);
  }

  GeometricBoundingBox blockBB = octree.blockBoundingBox( leafBlock);
  bool containsPt = blockBB.contains(queryPt);

  SLIC_INFO(
    fmt::format("\t(gridPt: {}; lev: {}) with bounds {} {} query point.",
                leafBlock.pt(), leafBlock.level(), blockBB,
                (containsPt ? " contains " : "does not contain ") ));
}

//------------------------------------------------------------------------------
int main( int argc, char** argv )
{
  axom::slic::UnitTestLogger logger;  // create & initialize logger
  // axom::slic::debug::checksAreErrors = true;

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

    stlFile =
      axom::utilities::filesystem::joinPath(defaultDir, defaultFileName);
  }

  stlFile = axom::slam::util::findFileInAncestorDirs(stlFile);
  SLIC_ASSERT( axom::utilities::filesystem::pathExists( stlFile));

  // STEP 2: read mesh file
  SLIC_INFO(fmt::format("\n\t{:*^80}"," Loading the mesh "));
  SLIC_INFO("Reading file: " << stlFile << "...");

  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( stlFile );
  reader->read();
  SLIC_INFO("done.");


  // STEP 3: create surface mesh
  axom::mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader->getMesh( static_cast<TriangleMesh*>( surface_mesh ) );
  // dump mesh info
  SLIC_INFO("Mesh has "
            << surface_mesh->getMeshNumberOfNodes() << " nodes and "
            << surface_mesh->getMeshNumberOfCells() << " cells.");

  // STEP 4: Delete the reader
  delete reader;
  reader = AXOM_NULLPTR;


  // STEP 5: Compute the bounding box and log some stats about the surface
  GeometricBoundingBox meshBB = compute_bounds( surface_mesh);
  SLIC_INFO( "Mesh bounding box: " << meshBB );
  print_surface_stats(surface_mesh);

  testIntersectionOnRegularGrid();

  axom::slic::setLoggingMsgLevel( axom::slic::message::Debug);


  // STEP 6: Create octree over mesh's bounding box and query a point in space
  SLIC_INFO(fmt::format("\n\t{:*^80}"," Generating the octree "));
  Octree3D octree(meshBB, surface_mesh);
  octree.generateIndex();


  print_surface_stats(surface_mesh);
  axom::mint::write_vtk(surface_mesh, "meldedTriMesh.vtk");

  SLIC_INFO(fmt::format("\n\t{:*^80}"," Querying the octree "));

  // Query on a slightly expanded bounding box
  GeometricBoundingBox queryBB
    = octree.boundingBox();

  // Bounding box for a problematic region on the plane_simp.stl model
  // = GeometricBoundingBox(  SpacePt::make_point(68.326,740.706,187.349),
  //                          SpacePt::make_point(68.5329,740.923,187.407));

  // We can scale the query region here
  //  queryBB.scale(1.1);


  // STEP 7: Query the mesh
  for(int i=1 ; i< MAX_CONTAINMENT_QUERY_LEVEL ; ++i)
    testContainmentOnRegularGrid( octree, queryBB, 1<<i);

  axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);

  // Other -- find leaf block of a given query point at various levels of
  // resolution
  SLIC_INFO(fmt::format("\n\t{:*^80}"," Other octree operations "));
  double alpha = 2./3.;
  SpacePt queryPt = SpacePt::lerp(meshBB.getMin(), meshBB.getMax(), alpha);

  SLIC_INFO("Finding associated grid point for query point: " << queryPt );
  for(int lev = 0 ; lev < octree.maxLeafLevel() ; ++lev)
  {
    GridPt gridPt = octree.findGridCellAtLevel(queryPt, lev);
    SLIC_INFO(
      fmt::format(
        "  {1} @ level {0}\n\t[max gridPt: {2}; spacing: {3};\n\t bounding box {4}]",
        lev,
        gridPt,
        octree.maxGridCellAtLevel(lev),
        octree.spacingAtLevel(lev),
        octree.blockBoundingBox(gridPt, lev)
        ));
  }


  SLIC_INFO("Recursively refining around query point: " << queryPt);
  refineAndPrint(octree, queryPt, false);
//  for(int i=0; i< octree.maxInternalLevel(); ++i)
//      refineAndPrint(octree, queryPt);

  // STEP 8: Reclaim memory
  if(surface_mesh != AXOM_NULLPTR)
  {
    delete surface_mesh;
    surface_mesh = AXOM_NULLPTR;
  }

  return 0;
}
