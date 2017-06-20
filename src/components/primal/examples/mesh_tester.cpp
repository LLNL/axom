/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



// Axom includes
#include "axom_utils/FileUtilities.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/intersection.hpp"
#include "primal/UniformGrid.hpp"

#include "quest/STLReader.hpp"

#include "mint/UniformMesh.hpp"
#include "mint/Mesh.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"


// C/C++ includes
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <cstdio>

using namespace axom;

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;


typedef struct Input
{
  std::string stlInput;
  std::string textOutput;
  int resolution;
  int errorCode;

  Input() : stlInput(""),
            textOutput("meshTestResults.txt"),
            resolution(0),
            errorCode(0)
  { };
} Input;

typedef struct TrianglePair
{
  int a, b;

  TrianglePair(const int na, const int nb) : a(na), b(nb) {};
} TrianglePair;

SpatialBoundingBox compute_bounds(mint::Mesh* mesh);
Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh);
inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);
bool checkTT(Triangle3& t1, Triangle3& t2);
std::vector<TrianglePair> naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
  std::vector<int> & degenerate);
void markSeen(const int a, const std::set<int> & bs, std::vector<TrianglePair> & list);
std::vector<TrianglePair> uGridIntersectionAlgorithm(mint::Mesh* surface_mesh,
  std::vector<int> & degenerate,
  int resolution);
void showhelp();
bool canOpenFile(const std::string & fname);
void init_params(Input& params, int argc, char ** argv);
bool writeCollisions(const std::vector<TrianglePair> & c,
  const std::vector<int> & d,
  const std::string & outfile);

SpatialBoundingBox compute_bounds(mint::Mesh* mesh)
{
  SLIC_ASSERT( mesh != AXOM_NULLPTR );

  SpatialBoundingBox meshBB;
  Point3 pt;

  for ( int i=0; i < mesh->getMeshNumberOfNodes(); ++i )
  {
    mesh->getMeshNode( i, pt.data() );
    meshBB.addPoint( pt );
  } // END for all nodes

  SLIC_ASSERT( meshBB.isValid() );

  return meshBB;
}

Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getMeshNumberOfCellNodes(i) == 3);
  primal::Point<int, 3> triCell;
  surface_mesh->getMeshCell( i, triCell.data());
  primal::Point< double,3 > A1;
  primal::Point< double,3 > B1;
  primal::Point< double,3 > C1;

  surface_mesh->getMeshNode(triCell[0], A1.data());
  surface_mesh->getMeshNode(triCell[1], B1.data());
  surface_mesh->getMeshNode(triCell[2], C1.data());
  Triangle3 t1 = Triangle3(A1,B1,C1);

  return t1;
}

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS=1.0e-9)
{
  return axom::utilities::isNearlyEqual(p1[0], p2[0], EPS) && \
    axom::utilities::isNearlyEqual(p1[1], p2[1], EPS) && \
    axom::utilities::isNearlyEqual(p1[2], p2[2], EPS);
}

bool checkTT(Triangle3& t1, Triangle3& t2)
{
  if (t2.degenerate()) return false;

  if (primal::intersect(t1, t2)) {
    return true;
  }
  return false;
}

std::vector<TrianglePair> naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
  std::vector<int> & degenerate)
{
  // For each triangle, check for intersection against
  // every other triangle with a greater index in the mesh, excluding
  // degenerate triangles.
  std::vector<TrianglePair> retval;

  const int ncells = surface_mesh->getMeshNumberOfCells();
  SLIC_INFO("Checking mesh with a total of "<< ncells<< " cells.");

  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();

  // For each triangle in the mesh
  for (int i = 0; i< ncells; i++) {
    t1 = getMeshTriangle(i, surface_mesh);

    // Skip if degenerate
    if (t1.degenerate()) {
      degenerate.push_back(i);
      continue;
    }

    // If the triangle is not degenerate, test against all other
    // triangles that this triangle has not been checked against
    for (int j = i + 1; j < ncells; j++) {
      t2 = getMeshTriangle(j, surface_mesh);
      if (checkTT(t1, t2)) {
        retval.push_back(TrianglePair(i, j));
      }
    }
  }

  return retval;
}

void markSeen(const int a, const std::set<int> & bs, std::vector<TrianglePair> & list)
{
  list.reserve(list.size() + bs.size());
  std::set<int>::const_iterator sit = bs.begin(), send = bs.end();
  for ( ; sit != send; ++sit) {
    list.push_back(TrianglePair(a, *sit));
  }
}

std::vector<TrianglePair> uGridIntersectionAlgorithm(mint::Mesh* surface_mesh,
  std::vector<int> & degenerate,
  int resolution)
{
  std::vector<TrianglePair> retval;
  std::set<int> seen, hit;

  int intersectingTriangleCount=0;
  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();
  SLIC_INFO("Running mesh_tester with virtual_grid");

  // Create a bounding box around mesh to find the minimum point
  SpatialBoundingBox meshBB  = compute_bounds(surface_mesh);
  const Point3 minBBPt= meshBB.getMin();
  const Point3 maxBBPt= meshBB.getMax();

  const int ncells = surface_mesh->getMeshNumberOfCells();
  const int reportInterval = (int)(ncells / 5.0);

  // find the specified resolution.  If we're passed a number less than one,
  // use the cube root of the number of triangles.
  if (resolution < 1) {
    resolution = (int)(1 + std::pow(ncells, 1/3.));
  }
  int resolutions[3]={resolution,resolution,resolution};

  std::cerr << "Building virtual grid..." << std::endl;
  UniformGrid3 ugrid(minBBPt.data(), maxBBPt.data(), resolutions);

  for (int i=0; i < ncells; i++) {
    if (ncells >= reportInterval && i % (ncells/reportInterval) == 0) {
      std::cerr << "Building grid is " << 100.0 * (double(i)/double(ncells)) <<
        " percent done" << std::endl;
    }
    SpatialBoundingBox triBB;
    t1=getMeshTriangle(i,  surface_mesh);
    triBB.addPoint(t1[0]);
    triBB.addPoint(t1[1]);
    triBB.addPoint(t1[2]);

    ugrid.insert(triBB, i);
  }


  // Iterate through triangle indices z from first index to last index.
  // Check against each other triangle with index greater than the index z
  // that also shares a virtual grid bin.
  SLIC_INFO("Checking mesh with a total of " << ncells << " cells.");
  for (size_t z=0; z< ncells; z++) {
    seen.clear();
    hit.clear();

    if (ncells >= reportInterval && z % (ncells/reportInterval) == 0) {
      std::cerr << "Querying grid is " << 100.0 * (double(z)/double(ncells)) <<
        " percent done" << std::endl;
    }

    // Retrieve the triangle at index z and construct a bounding box around it
    SpatialBoundingBox triBB2;
    t1 = getMeshTriangle(z,  surface_mesh);
 
    if (t1.degenerate()) { 
      degenerate.push_back(z);
      continue;
    }
    triBB2.addPoint(t1[0]);
    triBB2.addPoint(t1[1]);
    triBB2.addPoint(t1[2]);

    Point3 minBBPt2,maxBBPt2;

    minBBPt2 = triBB2.getMin();
    maxBBPt2 = triBB2.getMax();   

    // For each grid bins that this triangle will touch,
    std::vector<int> trianglesInBin;
    const std::vector<int> binsToCheck = ugrid.getBinsForBbox(triBB2);
    for (size_t curbin = 0; curbin < binsToCheck.size(); ++curbin) {
      //Step 7.5: retrieve triangles for the bin
      int bidx = binsToCheck[curbin];
      trianglesInBin = ugrid.getBinContents(bidx);
      int binSize = trianglesInBin.size();

      // For each triangle in the current bin,
      for (size_t l = 0; l < binSize; ++l) {
        int t2Index = trianglesInBin[l];
        SLIC_ASSERT(trianglesInBin[l] >= 0);
        if (t2Index <= z || seen.count(t2Index) > 0) {
          continue;  
        } else { 
          seen.insert(t2Index);

          t2 = getMeshTriangle(t2Index, surface_mesh);
          // Test for intersection.
          if (checkTT(t1, t2)) {
            hit.insert(t2Index);
          }
        }
      }  // end for triangles in bin
    }  // bins to check
    markSeen(z, hit, retval);
  }

  return retval;
}

void showhelp()
{
  std::cout << "Argument usage:" << std::endl <<
    "  --help           Show this help message." << std::endl <<
    "  --resolution N   Resolution of uniform grid.  Default N = 10.  " <<
    "       Set to 1 to run the naive algorithm, without the spatial index." << std::endl <<
    "       Set to less than 1 to use the spatial index with a resolution of the" << std::endl <<
    "         cube root of the number of triangles." << std::endl <<
    "  --infile fname   The STL input file (must be specified)." << std::endl <<
    "  --outfile fname  The text output file (defaults to meshTestResults.txt)." << std::endl <<
    std::endl;
}

bool canOpenFile(const std::string & fname)
{
  std::ifstream teststream(fname);
  return teststream.good();
}

void init_params(Input& params, int argc, char ** argv)
{
  if (argc < 2) {
    showhelp();
    params.errorCode = 1;
    return;
  } else {
    std::string help = argv[1];
    if(help == "--help") {
      showhelp();
      params.errorCode = 1;
      return;
    }
    for (int i = 1; i < argc; /* increment i in loop */){
      std::string arg = argv[i];
      if (arg == "--resolution"){
        params.resolution = atoi(argv[++i]);
      }
      else if (arg == "--infile"){
        params.stlInput = argv[++i];
      }
      else if (arg == "--outfile"){
        params.textOutput = argv[++i];
      }
      ++i;
    }
  }

  if (!canOpenFile(params.stlInput)) {
    params.errorCode = 2;
    return;
  }

  SLIC_INFO ("Using parameter values: " << std::endl << 
      "  resolution = " << params.resolution <<
         (params.resolution < 1? " (use cube root of triangle count)": "") <<
         std::endl <<
      "  infile = " << params.stlInput << std::endl << 
      "  outfile = " << params.textOutput << std::endl);
}

bool writeCollisions(const std::vector<TrianglePair> & c,
  const std::vector<int> & d,
  const std::string & outfile)
{
  std::ofstream outf(outfile);
  if (!outf) {
    return false;
  }

  outf << c.size() << " intersecting triangle pairs:" << std::endl;
  for (int i = 0; i < c.size(); ++i) {
    outf << c[i].a << " " << c[i].b << std::endl;
  }

  outf << d.size() << " degenerate triangles:" << std::endl;
  for (int i = 0; i < d.size(); ++i) {
    outf << d[i] << std::endl;
  }

  return true;
}


/***************************************************************************
 * The mesh tester checks a triangulated surface mesh for several problems.
 *
 * Currently the mesh tester checks for intersecting triangles.  This is
 * implemented in two ways.  First, a naive algorithm tests each triangle
 * against each other triangle.  This is easy to understand and verify,
 * but slow.  A second algorithm uses a UniformGrid to index the bounding
 * box of the mesh, and checks each triangle for intersection with the
 * triangles in all the bins the triangle's bounding box falls into.
 *
 * Currently, the mesh tester works only with Triangle meshes.
 ***************************************************************************
 */

int main( int argc, char** argv )
{
  int retval = EXIT_SUCCESS;

  // Initialize the SLIC logger
  slic::initialize();
  slic::setLoggingMsgLevel( axom::slic::message::Debug );

  // Customize logging levels and formatting
  std::string slicFormatStr = "[<LEVEL>] <MESSAGE> \n";
  slic::GenericOutputStream* defaultStream =
    new slic::GenericOutputStream(&std::cout);
  slic::GenericOutputStream* compactStream =
    new slic::GenericOutputStream(&std::cout, slicFormatStr);
  slic::addStreamToMsgLevel(defaultStream, axom::slic::message::Error);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Warning);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Info);
  slic::addStreamToMsgLevel(compactStream, axom::slic::message::Debug);

  // Initialize default parameters and update parameters with command line arguments:
  Input params;
  init_params(params, argc, argv);

  if (params.errorCode > 0) {
    if (params.errorCode == 1) {
      // user requested help message; don't print anything
      return EXIT_SUCCESS;
    } else if (params.errorCode == 2) {
      std::cerr << "Can't open STL file " << params.stlInput <<
        " for reading." << std::endl;
      return EXIT_FAILURE;
    } else {
      std::cerr << "Error " << params.errorCode <<
        " while parsing arguments." << std::endl;
    }
  }

  SLIC_ASSERT(params.resolution > 0);

  // Read file
  SLIC_INFO("Reading file: " <<  params.stlInput << "...\n");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( params.stlInput );
  reader->read();
  SLIC_INFO("done\n");

  // Get surface mesh
  mint::Mesh* surface_mesh = new TriangleMesh( 3 );
  reader->getMesh( static_cast<TriangleMesh*>( surface_mesh ) );

  // Delete the reader
  delete reader;
  reader = AXOM_NULLPTR;

  std::vector<TrianglePair> collisions;
  std::vector<int> degenerate;
  if (params.resolution == 1) {
    // Naive method
    collisions = naiveIntersectionAlgorithm(surface_mesh, degenerate);
  } else {
    // Use a spatial index
    collisions = uGridIntersectionAlgorithm(surface_mesh, degenerate, params.resolution);
  }
  if (!writeCollisions(collisions, degenerate, params.textOutput)) {
    SLIC_ERROR("Couldn't write results to " << params.textOutput);
  }

  return retval;
}

