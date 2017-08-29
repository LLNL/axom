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

#include "mint/Mesh.hpp"
#include "mint/UniformMesh.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/intersect.hpp"
#include "primal/Point.hpp"
#include "primal/Triangle.hpp"
#include "primal/UniformGrid.hpp"

#include "quest/MeshTester.hpp"
#include "quest/STLReader.hpp"

#include "slic/GenericOutputStream.hpp"
#include "slic/slic.hpp"


// C/C++ includes
#include <iostream>
#include <utility>
#include <vector>

using namespace axom;

typedef mint::UnstructuredMesh< MINT_TRIANGLE > TriangleMesh;
typedef primal::Triangle<double, 3> Triangle3;

typedef primal::Point<double, 3> Point3;
typedef primal::BoundingBox<double, 3> SpatialBoundingBox;
typedef primal::UniformGrid<int, 3> UniformGrid3;
typedef primal::Vector<double, 3> Vector3;
typedef primal::Segment<double, 3> Segment3;


enum InputStatus
{
  SUCCESS,
  SHOWHELP,
  CANTOPENFILE
};


struct Input
{
  std::string stlInput;
  std::string textOutput;
  int resolution;
  InputStatus errorCode;

  Input() : stlInput(""),
            textOutput("meshTestResults.txt"),
            resolution(0),
            errorCode(SUCCESS)
  { };

  Input(int argc, char ** argv);

  void showhelp()
  {
    std::cout << "Argument usage:" << std::endl <<
      "  --help           Show this help message." << std::endl <<
      "  --resolution N   Resolution of uniform grid.  Default N = 10.  " << std::endl <<
      "       Set to 1 to run the naive algorithm, without the spatial index." << std::endl <<
      "       Set to less than 1 to use the spatial index with a resolution of the" << std::endl <<
      "         cube root of the number of triangles." << std::endl <<
      "  --infile fname   The STL input file (must be specified)." << std::endl <<
      "  --outfile fname  The text output file (defaults to meshTestResults.txt)." << std::endl <<
      std::endl;
  };
};

Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh);
inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);
bool checkTT(Triangle3& t1, Triangle3& t2);
std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
                                                              std::vector<int> & degenerate);
bool canOpenFile(const std::string & fname);
bool writeCollisions(const std::vector< std::pair<int, int> > & c,
                     const std::vector<int> & d,
                     const std::string & outfile);

Input::Input(int argc, char ** argv) :
    stlInput(""),
    textOutput("meshTestResults.txt"),
    resolution(0),
    errorCode(SUCCESS)
{
    if (argc < 2) {
        errorCode = SHOWHELP;
        return;
    } else {
        std::string help = argv[1];
        if(help == "--help") {
            errorCode = SHOWHELP;
            return;
        }
        for (int i = 1; i < argc; /* increment i in loop */){
            std::string arg = argv[i];
            if (arg == "--resolution"){
                resolution = atoi(argv[++i]);
            }
            else if (arg == "--infile"){
                stlInput = argv[++i];
            }
            else if (arg == "--outfile"){
                textOutput = argv[++i];
            }
            ++i;
        }
    }

    if (!canOpenFile(stlInput)) {
        errorCode = CANTOPENFILE;
        return;
    }

    SLIC_INFO ("Using parameter values: " << std::endl <<
            "  resolution = " << resolution <<
            (resolution < 1? " (use cube root of triangle count)": "") <<
            std::endl <<
            "  infile = " << stlInput << std::endl <<
            "  outfile = " << textOutput << std::endl);
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

std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(mint::Mesh* surface_mesh,
  std::vector<int> & degenerate)
{
  // For each triangle, check for intersection against
  // every other triangle with a greater index in the mesh, excluding
  // degenerate triangles.
  std::vector< std::pair<int, int> > retval;

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
        retval.push_back(std::make_pair(i, j));
      }
    }
  }

  return retval;
}

bool canOpenFile(const std::string & fname)
{
  std::ifstream teststream(fname.c_str());
  return teststream.good();
}

bool writeCollisions(const std::vector< std::pair<int, int> > & c,
  const std::vector<int> & d,
  const std::string & outfile)
{
  std::ofstream outf(outfile.c_str());
  if (!outf) {
    return false;
  }

  outf << c.size() << " intersecting triangle pairs:" << std::endl;
  for (size_t i = 0; i < c.size(); ++i) {
    outf << c[i].first << " " << c[i].second << std::endl;
  }

  outf << d.size() << " degenerate triangles:" << std::endl;
  for (size_t i = 0; i < d.size(); ++i) {
    outf << d[i] << std::endl;
  }

  return true;
}


/*!
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

  // Initialize default parameters and update with command line arguments:
  Input params(argc, argv);

  if (params.errorCode != SUCCESS) {
    if (params.errorCode == SHOWHELP) {
      params.showhelp();
      return EXIT_SUCCESS;
    } else if (params.errorCode == CANTOPENFILE) {
      std::cerr << "Can't open STL file " << params.stlInput <<
        " for reading." << std::endl;
      return EXIT_FAILURE;
    } else {
      std::cerr << "Unknown error " << (int)params.errorCode <<
        " while parsing arguments." << std::endl;
      return EXIT_FAILURE;
    }
  }

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

  std::vector< std::pair<int, int> > collisions;
  std::vector<int> degenerate;
  int status = 0;
  if (params.resolution == 1) {
    // Naive method
    collisions = naiveIntersectionAlgorithm(surface_mesh, degenerate);
  } else {
    // Use a spatial index
    status = quest::findTriMeshIntersections(surface_mesh,
                                             collisions,
                                             degenerate,
                                             params.resolution);
  }
  if (!writeCollisions(collisions, degenerate, params.textOutput)) {
    SLIC_ERROR("Couldn't write results to " << params.textOutput);
  }

  return retval;
}

