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

// Axom includes
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/Timer.hpp"

#include "axom/mint/mesh/FieldVariable.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp" // for write_vtk

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/spatial_acceleration/UniformGrid.hpp"

#include "axom/quest/MeshTester.hpp"
#include "axom/quest/stl/STLReader.hpp"

#include "axom/slic/streams/GenericOutputStream.hpp"
#include "axom/slic/interface/slic.hpp"


// C/C++ includes
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

using namespace axom;

typedef mint::UnstructuredMesh< mint::SINGLE_SHAPE > UMesh;
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
  std::string vtkOutput;

  int resolution;
  double weldThreshold;
  bool skipWeld;
  InputStatus errorCode;

  Input() :
    stlInput(""),
    vtkOutput(""),
    resolution(0),
    weldThreshold(1e-6),
    skipWeld(false),
    errorCode(SUCCESS)
  { };

  Input(int argc, char** argv);

  void showhelp()
  {
    std::cout
      << "Argument usage:"
       "\n  --help           Show this help message."
       "\n  --resolution N   Resolution of uniform grid.  Default N = 0."
       "\n                   - Set to 1 to run the naive algorithm, without"
       "\n                   a spatial index."
       "\n                   - Set to less than 1 to use the spatial index"
       "\n                   with a resolution of the cube root of the"
       "\n                   number of triangles."
       "\n  --infile fname   The STL input file (must be specified)."
       "\n  --outfile fname  Output file name for collisions and welded mesh"
       "\n                   (defaults to a file in the CWD with the input file name)"
       "\n                   Collisions mesh will end with '.collisions.vtk' and"
       "\n                   welded mesh will end with '.welded.vtk'."
       "\n  --weldThresh eps Distance threshold for welding vertices. "
       "\n                   Default: eps = 1e-6"
       "\n  --skipWeld       Don't weld vertices (useful for testing,"
       "\n                   not helpful otherwise)."
      << std::endl << std::endl;
  };

  std::string collisionsMeshName() { return vtkOutput + ".collisions.vtk"; }
  std::string collisionsTextName() { return vtkOutput + ".collisions.txt"; }
  std::string weldMeshName() { return vtkOutput + ".weld.vtk"; }
};

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);
bool checkTT(Triangle3& t1, Triangle3& t2);
std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(
  mint::Mesh* surface_mesh,
  std::vector<int> & degenerate);
bool canOpenFile(const std::string & fname);
void announceMeshProblems(int triangleCount,
                          int intersectPairCount,
                          int degenerateCount);
void saveProblemFlagsToMesh(mint::Mesh* surface_mesh,
                            const std::vector< std::pair<int, int> > & c,
                            const std::vector<int> & d);
bool writeAnnotatedMesh(mint::Mesh* surface_mesh,
                        const std::string & outfile);
bool writeCollisions(const std::vector< std::pair<int, int> > & c,
                     const std::vector<int> & d,
                     std::string basename);

Input::Input(int argc, char** argv) :
  stlInput(""),
  vtkOutput(""),
  resolution(0),
  weldThreshold(1e-6),
  skipWeld(false),
  errorCode(SUCCESS)
{
  if (argc < 2)
  {
    errorCode = SHOWHELP;
    return;
  }
  else
  {
    for (int i = 1 ; i < argc ; /* increment i in loop */)
    {
      std::string arg = argv[i];
      if (arg == "--resolution")
      {
        resolution = atoi(argv[++i]);
      }
      else if (arg == "--infile")
      {
        stlInput = argv[++i];
      }
      else if (arg == "--outfile")
      {
        vtkOutput = argv[++i];
      }
      else if (arg == "--weldThresh")
      {
        weldThreshold = atof(argv[++i]);
      }
      else if (arg == "--skipWeld")
      {
        skipWeld = true;
      }
      else // help or unknown parameter
      {
        if(arg != "--help" && arg != "-h")
        {
          SLIC_WARNING("Unrecognized parameter: " << arg);
        }

        errorCode = SHOWHELP;
        return;
      }
      ++i;
    }
  }

  if (!canOpenFile(stlInput))
  {
    errorCode = CANTOPENFILE;
    return;
  }

  // Set the output file name
  {
    // Extract the stem of the input file, so can output files in the CWD
    std::string inFileDir;
    axom::utilities::filesystem::getDirName(inFileDir, stlInput);
    int separatorSkip = (int)(inFileDir.size() > 0);
    std::string inFileStem = stlInput.substr(inFileDir.size() + separatorSkip);
    std::string outFileBase = axom::utilities::filesystem::joinPath(
      axom::utilities::filesystem::getCWD(),inFileStem);

    // set output file name when not provided
    int sz = vtkOutput.size();
    if (sz < 1)
    {
      vtkOutput = outFileBase;
      sz = vtkOutput.size();
    }

    // ensure that output file does not end with '.vtk'
    if(sz > 4)
    {
      std::string ext = vtkOutput.substr(sz-4, 4);
      if( ext == ".vtk" || ext == ".stl")
      {
        vtkOutput = vtkOutput.substr(0, sz-ext.size());
      }
    }
  }

  SLIC_INFO (
    "Using parameter values: "
    <<"\n  resolution = " << resolution
    << (resolution < 1 ? " (use cube root of triangle count)" : "")
    <<"\n  weld threshold = " <<  weldThreshold
    <<"\n  " << (skipWeld ? "" : "not ") << "skipping weld"
    <<"\n  infile = " << stlInput
    <<"\n  collisions outfile = " << collisionsMeshName()
    <<"\n  weld outfile = " << weldMeshName()  );
}

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS=1.0e-9)
{
  return axom::utilities::isNearlyEqual(p1[0], p2[0], EPS) && \
         axom::utilities::isNearlyEqual(p1[1], p2[1], EPS) && \
         axom::utilities::isNearlyEqual(p1[2], p2[2], EPS);
}

bool checkTT(Triangle3& t1, Triangle3& t2)
{
  if (t2.degenerate())
    return false;

  if (primal::intersect(t1, t2))
  {
    return true;
  }
  return false;
}

inline Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getCellType( i ) == mint::TRIANGLE );
  primal::Point<mint::IndexType, 3> triCell;
  Triangle3 tri;
  surface_mesh->getCellNodeIDs(i, triCell.data());

  surface_mesh->getNode(triCell[0], tri[0].data());
  surface_mesh->getNode(triCell[1], tri[1].data());
  surface_mesh->getNode(triCell[2], tri[2].data());

  return tri;
}

std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(
  mint::Mesh* surface_mesh,
  std::vector<int> & degenerate)
{
  // For each triangle, check for intersection against
  // every other triangle with a greater index in the mesh, excluding
  // degenerate triangles.
  std::vector< std::pair<int, int> > retval;

  const int ncells = surface_mesh->getNumberOfCells();
  SLIC_INFO("Checking mesh with a total of "<< ncells<< " cells.");

  Triangle3 t1 = Triangle3();
  Triangle3 t2 = Triangle3();

  // For each triangle in the mesh
  for (int i = 0 ; i< ncells ; i++)
  {
    t1 = getMeshTriangle(i, surface_mesh);

    // Skip if degenerate
    if (t1.degenerate())
    {
      degenerate.push_back(i);
      continue;
    }

    // If the triangle is not degenerate, test against all other
    // triangles that this triangle has not been checked against
    for (int j = i + 1 ; j < ncells ; j++)
    {
      t2 = getMeshTriangle(j, surface_mesh);
      if (checkTT(t1, t2))
      {
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

void announceMeshProblems(int triangleCount,
                          int intersectPairCount,
                          int degenerateCount)
{
  std::cout << triangleCount << " triangles, with " << intersectPairCount
            << " intersecting tri pairs, " << degenerateCount
            <<  " degenerate tris." <<  std::endl;
}

void saveProblemFlagsToMesh(mint::Mesh* mesh,
                            const std::vector< std::pair<int, int> > & c,
                            const std::vector<int> & d)
{
  // Create new Field variables to hold degenerate and intersecting info
  const int num_cells = mesh->getNumberOfCells();
  int* intersectptr =
    mesh->createField< int >("nbr_intersection", mint::CELL_CENTERED );
  int* dgnptr =
    mesh->createField< int >("degenerate_triangles", mint::CELL_CENTERED );

  // Initialize everything to 0
  for (int i = 0 ; i < num_cells ; ++i)
  {
    intersectptr[i] = 0;
    dgnptr[i] = 0;
  }

  // Fill in intersect flag
  size_t csize = c.size();
  for (size_t i = 0 ; i < csize ; ++i)
  {
    std::pair<int, int> theC = c[i];
    intersectptr[theC.first] += 1;
    intersectptr[theC.second] += 1;
  }

  // Fill in degenerate flag
  size_t dsize = d.size();
  for (size_t i = 0 ; i < dsize ; ++i)
  {
    dgnptr[d[i]] = 1;
  }
}

bool writeAnnotatedMesh(mint::Mesh* surface_mesh,
                        const std::string & outfile)
{
  return write_vtk(surface_mesh, outfile) == 0;
}

bool writeCollisions(const std::vector< std::pair<int, int> > & c,
                     const std::vector<int> & d,
                     std::string filename)
{
  std::ofstream outf(filename.c_str());
  if (!outf)
  {
    return false;
  }

  outf << c.size() << " intersecting triangle pairs:" << std::endl;
  for (size_t i = 0 ; i < c.size() ; ++i)
  {
    outf << c[i].first << " " << c[i].second << std::endl;
  }

  outf << d.size() << " degenerate triangles:" << std::endl;
  for (size_t i = 0 ; i < d.size() ; ++i)
  {
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

  if (params.errorCode != SUCCESS)
  {
    if (params.errorCode == SHOWHELP)
    {
      params.showhelp();
      return EXIT_SUCCESS;
    }
    else if (params.errorCode == CANTOPENFILE)
    {
      std::cerr << "Can't open STL file " << params.stlInput <<
        " for reading." << std::endl;
      return EXIT_FAILURE;
    }
    else
    {
      std::cerr << "Unknown error " << (int)params.errorCode <<
        " while parsing arguments." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // Read file
  SLIC_INFO("Reading file: '" <<  params.stlInput << "'...\n");
  quest::STLReader* reader = new quest::STLReader();
  reader->setFileName( params.stlInput );
  reader->read();

  // Get surface mesh
  UMesh* surface_mesh = new UMesh( 3, mint::TRIANGLE );
  reader->getMesh( surface_mesh );

  // Delete the reader
  delete reader;
  reader = nullptr;

  SLIC_INFO(
    "Mesh has " << surface_mesh->getNumberOfNodes() << " vertices and "
                <<  surface_mesh->getNumberOfCells() << " triangles.");

  // Vertex welding
  if (!params.skipWeld)
  {
    axom::utilities::Timer timer(true);

    quest::weldTriMeshVertices(&surface_mesh, params.weldThreshold);

    timer.stop();
    SLIC_INFO("Vertex welding took "
              << timer.elapsedTimeInSec() << " seconds.");
    SLIC_INFO("After welding, mesh has "
              << surface_mesh->getNumberOfNodes() << " vertices and "
              <<  surface_mesh->getNumberOfCells() << " triangles.");

    mint::write_vtk(surface_mesh, params.weldMeshName() );
  }

  // Detect collisions
  {
    std::vector< std::pair<int, int> > collisions;
    std::vector<int> degenerate;

    axom::utilities::Timer timer(true);
    if (params.resolution == 1)
    {
      // Naive method
      collisions = naiveIntersectionAlgorithm(surface_mesh, degenerate);
    }
    else
    {
      // Use a spatial index
      quest::findTriMeshIntersections(surface_mesh,
                                      collisions,
                                      degenerate,
                                      params.resolution);
    }
    timer.stop();
    SLIC_INFO("Detecting intersecting triangles took "
              << timer.elapsedTimeInSec() << " seconds.");

    announceMeshProblems(surface_mesh->getNumberOfCells(),
                         collisions.size(), degenerate.size());

    saveProblemFlagsToMesh(surface_mesh, collisions, degenerate);

    if (!writeAnnotatedMesh(surface_mesh, params.collisionsMeshName()) )
    {
      SLIC_ERROR("Couldn't write results to " << params.collisionsMeshName());
    }

    if (!writeCollisions(collisions,degenerate,params.collisionsTextName()) )
    {
      SLIC_ERROR("Couldn't write results to "<< params.collisionsTextName());
    }
  }

  // Look for holes---must be welded
  if (params.skipWeld)
  {
    SLIC_INFO("Watertight check depends on vertex welding, which was skipped.");
  }
  else
  {
    SLIC_INFO("Checking for watertight mesh.");
    axom::utilities::Timer timer2(true);
    quest::WatertightStatus wtstat =
      quest::isSurfaceMeshWatertight(surface_mesh);
    timer2.stop();
    switch (wtstat)
    {
      case quest::WatertightStatus::WATERTIGHT:
        std::cout << "The mesh is watertight." << std::endl;
        break;
      case quest::WatertightStatus::NOT_WATERTIGHT:
        std::cout << "The mesh is not watertight: at least one " <<
          "boundary edge was detected." << std::endl;
        break;
      default:
        std::cout << "An error was encountered while checking." << std::endl <<
          "This may be due to a non-manifold mesh." << std::endl;
        break;
    }
    SLIC_INFO("Testing for watertightness took "
              << timer2.elapsedTimeInSec() << " seconds.");
  }

  // Delete the mesh
  delete surface_mesh;
  surface_mesh = nullptr;

  return retval;
}
