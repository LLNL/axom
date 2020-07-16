// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/core/Macros.hpp"
#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/Timer.hpp"

#include "axom/mint/mesh/FieldVariable.hpp"
// _read_stl_include2_start
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
// _read_stl_include2_end
#include "axom/mint/mesh/UniformMesh.hpp"
#include "axom/mint/utils/vtk_utils.hpp" // for write_vtk

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/operators/intersect.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/spin/UniformGrid.hpp"

// _read_stl_include1_start
#include "axom/quest/stl/STLReader.hpp"
// _read_stl_include1_end
// _check_repair_include_start
#include "axom/quest/MeshTester.hpp"
// _check_repair_include_end

#include "axom/slic/streams/GenericOutputStream.hpp"
#include "axom/slic/interface/slic.hpp"

// RAJA
#ifdef AXOM_USE_RAJA
  #include "RAJA/RAJA.hpp"
#endif

// RAJA policies
#include "axom/mint/execution/internal/structured_exec.hpp"

#if defined (AXOM_USE_RAJA)
  using seq_exec = axom::SEQ_EXEC;

  #if defined(AXOM_USE_OPENMP)
    using omp_exec = axom::OMP_EXEC;
  #else
    using omp_exec = seq_exec;
  #endif

  #if defined(AXOM_USE_CUDA)
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
  #else
    using cuda_exec = seq_exec;
  #endif
#endif

#include "CLI11/CLI11.hpp"

// C/C++ includes
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <iomanip>

// _read_stl_typedefs_start
using namespace axom;

using UMesh = mint::UnstructuredMesh< mint::SINGLE_SHAPE >;
// _read_stl_typedefs_end
using Triangle3 = primal::Triangle<double, 3>;

using Point3 = primal::Point<double, 3>;
using SpatialBoundingBox = primal::BoundingBox<double, 3>;
using UniformGrid3 = spin::UniformGrid<int, 3>;
using Vector3 = primal::Vector<double, 3>;
using Segment3 = primal::Segment<double, 3>;

enum RuntimePolicy { seq = 0,
                     raja_seq = 1,
                     raja_omp = 2,
                     raja_cuda = 3};

struct Input
{
  static const std::set<RuntimePolicy> s_validPolicies;

  std::string stlInput {""};
  std::string vtkOutput {""};
  RuntimePolicy policy {seq};

  int resolution {0};
  double weldThreshold {1e-6};
  double intersectionThreshold {1e-08};
  bool skipWeld {false};
  bool verboseOutput {false};

  Input() = default;

  void parse(int argc, char** argv, CLI::App& app);

  std::string collisionsMeshName() { return vtkOutput + ".collisions.vtk"; }
  std::string collisionsTextName() { return vtkOutput + ".collisions.txt"; }
  std::string weldMeshName() { return vtkOutput + ".weld.vtk"; }

private:
  void fixOutfilePath();
};

const std::set<RuntimePolicy> Input::s_validPolicies({
    seq
  #ifdef AXOM_USE_RAJA
    , raja_seq
    #ifdef AXOM_USE_OPENMP
    , raja_omp
    #endif
    #ifdef AXOM_USE_CUDA
    , raja_cuda
    #endif
  #endif
});

void Input::parse(int argc, char** argv, CLI::App& app)
{
  app.add_option("-r,--resolution", resolution,
                 "Resolution of uniform grid. \n"
                 "Set to 1 to run the naive algorithm (without a spatial index). \n"
                 "Set to less than 1 to use the spatial index with a resolution \n"
                 "of the cube root of the number of triangles.")
  ->capture_default_str();

  app.add_option("-p, --policy", policy,
                 "With \'-r 1\', set runtime policy. \n"
                 "Set to 0 to use the sequential algorithm (w/o RAJA). \n"
  #ifdef AXOM_USE_RAJA
                 "Set to 1 to use the RAJA sequential policy. \n"
    #ifdef AXOM_USE_OPENMP
                 "Set to 2 to use the RAJA OpenMP policy. \n"
    #endif
    #ifdef AXOM_USE_CUDA
                 "Set to 3 to use the RAJA CUDA policy."
    #endif
  #endif
  )
  ->capture_default_str()
  ->check(CLI::IsMember{Input::s_validPolicies});

  app.add_option("-i,--infile", stlInput,"The STL input file")
  ->required()
  ->check(CLI::ExistingFile);

  app.add_option("-o,--outfile", vtkOutput,
                 "Output file name for collisions and welded mesh. \n"
                 "Defaults to a file in the CWD with the input file name. \n"
                 "Collisions mesh will end with '.collisions.vtk' and \n"
                 "welded mesh will end with '.welded.vtk'.");

  app.add_option("--weldThresh", weldThreshold,
                 "Distance threshold for welding vertices.")
  ->capture_default_str();

  app.add_option("--intersectionThresh", intersectionThreshold,
                 "Tolerance threshold to use when testing for intersecting triangles")
  ->capture_default_str();

  app.add_flag("--skipWeld", skipWeld,
               "Don't weld vertices (useful for testing, not helpful otherwise).");

  app.add_flag("-v,--verbose", verboseOutput,
               "Increase logging verbosity.")
  ->capture_default_str();

  app.get_formatter()->column_width(35);

  // Could throw an exception
  app.parse(argc, argv);

  // Fix any options that need fixing
  fixOutfilePath();

  // Output parsed information
  SLIC_INFO (
    "Using parameter values: "
    <<"\n  resolution = " << resolution
    << (resolution < 1 ? " (use cube root of triangle count)" : "")
    << (resolution == 1 ? " (use naive algorithm)" : "")
    << (resolution == 1 ? "\n  policy = " : "")
    << (resolution == 1 ? std::to_string(policy) : "")
    << (resolution == 1 && policy == seq ? " (use sequential policy)" : "")
    << (resolution == 1 && policy == raja_seq ? " (use RAJA sequential policy)" : "")
    << (resolution == 1 && policy == raja_omp ? " (use RAJA OpenMP policy)" : "")
    << (resolution == 1 && policy == raja_cuda ? " (use RAJA CUDA policy)" : "")
    <<"\n  weld threshold = " <<  weldThreshold
    <<"\n  " << (skipWeld ? "" : "not ") << "skipping weld"
    <<"\n  intersection tolerance = " <<  intersectionThreshold
    <<"\n  infile = " << stlInput
    <<"\n  collisions outfile = " << collisionsMeshName()
    <<"\n  weld outfile = " << weldMeshName()  );
}

void Input::fixOutfilePath()
{
  namespace futil = axom::utilities::filesystem;

  // Extract the stem of the input file, to output files in the CWD
  std::string inFileDir;
  futil::getDirName(inFileDir, stlInput);
  int sepSkip = (int)(inFileDir.size() > 0);
  std::string inFileStem = stlInput.substr(inFileDir.size() + sepSkip);
  std::string outFileBase = futil::joinPath(futil::getCWD(),inFileStem);

  // set output file name when not provided
  int sz = static_cast< int >( vtkOutput.size() );
  if (sz < 1)
  {
    vtkOutput = outFileBase;
    sz = static_cast< int >( vtkOutput.size() );
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


inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS);

AXOM_HOST_DEVICE bool checkTT(Triangle3& t1, Triangle3& t2, double EPS);

std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(
  mint::Mesh* surface_mesh,
  std::vector<int> & degenerate);

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

inline bool pointIsNearlyEqual(Point3& p1, Point3& p2, double EPS=1.0e-9)
{
  return axom::utilities::isNearlyEqual(p1[0], p2[0], EPS) && \
         axom::utilities::isNearlyEqual(p1[1], p2[1], EPS) && \
         axom::utilities::isNearlyEqual(p1[2], p2[2], EPS);
}

AXOM_HOST_DEVICE
bool checkTT(Triangle3& t1, Triangle3& t2, double EPS)
{
  if (t2.degenerate())
    return false;

  const bool includeBoundaries = false; // only check for internal intersections
  if (primal::intersect(t1, t2, includeBoundaries, EPS))
  {
    return true;
  }
  return false;
}

inline Triangle3 getMeshTriangle(int i, mint::Mesh* surface_mesh)
{
  SLIC_ASSERT(surface_mesh->getCellType( i ) == mint::TRIANGLE );
  primal::Point<axom::IndexType, 3> triCell;
  Triangle3 tri;
  surface_mesh->getCellNodeIDs(i, triCell.data());

  surface_mesh->getNode(triCell[0], tri[0].data());
  surface_mesh->getNode(triCell[1], tri[1].data());
  surface_mesh->getNode(triCell[2], tri[2].data());

  return tri;
}

std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(
  mint::Mesh* surface_mesh,
  std::vector<int> & degenerate,
  double EPS)
{
  SLIC_INFO("Running naive intersection algorithm.");

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
      if (checkTT(t1, t2, EPS))
      {
        retval.push_back(std::make_pair(i, j));
      }
    }
  }

  return retval;
}

#if defined(AXOM_USE_RAJA)
template < typename ExecSpace >
std::vector< std::pair<int, int> > naiveIntersectionAlgorithm(
  mint::Mesh* surface_mesh,
  std::vector<int> & degenerate,
  double EPS)
{
  SLIC_INFO("Running naive intersection algorithm "
    << " in execution Space: "
    << axom::execution_space< ExecSpace >::name());

  // Get allocator
  int allocatorID = axom::execution_space< ExecSpace >::allocatorID();
  axom::setDefaultAllocator( allocatorID );

  std::vector< std::pair<int, int> > retval;

  const int ncells = surface_mesh->getNumberOfCells();
  SLIC_INFO("Checking mesh with a total of "<< ncells<< " cells.");

  Triangle3 * tris = axom::allocate <Triangle3> (ncells);

  // Get each triangle in the mesh and check for degeneracies
  for (int i = 0; i < ncells ; i++)
  {
    tris[i] = getMeshTriangle(i, surface_mesh);
    if (tris[i].degenerate())
    {
      degenerate.push_back(i);
    }
  }

  RAJA::RangeSegment row_range(0, ncells);
  RAJA::RangeSegment col_range(0, ncells);

  using KERNEL_POL =
    typename axom::mint::internal::structured_exec< ExecSpace >::loop2d_policy;
  using REDUCE_POL =
    typename axom::execution_space< ExecSpace >::reduce_policy;
  using ATOMIC_POL =
    typename axom::execution_space< ExecSpace >::atomic_policy;

  RAJA::ReduceSum< REDUCE_POL, int > numIntersect(0);

  // Compute the number of intersections
  RAJA::kernel<KERNEL_POL>( RAJA::make_tuple(col_range, row_range),
    AXOM_LAMBDA(int col, int row) {
    if (row > col)
    {
      if (checkTT (tris[row], tris[col], EPS))
      {
        numIntersect += 1;
      }
    }
  });

  // Allocation to hold intersection pairs and counter to know where to store
  int * intersections =
    axom::allocate <int> (numIntersect.get() * 2);
  int * counter = axom::allocate <int> (1);

  counter[0] = 0;

  // RAJA loop to populate with intersections
  RAJA::kernel<KERNEL_POL>( RAJA::make_tuple(col_range, row_range),
    AXOM_LAMBDA(int col, int row) {
    if (row > col)
    {
      if (checkTT (tris[row], tris[col], EPS))
      {
        auto idx = RAJA::atomicAdd<ATOMIC_POL>(counter, 2);
        intersections[idx] = row;
        intersections[idx+1] = col;
      }
    }
  });

  // Initialize pairs of clashes
  for (auto i = 0 ; i < numIntersect.get() * 2 ; i += 2)
  {
    retval.push_back(std::make_pair(intersections[i], intersections[i + 1]));
  }

  // Deallocate
  axom::deallocate(tris);
  axom::deallocate(intersections);
  axom::deallocate(counter);

  return retval;
}
#endif // defined(AXOM_USE_RAJA)

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

void initializeLogger()
{
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
}

/*!
 * The mesh tester checks a triangulated surface mesh for several problems.
 *
 * Currently the mesh tester checks for intersecting triangles.  This is
 * implemented in three ways.  First, a naive algorithm tests each triangle
 * against each other triangle.  This is easy to understand and verify,
 * but slow. Second, the same naive algorithm is run using raja. A third
 * algorithm uses a UniformGrid to index the bounding box of the mesh,
 * and checks each triangle for intersection with the triangles in all the
 * bins the triangle's bounding box falls into.
 *
 * Currently, the mesh tester works only with Triangle meshes.
 */

int main( int argc, char** argv )
{
  initializeLogger();

  // Parse the command line arguments
  Input params;
  CLI::App app {"MeshTester example"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch (const CLI::ParseError &e)
  {
    return app.exit(e);
  }

  // _read_stl_file_start
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
  // _read_stl_file_end

  // Vertex welding
  if (!params.skipWeld)
  {
    axom::utilities::Timer timer(true);

    // _check_repair_weld_start
    quest::weldTriMeshVertices(&surface_mesh, params.weldThreshold);
    // _check_repair_weld_end

    timer.stop();
    SLIC_INFO("Vertex welding took "
              << timer.elapsedTimeInSec() << " seconds.");
    SLIC_INFO("After welding, mesh has "
              << surface_mesh->getNumberOfNodes() << " vertices and "
              <<  surface_mesh->getNumberOfCells() << " triangles.");
  }

  // Detect collisions
  {
    // _check_repair_intersections_containers_start
    std::vector< std::pair<int, int> > collisions;
    std::vector<int> degenerate;
    // _check_repair_intersections_containers_end

    axom::utilities::Timer timer(true);
    if (params.resolution == 1)
    {
      switch (params.policy)
      {
      case seq:
        collisions = naiveIntersectionAlgorithm(surface_mesh, degenerate, params.intersectionThreshold);
        break;
    #ifdef AXOM_USE_RAJA
      case raja_seq:
        collisions = naiveIntersectionAlgorithm< seq_exec >(surface_mesh,
                                                            degenerate, params.intersectionThreshold);
        break;
     #ifdef AXOM_USE_OPENMP
      case raja_omp:
        collisions = naiveIntersectionAlgorithm< omp_exec >(surface_mesh,
                                                            degenerate, params.intersectionThreshold);
        break;
     #endif
     #ifdef AXOM_USE_CUDA
      case raja_cuda:
        collisions = naiveIntersectionAlgorithm< cuda_exec >(surface_mesh,
                                                             degenerate, params.intersectionThreshold);
        break;
     #endif
    #endif // AXOM_USE_RAJA

     default:
        SLIC_ERROR("Unhandled runtime policy case " << params.policy );
        break;
     }
    }
    else
    {
      // _check_repair_intersections_start
      // Use a spatial index
      quest::findTriMeshIntersections(surface_mesh,
                                      collisions,
                                      degenerate,
                                      params.resolution,
                                      params.intersectionThreshold);
      // _check_repair_intersections_end
    }
    timer.stop();
    SLIC_INFO("Detecting intersecting triangles took "
              << timer.elapsedTimeInSec() << " seconds.");

    announceMeshProblems(surface_mesh->getNumberOfCells(),
                         static_cast< int >( collisions.size() ),
                         static_cast< int >( degenerate.size() )
                         );

    saveProblemFlagsToMesh(surface_mesh, collisions, degenerate);

    if (!writeAnnotatedMesh(surface_mesh, params.collisionsMeshName()) )
    {
      SLIC_ERROR("Couldn't write results to " << params.collisionsMeshName());
    }

    if (!writeCollisions(collisions,degenerate,params.collisionsTextName()) )
    {
      SLIC_ERROR("Couldn't write results to "<< params.collisionsTextName());
    }

    if( params.verboseOutput && !collisions.empty() )
    {
      SLIC_INFO("Intersecting triangle pairs:");
      // Initialize pairs of clashes
      for (auto i : collisions)
      {
        auto t1 = getMeshTriangle(i.first, surface_mesh);
        auto t2 = getMeshTriangle(i.second, surface_mesh);

        SLIC_INFO( "  Triangle " << i.first << " -- "
                                 << std::setprecision(17) << t1);
        SLIC_INFO( "  Triangle " << i.second << " -- "
                                 << std::setprecision(17) << t2 << "\n");
      }
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
    // _check_watertight_start
    quest::WatertightStatus wtstat =
      quest::isSurfaceMeshWatertight(surface_mesh);
    // _check_watertight_end
    timer2.stop();
    // _report_watertight_start
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
    // _report_watertight_end
    SLIC_INFO("Testing for watertightness took "
              << timer2.elapsedTimeInSec() << " seconds.");

    mint::write_vtk(surface_mesh, params.weldMeshName() );
  }

  // Delete the mesh
  delete surface_mesh;
  surface_mesh = nullptr;

  return 0;
}
