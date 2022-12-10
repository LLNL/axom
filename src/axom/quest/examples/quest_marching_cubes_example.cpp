// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file marching_cubes_example.cpp
 * \brief Driver for a marching cubes iso-surface generation
 */

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/quest.hpp"
#include "axom/slam.hpp"
#include "axom/core/Types.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_blueprint_mpi.hpp"
#include "conduit_relay_io_blueprint.hpp"
#include "conduit_relay_mpi_io_blueprint.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#ifdef AXOM_USE_MPI
#include "mpi.h"
#endif

// C/C++ includes
#include <string>
#include <limits>
#include <map>
#include <vector>
#include <cmath>

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace mint = axom::mint;
namespace slam = axom::slam;
namespace spin = axom::spin;
namespace primal = axom::primal;
namespace mint = axom::mint;
namespace numerics = axom::numerics;

using RuntimePolicy = axom::quest::DistributedClosestPoint::RuntimePolicy;

// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str)
{
  return axom::fmt::format("{:=^80}", str);
}

/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;

  double circleRadius {1.0};
  std::vector<double> circleCenter {0.0, 0.0};
  // TODO: Ensure that circleCenter size matches dimensionality.
  int circlePoints {100};
  RuntimePolicy policy {RuntimePolicy::seq};

  bool checkResults {false};

private:
  bool _verboseOutput {false};

  // clang-format off
  const std::map<std::string, RuntimePolicy> s_validPolicies
  {
      {"seq", RuntimePolicy::seq}
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #ifdef AXOM_USE_OPENMP
    , {"omp", RuntimePolicy::omp}
  #endif
  #ifdef AXOM_USE_CUDA
    , {"cuda", RuntimePolicy::cuda}
  #endif
  #ifdef AXOM_USE_HIP
    , {"hip", RuntimePolicy::hip}
  #endif
#endif
  };
  // clang-format on

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-m,--mesh-file", meshFile)
      ->description(
        "Path to multidomain computational mesh following conduit blueprint "
        "convention.")
      ->check(axom::CLI::ExistingFile);

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    app.add_option("-r,--radius", circleRadius)
      ->description("Radius for circle")
      ->capture_default_str();

    auto* circle_options =
      app.add_option_group("circle",
                           "Options for setting up the circle of points");
    circle_options->add_option("--center", circleCenter)
      ->description("Center for object (x,y[,z])")
      ->expected(2, 3);

    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for marching cubes execution")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(s_validPolicies));

    app.add_flag("-c,--check-results,!--no-check-results", checkResults)
      ->description(
        "Enable/disable checking results against analytical solution")
      ->capture_default_str();

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug
                                             : slic::message::Info);
  }
};

/**
 \brief Generic computational mesh, to hold cell and node data.
*/
struct BlueprintStructuredMesh
{
public:
  using Point2D = primal::Point<double, 2>;
  using Point3D = primal::Point<double, 3>;
  using PointArray2D = axom::Array<Point2D>;
  using PointArray3D = axom::Array<Point3D>;

  explicit BlueprintStructuredMesh(
    const std::string& meshFile,
    const std::string& coordset = "coords",
    const std::string& topology = "mesh")
    : _coordsetPath("coordset/" + coordset)
    , _topologyPath("topology/" + topology)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_nranks);
    read_blueprint_mesh(meshFile);
  }

  /// Return the blueprint mesh in a conduit::Node
  conduit::Node& as_conduit_node()
  {
    return _mdMesh;
  }

  /// Get number of domains in the multidomain particle mesh
  axom::IndexType domain_count() const { return _domCount; }

  /// Get a domain group.
  conduit::Node& domain(axom::IndexType domainIdx)
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  /// Gets the MPI rank for this mesh
  int getRank() const { return _rank; }
  /// Gets the number of ranks in the problem
  int getNumRanks() const { return _nranks; }

  /// Returns the number of points in a particle mesh domain
  int point_count() const
  {
    SLIC_ASSERT_MSG(
      _mdMesh.fetch_existing(_coordsetPath + "/type").as_string() == "explicit",
      "Currently only supporting explicit coordinate types.");

    int rval = 0;
    for(const auto& dom : _mdMesh.children())
    {
      const auto& dims = dom.fetch_existing(_topologyPath + "/elements/dims");
      int n = 1;
      for(const auto& l : dims.children())
      {
        n *= 1 + l.as_int();
      }
      rval += n;
    }
    return rval;
  }

  int dimension() const { return _dimension; }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    {
      conduit::Node info;
      if(!conduit::blueprint::mpi::verify("mesh", _mdMesh, info, MPI_COMM_WORLD))
      {
        SLIC_INFO("Invalid blueprint for particle mesh: \n" << info.to_yaml());
        slic::flushStreams();
        return false;
      }
      // info.print();
    }
    return true;
  }

  /// Outputs the particle mesh to disk
  void saveMesh(const std::string& filename)
  {
    conduit::relay::mpi::io::blueprint::save_mesh(_mdMesh,
                                                  filename,
                                                  "hdf5",
                                                  MPI_COMM_WORLD);
  }

  void print_mesh_info() const
  {
    // Copy to conduit::Node.  It's output is easier to read, especially in parallel.
    _mdMesh.print();
  }

private:
  int _rank;
  int _nranks;
  int _dimension {-1};
  conduit::Node _mdMesh;
  axom::IndexType _domCount;
  const std::string _coordsetPath;
  const std::string _topologyPath;

  /*!
    @brief Read a blueprint mesh.
  */
  void read_blueprint_mesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    _mdMesh.reset();
    conduit::relay::mpi::io::blueprint::load_mesh(meshFilename,
                                                  _mdMesh,
                                                  MPI_COMM_WORLD);
    assert(conduit::blueprint::mesh::is_multi_domain(_mdMesh));
    _domCount = conduit::blueprint::mesh::number_of_domains(_mdMesh);

    if(_domCount > 0)
    {
      const conduit::Node coordsetNode =
        _mdMesh[0].fetch_existing(_coordsetPath);
      _dimension = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }
    MPI_Allreduce(MPI_IN_PLACE, &_dimension, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    SLIC_ASSERT(_dimension > 0);

    bool valid = isValid();
    SLIC_ASSERT(valid);
  }
};  // BlueprintStructuredMesh


void make_coords_contiguous(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_contiguous(oldValues, coordValues);
  }
}

void make_coords_interleaved(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(!isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_interleaved(oldValues, coordValues);
  }
}

/// Utility function to initialize the logger
void initializeLogger()
{
  // Initialize Logger
  slic::initialize();
  slic::setLoggingMsgLevel(slic::message::Info);

  slic::LogStream* logStream;

#ifdef AXOM_USE_MPI
  std::string fmt = "[<RANK>][<LEVEL>]: <MESSAGE>\n";
  #ifdef AXOM_USE_LUMBERJACK
  const int RLIMIT = 8;
  logStream = new slic::LumberjackStream(&std::cout, MPI_COMM_WORLD, RLIMIT, fmt);
  #else
  logStream = new slic::SynchronizedStream(&std::cout, MPI_COMM_WORLD, fmt);
  #endif
#else
  std::string fmt = "[<LEVEL>]: <MESSAGE>\n";
  logStream = new slic::GenericOutputStream(&std::cout, fmt);
#endif  // AXOM_USE_MPI

  slic::addStreamToAllMsgLevels(logStream);
}

/// Utility function to finalize the logger
void finalizeLogger()
{
  if(slic::isInitialized())
  {
    slic::flushStreams();
    slic::finalize();
  }
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  int myRank, numRanks;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

  initializeLogger();
  //slic::setAbortOnWarning(true);

  //---------------------------------------------------------------------------
  // Set up and parse command line arguments
  //---------------------------------------------------------------------------
  Input params;
  axom::CLI::App app {"Driver for marching cubes code"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    int retval = -1;
    if(myRank == 0)
    {
      retval = app.exit(e);
    }

    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();

    exit(retval);
  }

  constexpr int DIM = 2;

  using PointType = primal::Point<double, DIM>;
  using Circle = primal::Sphere<double, DIM>;

#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
  const std::string umpireResourceName =
    params.policy == RuntimePolicy::seq || params.policy == RuntimePolicy::omp
    ? "HOST"
    :
  #if defined(UMPIRE_ENABLE_DEVICE)
    "DEVICE"
  #elif defined(UMPIRE_ENABLE_UM)
    "UM"
  #elif defined(UMPIRE_ENABLE_PINNED)
    "PINNED"
  #else
    "HOST"
  #endif
    ;
  auto& rm = umpire::ResourceManager::getInstance();
  umpire::Allocator umpireAllocator = rm.getAllocator(umpireResourceName);
#endif

  //---------------------------------------------------------------------------
  // Load/generate object mesh
  //---------------------------------------------------------------------------
  const Circle circle(
    PointType(params.circleCenter.data(), params.circleCenter.size()),
    params.circleRadius);

  //---------------------------------------------------------------------------
  // Load computational mesh.
  //---------------------------------------------------------------------------
  BlueprintStructuredMesh computationalMesh(params.meshFile);
  // computationalMesh.print_mesh_info();

  SLIC_INFO_IF(
    params.isVerbose(),
    axom::fmt::format("Query mesh has {} points in {} domains locally",
                      computationalMesh.point_count(),
                      computationalMesh.domain_count()));
  slic::flushStreams();

  auto getIntMinMax = [](int inVal, int& minVal, int& maxVal, int& sumVal) {
    MPI_Allreduce(&inVal, &minVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &maxVal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &sumVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  };

  // Output some global mesh size stats
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.point_count(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Query mesh has {{min:{}, max:{}, sum:{}, avg:{}}} points",
      mn,
      mx,
      sum,
      (double)sum / numRanks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.domain_count(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Query mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
      mn,
      mx,
      sum,
      (double)sum / numRanks));
  }

  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Initialize spatial index for querying points, and run query
  //---------------------------------------------------------------------------

  auto init_str =
    banner(axom::fmt::format("Initializing BVH tree over {} points",
                             params.circlePoints));

  axom::utilities::Timer computeTimer(false);

  // To test with contiguous and interleaved coordinate storage,
  // make half them contiguous.
  for(int di = 0; di < computationalMesh.domain_count(); ++di)
  {
    auto& dom = computationalMesh.as_conduit_node().child(di);
    if((myRank + di) % 2 == 1)
    {
      make_coords_contiguous(dom.fetch_existing("coordsets/coords/values"));
    }
  }

  // Create marching cubes algorithm object and set some parameters
  quest::MarchingCubesAlgo mca(computationalMesh.as_conduit_node(), "fVal");

  // Build the spatial index over the object on each rank
  SLIC_INFO(init_str);
  slic::flushStreams();

  // Run the marching cubes algorithm on the computational mesh.
  // To test support for single-domain format, use single-domain when possible.
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> surfaceMesh(
    DIM,
    DIM == 2 ? mint::CellType::SEGMENT : mint::CellType::TRIANGLE);

  slic::flushStreams();
  computeTimer.start();
  mca.compute_iso_surface(0.0);
  computeTimer.stop();

  auto getDoubleMinMax =
    [](double inVal, double& minVal, double& maxVal, double& sumVal) {
      MPI_Allreduce(&inVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &sumVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    };

  // Output some timing stats
  {
    double minQuery, maxQuery, sumQuery;
    getDoubleMinMax(computeTimer.elapsedTimeInSec(), minQuery, maxQuery, sumQuery);

    SLIC_INFO(axom::fmt::format(
      "Query with policy {} took {{avg:{}, min:{}, max:{}}} seconds",
      params.policy,
      sumQuery / numRanks,
      minQuery,
      maxQuery));
  }
  slic::flushStreams();

  int errCount = 0;
  int localErrCount = 0;
  if(params.checkResults)
  {
    // TODO: Write a correctness checker.
  }
  MPI_Allreduce(&localErrCount, &errCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  // TODO: Write surface mesh to file for viz.

  if(errCount)
  {
    SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", errCount));
  }
  else
  {
    SLIC_INFO("Normal exit.");
  }

  finalizeLogger();
  MPI_Finalize();

  return errCount != 0;
}
