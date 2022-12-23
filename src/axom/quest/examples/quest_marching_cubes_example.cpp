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

#ifndef __WHERE
#define __STRINGIZE(x) __STRINGIZE2(x)
#define __STRINGIZE2(x) #x
//!@brief String literal for code location
#define __WHERE __FILE__ ":" __STRINGIZE(__LINE__) "(" + std::string(__func__) + ") "
#endif

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
  std::string fieldsFile{"fields"};

  double contourVal {1.0};
  std::vector<double> fcnCenter {0.0, 0.0};
  // TODO: Ensure that fcnCenter size matches dimensionality.
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

    app.add_option("-s,--fields-file", fieldsFile)
      ->description("Name of output mesh file with all its fields.")
      ->capture_default_str();

    app.add_flag("-v,--verbose,!--no-verbose", _verboseOutput)
      ->description("Enable/disable verbose output")
      ->capture_default_str();

    auto* fcn_options =
      app.add_option_group("fcn",
                           "Options for setting up the nodal distance function");

    fcn_options->add_option("--center", fcnCenter)
      ->description("Center for nodal distance function (x,y[,z])")
      ->expected(2, 3);

    app.add_option("--contourVal", contourVal)
      ->description("Contour value")
      ->capture_default_str();

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

  template<int DIM>
  axom::primal::Point<double, DIM> fcn_center() const
  {
    SLIC_ASSERT(fcnCenter.size() == DIM);
    return axom::primal::Point<double, DIM>(fcnCenter.data());
  }
};

Input params;

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
    : _coordsetPath("coordsets/" + coordset)
    , _topologyPath("topologies/" + topology)
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

  /// Get domain group.
  conduit::Node& domain(axom::IndexType domainIdx)
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  /// Get the number of cells in each direction of a blueprint single domain.
  template<int DIM>
  axom::StackArray<axom::IndexType, DIM> domain_lengths(const conduit::Node &dom) const
  {
    SLIC_ASSERT(_dimension == DIM);
    const conduit::Node& dimsNode = dom.fetch_existing(_topologyPath + "/elements/dims");
    axom::StackArray<axom::IndexType, DIM> cellCounts;
    for(int i=0; i<_dimension; ++i)
    {
      cellCounts[i] = dimsNode[i].as_int();
    }
    return cellCounts;
  }

  /// Gets the MPI rank for this mesh
  int get_mpi_rank() const { return _rank; }
  /// Gets the number of ranks in the problem
  int get_mpi_size() const { return _nranks; }

  /// Returns the number of points in a particle mesh domain
  int point_count() const
  {
    // _mdMesh.print();
    int rval = 0;
    for(const auto& dom : _mdMesh.children())
    {
      // dom.print();
      SLIC_ASSERT_MSG(
        dom.fetch_existing(_coordsetPath + "/type").as_string() == "explicit",
        axom::fmt::format("Currently only supporting explicit coordinate types."
                          "  '{}/type' is '{}'",
                          _coordsetPath, dom.fetch_existing(_coordsetPath + "/type").as_string()));

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
  void save_mesh(const std::string& filename)
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

    // Interleave coordinates if they aren't already.
    for(auto& dom : _mdMesh.children())
    {
      conduit::Node& coordsNode = dom.fetch_existing(_coordsetPath);
      conduit::Node& values = coordsNode.fetch_existing("values");
      bool isInterleaved =
        conduit::blueprint::mcarray::is_interleaved(values);
      if(!isInterleaved)
      {
        conduit::Node interleavedValues;
        conduit::blueprint::mcarray::to_interleaved(values, interleavedValues);
        values.set(interleavedValues);
      }
    }
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


template<int DIM>
struct DistanceFromPoint {
  DistanceFromPoint(const axom::primal::Point<double, DIM> &pt)
    : _center(pt)
    {}
  const axom::primal::Point<double, DIM> _center;
  void compute_nodal_distance(BlueprintStructuredMesh& bpMesh,
                              const std::string fieldName)
  {
    SLIC_ASSERT(bpMesh.dimension() == DIM);
    conduit::Node& bpNode = bpMesh.as_conduit_node();
    for(conduit::Node& dom : bpNode.children())
    {
      // Access the coordinates
      axom::StackArray<axom::IndexType, DIM> shape = bpMesh.domain_lengths<DIM>(dom);
      for(int i=0; i<DIM; ++i) { ++shape[i]; }
      conduit::index_t count = 1;
      for(auto i : shape) count *= i;
      conduit::Node& coordsValues = dom.fetch_existing("coordsets/coords/values");
      double *coordsPtrs[DIM];
      axom::ArrayView<double, DIM> coordsViews[DIM];
      for(int d=0; d<DIM; ++d)
      {
        coordsPtrs[d] = coordsValues[d].as_double_ptr();
        coordsViews[d] = axom::ArrayView<double, DIM>(coordsPtrs[d], shape);
      }
      axom::ArrayView<axom::primal::Point<double, DIM>, DIM>
        coordsView( (axom::primal::Point<double, DIM>*)coordsValues.data_ptr(),
                    shape );

      // Create the nodal function data.
      conduit::Node &fieldNode = dom["fields"][fieldName];
      fieldNode["association"] = "vertex";
      fieldNode["topology"] = "mesh";
      fieldNode["volume_dependent"] = "false";
      fieldNode["values"].set(conduit::DataType::float64(count));
      double* d = fieldNode["values"].value();
      axom::ArrayView<double, DIM> fieldView(d, shape);

      // Set the nodal data to the distance from center.
      axom::IndexType nodeCount = fieldView.size();
      for(int i=0; i<nodeCount; ++i)
      {
        axom::primal::Point<double, DIM> pt;
        for(int d=0; d<DIM; ++d)
        {
          pt[d] = coordsViews[d].flatIndex(i);
        };
        fieldView.flatIndex(i) = distance(pt);
      }
    }
  }

  double distance(axom::primal::Point<double, DIM> &pt) const
  {
    double dist = primal::squared_distance(_center, pt);
    dist = sqrt(dist);
    return dist;
  }

  /**
     Check for errors in the surface contour mesh.

     - surface triangle should be contained in their cells.
     - surface points should lie on computational mesh edges.
     - analytical scalar value at surface points should be
       contourVal, within machine zero.
  */
  int check_contour_surface(
    const axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &contourMesh,
    double contourVal,
    double tol = 1e-12)
  {
    int sumErrCount = 0;
    const axom::IndexType nodeCount = contourMesh.getNumberOfNodes();
    axom::primal::Point<double, DIM> pt;
    for(axom::IndexType i=0; i<nodeCount; ++i)
    {
      contourMesh.getNode(i, pt.data());
      double dist = distance(pt);
      double diff = std::abs(dist - contourVal);
      if(diff > tol)
      {
        ++sumErrCount;
        SLIC_INFO_IF(
          params.isVerbose(),
          axom::fmt::format("check_contour_surface: node {} has dist {}, off by {}",
                            i, dist, diff));
      }
    }
    SLIC_INFO_IF(
      params.isVerbose(),
      axom::fmt::format("check_contour_surface: found {} errors",
                        sumErrCount));
    return sumErrCount;
  }
};


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
  axom::CLI::App app {"Driver/test code for marching cubes algorithm"};

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

#if defined(AXOM_USE_UMPIRE)
  //---------------------------------------------------------------------------
  // Memory resource.  For testing, choose device memory if appropriate.
  //---------------------------------------------------------------------------
#if 0
  // Temporarily disable to prevent warnings.
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
#endif

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

  DistanceFromPoint<DIM> dfp(params.fcn_center<DIM>());

  //---------------------------------------------------------------------------
  // Initialize nodal scalar function.
  //---------------------------------------------------------------------------
  dfp.compute_nodal_distance(computationalMesh, "fVal");

  //---------------------------------------------------------------------------
  // Initialize spatial index for querying points, and run query
  //---------------------------------------------------------------------------

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
  slic::flushStreams();

  // Run the marching cubes algorithm on the computational mesh.
  // To test support for single-domain format, use single-domain when possible.
  axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> surfaceMesh(
    DIM,
    DIM == 2 ? mint::CellType::SEGMENT : mint::CellType::TRIANGLE);
  axom::Array<int> domainIds;

  mca.set_output_mesh(&surfaceMesh);
  mca.set_cell_id_field("zoneIds");
  mca.set_domain_id_field("domainIds");

  slic::flushStreams();
  computeTimer.start();
  mca.compute_iso_surface(params.contourVal);
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

  assert(computationalMesh.isValid());
  computationalMesh.save_mesh(params.fieldsFile);
  mint::write_vtk(&surfaceMesh, "surface_mesh.vtk");

  int errCount = 0;
  int localErrCount = 0;
  if(params.checkResults)
  {
    localErrCount = dfp.check_contour_surface(surfaceMesh, params.contourVal, 1e-3);
    MPI_Allreduce(&localErrCount, &errCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    if(errCount)
    {
      SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", errCount));
    }
    else
    {
      SLIC_INFO("Normal exit.");
    }
  }
  else
  {
    SLIC_INFO("Not checking results.");
  }

  // TODO: Write surface mesh to file for viz.


  finalizeLogger();
  MPI_Finalize();

  return errCount != 0;
}
