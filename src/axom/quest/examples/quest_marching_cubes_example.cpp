// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 \file marching_cubes_example.cpp
 \brief Driver and test for a marching cubes isocontour generation

  The test can generate planar and round contours.  Planar contours
  can be checked to machine-zero accuracy, but it doesn't test a great
  variety of contour-mesh intersection types.  Round contours can
  check more intersection types but requires a tolerance to allow
  for the function not varying linearly along mesh lines.
*/

// Axom includes
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/mint/execution/internal/structured_exec.hpp"
#include "axom/quest/ArrayIndexer.hpp"
#include "axom/quest/MarchingCubes.hpp"
#include "axom/quest/MeshViewUtil.hpp"
#include "axom/sidre.hpp"
#include "axom/core/Types.hpp"

#include "conduit_blueprint.hpp"
#include "conduit_relay_io_blueprint.hpp"
#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
  #include "conduit_relay_mpi_io_blueprint.hpp"
#endif

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
namespace sidre = axom::sidre;
namespace primal = axom::primal;
namespace mint = axom::mint;
namespace numerics = axom::numerics;

///////////////////////////////////////////////////////////////
// converts the input string into an 80 character string
// padded on both sides with '=' symbols
std::string banner(const std::string& str)
{
  return axom::fmt::format("{:=^80}", str);
}

///////////////////////////////////////////////////////////////
/// Struct to parse and store the input parameters
struct Input
{
public:
  std::string meshFile;
  std::string fieldsFile {"fields"};

  // Center of round contour function
  bool usingRound {false};
  std::vector<double> fcnCenter;

  // Parameters for planar contour function
  bool usingPlanar {false};
  std::vector<double> inPlane;
  std::vector<double> perpDir;

  std::size_t ndim {0};

  double contourVal {1.0};

  bool checkResults {false};

  quest::MarchingCubesRuntimePolicy policy {
    quest::MarchingCubesRuntimePolicy::seq};

private:
  bool _verboseOutput {false};

  // clang-format off
  const std::map<std::string, quest::MarchingCubesRuntimePolicy> s_validPolicies
  {
      {"seq", quest::MarchingCubesRuntimePolicy::seq}
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
    , {"omp", quest::MarchingCubesRuntimePolicy::omp}
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
    , {"cuda", quest::MarchingCubesRuntimePolicy::cuda}
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
    , {"hip", quest::MarchingCubesRuntimePolicy::hip}
  #endif
#endif
  };
  // clang-format on

public:
  bool isVerbose() const { return _verboseOutput; }

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    app.add_option("-p, --policy", policy)
      ->description("Set runtime policy for point query method")
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(s_validPolicies));

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

    auto* distFromPtOption = app.add_option_group(
      "distFromPtOption",
      "Options for setting up distance-from-point function");
    distFromPtOption->add_option("--center", fcnCenter)
      ->description("Center for distance-from-point function (x,y[,z])")
      ->expected(2, 3);

    auto* distFromPlaneOption = app.add_option_group(
      "distFromPlaneOption",
      "Options for setting up distance-from-plane function");
    distFromPlaneOption->add_option("--inPlane", inPlane)
      ->description("In-plane point for distance-from-plane function (x,y[,z])")
      ->expected(2, 3);
    distFromPlaneOption->add_option("--dir", perpDir)
      ->description(
        "Positive direction for distance-from-plane function (x,y[,z])")
      ->expected(2, 3);

    app.add_option("--contourVal", contourVal)
      ->description("Contour value")
      ->capture_default_str();

    app.add_flag("-c,--check-results,!--no-check-results", checkResults)
      ->description(
        "Enable/disable checking results against analytical solution")
      ->capture_default_str();

    app.get_formatter()->column_width(60);

    // could throw an exception
    app.parse(argc, argv);

    slic::setLoggingMsgLevel(_verboseOutput ? slic::message::Debug
                                            : slic::message::Info);

    ndim = std::max(ndim, fcnCenter.size());
    ndim = std::max(ndim, inPlane.size());
    ndim = std::max(ndim, perpDir.size());
    SLIC_ASSERT_MSG((fcnCenter.empty() || fcnCenter.size() == ndim) &&
                      (inPlane.empty() || inPlane.size() == ndim) &&
                      (perpDir.empty() || perpDir.size() == ndim),
                    "fcnCenter, inPlane and perpDir must have consistent sizes "
                    "if specified.");

    usingPlanar = !perpDir.empty();
    usingRound = !fcnCenter.empty();
    SLIC_ASSERT_MSG(
      usingPlanar || usingRound,
      "You must specify a planar scalar function or a round scalar"
      " function or both.");

    // inPlane defaults to origin if omitted.
    if(usingPlanar && inPlane.empty())
    {
      inPlane.insert(inPlane.begin(), ndim, 0.0);
    }
  }

  template <int DIM>
  axom::primal::Point<double, DIM> roundContourCenter() const
  {
    SLIC_ASSERT(fcnCenter.size() == DIM);
    return axom::primal::Point<double, DIM>(fcnCenter.data());
  }

  template <int DIM>
  axom::primal::Point<double, DIM> inplanePoint() const
  {
    SLIC_ASSERT(inPlane.size() == DIM);
    return axom::primal::Point<double, DIM>(inPlane.data());
  }

  template <int DIM>
  axom::primal::Vector<double, DIM> planeNormal() const
  {
    SLIC_ASSERT(perpDir.size() == DIM);
    return axom::primal::Vector<double, DIM>(perpDir.data());
  }
};

//!@brief Our allocator id, based on execution policy.
static int s_allocatorId = axom::INVALID_ALLOCATOR_ID;  // Set in main.
static int allocatorIdForPolicy(quest::MarchingCubesRuntimePolicy policy)
{
  //---------------------------------------------------------------------------
  // Set default allocator for possibly testing on devices
  //---------------------------------------------------------------------------
  int aid = axom::INVALID_ALLOCATOR_ID;

  // clang-format off
  if(policy == axom::quest::MarchingCubesRuntimePolicy::seq)
  {
    aid = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  }
#if defined(AXOM_USE_RAJA)
#ifdef _AXOM_MC_USE_OPENMP
  else if(policy == axom::quest::MarchingCubesRuntimePolicy::omp)
  {
    aid = axom::execution_space<axom::OMP_EXEC>::allocatorID();
  }
#endif
#ifdef _AXOM_MC_USE_CUDA
  else if(policy == axom::quest::MarchingCubesRuntimePolicy::cuda)
  {
    // aid = axom::execution_space<axom::CUDA_EXEC<256>>::allocatorID();
    aid = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  }
#endif
#ifdef _AXOM_MC_USE_HIP
  else if(policy == axom::quest::MarchingCubesRuntimePolicy::hip)
  {
    // aid = axom::execution_space<axom::HIP_EXEC<256>>::allocatorID();
    aid = axom::getUmpireResourceAllocatorID(umpire::resource::Device);
  }
#endif
#endif
  // clang-format on

  SLIC_ERROR_IF(
    aid == axom::INVALID_ALLOCATOR_ID,
    axom::fmt::format("Cannot find allocator id for policy '{}'", policy));
  return aid;
}

//!@brief Put a conduit::Node array data into the specified memory space.
template <typename T>
void moveConduitDataToNewMemorySpace(conduit::Node& node,
                                     const std::string& path,
                                     int allocId)
{
  conduit::Node& dataNode = node.fetch_existing(path);
  SLIC_ASSERT(!dataNode.dtype().is_empty() && !dataNode.dtype().is_object() &&
              !dataNode.dtype().is_list());

  std::size_t count = dataNode.dtype().number_of_elements();
  T* oldPtr = static_cast<T*>(dataNode.data_ptr());
  bool deleteOld = dataNode.is_data_external();
  T* newPtr = axom::allocate<T>(count, allocId);
  axom::copy(newPtr, oldPtr, count * sizeof(T));
  dataNode.set_external(newPtr, count);
  if(deleteOld)
  {
    axom::deallocate(oldPtr);
  }
}

Input params;

int myRank = -1, numRanks = -1;  // MPI stuff, set in main().

/**
 \brief Generic computational mesh, to hold cell and node data.
*/
struct BlueprintStructuredMesh
{
public:
  explicit BlueprintStructuredMesh(const std::string& meshFile,
                                   const std::string& topologyName)
    : _topologyPath("topologies/" + topologyName)
  {
    readBlueprintMesh(meshFile);
    for(int d = 0; d < _mdMesh.number_of_children(); ++d)
    {
      auto dl = domainLengths(d);
      SLIC_INFO(axom::fmt::format("dom[{}] size={}", d, dl));
    }
    _maxSpacing = maxSpacing();
  }

  /// Return the blueprint mesh in a conduit::Node
  conduit::Node& asConduitNode() { return _mdMesh; }

  /// Get number of domains in the multidomain mesh
  axom::IndexType domainCount() const { return _domCount; }

  //// Whether mesh is empty locally.
  bool empty() const { return _domCount == 0; }

  /// Get domain group.
  conduit::Node& domain(axom::IndexType domainIdx)
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  const conduit::Node& domain(axom::IndexType domainIdx) const
  {
    SLIC_ASSERT(domainIdx >= 0 && domainIdx < _domCount);
    return _mdMesh.child(domainIdx);
  }

  /*!
    @brief Get the number of cells in each direction of a blueprint single domain.

    @param domId Index of domain
    @lengths Space for dimension() numbers.
  */
  void domainLengths(axom::IndexType domId, axom::IndexType* lengths) const
  {
    const conduit::Node& dom = domain(domId);
    SLIC_ASSERT_MSG(
      dom.fetch_existing(_coordsetPath + "/type").as_string() == "explicit",
      axom::fmt::format("Currently only supporting explicit coordinate types."
                        "  '{}/type' is '{}'",
                        _coordsetPath,
                        dom.fetch_existing(_coordsetPath + "/type").as_string()));
    const conduit::Node& dimsNode =
      dom.fetch_existing(_topologyPath + "/elements/dims");
    for(int i = 0; i < _ndims; ++i)
    {
      lengths[i] = dimsNode[i].as_int();
    }
  }

  axom::Array<axom::IndexType> domainLengths(axom::IndexType domainId) const
  {
    axom::Array<axom::IndexType> rval(_ndims, _ndims);
    domainLengths(domainId, rval.data());
    return rval;
  }

  /// Returns the number of cells in a domain
  int cellCount(axom::IndexType domId) const
  {
    auto shape = domainLengths(domId);
    int rval = 1;
    for(const auto& l : shape)
    {
      rval *= l;
    }
    return rval;
  }

  /// Returns the number of cells in all mesh domains
  int cellCount() const
  {
    int rval = 0;
    for(int domId = 0; domId < _mdMesh.number_of_children(); ++domId)
    {
      rval += cellCount(domId);
    }
    return rval;
  }

  /// Returns the number of nodes in a domain
  int nodeCount(axom::IndexType domId) const
  {
    auto shape = domainLengths(domId);
    int rval = 1;
    for(const auto& l : shape)
    {
      rval *= 1 + l;
    }
    return rval;
  }

  /// Returns the number of nodes in all mesh domains
  int nodeCount() const
  {
    int rval = 0;
    for(int domId = 0; domId < _mdMesh.number_of_children(); ++domId)
    {
      rval += nodeCount(domId);
    }
    return rval;
  }

  int dimension() const { return _ndims; }

  const std::string& coordsetPath() const { return _coordsetPath; }

  /*!
    @return largest mesh spacing.

    Compute only once, because after that, coordinates data may be
    moved to devices.
  */
  double maxSpacing() const
  {
    if(_maxSpacing >= 0)
    {
      // Max spacing has been computed and cached.
      return _maxSpacing;
    }

    double localRval = 0.0;
    for(axom::IndexType domId = 0; domId < domainCount(); ++domId)
    {
      localRval = std::max(localRval, maxSpacing1(domId));
    }

    double rval = localRval;
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&localRval, &rval, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

    return rval;
  }

  /*!
    @return largest mesh spacing in a domain.

    This method takes shortcuts by assuming
    the mesh is structured and cartesian, with explicit coordinates.
  */
  double maxSpacing1(axom::IndexType domId) const
  {
    const conduit::Node& dom = domain(domId);
    const conduit::Node& dimsNode =
      dom.fetch_existing("topologies/mesh/elements/dims");
    axom::Array<axom::IndexType> ls(_ndims);
    for(int d = 0; d < _ndims; ++d)
    {
      ls[d] = 1 + dimsNode[d].as_int();
    }

    double rval = 0.0;

    const conduit::Node& cVals = dom.fetch_existing("coordsets/coords/values");
    if(_ndims == 2)
    {
      axom::ArrayView<const double, 2> xs(cVals["x"].as_double_ptr(),
                                          ls[1],
                                          ls[0]);
      axom::ArrayView<const double, 2> ys(cVals["y"].as_double_ptr(),
                                          ls[1],
                                          ls[0]);
      rval = std::max(rval, std::abs(xs(0, 0) - xs(0, 1)));
      rval = std::max(rval, std::abs(ys(0, 0) - ys(1, 0)));
    }
    else
    {
      axom::ArrayView<const double, 3> xs(cVals["x"].as_double_ptr(),
                                          ls[2],
                                          ls[1],
                                          ls[0]);
      axom::ArrayView<const double, 3> ys(cVals["y"].as_double_ptr(),
                                          ls[2],
                                          ls[1],
                                          ls[0]);
      axom::ArrayView<const double, 3> zs(cVals["z"].as_double_ptr(),
                                          ls[2],
                                          ls[1],
                                          ls[0]);
      rval = std::max(rval, std::abs(xs(0, 0, 0) - xs(0, 0, 1)));
      rval = std::max(rval, std::abs(ys(0, 0, 0) - ys(0, 1, 0)));
      rval = std::max(rval, std::abs(zs(0, 0, 0) - zs(1, 0, 0)));
    }
    return rval;
  }

  /// Checks whether the blueprint is valid and prints diagnostics
  bool isValid() const
  {
    conduit::Node info;
#ifdef AXOM_USE_MPI
    if(!conduit::blueprint::mpi::verify("mesh", _mdMesh, info, MPI_COMM_WORLD))
#else
    if(!conduit::blueprint::verify("mesh", _mdMesh, info))
#endif
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      return false;
    }
    return true;
  }

  void printMeshInfo() const { _mdMesh.print(); }

  /*!
    @param[in] path Path to existing data in the blueprint mesh,
    relative to each domain in the mesh.
    @param[in] allocId Allocator id for the new memory space.
    @tparam Type of data being moved.  Should be something Conduit
    supports, i.e., not custom user data.
  */
  template <typename T>
  void moveMeshDataToNewMemorySpace(const std::string& path, int allocId)
  {
    for(auto& dom : _mdMesh.children())
    {
      moveConduitDataToNewMemorySpace<T>(dom, path, allocId);
    }
  }

private:
  int _ndims {-1};
  conduit::Node _mdMesh;
  axom::IndexType _domCount;
  bool _coordsAreStrided = false;
  const std::string _topologyPath;
  std::string _coordsetPath;
  double _maxSpacing = -1.0;

  /*!
    @brief Read a blueprint mesh into conduit::Node _mdMesh.
  */
  void readBlueprintMesh(const std::string& meshFilename)
  {
    SLIC_ASSERT(!meshFilename.empty());

    _mdMesh.reset();
#ifdef AXOM_USE_MPI
    conduit::relay::mpi::io::blueprint::load_mesh(meshFilename,
                                                  _mdMesh,
                                                  MPI_COMM_WORLD);
#else
    conduit::relay::io::blueprint::load_mesh(meshFilename, _mdMesh);
#endif
    SLIC_ASSERT(conduit::blueprint::mesh::is_multi_domain(_mdMesh));
    _domCount = conduit::blueprint::mesh::number_of_domains(_mdMesh);

    if(_domCount > 0)
    {
      SLIC_ASSERT(_mdMesh[0].has_path(_topologyPath));
      auto coordsetName =
        _mdMesh[0].fetch_existing(_topologyPath + "/coordset").as_string();
      _coordsetPath = axom::fmt::format("coordsets/{}/", coordsetName);
      SLIC_ASSERT(_mdMesh[0].has_path(_coordsetPath));

      _coordsAreStrided = _mdMesh[0]
                            .fetch_existing(_topologyPath + "/elements/dims")
                            .has_child("strides");
      if(_coordsAreStrided)
      {
        SLIC_WARNING(axom::fmt::format(
          "Mesh '{}' is strided.  Stride support is under development.",
          meshFilename));
      }
      const conduit::Node coordsetNode = _mdMesh[0].fetch_existing(_coordsetPath);
      _ndims = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }
#ifdef AXOM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &_ndims, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    SLIC_ASSERT(_ndims > 0);

    SLIC_ASSERT(isValid());
  }
};  // BlueprintStructuredMesh

/// Output some timing stats
void printTimingStats(axom::utilities::Timer& t, const std::string& description)
{
  auto getDoubleMinMax =
    [](double inVal, double& minVal, double& maxVal, double& sumVal) {
#ifdef AXOM_USE_MPI
      MPI_Allreduce(&inVal, &minVal, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &maxVal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&inVal, &sumVal, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      minVal = inVal;
      maxVal = inVal;
      sumVal = inVal;
#endif
    };

  {
    double minCompute, maxCompute, sumCompute;
    getDoubleMinMax(t.elapsedTimeInSec(), minCompute, maxCompute, sumCompute);

    SLIC_INFO(axom::fmt::format("'{}' took {{avg:{}, min:{}, max:{}}} seconds",
                                description,
                                sumCompute / numRanks,
                                minCompute,
                                maxCompute));
  }
}

/// Write blueprint mesh to disk
void saveMesh(const conduit::Node& mesh, const std::string& filename)
{
#ifdef AXOM_USE_MPI
  conduit::relay::mpi::io::blueprint::save_mesh(mesh,
                                                filename,
                                                "hdf5",
                                                MPI_COMM_WORLD);
#else
  conduit::relay::io::blueprint::save_mesh(mesh, filename, "hdf5");
#endif
}

/// Write blueprint mesh to disk
void saveMesh(const sidre::Group& mesh, const std::string& filename)
{
  conduit::Node tmpMesh;
  mesh.createNativeLayout(tmpMesh);
  {
    conduit::Node info;
#ifdef AXOM_USE_MPI
    if(!conduit::blueprint::mpi::verify("mesh", tmpMesh, info, MPI_COMM_WORLD))
#else
    if(!conduit::blueprint::verify("mesh", tmpMesh, info))
#endif
    {
      SLIC_INFO("Invalid blueprint for mesh: \n" << info.to_yaml());
      slic::flushStreams();
      assert(false);
    }
    // info.print();
  }
  saveMesh(tmpMesh, filename);
}

template <typename T, int DIM>
T product(const axom::StackArray<T, DIM>& a)
{
  T rval = a[0];
  for(int d = 1; d < DIM; ++d)
  {
    rval *= a[d];
  }
  return rval;
}

template <typename T, int DIM, typename U>
static void addToStackArray(axom::StackArray<T, DIM>& a, U b)
{
  for(int d = 0; d < DIM; ++d)
  {
    a[d] += b;
  }
}

template <int DIM, typename ExecSpace>
struct ContourTestBase
{
  static constexpr auto MemorySpace =
    axom::execution_space<ExecSpace>::memory_space;
  using PointType = axom::primal::Point<double, DIM>;
  ContourTestBase()
    : m_parentCellIdField("parentCellIds")
    , m_domainIdField("domainIdField")
  { }
  virtual ~ContourTestBase() { }

  //!@brief Return field name for storing nodal function.
  virtual std::string name() const = 0;

  //!@brief Return field name for storing nodal function.
  virtual std::string functionName() const = 0;

  //!@brief Return function value at a point.
  virtual AXOM_HOST_DEVICE double value(const PointType& pt) const = 0;

  //!@brief Return error tolerance for contour mesh accuracy check.
  virtual double errorTolerance() const = 0;

  const std::string m_parentCellIdField;
  const std::string m_domainIdField;

  int runTest(BlueprintStructuredMesh& computationalMesh, quest::MarchingCubes& mc)
  {
    SLIC_INFO(banner(axom::fmt::format("Testing {} contour.", name())));

    // Conduit data is in host memory, move to devices for testing.
    if(s_allocatorId != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
    {
      const std::string axes[3] = {"x", "y", "z"};
      for(int d = 0; d < DIM; ++d)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          computationalMesh.coordsetPath() + "/values/" + axes[d],
          s_allocatorId);
      }
      computationalMesh.moveMeshDataToNewMemorySpace<double>(
        axom::fmt::format("fields/{}/values", functionName()),
        s_allocatorId);
    }

#if defined(AXOM_USE_UMPIRE)
    /*
      Make sure data is correctly on host or device.
      We don't test with Unified memory because it's too forgiving.
    */
    if(!computationalMesh.empty())
    {
      std::string resourceName = "HOST";
      umpire::ResourceManager& rm = umpire::ResourceManager::getInstance();
      const std::string dataPath =
        axom::fmt::format("fields/{}/values", functionName());
      void* dataPtr =
        computationalMesh.domain(0).fetch_existing(dataPath).data_ptr();
      bool dataFromUmpire = rm.hasAllocator(dataPtr);
      if(dataFromUmpire)
      {
        umpire::Allocator allocator = rm.getAllocator(dataPtr);
        resourceName = allocator.getName();
      }
      SLIC_INFO(
        axom::fmt::format("Testing with policy {} and function data on {}",
                          params.policy,
                          resourceName));
      if(params.policy == quest::MarchingCubesRuntimePolicy::seq ||
         params.policy == quest::MarchingCubesRuntimePolicy::omp)
      {
        SLIC_ASSERT(resourceName == "HOST");
      }
      else
      {
        SLIC_ASSERT(resourceName == "DEVICE");
      }
    }
#endif

    mc.setFunctionField(functionName());

    axom::utilities::Timer computeTimer(false);
#ifdef AXOM_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    computeTimer.start();
    mc.computeIsocontour(params.contourVal);
    computeTimer.stop();
    printTimingStats(computeTimer, name() + " contour");

    SLIC_INFO(axom::fmt::format("Surface mesh has locally {} cells, {} nodes.",
                                mc.getContourCellCount(),
                                mc.getContourNodeCount()));

    // Return conduit data to host memory.
    if(s_allocatorId != axom::execution_space<axom::SEQ_EXEC>::allocatorID())
    {
      const std::string axes[3] = {"x", "y", "z"};
      for(int d = 0; d < DIM; ++d)
      {
        computationalMesh.moveMeshDataToNewMemorySpace<double>(
          computationalMesh.coordsetPath() + "/values/" + axes[d],
          axom::execution_space<axom::SEQ_EXEC>::allocatorID());
      }
      computationalMesh.moveMeshDataToNewMemorySpace<double>(
        axom::fmt::format("fields/{}/values", functionName()),
        axom::execution_space<axom::SEQ_EXEC>::allocatorID());
    }

    // Put mesh mesh in a mint object for error checking and output.
    std::string sidreGroupName = name() + "_mesh";
    sidre::DataStore objectDS;
    sidre::Group* meshGroup = objectDS.getRoot()->createGroup(sidreGroupName);
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> contourMesh(
      DIM,
      DIM == 2 ? mint::CellType::SEGMENT : mint::CellType::TRIANGLE,
      meshGroup,
      mc.getContourNodeCount(),
      mc.getContourCellCount());
    mc.populateContourMesh(contourMesh, m_parentCellIdField, m_domainIdField);

    int localErrCount = 0;
    if(params.checkResults)
    {
      localErrCount +=
        checkContourSurface(contourMesh, params.contourVal, "diff");

      localErrCount += checkContourCellLimits(computationalMesh, contourMesh);

      localErrCount +=
        checkCellsContainingContour(computationalMesh, contourMesh);
    }

    // Write contour mesh to file.
    addRankOffsetToContourMeshDomainIds(computationalMesh.domainCount(),
                                        contourMesh);
    std::string outputName = name() + "_contour_mesh";
    saveMesh(*meshGroup, outputName);
    SLIC_INFO(axom::fmt::format("Wrote {} contour in {}", name(), outputName));

    objectDS.getRoot()->destroyGroupAndData(sidreGroupName);

    return localErrCount;
  }

  void computeNodalDistance(BlueprintStructuredMesh& bpMesh)
  {
    SLIC_ASSERT(bpMesh.dimension() == DIM);
    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      conduit::Node& dom = bpMesh.domain(domId);
      axom::quest::MeshViewUtil<DIM, axom::execution_space<ExecSpace>::memory_space>
        mvu(dom, "mesh");

      // Create nodal function data with ghosts like node coords.
      mvu.createField(functionName(),
                      "vertex",
                      conduit::DataType::float64(mvu.getCoordsCountWithGhosts()),
                      mvu.getCoordsStrides(),
                      mvu.getCoordsOffsets());
    }

    for(int domId = 0; domId < bpMesh.domainCount(); ++domId)
    {
      conduit::Node& dom = bpMesh.domain(domId);
      axom::quest::MeshViewUtil<DIM, axom::execution_space<ExecSpace>::memory_space>
        mvu(dom, "mesh");
      const auto coordsViews = mvu.getConstCoordsViews(false);
      axom::ArrayView<double, DIM, MemorySpace> fieldView =
        mvu.template getFieldView<double>(functionName(), false);
      const auto& fieldShape = fieldView.shape();
      for(int d = 0; d < DIM; ++d)
      {
        SLIC_ASSERT(coordsViews[d].shape() == fieldShape);
      }
      populateNodalDistance(coordsViews, fieldView);
    }
  }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2>::type populateNodalDistance(
    const axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>&
      coordsViews,
    axom::ArrayView<double, DIM, MemorySpace>& fieldView)
  {
    const auto& fieldShape = fieldView.shape();
    for(int d = 0; d < DIM; ++d)
    {
      SLIC_ASSERT(coordsViews[d].shape() == fieldShape);
    }

#if defined(AXOM_USE_RAJA)
    RAJA::RangeSegment iRange(0, fieldShape[0]);
    RAJA::RangeSegment jRange(0, fieldShape[1]);
    using EXEC_POL =
      typename axom::mint::internal::structured_exec<axom::SEQ_EXEC>::loop2d_policy;
    RAJA::kernel<EXEC_POL>(
      RAJA::make_tuple(iRange, jRange),
      AXOM_LAMBDA(axom::IndexType i, axom::IndexType j) {
        PointType pt;
        for(int d = 0; d < DIM; ++d)
        {
          pt[d] = coordsViews[d](i, j);
        }
        // auto v = value(pt);
        fieldView(i, j) = 0.0;  // value(pt);
      });
#else
    for(axom::IndexType j = 0; j < fieldShape[1]; ++j)
    {
      for(axom::IndexType i = 0; i < fieldShape[0]; ++i)
      {
        PointType pt;
        for(int d = 0; d < DIM; ++d)
        {
          pt[d] = coordsViews[d](i, j);
        }
        fieldView(i, j) = value(pt);
      }
    }
#endif
  }
  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3>::type populateNodalDistance(
    const axom::StackArray<axom::ArrayView<const double, DIM, MemorySpace>, DIM>&
      coordsViews,
    axom::ArrayView<double, DIM, MemorySpace>& fieldView)
  {
    const auto& fieldShape = fieldView.shape();
    for(int d = 0; d < DIM; ++d)
    {
      SLIC_ASSERT(coordsViews[d].shape() == fieldShape);
    }

    for(axom::IndexType k = 0; k < fieldShape[2]; ++k)
    {
      for(axom::IndexType j = 0; j < fieldShape[1]; ++j)
      {
        for(axom::IndexType i = 0; i < fieldShape[0]; ++i)
        {
          PointType pt;
          for(int d = 0; d < DIM; ++d)
          {
            pt[d] = coordsViews[d](i, j, k);
          }
          fieldView(i, j, k) = value(pt);
        }
      }
    }
  }

  /*
    TODO: Additional tests:
     - surface points should lie on computational mesh edges.
     - computational cells stradling contourVal should be
       parents to at least 1 contour mesh cell.
  */
  /**
     Check for errors in the surface contour mesh.
     - analytical scalar value at surface points should be
       contourVal, within tolerance zero.
  */
  int checkContourSurface(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh,
    double contourVal,
    const std::string& diffField = {})
  {
    double* diffPtr = nullptr;
    if(!diffField.empty())
    {
      diffPtr =
        contourMesh.createField<double>(diffField, axom::mint::NODE_CENTERED);
    }

    double tol = errorTolerance();
    int errCount = 0;
    const axom::IndexType nodeCount = contourMesh.getNumberOfNodes();
    PointType pt;
    for(axom::IndexType i = 0; i < nodeCount; ++i)
    {
      contourMesh.getNode(i, pt.data());
      double analyticalVal = value(pt);
      double diff = std::abs(analyticalVal - contourVal);
      if(diffPtr)
      {
        diffPtr[i] = diff;
      }
      if(diff > tol)
      {
        ++errCount;
        SLIC_INFO_IF(
          params.isVerbose(),
          axom::fmt::format(
            "checkContourSurface: node {} at {} has dist {}, off by {}",
            i,
            pt,
            analyticalVal,
            diff));
      }
    }
    SLIC_INFO_IF(
      params.isVerbose(),
      axom::fmt::format(
        "checkContourSurface: found {} errors outside tolerance of {}",
        errCount,
        tol));
    return errCount;
  }

  //!@brief Get view of output domain id data.
  axom::ArrayView<const axom::IndexType> getDomainIdView(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh) const
  {
    const auto* ptr =
      contourMesh.getFieldPtr<axom::IndexType>(m_domainIdField,
                                               axom::mint::CELL_CENTERED);
    axom::ArrayView<const axom::IndexType> view(ptr,
                                                contourMesh.getNumberOfCells());
    return view;
  }

  //!@brief Get view of output parent cell idx data.
  axom::ArrayView<const axom::StackArray<axom::IndexType, DIM>> get_parent_cell_idx_view(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh) const
  {
    axom::IndexType numIdxComponents = -1;
    axom::IndexType* ptr =
      contourMesh.getFieldPtr<axom::IndexType>(m_parentCellIdField,
                                               axom::mint::CELL_CENTERED,
                                               numIdxComponents);

    SLIC_ASSERT(numIdxComponents == DIM);

    axom::ArrayView<const axom::StackArray<axom::IndexType, DIM>> view(
      (axom::StackArray<axom::IndexType, DIM>*)ptr,
      contourMesh.getNumberOfCells());
    return view;
  }

  /**
     Check that generated cells fall within their parents.
  */
  int checkContourCellLimits(
    BlueprintStructuredMesh& computationalMesh,
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh)
  {
    int errCount = 0;
    const axom::IndexType cellCount = contourMesh.getNumberOfCells();

    auto parentCellIdxView = get_parent_cell_idx_view(contourMesh);
    auto domainIdView = getDomainIdView(contourMesh);

    const axom::IndexType domainCount = computationalMesh.domainCount();
    axom::Array<typename axom::quest::MeshViewUtil<DIM, MemorySpace>::ConstCoordsViewsType>
      allCoordsViews(domainCount);
    for(int n = 0; n < domainCount; ++n)
    {
      const auto& dom = computationalMesh.domain(n);
      axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(dom, "mesh");
      allCoordsViews[n] = mvu.getConstCoordsViews(false);
    }

    for(axom::IndexType contourCellNum = 0; contourCellNum < cellCount;
        ++contourCellNum)
    {
      axom::IndexType domainId = domainIdView[contourCellNum];
      typename axom::quest::MeshViewUtil<DIM, MemorySpace>::ConstCoordsViewsType&
        coordsViews = allCoordsViews[domainId];

      axom::StackArray<axom::IndexType, DIM> parentCellIdx =
        parentCellIdxView[contourCellNum];
      axom::StackArray<axom::IndexType, DIM> upperIdx = parentCellIdx;
      addToStackArray(upperIdx, 1);

      PointType lower, upper;
      for(int d = 0; d < DIM; ++d)
      {
        lower[d] = coordsViews[d][parentCellIdx];
        upper[d] = coordsViews[d][upperIdx];
      }
      axom::primal::BoundingBox<double, DIM> parentCellBox(lower, upper);
      double tol = errorTolerance();
      axom::primal::BoundingBox<double, DIM> big(parentCellBox);
      axom::primal::BoundingBox<double, DIM> small(parentCellBox);
      big.expand(tol);
      small.expand(-tol);

      axom::IndexType* cellNodeIds = contourMesh.getCellNodeIDs(contourCellNum);
      const axom::IndexType cellNodeCount =
        contourMesh.getNumberOfCellNodes(contourCellNum);

      for(axom::IndexType nn = 0; nn < cellNodeCount; ++nn)
      {
        PointType nodeCoords;
        contourMesh.getNode(cellNodeIds[nn], nodeCoords.data());

        if(!big.contains(nodeCoords) || small.contains(nodeCoords))
        {
          ++errCount;
          SLIC_INFO_IF(
            params.isVerbose(),
            axom::fmt::format("checkContourCellLimits: node {} at {} is not "
                              "on parent cell boundary.",
                              cellNodeIds[nn],
                              nodeCoords));
        }
      }
    }

    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("checkContourCellLimits: found {} nodes "
                                   "not on parent cell boundary.",
                                   errCount));
    return errCount;
  }

  /*!
    Check that computational cells that contain the contour value
    have at least one contour mesh cell.
  */
  int checkCellsContainingContour(
    BlueprintStructuredMesh& computationalMesh,
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh)
  {
    int errCount = 0;
    const axom::IndexType cellCount = contourMesh.getNumberOfCells();

    auto parentCellIdxView = get_parent_cell_idx_view(contourMesh);
    auto domainIdView = getDomainIdView(contourMesh);

    const axom::IndexType domainCount = computationalMesh.domainCount();

    /*
      Space to store function views and whether a computational cell
      contains the contour.  We set these up for all domains ahead
      of time for accessing as needed later.
    */
    axom::Array<axom::ArrayView<const double, DIM, MemorySpace>> fcnViews(
      domainCount);
    axom::Array<axom::Array<bool, DIM>> hasContours(domainCount);
    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      const auto& dom = computationalMesh.domain(domId);
      axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(dom, "mesh");

      const axom::StackArray<axom::IndexType, DIM> domLengths =
        mvu.getRealExtents("element");

      axom::Array<bool, DIM>& hasContour = hasContours[domId];
      hasContour.resize(domLengths, false);

      fcnViews[domId] =
        mvu.template getConstFieldView<double>(functionName(), false);
    }
    for(axom::IndexType contourCellNum = 0; contourCellNum < cellCount;
        ++contourCellNum)
    {
      axom::IndexType domainId = domainIdView[contourCellNum];
      const axom::StackArray<axom::IndexType, DIM>& parentCellIdx =
        parentCellIdxView[contourCellNum];
      hasContours[domainId][parentCellIdx] = true;
    }

    // Verify that cells marked by hasContours touches the contour
    // and other cells don't.
    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      const auto& dom = computationalMesh.domain(domId);
      axom::quest::MeshViewUtil<DIM, MemorySpace> mvu(dom, "mesh");

      axom::StackArray<axom::IndexType, DIM> domLengths;
      computationalMesh.domainLengths(domId, domLengths);
      assert(domLengths == mvu.getRealExtents("element"));

      const auto& fcnView = fcnViews[domId];

      axom::ArrayIndexer<axom::IndexType, DIM> rowMajor(domLengths, 'r');
      const axom::IndexType cellCount = product(domLengths);
      for(axom::IndexType cellId = 0; cellId < cellCount; ++cellId)
      {
        axom::StackArray<axom::IndexType, DIM> cellIdx =
          rowMajor.toMultiIndex(cellId);

        // Compute min and max function value in the cell.
        double minFcnValue = std::numeric_limits<double>::max();
        double maxFcnValue = std::numeric_limits<double>::min();
        constexpr short int cornerCount =
          (1 << DIM);  // Number of nodes in a cell.
        for(short int cornerId = 0; cornerId < cornerCount; ++cornerId)
        {
          axom::StackArray<axom::IndexType, DIM> cornerIdx = cellIdx;
          for(int d = 0; d < DIM; ++d)
          {
            if(cornerId & (1 << d))
            {
              ++cornerIdx[d];
            }
          }

          double fcnValue = fcnView[cornerIdx];
          minFcnValue = std::min(minFcnValue, fcnValue);
          maxFcnValue = std::max(maxFcnValue, fcnValue);
        }

        const bool touchesContour =
          (minFcnValue <= params.contourVal && maxFcnValue >= params.contourVal);
        const bool hasCont = hasContours[domId][cellIdx];
        if(touchesContour != hasCont)
        {
          ++errCount;
          SLIC_INFO_IF(params.isVerbose(),
                       axom::fmt::format(
                         "checkCellsContainingContour: cell {}: hasContour "
                         "({}) and touchesContour ({}) don't agree.",
                         cellIdx,
                         hasCont,
                         touchesContour));
        }
      }
    }

    SLIC_INFO_IF(params.isVerbose(),
                 axom::fmt::format("checkCellsContainingContour: found {} "
                                   "misrepresented computational cells.",
                                   errCount));
    return errCount;
  }

  /**
   * Change cp_domain data from a local index to a global domain index
   * by adding rank offsets.
   * This is an optional step to make domain ids globally unique.
   */
  void addRankOffsetToContourMeshDomainIds(
    int localDomainCount,
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>& contourMesh)
  {
#ifdef AXOM_USE_MPI
    axom::Array<axom::IndexType> starts(numRanks, numRanks);
    {
      axom::Array<axom::IndexType> indivDomainCounts(numRanks, numRanks);
      indivDomainCounts.fill(-1);
      MPI_Allgather(&localDomainCount,
                    1,
                    MPI_INT,
                    indivDomainCounts.data(),
                    1,
                    MPI_INT,
                    MPI_COMM_WORLD);
      starts[0] = 0;
      for(int i = 1; i < numRanks; ++i)
      {
        starts[i] = starts[i - 1] + indivDomainCounts[i - 1];
      }
    }

    const std::string domainIdField = m_domainIdField;
    auto* domainIdPtr =
      contourMesh.getFieldPtr<axom::IndexType>(domainIdField,
                                               axom::mint::CELL_CENTERED);
    int cellCount = contourMesh.getNumberOfCells();

    for(int i = 0; i < cellCount; ++i)
    {
      domainIdPtr[i] += starts[myRank];
    }
#endif
  }
};

/*!
  @brief Function providing distance from a point.
*/
template <int DIM, typename ExecSpace>
struct RoundContourTest : public ContourTestBase<DIM, ExecSpace>
{
  static constexpr auto MemorySpace =
    axom::execution_space<ExecSpace>::memory_space;
  using PointType = axom::primal::Point<double, DIM>;
  RoundContourTest(const PointType& center)
    : ContourTestBase<DIM, ExecSpace>()
    , _sphere(center, 0.0)
    , _errTol(1e-3)
  { }
  virtual ~RoundContourTest() { }
  const axom::primal::Sphere<double, DIM> _sphere;
  double _errTol;

  virtual std::string name() const override { return std::string("round"); }

  virtual std::string functionName() const override
  {
    return std::string("dist_to_center");
  }

  AXOM_HOST_DEVICE double value(const PointType& pt) const override
  {
    return _sphere.computeSignedDistance(pt);
  }

  double errorTolerance() const override { return _errTol; }

  void setToleranceByLongestEdge(const BlueprintStructuredMesh& bsm)
  {
    double maxSpacing = bsm.maxSpacing();
    _errTol = 0.1 * maxSpacing;
  }
};

/*!
  @brief Function providing signed distance from a plane.
*/
template <int DIM, typename ExecSpace>
struct PlanarContourTest : public ContourTestBase<DIM, ExecSpace>
{
  static constexpr auto MemorySpace =
    axom::execution_space<ExecSpace>::memory_space;
  using PointType = axom::primal::Point<double, DIM>;
  /*!
    @brief Constructor.

    @param inPlane [in] A point in the plane.
    @param perpDir [in] Perpendicular direction on positive side.
  */
  PlanarContourTest(const PointType& inPlane,
                    const axom::primal::Vector<double, DIM>& perpDir)
    : ContourTestBase<DIM, ExecSpace>()
    , _plane(perpDir.unitVector(), inPlane)
  { }
  virtual ~PlanarContourTest() { }
  const axom::primal::Plane<double, DIM> _plane;

  virtual std::string name() const override { return std::string("planar"); }

  virtual std::string functionName() const override
  {
    return std::string("dist_to_plane");
  }

  AXOM_HOST_DEVICE double value(const PointType& pt) const override
  {
    return _plane.signedDistance(pt);
  }

  double errorTolerance() const override { return 1e-15; }
};

/// Utility function to transform blueprint node storage.
void makeCoordsContiguous(conduit::Node& coordValues)
{
  bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
  if(isInterleaved)
  {
    conduit::Node oldValues = coordValues;
    conduit::blueprint::mcarray::to_contiguous(oldValues, coordValues);
  }
}

/// Utility function to transform blueprint node storage.
void makeCoordsInterleaved(conduit::Node& coordValues)
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

  conduit::utils::set_error_handler([](auto& msg, auto& file, int line) {
    slic::logErrorMessage(msg, file, line);
  });
  conduit::utils::set_warning_handler([](auto& msg, auto& file, int line) {
    slic::logWarningMessage(msg, file, line);
  });
  conduit::utils::set_info_handler([](auto& msg, auto& file, int line) {
    slic::logMessage(slic::message::Info, msg, file, line);
  });
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

/*!
  All the test code that depends on DIM to instantiate.
*/
template <int DIM, typename ExecSpace>
int testNdimInstance(BlueprintStructuredMesh& computationalMesh)
{
  // Create marching cubes algorithm object and set some parameters
  quest::MarchingCubes mc(params.policy,
                          computationalMesh.asConduitNode(),
                          "mesh");

  //---------------------------------------------------------------------------
  // params specify which tests to run.
  //---------------------------------------------------------------------------

  std::shared_ptr<RoundContourTest<DIM, ExecSpace>> roundTest;
  std::shared_ptr<PlanarContourTest<DIM, ExecSpace>> planarTest;

  if(params.usingRound)
  {
    roundTest = std::make_shared<RoundContourTest<DIM, ExecSpace>>(
      params.roundContourCenter<DIM>());
    roundTest->setToleranceByLongestEdge(computationalMesh);
    roundTest->computeNodalDistance(computationalMesh);
  }
  if(params.usingPlanar)
  {
    planarTest = std::make_shared<PlanarContourTest<DIM, ExecSpace>>(
      params.inplanePoint<DIM>(),
      params.planeNormal<DIM>());
    planarTest->computeNodalDistance(computationalMesh);
  }
  if(params.isVerbose())
  {
    computationalMesh.printMeshInfo();
  }

  // Write computational mesh with contour functions.
  saveMesh(computationalMesh.asConduitNode(), params.fieldsFile);

  int localErrCount = 0;

  if(planarTest)
  {
    localErrCount += planarTest->runTest(computationalMesh, mc);
  }
  slic::flushStreams();

  if(roundTest)
  {
    localErrCount += roundTest->runTest(computationalMesh, mc);
  }
  slic::flushStreams();

  // Check results

  int errCount = 0;
  if(params.checkResults)
  {
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&localErrCount, &errCount, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    errCount = localErrCount;
#endif

    if(errCount)
    {
      SLIC_INFO(axom::fmt::format(" Error exit: {} errors found.", errCount));
    }
    else
    {
      SLIC_INFO(banner("Normal exit."));
    }
  }
  else
  {
    SLIC_INFO("Results not checked.");
  }

  return errCount;
}

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
#ifdef AXOM_USE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);
#else
  numRanks = 1;
  myRank = 0;
#endif

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

#ifdef AXOM_USE_MPI
    MPI_Bcast(&retval, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Finalize();
#endif

    exit(retval);
  }

  s_allocatorId = allocatorIdForPolicy(params.policy);

  //---------------------------------------------------------------------------
  // Load computational mesh.
  //---------------------------------------------------------------------------
  BlueprintStructuredMesh computationalMesh(params.meshFile, "mesh");

  SLIC_INFO_IF(
    params.isVerbose(),
    axom::fmt::format("Computational mesh has {} cells in {} domains locally",
                      computationalMesh.cellCount(),
                      computationalMesh.domainCount()));
  slic::flushStreams();

  auto getIntMinMax = [](int inVal, int& minVal, int& maxVal, int& sumVal) {
#ifdef AXOM_USE_MPI
    MPI_Allreduce(&inVal, &minVal, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &maxVal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&inVal, &sumVal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
    minVal = inVal;
    maxVal = inVal;
    sumVal = inVal;
#endif
  };

  // Output some global mesh size stats
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.cellCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} cells",
      mn,
      mx,
      sum,
      (double)sum / numRanks));
  }
  {
    int mn, mx, sum;
    getIntMinMax(computationalMesh.domainCount(), mn, mx, sum);
    SLIC_INFO(axom::fmt::format(
      "Computational mesh has {{min:{}, max:{}, sum:{}, avg:{}}} domains",
      mn,
      mx,
      sum,
      (double)sum / numRanks));
  }

  slic::flushStreams();

  //---------------------------------------------------------------------------
  // Run test in the execution space set by command line.
  //---------------------------------------------------------------------------
  int errCount = 0;
  if(params.policy == quest::MarchingCubesRuntimePolicy::seq)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::SEQ_EXEC>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::SEQ_EXEC>(computationalMesh);
    }
  }
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
  else if(params.policy == quest::MarchingCubesRuntimePolicy::omp)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::OMP_EXEC>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::OMP_EXEC>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == quest::MarchingCubesRuntimePolicy::cuda)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::CUDA_EXEC<256>>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::CUDA_EXEC<256>>(computationalMesh);
    }
  }
  #endif
  #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
  else if(params.policy == quest::MarchingCubesRuntimePolicy::hip)
  {
    if(params.ndim == 2)
    {
      errCount = testNdimInstance<2, axom::HIP_EXEC<256>>(computationalMesh);
    }
    else if(params.ndim == 3)
    {
      errCount = testNdimInstance<3, axom::HIP_EXEC<256>>(computationalMesh);
    }
  }
  #endif
#endif

  finalizeLogger();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return errCount != 0;
}
