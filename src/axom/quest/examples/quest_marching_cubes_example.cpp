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
#include "axom/mint.hpp"
#include "axom/quest.hpp"
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

#ifndef __WHERE
  #define __STRINGIZE(x) __STRINGIZE2(x)
  #define __STRINGIZE2(x) #x
  //!@brief String literal for code location
  #define __WHERE \
    __FILE__ ":" __STRINGIZE(__LINE__) "(" + std::string(__func__) + ") "
#endif

namespace quest = axom::quest;
namespace slic = axom::slic;
namespace mint = axom::mint;
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

  size_t ndim {0};

  // TODO: Ensure that fcnCenter, inPlane and perpDir sizes match dimensionality.

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

Input params;

int myRank = -1, numRanks = -1;  // MPI stuff, set in main().

/**
 \brief Generic computational mesh, to hold cell and node data.
*/
struct BlueprintStructuredMesh
{
public:
  explicit BlueprintStructuredMesh(const std::string& meshFile,
                                   const std::string& coordset = "coords",
                                   const std::string& topology = "mesh")
    : _coordsetPath("coordsets/" + coordset)
    , _topologyPath("topologies/" + topology)
  {
    readBlueprintMesh(meshFile);
    for(int d = 0; d < _mdMesh.number_of_children(); ++d)
    {
      auto dl = domainLengths(d);
      SLIC_INFO(axom::fmt::format("dom[{}] size={}", d, dl));
    }
    if(_ndims == 2)
    {
      checkMeshStorage<2>();
    }
    if(_ndims == 3)
    {
      checkMeshStorage<3>();
    }
  }

  /// Return the blueprint mesh in a conduit::Node
  conduit::Node& asConduitNode() { return _mdMesh; }

  /// Get number of domains in the multidomain mesh
  axom::IndexType domainCount() const { return _domCount; }

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
  */
  double maxSpacing() const
  {
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

    This method takes shorcuts by assuming
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
    // info.print();
    return true;
  }

  void printMeshInfo() const { _mdMesh.print(); }

  template <int DIM>
  axom::StackArray<axom::IndexType, DIM> getCellsShape(int domainNum)
  {
    auto domLengths = domainLengths(domainNum);
    axom::StackArray<axom::IndexType, DIM> shape;
    for(int i = 0; i < DIM; ++i)
    {
      shape[i] = domLengths[i];
    }
    return shape;
  }

  template <int DIM>
  axom::StackArray<axom::IndexType, DIM> getNodesShape(int domainNum)
  {
    auto domLengths = domainLengths(domainNum);
    axom::StackArray<axom::IndexType, DIM> shape;
    for(int i = 0; i < DIM; ++i)
    {
      shape[i] = 1 + domLengths[i];
    }
    return shape;
  }

  template <int DIM>
  axom::Array<axom::ArrayView<double, DIM>> get_coords_view(int domainNum)
  {
    axom::StackArray<axom::IndexType, DIM> shape = getNodesShape<DIM>(domainNum);
    for(int d = 0; d < DIM / 2; ++d) std::swap(shape[d], shape[DIM - 1 - d]);

    conduit::Node& dom = _mdMesh[domainNum];
    conduit::Node& coordsValues = dom.fetch_existing("coordsets/coords/values");
    axom::Array<axom::ArrayView<double, DIM>> coordsViews(DIM);
    for(int d = 0; d < DIM; ++d)
    {
      double* ptr = coordsValues[d].as_double_ptr();
      coordsViews[d] = axom::ArrayView<double, DIM>(ptr, shape);
    }
    return coordsViews;
  }

  template <int DIM>
  typename std::enable_if<DIM == 2>::type checkMeshStorage()
  {
    SLIC_ASSERT(dimension() == DIM);
    for(int d = 0; d < _mdMesh.number_of_children(); ++d)
    {
      axom::Array<axom::ArrayView<double, DIM>> coordsViews =
        get_coords_view<DIM>(d);
      const axom::StackArray<axom::IndexType, DIM>& shape =
        coordsViews[0].shape();

      // Verify that i is slowest in m_coordsViews.
      // It appears conduit stores column major and ArrayView computes offsets
      // assuming row major.
      int n = 0, errCount = 0;
      for(int j = 0; j < shape[0]; ++j)
      {
        for(int i = 0; i < shape[1]; ++i)
        {
          errCount += (&coordsViews[0].data()[n++] != &coordsViews[0](j, i));
        }
      }
      assert(errCount == 0);
    }
  }
  template <int DIM>
  typename std::enable_if<DIM == 3>::type checkMeshStorage()
  {
    SLIC_ASSERT(dimension() == DIM);
    for(int d = 0; d < _mdMesh.number_of_children(); ++d)
    {
      axom::Array<axom::ArrayView<double, DIM>> coordsViews =
        get_coords_view<DIM>(d);
      const axom::StackArray<axom::IndexType, DIM>& shape =
        coordsViews[0].shape();

      // Verify that i is slowest in m_coordsViews.
      // It appears conduit stores column major and ArrayView computes offsets
      // assuming row major.
      int n = 0, errCount = 0;
      for(int k = 0; k < shape[0]; ++k)
      {
        for(int j = 0; j < shape[1]; ++j)
        {
          for(int i = 0; i < shape[2]; ++i)
          {
            errCount += (&coordsViews[0].data()[n++] != &coordsViews[0](k, j, i));
          }
        }
      }
      assert(errCount == 0);
    }
  }

private:
  int _ndims {-1};
  conduit::Node _mdMesh;
  axom::IndexType _domCount;
  const std::string _coordsetPath;
  const std::string _topologyPath;

  /*!
    @brief Read a blueprint mesh.
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
    assert(conduit::blueprint::mesh::is_multi_domain(_mdMesh));
    _domCount = conduit::blueprint::mesh::number_of_domains(_mdMesh);

    if(_domCount > 0)
    {
      const conduit::Node coordsetNode = _mdMesh[0].fetch_existing(_coordsetPath);
      _ndims = conduit::blueprint::mesh::coordset::dims(coordsetNode);
    }
#ifdef AXOM_USE_MPI
    MPI_Allreduce(MPI_IN_PLACE, &_ndims, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
#endif
    SLIC_ASSERT(_ndims > 0);

    bool valid = isValid();
    SLIC_ASSERT(valid);
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
  conduit::relay::io::blueprint::saveMesh(mesh, filename, "hdf5");
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
#ifdef AXOM_USE_MPI
  conduit::relay::mpi::io::blueprint::save_mesh(tmpMesh,
                                                filename,
                                                "hdf5",
                                                MPI_COMM_WORLD);
#else
  conduit::relay::io::blueprint::save_mesh(tmpMesh, filename, "hdf5");
#endif
}

//!@brief Reverse the order of a StackArray.
template <typename T, int DIM>
void reverse(axom::StackArray<T, DIM>& a)
{
  for(int d = 0; d < DIM / 2; ++d)
  {
    std::swap(a[d], a[DIM - 1 - d]);
  }
  return;
}

template <typename T, int DIM>
T product(const axom::StackArray<T, DIM>& a)
{
  T rval = a[0];
  for(int d = 1; d < DIM; ++d) rval *= a[d];
  return rval;
}

template <typename T, int DIM, typename U>
static void addToStackArray(axom::StackArray<T, DIM>& a, U b)
{
  for(int d = 0; d < DIM; ++d) a[d] += b;
}

/*!
  @brief Compute multidmensional index corresponding to flat index
  in rectangular index space.

  The flatId corresponds to the multidimensional index through
  an order that advances the first index fastest.
*/
template <int DIM>
axom::StackArray<axom::IndexType, DIM> flatToMultidimIndex(
  axom::IndexType flatId,
  const axom::StackArray<axom::IndexType, DIM>& sizes)
{
  axom::IndexType strides[DIM] = {1};
  for(int d = 1; d < DIM; ++d) strides[d] = strides[d - 1] * sizes[d - 1];
  if(flatId >= strides[DIM - 1] * sizes[DIM - 1])
  {
    SLIC_ERROR("flatId is too big.");
  }

  axom::StackArray<axom::IndexType, DIM> rval;
  for(int d = DIM - 1; d >= 0; --d)
  {
    rval[d] = flatId / strides[d];
    flatId -= rval[d] * strides[d];
  }
  return rval;
}

template <int DIM>
struct ContourTestBase
{
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
  virtual double value(const axom::primal::Point<double, DIM>& pt) const = 0;

  //!@brief Return error tolerance for contour mesh accuracy check.
  virtual double errorTolerance() const = 0;

  const std::string m_parentCellIdField;
  const std::string m_domainIdField;

  int runTest(BlueprintStructuredMesh& computationalMesh, quest::MarchingCubes& mc)
  {
    SLIC_INFO(banner(axom::fmt::format("Testing {} contour.", name())));

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
      // Access the coordinates
      auto cellCounts = bpMesh.domainLengths(domId);
      axom::StackArray<axom::IndexType, DIM> shape;
      for(int i = 0; i < DIM; ++i)
      {
        shape[i] = 1 + cellCounts[i];
      }
      conduit::index_t pointCount = bpMesh.nodeCount(domId);
      conduit::Node& coordsValues =
        dom.fetch_existing("coordsets/coords/values");
      double* coordsPtrs[DIM];
      axom::ArrayView<double, DIM> coordsViews[DIM];
      for(int d = 0; d < DIM; ++d)
      {
        coordsPtrs[d] = coordsValues[d].as_double_ptr();
        coordsViews[d] = axom::ArrayView<double, DIM>(coordsPtrs[d], shape);
      }
      axom::ArrayView<axom::primal::Point<double, DIM>, DIM> coordsView(
        (axom::primal::Point<double, DIM>*)coordsValues.data_ptr(),
        shape);

      // Create the nodal function data.
      conduit::Node& fieldNode = dom["fields"][functionName()];
      fieldNode["association"] = "vertex";
      fieldNode["topology"] = "mesh";
      fieldNode["volume_dependent"] = "false";
      fieldNode["values"].set(conduit::DataType::float64(pointCount));
      double* d = fieldNode["values"].value();
      axom::ArrayView<double, DIM> fieldView(d, shape);

      // Set the function value at the nodes.
      // value(pt) is the virtual function defining the
      // distance in a derived class.
      for(int n = 0; n < pointCount; ++n)
      {
        axom::primal::Point<double, DIM> pt;
        for(int d = 0; d < DIM; ++d)
        {
          pt[d] = coordsViews[d].flatIndex(n);
        }
        fieldView.flatIndex(n) = value(pt);
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
    axom::primal::Point<double, DIM> pt;
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
    axom::Array<axom::ArrayView<const double, DIM>> coordsViews(domainCount * DIM);

    // Get info about the computational domains available for look-up.
    axom::Array<axom::StackArray<axom::IndexType, DIM>> domainLengths(domainCount);
    for(axom::IndexType n = 0; n < domainCount; ++n)
    {
      auto& domain = computationalMesh.domain(n);
      axom::ArrayView<const double, DIM>* domainCoordsView =
        &coordsViews[DIM * n];
      getCoordsViews(domain, computationalMesh.coordsetPath(), domainCoordsView);

      axom::Array<axom::IndexType> domLengths =
        computationalMesh.domainLengths(n);
      for(int d = 0; d < DIM; ++d)
      {
        domainLengths[n][d] = domLengths[d];
      }
    }

    for(axom::IndexType cn = 0; cn < cellCount; ++cn)
    {
      axom::IndexType domainId = domainIdView[cn];

      axom::StackArray<axom::IndexType, DIM> parentCellIdx =
        parentCellIdxView[cn];
      reverse(parentCellIdx);  // ArrayView expects indices in reverse order.
      // This should work but breaks gcc11 on 64-bit linux:
      // axom::StackArray<axom::IndexType, DIM> upperIdx = parentCellIdx + 1;
      axom::StackArray<axom::IndexType, DIM> upperIdx = parentCellIdx;
      addToStackArray(upperIdx, 1);

      axom::ArrayView<const double, DIM>* domainCoordsView =
        &coordsViews[DIM * domainId];
      axom::primal::Point<double, DIM> lower, upper;
      for(int d = 0; d < DIM; ++d)
      {
        lower[d] = domainCoordsView[d][parentCellIdx];
        upper[d] = domainCoordsView[d][upperIdx];
      }
      axom::primal::BoundingBox<double, DIM> parentCellBox(lower, upper);
      double tol = errorTolerance();
      axom::primal::BoundingBox<double, DIM> big(parentCellBox);
      axom::primal::BoundingBox<double, DIM> small(parentCellBox);
      big.expand(tol);
      small.expand(-tol);

      // WRONG: the node ids should increased by the number of nodes in all previous domains.
      axom::IndexType* cellNodeIds = contourMesh.getCellNodeIDs(cn);
      const axom::IndexType cellNodeCount = contourMesh.getNumberOfCellNodes(cn);

      for(axom::IndexType nn = 0; nn < cellNodeCount; ++nn)
      {
        axom::primal::Point<double, DIM> nodeCoords;
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

    // Nodal values of functions for each domain.
    axom::Array<axom::ArrayView<const double, DIM>> fcnViews(domainCount);
    // Whether computational cells have parts of the contour mesh.
    axom::Array<axom::Array<bool, DIM>> hasContours(domainCount);

    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      axom::StackArray<axom::IndexType, DIM> domLengths;
      computationalMesh.domainLengths(domId, domLengths);

      axom::Array<bool, DIM>& hasContour = hasContours[domId];
      hasContour.resize(domLengths, false);

      // Add 1 to count nodes.  Reverse to match Conduit data.
      reverse(domLengths);
      // This should work but breaks gcc11 on 64-bit linux:
      // domLengths = domLengths + 1;
      addToStackArray(domLengths, 1);

      conduit::Node& dom = computationalMesh.domain(domId);
      double* fcnPtr =
        dom.fetch_existing("fields/" + functionName() + "/values").as_double_ptr();
      fcnViews[domId] = axom::ArrayView<const double, DIM>(fcnPtr, domLengths);
    }
    for(axom::IndexType cn = 0; cn < cellCount; ++cn)
    {
      axom::IndexType domainId = domainIdView[cn];
      const axom::StackArray<axom::IndexType, DIM>& parentCellIdx =
        parentCellIdxView[cn];
      hasContours[domainId][parentCellIdx] = true;
    }

    // Verify that marked cells contain the contour value
    // unmarked ones don't.
    for(axom::IndexType domId = 0; domId < domainCount; ++domId)
    {
      axom::StackArray<axom::IndexType, DIM> domLengths;
      computationalMesh.domainLengths(domId, domLengths);

      axom::ArrayView<const double, DIM>& fcnView = fcnViews[domId];

      const axom::IndexType cellCount = product(domLengths);
      for(axom::IndexType cellId = 0; cellId < cellCount; ++cellId)
      {
        axom::StackArray<axom::IndexType, DIM> cellIdx =
          flatToMultidimIndex(cellId, domLengths);

        // Compute min and max function value in the cell.
        double minFcnValue = std::numeric_limits<double>::max();
        double maxFcnValue = std::numeric_limits<double>::min();
        constexpr short int cornerCount =
          (1 << DIM);  // Number of nodes in a cell.
        for(short int cornerId = 0; cornerId < cornerCount; ++cornerId)
        {
          axom::StackArray<axom::IndexType, DIM> cornerIdx = cellIdx;
          for(int d = 0; d < DIM; ++d)
            if(cornerId & (1 << d)) ++cornerIdx[d];

          reverse(cornerIdx);
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

  void getCoordsViews(conduit::Node& domain,
                      const std::string& coordsetPath,
                      axom::ArrayView<const double, DIM> coordsViews[DIM])
  {
    const conduit::Node& dimsNode =
      domain.fetch_existing("topologies/mesh/elements/dims");
    axom::StackArray<axom::IndexType, DIM> coordsViewShape;
    for(int d = 0; d < DIM; ++d)
      coordsViewShape[d] = 1 + dimsNode[DIM - 1 - d].as_int();

    const conduit::Node& coordValues =
      domain.fetch_existing(coordsetPath + "/values");
    bool isInterleaved = conduit::blueprint::mcarray::is_interleaved(coordValues);
    const int coordSp = isInterleaved ? DIM : 1;

    for(int d = 0; d < DIM; ++d)
    {
      auto* coordsPtr = coordValues[d].as_double_ptr();
      coordsViews[d] =
        axom::ArrayView<const double, DIM>(coordsPtr, coordsViewShape, coordSp);
    }
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
template <int DIM>
struct RoundContourTest : public ContourTestBase<DIM>
{
  RoundContourTest(const axom::primal::Point<double, DIM>& pt)
    : ContourTestBase<DIM>()
    , _center(pt)
    , _errTol(1e-3)
  { }
  virtual ~RoundContourTest() { }
  const axom::primal::Point<double, DIM> _center;
  double _errTol;

  virtual std::string name() const override { return std::string("round"); }

  virtual std::string functionName() const override
  {
    return std::string("dist_to_center");
  }

  double value(const axom::primal::Point<double, DIM>& pt) const override
  {
    double dist = primal::squared_distance(_center, pt);
    dist = sqrt(dist);
    return dist;
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
template <int DIM>
struct PlanarContourTest : public ContourTestBase<DIM>
{
  /*!
    @brief Constructor.

    @param inPlane [in] A point in the plane.
    @param perpDir [in] Perpendicular direction on positive side.
  */
  PlanarContourTest(const axom::primal::Point<double, DIM>& inPlane,
                    const axom::primal::Vector<double, DIM>& perpDir)
    : ContourTestBase<DIM>()
    , _inPlane(inPlane)
    , _normal(perpDir.unitVector())
  { }
  virtual ~PlanarContourTest() { }
  const axom::primal::Point<double, DIM> _inPlane;
  const axom::primal::Vector<double, DIM> _normal;

  virtual std::string name() const override { return std::string("planar"); }

  virtual std::string functionName() const override
  {
    return std::string("dist_to_plane");
  }

  double value(const axom::primal::Point<double, DIM>& pt) const override
  {
    axom::primal::Vector<double, DIM> r(_inPlane, pt);
    double dist = r.dot(_normal);
    return dist;
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
template <int DIM>
int testNdimInstance(BlueprintStructuredMesh& computationalMesh)
{
  // Create marching cubes algorithm object and set some parameters
  quest::MarchingCubes mc(params.policy,
                          computationalMesh.asConduitNode(),
                          "coords");

  //---------------------------------------------------------------------------
  // params specify which tests to run.
  //---------------------------------------------------------------------------

  std::shared_ptr<RoundContourTest<DIM>> roundTest;
  std::shared_ptr<PlanarContourTest<DIM>> planarTest;

  if(params.usingRound)
  {
    roundTest =
      std::make_shared<RoundContourTest<DIM>>(params.roundContourCenter<DIM>());
    roundTest->setToleranceByLongestEdge(computationalMesh);
    roundTest->computeNodalDistance(computationalMesh);
  }
  if(params.usingPlanar)
  {
    planarTest =
      std::make_shared<PlanarContourTest<DIM>>(params.inplanePoint<DIM>(),
                                               params.planeNormal<DIM>());
    planarTest->computeNodalDistance(computationalMesh);
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

  //---------------------------------------------------------------------------
  // Load computational mesh.
  //---------------------------------------------------------------------------
  BlueprintStructuredMesh computationalMesh(params.meshFile);
  // computationalMesh.printMeshInfo();

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

  int errCount = 0;
  if(params.ndim == 2)
  {
    errCount = testNdimInstance<2>(computationalMesh);
  }
  else if(params.ndim == 3)
  {
    errCount = testNdimInstance<3>(computationalMesh);
  }

  finalizeLogger();
#ifdef AXOM_USE_MPI
  MPI_Finalize();
#endif

  return errCount != 0;
}
