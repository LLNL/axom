// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/mir.hpp"

#include <string>

// namespace aliases
namespace numerics = axom::numerics;
namespace slam = axom::slam;
namespace mir = axom::mir;
namespace fs = axom::utilities::filesystem;
namespace bputils = axom::mir::utilities::blueprint;

using RuntimePolicy = axom::runtime_policy::Policy;

//--------------------------------------------------------------------------------
/// Contain program options.
struct Input
{
  int m_test_case{1};  // valid values 1,2,3,4,5
  bool m_should_iterate{false};
  int m_iter_count{0};
  double m_iter_percent{0.};
  bool m_verbose{false};
  std::string m_output_dir{};
  RuntimePolicy m_policy {RuntimePolicy::seq};
  std::string m_annotationMode{"report"};
  axom::CLI::App m_app{};

  /// Parse command line.
  int parse(int argc, char** argv)
  {
    m_app.add_option("--test-case", m_test_case)
      ->description("Select the test case.");

    m_app.add_option("--output-dir", m_output_dir)
      ->description("The directory for output files");

    m_app.add_option("--iter-count", m_iter_count)
      ->description("The number of iterations for MIR");

    m_app.add_option("--iter-percent", m_iter_percent)
      ->description("The percent error for iterative MIR");

    m_app.add_flag("--verbose", m_verbose)
      ->description("Verbose output");

#if defined(AXOM_USE_CALIPER)
    m_app.add_option("--caliper", m_annotationMode)
      ->description(
        "caliper annotation mode. Valid options include 'none' and 'report'. "
        )
      ->capture_default_str()
      ->check(axom::utilities::ValidCaliperMode);
#endif

    std::stringstream pol_sstr;
    pol_sstr << "Set MIR runtime policy.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    pol_sstr << "\nSet to 'seq' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
    pol_sstr << "\nSet to 'omp' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
    pol_sstr << "\nSet to 'cuda' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
    pol_sstr << "\nSet to 'hip' or 3 to use the RAJA HIP policy.";
  #endif
#endif
    m_app.add_option("-p, --policy", m_policy, pol_sstr.str())
      ->capture_default_str()
      ->transform(
        axom::CLI::CheckedTransformer(axom::runtime_policy::s_nameToPolicy));

    // Parse command line options.
    try
    {
      m_app.parse(argc, argv);
    }
    catch (const axom::CLI::ParseError &e)
    {
      return m_app.exit(e);
    }

    int retval = 0;
    checkTestCase(retval);
    checkOutputDir();
    checkIterationParams(retval);
    return retval;
  }

  bool shouldIterate() const { return m_should_iterate; }
  int numIterations() const { return m_iter_count; }
  int iterPercentage() const { return m_iter_percent; }

private:
  void checkTestCase(int &retval)
  {
    if(m_test_case < 1 || m_test_case > 6)
    {
      retval = -1;
      SLIC_ERROR("Invalid test case " << m_test_case);
    }
  }

  void checkOutputDir()
  {
    if(!fs::pathExists(m_output_dir))
    {
      fs::makeDirsForPath(m_output_dir);
    }
  }

  void checkIterationParams(int &retval)
  {
    if(m_should_iterate)
    {
      if(m_iter_count < 1)
      {
        retval = -2;
        SLIC_ERROR("Invalid iteration count " << m_iter_count);
      }

      if(m_iter_percent <= 0. || m_iter_percent > 1.)
      {
        retval = -3;
        SLIC_ERROR("Invalid iteration percentage " << m_iter_percent);
      }
    }
  }
};

//--------------------------------------------------------------------------------
/// Print a Conduit node.
void printNode(const conduit::Node &n)
{
  conduit::Node options;
  options["num_children_threshold"] = 10000;
  options["num_elements_threshold"] = 10000;
  n.to_summary_string_stream(std::cout, options);
}

//--------------------------------------------------------------------------------
/*!
 * \brief Run MIR on the input mesh.
 *
 * \tparam ExecSpace The execution space where the algorithm will run.
 *
 * \param hostMesh A conduit node that contains the test mesh.
 * \param options A conduit node that contains the test mesh.
 * \param hostResult A conduit node that will contain the MIR results.
 */
template <typename ExecSpace>
int runMIR(const conduit::Node &hostMesh, const conduit::Node &options, conduit::Node &hostResult)
{
  std::string shape = hostMesh["topologies/mesh/elements/shape"].as_string();
  SLIC_INFO(axom::fmt::format("Using policy {}",
                              axom::execution_space<ExecSpace>::name()));

  // host->device
  conduit::Node deviceMesh;
  bputils::copy<ExecSpace>(deviceMesh, hostMesh);

  conduit::Node &n_coordset = deviceMesh["coordsets/coords"];
  conduit::Node &n_topo = deviceMesh["topologies/mesh"];
  conduit::Node &n_matset = deviceMesh["matsets/mat"];
  auto connView = bputils::make_array_view<int>(n_topo["elements/connectivity"]);

  // Make matset view. (There's often 1 more material so add 1)
  constexpr int MAXMATERIALS = 12;
  using MatsetView = axom::mir::views::UnibufferMaterialView<int, float, MAXMATERIALS + 1>;
  MatsetView matsetView;
  matsetView.set(
    bputils::make_array_view<int>(n_matset["material_ids"]),
    bputils::make_array_view<float>(n_matset["volume_fractions"]),
    bputils::make_array_view<int>(n_matset["sizes"]),
    bputils::make_array_view<int>(n_matset["offsets"]),
    bputils::make_array_view<int>(n_matset["indices"]));

  // Coord/Topo views differ.
  conduit::Node deviceResult;
  if(shape == "tri")
  {
    auto coordsetView = axom::mir::views::make_explicit_coordset<float, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::TriShape<int>>;
    TopologyView topologyView(connView);

    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceResult);   
  }
  else if(shape == "quad")
  {
    auto coordsetView = axom::mir::views::make_explicit_coordset<float, 2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::QuadShape<int>>;
    TopologyView topologyView(connView);

    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceResult);   
  }
  else if(shape == "hex")
  {
    auto coordsetView = axom::mir::views::make_explicit_coordset<float, 3>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::HexShape<int>>;
    TopologyView topologyView(connView);

    using MIR = axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    m.execute(deviceMesh, options, deviceResult);   
  }

  // device->host
  bputils::copy<axom::SEQ_EXEC>(hostResult, deviceResult);

  return 0;
}

//--------------------------------------------------------------------------------
/*!
 * \brief Tutorial main showing how to initialize test cases and perform mir.
 */
int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;  // create & initialize test logger
  axom::slic::setLoggingMsgLevel(axom::slic::message::Info);

  // Parse arguments
  Input params;
  int retval = params.parse(argc, argv);
  if(retval != 0)
  {
    return retval;
  }
#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    params.m_annotationMode);
#endif

  // Make the mesh
  conduit::Node mesh;
  mir::MeshTester tester;
  auto timer = axom::utilities::Timer(true);
  switch(params.m_test_case)
  {
  case 1:
    tester.initTestCaseOne(mesh);
    break;
  case 2:
    tester.initTestCaseTwo(mesh);
    break;
  case 3:
    tester.initTestCaseThree(mesh);
    break;
  case 4:
    tester.initTestCaseFour(mesh);
    break;
  case 5:
    {
      constexpr int GRIDSIZE = 25;
      constexpr int MAXMATERIALS = 12;
      tester.initTestCaseFive(GRIDSIZE, MAXMATERIALS, mesh);
    }
    break;
  case 6:
    {
      constexpr int GRIDSIZE = 15;
      constexpr int MAXMATERIALS = 3;
      tester.initTestCaseSix(GRIDSIZE, MAXMATERIALS, mesh);
    }
    break;
  }
  timer.stop();
  SLIC_INFO("Mesh init time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Save input mesh
  std::string filepath, filename("inputMesh");
  if(params.m_output_dir.empty())
    filepath = filename;
  else
    filepath = axom::utilities::filesystem::joinPath(params.m_output_dir, filename);
  conduit::relay::io::blueprint::save_mesh(mesh, filepath, "hdf5");

  if(params.m_verbose)
  {
    SLIC_INFO("Initial mesh:");
    printNode(mesh);
  }

  // Begin material interface reconstruction
  timer.start();

  // Set up options.
  conduit::Node options;
  options["matset"] = "mat";
  // Future options
  options["iterate"] = params.shouldIterate() ? 1 : 0;
  options["iterate_percentage"] = params.iterPercentage();

  // Run MIR
  conduit::Node resultMesh;
  if(params.m_policy == RuntimePolicy::seq)
  {
    retval = runMIR<axom::SEQ_EXEC>(mesh, options, resultMesh);
  }
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #if defined(AXOM_USE_OPENMP)
  else if(params.m_policy == RuntimePolicy::omp)
  {
    retval = runMIR<axom::OMP_EXEC>(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_CUDA)
  else if(params.m_policy == RuntimePolicy::cuda)
  {
    constexpr int CUDA_BLOCK_SIZE = 256;
    using cuda_exec = axom::CUDA_EXEC<CUDA_BLOCK_SIZE>;
    retval = runMIR<cuda_exec>(mesh, options, resultMesh);
  }
  #endif
  #if defined(AXOM_USE_HIP)
  else if(params.m_policy == RuntimePolicy::hip)
  {
    constexpr int HIP_BLOCK_SIZE = 64;
    using hip_exec = axom::HIP_EXEC<HIP_BLOCK_SIZE>;
    retval = runMIR<hip_exec>(mesh, options, resultMesh);
  }
  #endif
#endif
  else
  {
    retval = -1;
    SLIC_ERROR("Unhandled policy.");
  }
  timer.stop();
  SLIC_INFO("Reconstruction time: " << timer.elapsedTimeInMilliSec() << " ms.");

  // Save output.
  if(retval == 0)
  {
    std::string filepath, filename("processedMesh");
    if(params.m_output_dir.empty())
      filepath = filename;
    else
      filepath = axom::utilities::filesystem::joinPath(params.m_output_dir, filename);
    conduit::relay::io::blueprint::save_mesh(resultMesh, filepath, "hdf5");
  }

  if(params.m_verbose)
  {
    SLIC_INFO("Final mesh:");
    printNode(resultMesh);
  }

  return retval;
}

//--------------------------------------------------------------------------------
