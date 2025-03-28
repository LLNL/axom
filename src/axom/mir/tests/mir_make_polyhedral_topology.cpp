// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

//------------------------------------------------------------------------------

// Uncomment to generate baselines
// #define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
// #define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(dataDirectory(), "mir", "regression",
               "mir_make_polyhedral_topology");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct make_polyhedral
{
  static void initialize(const std::string &type, conduit::Node &n_mesh)
  {
    axom::StackArray<axom::IndexType, 3> dims {4, 4, 4};
    axom::mir::testing::data::braid(type, dims, n_mesh);
  }

  static void test(const std::string &type, const std::string &name)
  {
    namespace bputils = axom::mir::utilities::blueprint;
    namespace views = axom::mir::views;

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(type, hostMesh);
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
    conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
    conduit::relay::io::save(hostMesh, name + "_orig.yaml", "yaml");
#endif

    //_mir_utilities_makepolyhedraltopology_begin
    // Run the algorithm
    const conduit::Node &n_input = deviceMesh["topologies/mesh"];
    conduit::Node &n_output = deviceMesh["topologies/polymesh"];
    if(type == "uniform")
    {
      auto topologyView = views::make_uniform<3>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bputils::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bputils::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(
        n_output);
    }
    //_mir_utilities_makepolyhedraltopology_end
    else if(type == "tets")
    {
      auto topologyView =
        views::make_unstructured_single_shape<views::TetShape<int>>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bputils::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bputils::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(
        n_output);
    }
    else if(type == "pyramids")
    {
      auto topologyView =
        views::make_unstructured_single_shape<views::PyramidShape<int>>::view(
          n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bputils::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bputils::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(
        n_output);
    }
    else if(type == "wedges")
    {
      auto topologyView =
        views::make_unstructured_single_shape<views::WedgeShape<int>>::view(
          n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bputils::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bputils::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(
        n_output);
    }
    else if(type == "hexs")
    {
      auto topologyView =
        views::make_unstructured_single_shape<views::HexShape<int>>::view(n_input);
      using TopologyView = decltype(topologyView);
      using ConnectivityType = typename TopologyView::ConnectivityType;

      bputils::MakePolyhedralTopology<ExecSpace, TopologyView> mp(topologyView);
      mp.execute(n_input, n_output);
      bputils::MergePolyhedralFaces<ExecSpace, ConnectivityType>::execute(
        n_output);
    }
    else
    {
      SLIC_ERROR(axom::fmt::format("Unsupported shape in the test: {}", type));
    }

    // device->host
    conduit::Node hostOutputMesh;
    axom::mir::utilities::blueprint::copy<seq_exec>(hostOutputMesh, deviceMesh);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  #if defined(AXOM_USE_HDF5)
    conduit::relay::io::blueprint::save_mesh(hostOutputMesh, name, "hdf5");
  #endif
    conduit::relay::io::save(hostOutputMesh, name + ".yaml", "yaml");
#endif
    // Handle baseline comparison.
    {
      std::string baselineName(yamlRoot(name));
      const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
      saveBaseline(paths, baselineName, hostOutputMesh);
#else
      constexpr double tolerance = 2.6e-06;
      EXPECT_TRUE(compareBaseline(paths, baselineName, hostOutputMesh, tolerance));
#endif
    }
  }
};

//------------------------------------------------------------------------------
TEST(mir_make_polyhedral_topology, uniform_seq)
{
  AXOM_ANNOTATE_SCOPE("uniform_seq");
  make_polyhedral<seq_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(mir_make_polyhedral_topology, tets_seq)
{
  AXOM_ANNOTATE_SCOPE("tets_seq");
  make_polyhedral<seq_exec>::test("tets", "make_polyhedral_tets");
}

TEST(mir_make_polyhedral_topology, pyramids_seq)
{
  AXOM_ANNOTATE_SCOPE("pyramids_seq");
  make_polyhedral<seq_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(mir_make_polyhedral_topology, wedges_seq)
{
  AXOM_ANNOTATE_SCOPE("wedges_seq");
  make_polyhedral<seq_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(mir_make_polyhedral_topology, hexs_seq)
{
  AXOM_ANNOTATE_SCOPE("hexs_seq");
  make_polyhedral<seq_exec>::test("hexs", "make_polyhedral_hexs");
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_make_polyhedral_topology, uniform_omp)
{
  AXOM_ANNOTATE_SCOPE("uniform_omp");
  make_polyhedral<omp_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(mir_make_polyhedral_topology, tets_omp)
{
  AXOM_ANNOTATE_SCOPE("tets_omp");
  make_polyhedral<omp_exec>::test("tets", "make_polyhedral_tets");
}

TEST(mir_make_polyhedral_topology, pyramids_omp)
{
  AXOM_ANNOTATE_SCOPE("pyramids_omp");
  make_polyhedral<omp_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(mir_make_polyhedral_topology, wedges_omp)
{
  AXOM_ANNOTATE_SCOPE("wedges_omp");
  make_polyhedral<omp_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(mir_make_polyhedral_topology, hexs_omp)
{
  AXOM_ANNOTATE_SCOPE("hexs_omp");
  make_polyhedral<omp_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_make_polyhedral_topology, uniform_cuda)
{
  AXOM_ANNOTATE_SCOPE("uniform_cuda");
  make_polyhedral<cuda_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(mir_make_polyhedral_topology, tets_cuda)
{
  AXOM_ANNOTATE_SCOPE("tets_cuda");
  make_polyhedral<cuda_exec>::test("tets", "make_polyhedral_tets");
}

TEST(mir_make_polyhedral_topology, pyramids_cuda)
{
  AXOM_ANNOTATE_SCOPE("pyramids_cuda");
  make_polyhedral<cuda_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(mir_make_polyhedral_topology, wedges_cuda)
{
  AXOM_ANNOTATE_SCOPE("wedges_cuda");
  make_polyhedral<cuda_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(mir_make_polyhedral_topology, hexs_cuda)
{
  AXOM_ANNOTATE_SCOPE("hexs_cuda");
  make_polyhedral<cuda_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_make_polyhedral_topology, uniform_hip)
{
  AXOM_ANNOTATE_SCOPE("uniform_hip");
  make_polyhedral<hip_exec>::test("uniform", "make_polyhedral_uniform");
}

TEST(mir_make_polyhedral_topology, tets_hip)
{
  AXOM_ANNOTATE_SCOPE("tets_hip");
  make_polyhedral<hip_exec>::test("tets", "make_polyhedral_tets");
}

TEST(mir_make_polyhedral_topology, pyramids_hip)
{
  AXOM_ANNOTATE_SCOPE("pyramids_hip");
  make_polyhedral<hip_exec>::test("pyramids", "make_polyhedral_pyramids");
}

TEST(mir_make_polyhedral_topology, wedges_hip)
{
  AXOM_ANNOTATE_SCOPE("wedges_hip");
  make_polyhedral<hip_exec>::test("wedges", "make_polyhedral_wedges");
}

TEST(mir_make_polyhedral_topology, hexs_hip)
{
  AXOM_ANNOTATE_SCOPE("hexs_hip");
  make_polyhedral<hip_exec>::test("hexs", "make_polyhedral_hexs");
}
#endif

//------------------------------------------------------------------------------
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  // Define command line options.
  axom::CLI::App app;
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  bool handlerEnabled = false;
  app.add_flag("--handler", handlerEnabled, "Enable Conduit handler.");

  // Parse command line options.
  try
  {
    app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
    axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
      annotationMode);
#endif

    axom::slic::SimpleLogger logger;  // create & initialize test logger,
    if(handlerEnabled)
    {
      conduit::utils::set_error_handler(conduit_debug_err_handler);
    }

    result = RUN_ALL_TESTS();
  }
  catch(axom::CLI::CallForHelp &e)
  {
    std::cout << app.help() << std::endl;
    result = 0;
  }
  catch(axom::CLI::ParseError &e)
  {
    // Handle other parsing errors
    std::cerr << e.what() << std::endl;
    result = app.exit(e);
  }

  return result;
}
