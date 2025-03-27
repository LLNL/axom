// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"
#include "axom/mir/MeshTester.hpp"

//------------------------------------------------------------------------------

// Uncomment to generate baselines
#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "mir"), "regression"),
               "mir_elvira3d");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_Elvira3D
{
  static const int gridSize = 20;   //10;
  static const int numSpheres = 2;  //4;

  static void initialize(conduit::Node &n_mesh)
  {
    AXOM_ANNOTATE_SCOPE("initialize");
    axom::mir::MeshTester M;
    M.setStructured(true);
    M.initTestCaseSix(gridSize, numSpheres, n_mesh);
  }

  // Select a chunk of zones.
  static void selectZones(conduit::Node &n_options)
  {
    std::vector<axom::IndexType> selected;
    for(int k = 0; k < gridSize; k++)
    {
      for(int j = 0; j < gridSize; j++)
      {
        for(int i = 0; i < gridSize; i++)
        {
          // Save all but an octant so we can see inside.
          if(i > gridSize / 2 && j > gridSize / 2 && k > gridSize / 2)
          {
            continue;
          }
          auto idx = k * gridSize * gridSize + j * gridSize + i;
          selected.push_back(idx);
        }
      }
    }
    n_options["selectedZones"].set(selected);
  }

  static void test(const std::string &name, bool selectedZones = false)
  {
    namespace bputils = axom::mir::utilities::blueprint;

    // Create the data
    conduit::Node hostMesh, deviceMesh;
    initialize(hostMesh);
    {
      AXOM_ANNOTATE_SCOPE("host_to_device");
      axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
    }
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    {
      AXOM_ANNOTATE_SCOPE("save_original");
  #if defined(AXOM_USE_HDF5)
      conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  #endif
      conduit::relay::io::save(hostMesh, name + "_orig.yaml", "yaml");
    }
#endif

    // Make views.
    auto coordsetView = axom::mir::views::make_explicit_coordset<float, 3>::view(
      deviceMesh["coordsets/coords"]);

    auto topologyView =
      axom::mir::views::make_structured<3>::view(deviceMesh["topologies/mesh"]);
    using CoordsetView = decltype(coordsetView);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView =
      axom::mir::views::make_unibuffer_matset<int, float, 3>::view(
        deviceMesh["matsets/mat"]);
    using MatsetView = decltype(matsetView);

    // Do MIR
    using MIR =
      axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node deviceMIRMesh;
    conduit::Node options;
    options["matset"] = "mat";
    if(selectedZones)
    {
      selectZones(options);
    }
    m.execute(deviceMesh, options, deviceMIRMesh);

    // device->host
    conduit::Node hostMIRMesh;
    {
      AXOM_ANNOTATE_SCOPE("device_to_host");
      axom::mir::utilities::blueprint::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);
    }

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
    {
      AXOM_ANNOTATE_SCOPE("save_baseline");
  #if defined(AXOM_USE_HDF5)
      conduit::relay::io::blueprint::save_mesh(hostMIRMesh, name, "hdf5");
  #endif
      conduit::relay::io::save(hostMIRMesh, name + ".yaml", "yaml");
    }
#endif
    // Handle baseline comparison.
    {
      AXOM_ANNOTATE_SCOPE("compare_baseline");
      std::string baselineName(yamlRoot(name));
      const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
      saveBaseline(paths, baselineName, hostMIRMesh);
#else
      constexpr double tolerance = 2.6e-06;
      EXPECT_TRUE(compareBaseline(paths, baselineName, hostMIRMesh, tolerance));
#endif
    }
  }
};

//------------------------------------------------------------------------------
TEST(mir_elvira3d, elvira3d_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_seq");
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer");
}
#if 0
TEST(mir_elvira3d, elvira3d_unibuffer_sel_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_seq");
  constexpr bool selectZones = true;
  test_Elvira3D<seq_exec>::test("elvira3d_unibuffer_sel",
                                selectZones);
}

  #if defined(AXOM_USE_OPENMP)
TEST(mir_elvira3d, elvira3d_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_omp");
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_omp)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_omp");
  constexpr bool selectZones = true;
  test_Elvira3D<omp_exec>::test("elvira3d_unibuffer_sel",
                                selectZones);
}
  #endif

  #if defined(AXOM_USE_CUDA)
TEST(mir_elvira3d, elvira3d_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_cuda");
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_cuda");
  constexpr bool selectZones = true;
  test_Elvira3D<cuda_exec>::test("elvira3d_unibuffer_sel",
                                 selectZones);
}
  #endif

  #if defined(AXOM_USE_HIP)
TEST(mir_elvira3d, elvira3d_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_hip");
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer");
}

TEST(mir_elvira3d, elvira3d_unibuffer_sel_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira3d_unibuffer_sel_hip");
  constexpr bool selectZones = true;
  test_Elvira3D<hip_exec>::test("elvira3d_unibuffer_sel",
                                selectZones);
}
  #endif
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
