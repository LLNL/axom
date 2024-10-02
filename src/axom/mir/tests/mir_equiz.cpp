// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/primal.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

//------------------------------------------------------------------------------

//#define DEBUGGING_TEST_CASES

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "mir"), "regression"), "mir_equiz");
}

//------------------------------------------------------------------------------
TEST(mir_equiz, miralgorithm)
{
  axom::mir::MIRAlgorithm *m = nullptr;
  EXPECT_EQ(m, nullptr);
}

//------------------------------------------------------------------------------
TEST(mir_equiz, materialinformation)
{
  conduit::Node matset;
  matset["material_map/a"] = 1;
  matset["material_map/b"] = 2;
  matset["material_map/c"] = 0;

  auto mi = axom::mir::views::materials(matset);
  EXPECT_EQ(mi.size(), 3);
  EXPECT_EQ(mi[0].number, 1);
  EXPECT_EQ(mi[0].name, "a");

  EXPECT_EQ(mi[1].number, 2);
  EXPECT_EQ(mi[1].name, "b");

  EXPECT_EQ(mi[2].number, 0);
  EXPECT_EQ(mi[2].name, "c");
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid2d_mat_test(const std::string &type,
                      const std::string &mattype,
                      const std::string &name)
{
  namespace bputils = axom::mir::utilities::blueprint;

  axom::StackArray<axom::IndexType, 2> dims {10, 10};
  axom::StackArray<axom::IndexType, 2> zoneDims {dims[0] - 1, dims[1] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  axom::mir::testing::data::make_matset(mattype, "mesh", zoneDims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Make views.
  auto coordsetView = axom::mir::views::make_uniform_coordset<2>::view(
    deviceMesh["coordsets/coords"]);
  auto topologyView =
    axom::mir::views::make_uniform<2>::view(deviceMesh["topologies/mesh"]);
  using CoordsetView = decltype(coordsetView);
  using TopologyView = decltype(topologyView);

  conduit::Node deviceMIRMesh;
  if(mattype == "unibuffer")
  {
    // clang-format off
    using MatsetView = axom::mir::views::UnibufferMaterialView<int, float, 3>;
    MatsetView matsetView;
    matsetView.set(bputils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                   bputils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));
    // clang-format on

    using MIR =
      axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "mat";
    m.execute(deviceMesh, options, deviceMIRMesh);
  }

  // device->host
  conduit::Node hostMIRMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMIRMesh, name, "hdf5");
#endif
  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostMIRMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostMIRMesh));
#endif
  }
}

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_seq");
  braid2d_mat_test<seq_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_uniform_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_omp");
  braid2d_mat_test<omp_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_cuda");
  braid2d_mat_test<cuda_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_uniform_unibuffer_hip");
  braid2d_mat_test<hip_exec>("uniform", "unibuffer", "equiz_uniform_unibuffer");
}
#endif

#if 0 // FIXME
//------------------------------------------------------------------------------
template <typename ExecSpace>
void braid3d_mat_test(const std::string &type,
                      const std::string &mattype,
                      const std::string &name)
{
  namespace bputils = axom::mir::utilities::blueprint;

  axom::StackArray<axom::IndexType, 3> dims {10, 10, 10};
  axom::StackArray<axom::IndexType, 3> zoneDims {dims[0] - 1, dims[1] - 1, dims[2] - 1};

  // Create the data
  conduit::Node hostMesh, deviceMesh;
  axom::mir::testing::data::braid(type, dims, hostMesh);
  axom::mir::testing::data::make_matset(mattype, "mesh", zoneDims, hostMesh);
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
#endif

  // Make views.
  auto coordsetView = axom::mir::views::make_explicit_coordset<float, 3>::view(
    deviceMesh["coordsets/coords"]);
  using CoordsetView = decltype(coordsetView);

  using ShapeType = axom::mir::views::HexShape<int>;
  using TopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<ShapeType>;
  TopologyView topologyView(deviceMesh["topologies/mesh"]);

  conduit::Node deviceMIRMesh;
  if(mattype == "unibuffer")
  {
    // clang-format off
    using MatsetView = axom::mir::views::UnibufferMaterialView<int, float, 3>;
    MatsetView matsetView;
    matsetView.set(bputils::make_array_view<int>(deviceMesh["matsets/mat/material_ids"]),
                   bputils::make_array_view<float>(deviceMesh["matsets/mat/volume_fractions"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/sizes"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/offsets"]),
                   bputils::make_array_view<int>(deviceMesh["matsets/mat/indices"]));
    // clang-format on

    using MIR =
      axom::mir::EquiZAlgorithm<ExecSpace, TopologyView, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "mat";
    m.execute(deviceMesh, options, deviceMIRMesh);
  }

  // device->host
  conduit::Node hostMIRMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
  conduit::relay::io::blueprint::save_mesh(hostMIRMesh, name, "hdf5");
#endif
  // Handle baseline comparison.
  {
    std::string baselineName(yamlRoot(name));
    const auto paths = baselinePaths<ExecSpace>();
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostMIRMesh);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostMIRMesh));
#endif
  }
}

//------------------------------------------------------------------------------
TEST(mir_equiz, equiz_hex_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("equiz_explicit_hex_seq");
  braid3d_mat_test<seq_exec>("hex", "unibuffer", "equiz_hex_unibuffer");
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_equiz, equiz_hex_unibuffer_omp)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_omp");
  braid3d_mat_test<omp_exec>("hex", "unibuffer", "equiz_hex_unibuffer");
}
#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_equiz, equiz_hex_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_cuda");
  braid3d_mat_test<cuda_exec>("hex", "unibuffer", "equiz_hex_unibuffer");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_equiz, equiz_hex_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("equiz_hex_unibuffer_hip");
  braid3d_mat_test<hip_exec>("hex", "unibuffer", "equiz_hex_unibuffer");
}
#endif
#endif

//------------------------------------------------------------------------------
#if defined(DEBUGGING_TEST_CASES)
void conduit_debug_err_handler(const std::string &s1, const std::string &s2, int i1)
{
  std::cout << "s1=" << s1 << ", s2=" << s2 << ", i1=" << i1 << std::endl;
  // This is on purpose.
  while(1)
    ;
}
#endif
//------------------------------------------------------------------------------

int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  // Define command line options.
  bool handler = true;
  axom::CLI::App app;
  app.add_option("--handler", handler)
    ->description("Install a custom error handler that loops forever.");
#if defined(AXOM_USE_CALIPER)
  std::string annotationMode("report");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif
  // Parse command line options.
  app.parse(argc, argv);

#if defined(AXOM_USE_CALIPER)
  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    annotationMode);
#endif

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
#if defined(DEBUGGING_TEST_CASES)
  if(handler)
  {
    conduit::utils::set_error_handler(conduit_debug_err_handler);
  }
#endif
  result = RUN_ALL_TESTS();
  return result;
}
