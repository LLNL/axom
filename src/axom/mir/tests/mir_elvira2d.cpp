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
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "mir"), "regression"), "mir_elvira");
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
#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMesh, name + "_orig", "hdf5");
  conduit::relay::io::save(hostMesh, name + "_orig.yaml", "yaml");
#endif

  // _elvira_mir_start
  // Make views.
  auto coordsetView = axom::mir::views::make_uniform_coordset<2>::view(
    deviceMesh["coordsets/coords"]);
  auto topologyView =
    axom::mir::views::make_uniform<2>::view(deviceMesh["topologies/mesh"]);
  using CoordsetView = decltype(coordsetView);
  using TopologyView = decltype(topologyView);
  using IndexingPolicy = typename TopologyView::IndexingPolicy;

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
      axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "mat";
    m.execute(deviceMesh, options, deviceMIRMesh);
  }
  // _elvira_mir_end

  // device->host
  conduit::Node hostMIRMesh;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostMIRMesh, deviceMIRMesh);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION) && defined(AXOM_USE_HDF5)
  conduit::relay::io::blueprint::save_mesh(hostMIRMesh, name, "hdf5");
#endif
  // Handle baseline comparison.
  {
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

template <typename ExecSpace>
void test_last_clip()
{
  constexpr int MAX_VERTS = 6;
  using value_type = double;
  using PointType = axom::primal::Point<value_type, 2>;
  using VectorType = axom::primal::Vector<value_type, 2>;
  using PlaneType = axom::primal::Plane<value_type, 2>;
  using PolygonType =
    axom::primal::Polygon<value_type, 2, axom::primal::PolygonArray::Static, MAX_VERTS>;

  // Make the zone to clip.
  PolygonType shape;
  shape.addVertex(PointType({-7.7777777777, -5.5555555555}));
  shape.addVertex(PointType({10., -5.5555555555}));
  shape.addVertex(PointType({10., -3.3333333333}));
  shape.addVertex(PointType({-7.7777777777, -3.3333333333}));

  // Make return data.
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
  axom::Array<PolygonType> polygons(2, 2, allocatorID);
  auto polygonsView = polygons.view();

  // Clip the polygons in the kernel
  axom::for_all<ExecSpace>(
    1,
    AXOM_LAMBDA(axom::IndexType index) {
      const VectorType normal({0.371391, 0.928477});
      const PointType pt({9.24395, -3.55679});

      // Clip one side.
      auto P = PlaneType(-normal, pt, false);
      polygonsView[0] = axom::primal::clip(shape, P);

      // Clip the other side. This is the one that I was having problems with on CUDA.
      P = PlaneType(normal, pt, false);
      polygonsView[1] = axom::primal::clip(shape, P);
    });

  // Get polygons back on host.
  axom::Array<PolygonType> hostPolygons;
  hostPolygons.resize(2);
  axom::copy(hostPolygons.data(), polygons.data(), 2 * sizeof(PolygonType));

  for(int i = 0; i < 2; i++)
  {
    std::cout << "hostPolygons[" << i << "]=" << hostPolygons[i] << std::endl;
  }
}

//------------------------------------------------------------------------------
TEST(mir_elvira, elvira_uniform_unibuffer_seq)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_seq");
  braid2d_mat_test<seq_exec>("uniform", "unibuffer", "elvira_uniform_unibuffer");
}

//#if defined(AXOM_USE_OPENMP)
//TEST(mir_elvira, elvira_uniform_unibuffer_omp)
//{
//  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_omp");
//  braid2d_mat_test<omp_exec>("uniform", "unibuffer", "elvira_uniform_unibuffer");
//}
//#endif

#if defined(AXOM_USE_CUDA)
TEST(mir_elvira, elvira_uniform_unibuffer_cuda)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_cuda");
  //  test_last_clip<cuda_exec>();
  braid2d_mat_test<cuda_exec>("uniform",
                              "unibuffer",
                              "elvira_uniform_unibuffer");
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_elvira, elvira_uniform_unibuffer_hip)
{
  AXOM_ANNOTATE_SCOPE("elvira_uniform_unibuffer_hip");
  braid2d_mat_test<hip_exec>("uniform", "unibuffer", "elvira_uniform_unibuffer");
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
  return result;
}
