// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <conduit/conduit_relay_io_blueprint.hpp>
#include <cmath>
#include <cstdlib>

namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------

// Uncomment to generate baselines
//#define AXOM_TESTING_GENERATE_BASELINES

// Uncomment to save visualization files for debugging (when making baselines)
//#define AXOM_TESTING_SAVE_VISUALIZATION

#include "axom/mir/tests/mir_testing_helpers.hpp"

std::string baselineDirectory()
{
  return pjoin(pjoin(pjoin(dataDirectory(), "mir"), "regression"),
               "mir_coupled");
}

//------------------------------------------------------------------------------
/*!

coarse

%42-------%43-------%44-------%45-------%46-------%47-------%48
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|   (3)   |   (3)   |     (3) |    (3)  |    (3)  |   (3)   |
%35-------%36-------%37#######%38#######%39#######%40-------%41
|         |         #\        #         #         #         |
|         |         #  \  (3) #         #         #         |
|         |         #    \    #         #         #         |
|         |         #      \  #         #         #         |
|   (2)   |   (2)   #  (2)   \#    (3)  #    (3)  #   (3)   |
%28-------%29------%30#######%31#######%32#######%33-------%34
|         |         #\ (2)    |\   (3)  |    (3)  #   (3)   |
|         |         #  \      |  \      |         #         |
|         |         #    \    |    - - -|- - - - -#- - - - -|
|         |         #      \  |         |         #         |
|   (1)   |   (1)   #  (1)   \|    (2)  |    (2)  #   (2)   |
%21-------%22-------%23#######%24#######%25#######%26-------%27
|         |   (0)   #\ (1)    |\   (2)  |    (2)  #   (2)   |
|         |         #  \      |(1)\     |         #         |
|         |         #    \- - |- - - - -|- - - - -#- - - - -|
|         |         #         |   (0)   |    (0)  #   (0)   |
|   (0)   |         #  (0)    |         |         #         |
%14-------%15-------%16#######%17#######%18#######%19-------%20
|   (0)   |   (0)   |  (0)    |   (0)   |    (0)  |   (0)   |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
%7--------%8--------%9--------%10-------%11-------%12-------%13
|   (0)   |   (0)   |  (0)    |   (0)   |    (0)  |   (0)   |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
|         |         |         |         |         |         |
%---------%---------%---------%---------%---------%---------%
0         1         2         3         4         5         6


fine - refines coarse with equal sized quads.

If fine has 2x2 refinement, it looks like this:

%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
|    |    |    |    |    |    |    |    |    |    |    |    |
|    |    |    |    |    |    |    |    |    |    |    |    |
%----%----%----%----%----%----%----%----%----%----%----%----%
0    1    2    3    4    5    6    7    8    9    10   11   12

*/
const char *yaml = R"(
coordsets:
  coarse_coords:
    type: uniform
    dims:
      i: 7
      j: 7
    origin:
      x: 0
      y: 0
    spacing:
      x: 1
      y: 1
  fine_coords:
    type: uniform
    dims:
      i: 13
      j: 13
    origin:
      x: 0
      y: 0
    spacing:
      x: 0.5
      y: 0.5
topologies:
  coarse_strided:
    type: structured
    coordset: coarse_coords
    elements:
      dims:
       i: 3
       j: 3
       offsets: [2, 2]
       strides: [1, 7]
  coarse:
    type: structured
    coordset: coarse_coords
    elements:
      dims:
       i: 6
       j: 6
  fine:
    type: structured
    coordset: fine_coords
    elements:
      dims:
       i: 12
       j: 12
matsets:
  coarse_matset:
   topology: coarse
   material_map:
     mat0: 0
     mat1: 1
     mat2: 2
     mat3: 3
   material_ids: [0, 0, 0, 0, 0, 0,               0, 0, 0, 0, 0, 0,           0, 0, 0,1, 0,1,2, 0,2, 0,2,                                 1, 1, 1,2, 2,3, 2,3, 2,3,                           2, 2, 2,3, 3, 3, 3,             3, 3, 3, 3, 3, 3]
   volume_fractions: [1., 1., 1., 1., 1., 1.,     1., 1., 1., 1., 1., 1.,     1., 1., 0.625,0.375, 0.5,0.125,0.375, 0.5,0.5, 0.5,0.5,     1., 1., 0.5,0.5, 0.625,0.375, 0.5,0.5, 0.5,0.5,     1., 1., 0.5,0.5, 1., 1., 1.,    1., 1., 1., 1., 1., 1.]
   indices: [0, 1, 2, 3, 4, 5,                    6, 7, 8, 9, 10, 11,         12, 13, 14,15, 16,17,18, 19,20, 21,22,                      23, 24, 25,26, 27,28, 29,30, 31,32,                 33, 34, 35,36, 37, 38, 39,      40, 41, 42, 43, 44, 45]
   sizes: [1, 1, 1, 1, 1, 1,     1, 1, 1, 1, 1, 1,     1, 1, 2, 3, 2, 2,    1, 1, 2, 2, 2, 2,    1, 1, 2, 1, 1, 1,    1, 1, 1, 1, 1, 1]
   offsets: [0, 1, 2, 3, 4, 5,     6, 7, 8, 9, 10, 11,     12, 13, 14, 16, 19, 21,     23, 24, 25, 27, 29, 31,     33, 34, 35, 37, 38, 39,     40, 41, 42, 43, 44, 45]
)";


//------------------------------------------------------------------------------
template <typename ExecSpace>
class test_coupling
{
public:
  static void test2D()
  {
    // Make the 2D input mesh.
    conduit::Node n_mesh;
    initialize(n_mesh);

    // host->device
    conduit::Node n_dev;
    axom::mir::utilities::blueprint::copy<ExecSpace>(n_dev, n_mesh);

    conduit::Node n_mir, n_map;
    mir2D(n_dev, n_mir);
    mapping2D(n_mir, n_map);

    // device->host
    conduit::Node hostResult;
    bputils::copy<seq_exec>(hostResult, n_map);

#if defined(AXOM_TESTING_SAVE_VISUALIZATION)
#if defined(AXOM_USE_HDF5)
    conduit::relay::io::blueprint::save_mesh(hostResult, "test2D", "hdf5");
#endif
    conduit::relay::io::save(hostResult, "test2D.yaml", "yaml");
#endif

    // Handle baseline comparison.
    const auto paths = baselinePaths<ExecSpace>();
    std::string baselineName(yamlRoot("test2D"));
#if defined(AXOM_TESTING_GENERATE_BASELINES)
    saveBaseline(paths, baselineName, hostResult);
#else
    EXPECT_TRUE(compareBaseline(paths, baselineName, hostResult));
#endif
  }

private:
  static constexpr int refinement = 2;

  static void initialize(conduit::Node &n_mesh)
  {
    // Make the 2D input mesh.
    n_mesh.parse(yaml);
  }

  static void mir2D(conduit::Node &n_input, conduit::Node &n_output)
  {
    const conduit::Node &n_coordset = n_input["coordsets/coarse_coords"];
    const conduit::Node &n_topology = n_input["topologies/coarse"];
    const conduit::Node &n_matset = n_input["matsets/coarse_matset"];

    // Wrap the coarse mesh in views.
    auto coordsetView = axom::mir::views::make_uniform_coordset<2>::view(n_coordset);
    using CoordsetView = decltype(coordsetView);

//    auto topologyView = axom::mir::views::make_strided_structured<2>::view(n_topology);
    auto topologyView = axom::mir::views::make_structured<2>::view(n_topology);
    using TopologyView = decltype(topologyView);
    using IndexingPolicy = typename TopologyView::IndexingPolicy;

    auto matsetView = axom::mir::views::make_unibuffer_matset<std::int64_t, double, 3>::view(n_matset);
    using MatsetView = decltype(matsetView);

    // Do MIR on the mesh.
    using MIR = axom::mir::ElviraAlgorithm<ExecSpace, IndexingPolicy, CoordsetView, MatsetView>;
    MIR m(topologyView, coordsetView, matsetView);
    conduit::Node options;
    options["matset"] = "coarse_matset";
    options["topologyName"] = "postmir";
    options["coordsetName"] = "postmir_coords";
    options["matsetName"] = "postmir_matset";
    m.execute(n_input, options, n_output);

    std::cout << "------ MIR output ------\n";
    printNode(n_output);
  }

  static void mapping2D(conduit::Node &n_input, conduit::Node &n_output)
  {
#if 0
    const conduit::Node &n_coordset = n_input["coordsets/coarse_coordset"];
    const conduit::Node &n_topology = n_input["topologies/coarse"];
    const conduit::Node &n_matset = n_input["matsets/coarse_matset"];

    // Wrap the coarse mesh in views.
    make_explicit_coordset<DataType, 3>::view(n_coordset);

    // Wrap coarse/post_mir mesh in views.
    using SrcCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using SrcTopologyView =
      axom::mir::views::UnstructuredTopologyMixedShapeView<conduit::index_t>;
    SrcCoordsetView srcCoordset(bputils::make_array_view<double>(
                                  n_dev["coordsets/postmir_coords/values/x"]),
                                bputils::make_array_view<double>(
                                  n_dev["coordsets/postmir_coords/values/y"]));

    axom::Array<axom::IndexType> shapeValues, shapeIds;
    const conduit::Node &n_srcTopo = n_dev["topologies/postmir"];
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();
    auto shapeMap = axom::mir::views::buildShapeMap(n_srcTopo,
                                                    shapeValues,
                                                    shapeIds,
                                                    allocatorID);
    SrcTopologyView srcTopo(
      bputils::make_array_view<conduit::index_t>(
        n_srcTopo["elements/connectivity"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/shapes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/sizes"]),
      bputils::make_array_view<conduit::index_t>(n_srcTopo["elements/offsets"]),
      shapeMap);

    // Wrap fine mesh in views.
    using TargetCoordsetView = axom::mir::views::ExplicitCoordsetView<double, 2>;
    using TargetShapeType = axom::mir::views::QuadShape<int>;
    using TargetTopologyView =
      axom::mir::views::UnstructuredTopologySingleShapeView<TargetShapeType>;

    TargetCoordsetView targetCoordset(
      bputils::make_array_view<double>(n_dev["coordsets/fine_coords/values/x"]),
      bputils::make_array_view<double>(
        n_dev["coordsets/fine_coords/values/y"]));
    const conduit::Node &n_targetTopo = n_dev["topologies/fine"];
    TargetTopologyView targetTopo(
      bputils::make_array_view<int>(n_targetTopo["elements/connectivity"]),
      bputils::make_array_view<int>(n_targetTopo["elements/sizes"]),
      bputils::make_array_view<int>(n_targetTopo["elements/offsets"]));

    // Make new VFs via mapper.
    using Mapper = bputils::TopologyMapper<ExecSpace,
                                           SrcTopologyView,
                                           SrcCoordsetView,
                                           TargetTopologyView,
                                           TargetCoordsetView>;
    Mapper mapper(srcTopo, srcCoordset, targetTopo, targetCoordset);
    conduit::Node n_opts;
    n_opts["source/matsetName"] = "postmir_matset";
    n_opts["target/topologyName"] = "fine";
    n_opts["target/matsetName"] = "fine_matset";
    mapper.execute(n_dev, n_opts, n_dev);
#endif
  }
};

//------------------------------------------------------------------------------
TEST(mir_coupling, coupling_2D_seq)
{
  test_coupling<seq_exec>::test2D();
}
/*#if defined(AXOM_USE_OPENMP)
TEST(mir_coupling, coupling_2D_omp)
{
  test_coupling<omp_exec>::test2D();
}
#endif
#if defined(AXOM_USE_CUDA)
TEST(mir_coupling, coupling_2D_cuda)
{
  test_coupling<cuda_exec>::test2D();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_coupling, coupling_2D_hip)
{
  test_coupling<hip_exec>::test2D();
}
#endif
*/

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,
#if defined(AXOM_USE_CALIPER)
  axom::CLI::App app;
  std::string annotationMode("none");
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);

  // Parse command line options.
  app.parse(argc, argv);

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    annotationMode);
#endif

  result = RUN_ALL_TESTS();
  return result;
}
