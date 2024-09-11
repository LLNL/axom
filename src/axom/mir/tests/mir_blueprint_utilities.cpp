// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/mir.hpp"
//#include "axom/mir/blueprint_utilities.hpp"

#include "axom/mir/tests/mir_testing_helpers.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <iostream>
#include <algorithm>

namespace mir = axom::mir;
namespace bputils = axom::mir::utilities::blueprint;

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_conduit_allocate
{
  static void test()
  {
    axom::mir::utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
    EXPECT_TRUE(c2a.getConduitAllocatorID() > 0);

    constexpr int nValues = 100;
    conduit::Node n;
    n.set_allocator(c2a.getConduitAllocatorID());
    n.set(conduit::DataType::int32(nValues));

    // Make sure we can store some values into the data that were allocated.
    auto nview = bputils::make_array_view<int>(n);
    axom::for_all<ExecSpace>(
      nValues,
      AXOM_LAMBDA(auto index) { nview[index] = index; });

    EXPECT_EQ(n.dtype().number_of_elements(), nValues);

    // Get the values back to the host.
    std::vector<int> hostValues(nValues);
    axom::copy(hostValues.data(), n.data_ptr(), sizeof(int) * nValues);

    // Check that the values were set.
    for(int i = 0; i < nValues; i++)
    {
      EXPECT_EQ(hostValues[i], i);
    }
  }
};

TEST(mir_blueprint_utilities, allocate_seq)
{
  test_conduit_allocate<seq_exec>::test();
}
#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, allocate_omp)
{
  test_conduit_allocate<omp_exec>::test();
}
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, allocate_cuda)
{
  test_conduit_allocate<cuda_exec>::test();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, allocate_hip)
{
  test_conduit_allocate<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_copy_braid
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // Copy the mesh to device.
    conduit::Node deviceMesh;
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);

    // Run some minmax operations on device (proves that the data was in the right place) and check the results.

    constexpr double eps = 1.e-7;

    auto x = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
      deviceMesh["coordsets/coords/values/x"]);
    //std::cout << std::setw(16) << "x={" << x.first << ", " << x.second << "}\n";
    EXPECT_NEAR(x.first, -10., eps);
    EXPECT_NEAR(x.second, 10., eps);

    auto y = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
      deviceMesh["coordsets/coords/values/y"]);
    //std::cout << std::setw(16) << "y={" << y.first << ", " << y.second << "}\n";
    EXPECT_NEAR(y.first, -10., eps);
    EXPECT_NEAR(y.second, 10., eps);

    auto c = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
      deviceMesh["topologies/mesh/elements/connectivity"]);
    //std::cout << std::setw(16) << "conn={" << c.first << ", " << c.second << "}\n";
    EXPECT_NEAR(c.first, 0., eps);
    EXPECT_NEAR(c.second, 999., eps);

    auto r = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
      deviceMesh["fields/radial/values"]);
    //std::cout << std::setw(16) << "radial={" << r.first << ", " << r.second << "}\n";
    EXPECT_NEAR(r.first, 19.2450089729875, eps);
    EXPECT_NEAR(r.second, 173.205080756888, eps);
  }

  static void create(conduit::Node &mesh)
  {
    const int d[3] = {10, 10, 10};
    conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);
  }
};

TEST(mir_blueprint_utilities, copy_seq)
{
  test_copy_braid<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, copy_omp)
{
  test_copy_braid<omp_exec>::test();
}
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, copy_cuda)
{
  test_copy_braid<cuda_exec>::test();
}
#endif
#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, copy_hip)
{
  test_copy_braid<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
TEST(mir_blueprint_utilities, to_unstructured)
{
  // TODO: to_unstructured
}

//------------------------------------------------------------------------------

template <typename IndexT = int>
void compareRelation(const conduit::Node &hostRelation,
                     const axom::ArrayView<const int> &zones,
                     const axom::ArrayView<const int> &sizes,
                     const axom::ArrayView<const int> &offsets)
{
  const auto zonesView =
    bputils::make_array_view<IndexT>(hostRelation["zones"]);
  const auto sizesView =
    bputils::make_array_view<IndexT>(hostRelation["sizes"]);
  const auto offsetsView =
    bputils::make_array_view<IndexT>(hostRelation["offsets"]);
  EXPECT_EQ(sizesView.size(), sizes.size());
  EXPECT_EQ(offsetsView.size(), offsets.size());
  for(axom::IndexType i = 0; i < sizesView.size(); i++)
  {
    EXPECT_EQ(sizes[i], sizesView[i]);
    EXPECT_EQ(offsets[i], offsetsView[i]);
  }
  for(axom::IndexType i = 0; i < sizesView.size(); i++)
  {
    // Sort the result so we can compare to the expected answer.
    IndexT *begin = zonesView.data() + offsetsView[i];
    IndexT *end = zonesView.data() + offsetsView[i] + sizesView[i];
    std::sort(begin, end);

    for(int j = 0; j < sizesView[i]; j++)
    {
      EXPECT_EQ(zones[offsets[i] + j], zonesView[offsetsView[i] + j]);
    }
  }
}

template <typename ExecSpace, typename IndexT = int>
void test_node_to_zone_relation_builder(const conduit::Node &hostMesh)
{
  // host -> device
  conduit::Node deviceMesh;
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
  const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
  const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

  // Run the algorithm on the device
  conduit::Node deviceRelation;
  axom::mir::utilities::blueprint::NodeToZoneRelationBuilder<ExecSpace> n2z;
  n2z.execute(deviceTopo, deviceCoordset, deviceRelation);

  // device -> host
  conduit::Node hostRelation;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostRelation, deviceRelation);
#if 0
  // Print the results.
  printNode(hostRelation);
#endif
  // Expected answers
  // clang-format off
  const int zones[] = {
   0,
   0, 1,
   1, 2,
   2,
   0, 3,        // NOTE: these are sorted here
   0, 1, 3, 4,
   1, 2, 4, 5,
   2, 5,
   3,
   3, 4,
   4, 5,
   5
  };
  const int sizes[] = {1, 2, 2, 1, 2, 4, 4, 2, 1, 2, 2, 1};
  const int offsets[] = {0, 1, 3, 5, 6, 8, 12, 16, 18, 19, 21, 23};
  // clang-format on

  // Compare answers.
  compareRelation<IndexT>(
    hostRelation,
    axom::ArrayView<const int>(zones, sizeof(zones) / sizeof(int)),
    axom::ArrayView<const int>(sizes, sizeof(sizes) / sizeof(int)),
    axom::ArrayView<const int>(offsets, sizeof(offsets) / sizeof(int)));
}

TEST(mir_blueprint_utilities, node_to_zone_relation_builder_unstructured)
{
  /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::mir::testing::data::braid("quads", dims, mesh);

  test_node_to_zone_relation_builder<seq_exec>(mesh);
#if defined(AXOM_USE_OPENMP)
  test_node_to_zone_relation_builder<omp_exec>(mesh);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_node_to_zone_relation_builder<cuda_exec>(mesh);
#endif

#if defined(AXOM_USE_HIP)
  test_node_to_zone_relation_builder<hip_exec>(mesh);
#endif
}

TEST(mir_blueprint_utilities, node_to_zone_relation_builder_rectilinear)
{
  /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims {{4, 3}};
  axom::mir::testing::data::braid("rectilinear", dims, mesh);
  //mesh.print();

  test_node_to_zone_relation_builder<seq_exec, conduit::index_t>(mesh);
#if defined(AXOM_USE_OPENMP)
  test_node_to_zone_relation_builder<omp_exec, conduit::index_t>(mesh);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_node_to_zone_relation_builder<cuda_exec, conduit::index_t>(mesh);
#endif

#if defined(AXOM_USE_HIP)
  test_node_to_zone_relation_builder<hip_exec, conduit::index_t>(mesh);
#endif
}

template <typename ExecSpace, typename IndexT = int>
void test_node_to_zone_relation_builder_polyhedral(const conduit::Node &hostMesh)
{
  // host -> device
  conduit::Node deviceMesh;
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
  const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
  const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

  // Run the algorithm on the device
  conduit::Node deviceRelation;
  axom::mir::utilities::blueprint::NodeToZoneRelationBuilder<ExecSpace> n2z;
  n2z.execute(deviceTopo, deviceCoordset, deviceRelation);

  // device -> host
  conduit::Node hostRelation;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostRelation, deviceRelation);
#if 0
  // Print the results.
  printNode(hostRelation);
#endif
  // Expected answers
  // clang-format off
  const int zones[] = {
   /*node 0*/ 0,
   /*node 1*/ 0, 1,
   /*node 2*/ 1,
   /*node 3*/ 0, 2,
   /*node 4*/ 0, 1, 2, 3,
   /*node 5*/ 1, 3,
   /*node 6*/ 2,
   /*node 7*/ 2, 3,
   /*node 8*/ 3,

   /*node 9*/ 0, 4,
   /*node 10*/ 0, 1, 4, 5,
   /*node 11*/ 1, 5,
   /*node 12*/ 0, 2, 4, 6,
   /*node 13*/ 0, 1, 2, 3, 4, 5, 6, 7,
   /*node 14*/ 1, 3, 5, 7,
   /*node 15*/ 2, 6,
   /*node 16*/ 2, 3, 6, 7,
   /*node 17*/ 3, 7,

   /*node 18*/ 4,
   /*node 19*/ 4, 5,
   /*node 20*/ 5,
   /*node 21*/ 4, 6,
   /*node 22*/ 4, 5, 6, 7,
   /*node 23*/ 5, 7,
   /*node 24*/ 6,
   /*node 25*/ 6, 7,
   /*node 26*/ 7
  };
  const int sizes[] = {1, 2, 1, 2, 4, 2, 1, 2, 1,
                       2, 4, 2, 4, 8, 4, 2, 4, 2,
                       1, 2, 1, 2, 4, 2, 1, 2, 1
                      };
  const int offsets[] = {0, 1, 3, 4, 6, 10, 12, 13, 15,
                         16, 18, 22, 24, 28, 36, 40, 42, 46,
                         48, 49, 51, 52, 54, 58, 60, 61, 63};
  // clang-format on

  // Compare answers.
  compareRelation<IndexT>(
    hostRelation,
    axom::ArrayView<const int>(zones, sizeof(zones) / sizeof(int)),
    axom::ArrayView<const int>(sizes, sizeof(sizes) / sizeof(int)),
    axom::ArrayView<const int>(offsets, sizeof(offsets) / sizeof(int)));
}

TEST(mir_blueprint_utilities, node_to_zone_relation_builder_polyhedral)
{
  conduit::Node mesh;
  conduit::blueprint::mesh::examples::basic("polyhedra", 3, 3, 3, mesh);
  // Make sure all the types are the same.
  conduit::blueprint::mesh::utils::convert(
    mesh,
    conduit::DataType::int32(),
    std::vector<std::string> {{"topologies/mesh/elements/connectivity",
                               "topologies/mesh/elements/sizes",
                               "topologies/mesh/elements/offsets",
                               "topologies/mesh/subelements/connectivity",
                               "topologies/mesh/subelements/sizes",
                               "topologies/mesh/subelements/offsets"}});
  //printNode(mesh);

  test_node_to_zone_relation_builder_polyhedral<seq_exec, conduit::int32>(mesh);

#if defined(AXOM_USE_OPENMP)
  test_node_to_zone_relation_builder_polyhedral<omp_exec, conduit::int32>(mesh);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_node_to_zone_relation_builder_polyhedral<cuda_exec, conduit::int32>(mesh);
#endif

#if defined(AXOM_USE_HIP)
  test_node_to_zone_relation_builder_polyhedral<hip_exec, conduit::int32>(mesh);
#endif
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_recenter_field
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host -> device
    conduit::Node deviceMesh;
    axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
    const conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];
    const conduit::Node &deviceCoordset = deviceMesh["coordsets/coords"];

    // Make a node to zone relation on the device.
    conduit::Node deviceRelation;
    axom::mir::utilities::blueprint::NodeToZoneRelationBuilder<ExecSpace> n2z;
    n2z.execute(deviceTopo, deviceCoordset, deviceRelation);

    // Recenter a field zonal->nodal on the device
    axom::mir::utilities::blueprint::RecenterField<ExecSpace> r;
    r.execute(deviceMesh["fields/easy_zonal"],
              deviceRelation,
              deviceMesh["fields/z2n"]);

    // Recenter a field nodal->zonal on the device. (The elements are an o2m relation)
    r.execute(deviceMesh["fields/z2n"],
              deviceMesh["topologies/mesh/elements"],
              deviceMesh["fields/n2z"]);

    // device -> host
    conduit::Node hostResultMesh;
    axom::mir::utilities::blueprint::copy<seq_exec>(hostResultMesh, deviceMesh);
#if 0
    // Print the results.
    printNode(hostResultMesh);
#endif
    const float n2z_result[] = {1., 2., 4., 5., 4., 5., 7., 8., 7., 8., 10., 11.};
    for(size_t i = 0; i < (sizeof(n2z_result) / sizeof(float)); i++)
    {
      EXPECT_EQ(n2z_result[i],
                hostResultMesh["fields/z2n/values"].as_float_accessor()[i]);
    }
    const float z2n_result[] = {3., 4.5, 6., 6., 7.5, 9.};
    for(size_t i = 0; i < (sizeof(z2n_result) / sizeof(float)); i++)
    {
      EXPECT_EQ(z2n_result[i],
                hostResultMesh["fields/n2z/values"].as_float_accessor()[i]);
    }
  }

  static void create(conduit::Node &mesh)
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    axom::StackArray<int, 2> dims {{4, 3}};
    axom::mir::testing::data::braid("quads", dims, mesh);
    mesh["topologies/mesh/elements/sizes"].set(
      std::vector<int> {{4, 4, 4, 4, 4, 4}});
    mesh["topologies/mesh/elements/offsets"].set(
      std::vector<int> {{0, 4, 8, 12, 16, 20}});
    mesh["fields/easy_zonal/topology"] = "mesh";
    mesh["fields/easy_zonal/association"] = "element";
    mesh["fields/easy_zonal/values"].set(std::vector<float> {{1, 3, 5, 7, 9, 11}});
  }
};

TEST(mir_blueprint_utilities, recenterfield_seq)
{
  test_recenter_field<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, recenterfield_omp)
{
  test_recenter_field<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, recenterfield_cuda)
{
  test_recenter_field<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, recenterfield_hip)
{
  test_recenter_field<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename Container1, typename Container2>
bool compare_views(const Container1 &a, const Container2 &b)
{
  bool eq = a.size() == b.size();
  for(axom::IndexType i = 0; i < a.size() && eq; i++)
  {
    eq &= a[i] == b[i];
  }
  return eq;
} 

template <typename ExecSpace>
struct test_matset_slice
{
  static void test()
  {
    conduit::Node hostMatset;
    create(hostMatset);

    // host->device
    conduit::Node deviceMatset;
    bputils::copy<ExecSpace>(deviceMatset, hostMatset);

    axom::Array<int> ids{{1,3,5}};
    axom::Array<int> selectedZones(3, 3, axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), 3 * sizeof(int));

    using MatsetView = axom::mir::views::UnibufferMaterialView<conduit::int64, conduit::float64, 3>;
    MatsetView matsetView;
    matsetView.set(
      bputils::make_array_view<conduit::int64>(deviceMatset["material_ids"]),
      bputils::make_array_view<conduit::float64>(deviceMatset["volume_fractions"]),
      bputils::make_array_view<conduit::int64>(deviceMatset["sizes"]),
      bputils::make_array_view<conduit::int64>(deviceMatset["offsets"]),
      bputils::make_array_view<conduit::int64>(deviceMatset["indices"]));

    // Slice it.
    bputils::MatsetSlicer<ExecSpace, MatsetView> slicer;
    conduit::Node newDeviceMatset;
    slicer.execute(matsetView, selectedZones.view(), deviceMatset, newDeviceMatset);

    // device->host
    conduit::Node newHostMatset;
    bputils::copy<ExecSpace>(newHostMatset, newDeviceMatset);

    // Expected answers.
    const axom::Array<conduit::int64> sizes{{2, 1, 2}};
    const axom::Array<conduit::int64> offsets{{0, 2, 3}};
    const axom::Array<conduit::int64> indices{{0, 1, 2, 3, 4}};
    const axom::Array<conduit::int64> material_ids{{1, 2, 2, 2, 3}};
    const axom::Array<conduit::float64> volume_fractions{{0.5, 0.5, 1.0, 0.8, 0.2}};

    EXPECT_TRUE(compare_views(sizes.view(), bputils::make_array_view<conduit::int64>(newHostMatset["sizes"])));
    EXPECT_TRUE(compare_views(offsets.view(), bputils::make_array_view<conduit::int64>(newHostMatset["offsets"])));
    EXPECT_TRUE(compare_views(indices.view(), bputils::make_array_view<conduit::int64>(newHostMatset["indices"])));
    EXPECT_TRUE(compare_views(material_ids.view(), bputils::make_array_view<conduit::int64>(newHostMatset["material_ids"])));
    EXPECT_TRUE(compare_views(volume_fractions.view(), bputils::make_array_view<conduit::float64>(newHostMatset["volume_fractions"])));
  }

  static void create(conduit::Node &matset)
  {
    /*
      8-------9------10------11
      |  2/1  | 1/0.1 | 2/0.8 |
      |       | 2/0.5 | 3/0.2 |
      |       | 3/0.4 |       |
      4-------5-------6-------7
      |       | 1/0.5 | 1/0.2 |
      |  1/1  | 2/0.5 | 2/0.8 |
      |       |       |       |
      0-------1-------2-------3
      */
    const char *yaml = R"xx(
topology: mesh
material_map:
  a: 1
  b: 2
  c: 3
material_ids: [1, 1,2, 1,2, 2, 1,2,3, 2,3]
volume_fractions: [1., 0.5,0.5, 0.2,0.8, 1., 0.1,0.5,0.4, 0.8,0.2]
indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
sizes: [1, 2, 2, 1, 3, 2]
offsets: [0, 1, 3, 5, 6, 9]
)xx";

    matset.parse(yaml); 
  }
};

TEST(mir_blueprint_utilities, matsetslice_seq)
{
  test_matset_slice<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, matsetslice_omp)
{
  test_matset_slice<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, matsetslice_cuda)
{
  test_matset_slice<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, matsetslice_hip)
{
  test_matset_slice<hip_exec>::test();
}
#endif


//------------------------------------------------------------------------------
template <typename ExecSpace, typename Func>
void test_coordsetslicer(const conduit::Node &hostCoordset, Func &&makeView)
{
  axom::Array<axom::IndexType> ids{{0,1,2,4,5,6}};

  const auto nnodes = ids.size();
  axom::Array<axom::IndexType> selectedNodes(nnodes, nnodes, axom::execution_space<ExecSpace>::allocatorID());
  axom::copy(selectedNodes.data(), ids.data(), nnodes * sizeof(axom::IndexType));

  bputils::SliceData slice;
  slice.m_indicesView = selectedNodes.view();

  // host->device
  conduit::Node deviceCoordset;
  bputils::copy<ExecSpace>(deviceCoordset, hostCoordset);

  // Make a view.
  auto coordsetView = makeView(deviceCoordset);
  using CoordsetView = decltype(coordsetView);

  // Pull out selected nodes
  bputils::CoordsetSlicer<ExecSpace, CoordsetView> slicer(coordsetView);
  conduit::Node newDeviceCoordset;
  slicer.execute(slice, deviceCoordset, newDeviceCoordset);

  // device->host
  conduit::Node newHostCoordset;
  bputils::copy<ExecSpace>(newHostCoordset, newDeviceCoordset);

  // We get an explicit coordset out of the slicer.
  const axom::Array<conduit::float64> x{{0., 1., 2., 0., 1., 2.}};
  const axom::Array<conduit::float64> y{{0., 0., 0., 1., 1., 1.}};
  EXPECT_TRUE(compare_views(x.view(), bputils::make_array_view<conduit::float64>(newHostCoordset["values/x"])));
  EXPECT_TRUE(compare_views(y.view(), bputils::make_array_view<conduit::float64>(newHostCoordset["values/y"])));
}

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_explicit
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: explicit
values:
  x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
  y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.]
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset)
    {
      return axom::mir::views::make_explicit_coordset<conduit::float64, 2>::view(deviceCoordset);
    };

    test_coordsetslicer<seq_exec>(coordset, makeView);
  }
};

TEST(mir_blueprint_utilities, coordsetslicer_explicit_seq)
{
  coordsetslicer_explicit<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, coordsetslicer_explicit_omp)
{
  coordsetslicer_explicit<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, coordsetslicer_explicit_cuda)
{
  coordsetslicer_explicit<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, coordsetslicer_explicit_hip)
{
  coordsetslicer_explicit<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_rectilinear
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: rectilinear
values:
  x: [0., 1., 2., 3.]
  y: [0., 1., 2.]
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset)
    {
      return axom::mir::views::make_rectilinear_coordset<conduit::float64, 2>::view(deviceCoordset);
    };
    test_coordsetslicer<ExecSpace>(coordset, makeView);
  }
};

TEST(mir_blueprint_utilities, coordsetslicer_rectilinear_seq)
{
  coordsetslicer_rectilinear<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, coordsetslicer_rectilinear_omp)
{
  coordsetslicer_rectilinear<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, coordsetslicer_rectilinear_cuda)
{
  coordsetslicer_rectilinear<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, coordsetslicer_rectilinear_hip)
{
  coordsetslicer_rectilinear<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct coordsetslicer_uniform
{
  static void test()
  {
    /*
      8---9--10--11
      |   |   |   |
      4---5---6---7
      |   |   |   |
      0---1---2---3
      */
    const char *yaml = R"xx(
type: uniform
dims:
  i: 4
  j: 3
)xx";

    conduit::Node coordset;
    coordset.parse(yaml);

    auto makeView = [](const conduit::Node &deviceCoordset)
    {
      return axom::mir::views::make_uniform_coordset<2>::view(deviceCoordset);
    };
    test_coordsetslicer<ExecSpace>(coordset, makeView);
  }
};

TEST(mir_blueprint_utilities, coordsetslicer_uniform_seq)
{
  coordsetslicer_uniform<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, coordsetslicer_uniform_omp)
{
  coordsetslicer_uniform<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, coordsetslicer_uniform_cuda)
{
  coordsetslicer_uniform<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, coordsetslicer_uniform_hip)
{
  coordsetslicer_uniform<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
template <typename ExecSpace>
struct test_extractzones
{
  static void test()
  {
    conduit::Node hostMesh;
    create(hostMesh);

    // host->device
    conduit::Node deviceMesh;
    bputils::copy<ExecSpace>(deviceMesh, hostMesh);

    axom::Array<axom::IndexType> ids{{1,3,5}};
    const auto nzones = ids.size();
    axom::Array<axom::IndexType> selectedZones(nzones, nzones, axom::execution_space<ExecSpace>::allocatorID());
    axom::copy(selectedZones.data(), ids.data(), nzones * sizeof(axom::IndexType));

    // Wrap the data in views.
    auto coordsetView = axom::mir::views::make_explicit_coordset<conduit::float64, 2>::view(deviceMesh["coordsets/coords"]);
    using CoordsetView = decltype(coordsetView);

    using TopologyView = axom::mir::views::UnstructuredTopologySingleShapeView<axom::mir::views::QuadShape<conduit::int64>>;
    TopologyView topoView(bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/connectivity"]),
                          bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/sizes"]),
                          bputils::make_array_view<conduit::int64>(deviceMesh["topologies/mesh/elements/offsets"]));

    // Pull out selected zones
    bputils::ExtractZones<ExecSpace, TopologyView, CoordsetView> extract(topoView, coordsetView);
    conduit::Node options, newDeviceMesh;
    options["topology"] = "mesh";
    extract.execute(selectedZones.view(), deviceMesh, options, newDeviceMesh);

    // device->host
    conduit::Node newHostMesh;
    bputils::copy<ExecSpace>(newHostMesh, newDeviceMesh);

    //printNode(newHostMesh);

    // Check some of the key arrays
    const axom::Array<conduit::int64> connectivity{{0, 1, 4, 3, 2, 3, 7, 6, 4, 5, 9, 8}};
    const axom::Array<conduit::int64> sizes{{4, 4, 4}};
    const axom::Array<conduit::int64> offsets{{0, 4, 8}};
    const axom::Array<conduit::float64> x{{1.0, 2.0, 0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.0, 3.0}};
    const axom::Array<conduit::float64> y{{0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0}};
    const axom::Array<conduit::float64> zonal{{1.0, 3.0, 5.0}};
    const axom::Array<conduit::float64> nodal{{1.0, 2.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0}};

    EXPECT_TRUE(compare_views(connectivity.view(), bputils::make_array_view<conduit::int64>(newHostMesh["topologies/mesh/elements/connectivity"])));
    EXPECT_TRUE(compare_views(sizes.view(), bputils::make_array_view<conduit::int64>(newHostMesh["topologies/mesh/elements/sizes"])));
    EXPECT_TRUE(compare_views(offsets.view(), bputils::make_array_view<conduit::int64>(newHostMesh["topologies/mesh/elements/offsets"])));
    EXPECT_TRUE(compare_views(x.view(), bputils::make_array_view<conduit::float64>(newHostMesh["coordsets/coords/values/x"])));
    EXPECT_TRUE(compare_views(y.view(), bputils::make_array_view<conduit::float64>(newHostMesh["coordsets/coords/values/y"])));
    EXPECT_TRUE(compare_views(zonal.view(), bputils::make_array_view<conduit::float64>(newHostMesh["fields/zonal/values"])));
    EXPECT_TRUE(compare_views(nodal.view(), bputils::make_array_view<conduit::float64>(newHostMesh["fields/nodal/values"])));

    // Do the material too.
    using MatsetView = axom::mir::views::UnibufferMaterialView<conduit::int64, conduit::float64, 3>;
    MatsetView matsetView;
    matsetView.set(
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/material_ids"]),
      bputils::make_array_view<conduit::float64>(deviceMesh["matsets/mat1/volume_fractions"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/sizes"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/offsets"]),
      bputils::make_array_view<conduit::int64>(deviceMesh["matsets/mat1/indices"]));

    // Pull out selected zones
    bputils::ExtractZonesAndMatset<ExecSpace, TopologyView, CoordsetView, MatsetView> extractM(topoView, coordsetView, matsetView);
    newDeviceMesh.reset();
    extractM.execute(selectedZones.view(), deviceMesh, options, newDeviceMesh);

    // device->host
    newHostMesh.reset();
    bputils::copy<ExecSpace>(newHostMesh, newDeviceMesh);

    // Check some of the key arrays in the sliced material
    const axom::Array<conduit::int64> mat_sizes{{2, 1, 2}};
    const axom::Array<conduit::int64> mat_offsets{{0, 2, 3}};
    const axom::Array<conduit::int64> mat_indices{{0, 1, 2, 3, 4}};
    const axom::Array<conduit::int64> mat_material_ids{{1, 2, 2, 2, 3}};
    const axom::Array<conduit::float64> mat_volume_fractions{{0.5, 0.5, 1.0, 0.8, 0.2}};

    EXPECT_TRUE(compare_views(mat_sizes.view(), bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/sizes"])));
    EXPECT_TRUE(compare_views(mat_offsets.view(), bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/offsets"])));
    EXPECT_TRUE(compare_views(mat_indices.view(), bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/indices"])));
    EXPECT_TRUE(compare_views(mat_material_ids.view(), bputils::make_array_view<conduit::int64>(newHostMesh["matsets/mat1/material_ids"])));
    EXPECT_TRUE(compare_views(mat_volume_fractions.view(), bputils::make_array_view<conduit::float64>(newHostMesh["matsets/mat1/volume_fractions"])));
  }

  static void create(conduit::Node &hostMesh)
  {
    /*
      8-------9------10------11
      |  2/1  | 1/0.1 | 2/0.8 |
      |       | 2/0.5 | 3/0.2 |
      |       | 3/0.4 |       |
      4-------5-------6-------7
      |       | 1/0.5 | 1/0.2 |
      |  1/1  | 2/0.5 | 2/0.8 |
      |       |       |       |
      0-------1-------2-------3
      */
    const char *yaml = R"xx(
coordsets:
  coords:
    type: explicit
    values:
      x: [0., 1., 2., 3., 0., 1., 2., 3., 0., 1., 2., 3.]
      y: [0., 0., 0., 0., 1., 1., 1., 1., 2., 2., 2., 2.]
topologies:
  mesh:
    coordset: coords
    type: unstructured
    elements:
      shape: quad
      connectivity: [0,1,5,4, 1,2,6,5, 2,3,7,6, 4,5,9,8, 5,6,10,9, 6,7,11,10]
      sizes: [4,4,4,4,4,4]
      offsets: [0,4,8,12,16,20]
fields:
  zonal:
    topology: mesh
    association: element
    values: [0.,1.,2.,3.,4.,5.]
  nodal:
    topology: mesh
    association: vertex
    values: [0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.]
matsets:
  mat1:
    topology: mesh
    material_map:
      a: 1
      b: 2
      c: 3
    material_ids: [1, 1,2, 1,2, 2, 1,2,3, 2,3]
    volume_fractions: [1., 0.5,0.5, 0.2,0.8, 1., 0.1,0.5,0.4, 0.8,0.2]
    indices: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    sizes: [1, 2, 2, 1, 3, 2]
    offsets: [0, 1, 3, 5, 6, 9]
)xx";

    hostMesh.parse(yaml); 
  }
};

TEST(mir_blueprint_utilities, extractzones_seq)
{
  test_extractzones<seq_exec>::test();
}

#if defined(AXOM_USE_OPENMP)
TEST(mir_blueprint_utilities, extractzones_omp)
{
  test_extractzones<omp_exec>::test();
}
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
TEST(mir_blueprint_utilities, extractzones_cuda)
{
  test_extractzones<cuda_exec>::test();
}
#endif

#if defined(AXOM_USE_HIP)
TEST(mir_blueprint_utilities, extractzones_hip)
{
  test_extractzones<hip_exec>::test();
}
#endif

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
