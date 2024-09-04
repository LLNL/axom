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

template <typename ExecSpace>
void test_conduit_allocate()
{
  axom::mir::utilities::blueprint::ConduitAllocateThroughAxom<ExecSpace> c2a;
  EXPECT_TRUE(c2a.getConduitAllocatorID() > 0);

  constexpr int nValues = 100;
  conduit::Node n;
  n.set_allocator(c2a.getConduitAllocatorID());
  n.set(conduit::DataType::int32(nValues));

  // Make sure we can store some values into the data that were allocated.
  int *ptr = static_cast<int *>(n.data_ptr());
  axom::for_all<ExecSpace>(
    nValues,
    AXOM_LAMBDA(auto index) { ptr[index] = index; });

  EXPECT_EQ(n.dtype().number_of_elements(), nValues);
}

TEST(mir_blueprint_utilities, allocate)
{
  test_conduit_allocate<seq_exec>();
#if defined(AXOM_USE_OPENMP)
  test_conduit_allocate<omp_exec>();
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_conduit_allocate<cuda_exec>();
#endif
#if defined(AXOM_USE_HIP)
  test_conduit_allocate<hip_exec>();
#endif
}

template <typename ExecSpace>
void test_copy_braid(const conduit::Node &mesh)
{
  // Copy the mesh to device.
  conduit::Node mesh_dev;
  axom::mir::utilities::blueprint::copy<ExecSpace>(mesh_dev, mesh);

  // Run some minmax operations on device (proves that the data was in the right place) and check the results.

  constexpr double eps = 1.e-7;

  auto x = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
    mesh["coordsets/coords/values/x"]);
  //std::cout << std::setw(16) << "x={" << x.first << ", " << x.second << "}\n";
  EXPECT_NEAR(x.first, -10., eps);
  EXPECT_NEAR(x.second, 10., eps);

  auto y = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
    mesh["coordsets/coords/values/y"]);
  //std::cout << std::setw(16) << "y={" << y.first << ", " << y.second << "}\n";
  EXPECT_NEAR(y.first, -10., eps);
  EXPECT_NEAR(y.second, 10., eps);

  auto c = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
    mesh["topologies/mesh/elements/connectivity"]);
  //std::cout << std::setw(16) << "conn={" << c.first << ", " << c.second << "}\n";
  EXPECT_NEAR(c.first, 0., eps);
  EXPECT_NEAR(c.second, 999., eps);

  auto r = axom::mir::utilities::blueprint::minmax<ExecSpace, double>(
    mesh["fields/radial/values"]);
  //std::cout << std::setw(16) << "radial={" << r.first << ", " << r.second << "}\n";
  EXPECT_NEAR(r.first, 19.2450089729875, eps);
  EXPECT_NEAR(r.second, 173.205080756888, eps);
}

TEST(mir_blueprint_utilities, copy)
{
  const int d[3] = {10, 10, 10};
  conduit::Node mesh;
  conduit::blueprint::mesh::examples::braid("hexs", d[0], d[1], d[2], mesh);

  test_copy_braid<seq_exec>(mesh);
#if defined(AXOM_USE_OPENMP)
  test_copy_braid<omp_exec>(mesh);
#endif
#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_copy_braid<cuda_exec>(mesh);
#endif
#if defined(AXOM_USE_HIP)
  test_copy_braid<hip_exec>(mesh);
#endif
}

TEST(mir_blueprint_utilities, to_unstructured)
{
  // TODO: to_unstructured
}

//----------------------------------------------------------------------

template <typename IndexT = int>
void compareRelation(const conduit::Node &hostRelation,
                     const axom::ArrayView<const int> &zones,
                     const axom::ArrayView<const int> &sizes,
                     const axom::ArrayView<const int> &offsets)
{
  namespace bputils = axom::mir::utilities::blueprint;

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
void test_recenter_field(const conduit::Node &hostMesh)
{
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

TEST(mir_blueprint_utilities, recenterfield)
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
  mesh["topologies/mesh/elements/sizes"].set(
    std::vector<int> {{4, 4, 4, 4, 4, 4}});
  mesh["topologies/mesh/elements/offsets"].set(
    std::vector<int> {{0, 4, 8, 12, 16, 20}});
  mesh["fields/easy_zonal/topology"] = "mesh";
  mesh["fields/easy_zonal/association"] = "element";
  mesh["fields/easy_zonal/values"].set(std::vector<float> {{1, 3, 5, 7, 9, 11}});

  test_recenter_field<seq_exec>(mesh);
#if defined(AXOM_USE_OPENMP)
  test_recenter_field<omp_exec>(mesh);
#endif

#if defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  test_recenter_field<cuda_exec>(mesh);
#endif

#if defined(AXOM_USE_HIP)
  test_recenter_field<hip_exec>(mesh);
#endif
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}
