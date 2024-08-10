// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MIR_UTILITIES_TEST_H_
#define MIR_UTILITIES_TEST_H_

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/mir.hpp"

#include "axom/mir/tests/mir_testing_helpers.hpp"
#include "axom/mir/tests/mir_testing_data_helpers.hpp"

#include <algorithm>

namespace mir = axom::mir;

//----------------------------------------------------------------------

TEST(mir_shape_tests, shape_dimesionality)
{
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Triangle), false);
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Quad), false);
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Tetrahedron),
            true);
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Pyramid), true);
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Triangular_Prism),
            true);
  EXPECT_EQ(mir::utilities::isShapeThreeDimensional(mir::Shape::Hexahedron),
            true);
}

//----------------------------------------------------------------------

TEST(mir_shape_tests, check_shape_central_vertex)
{
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Triangle, 6), false);

  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Quad, 8), false);

  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 9), false);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 10), true);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Tetrahedron, 11), false);

  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Pyramid, 12), false);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Pyramid, 13), true);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Pyramid, 14), false);

  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 14),
            false);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 15),
            true);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Triangular_Prism, 16),
            false);

  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 19), false);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 20), true);
  EXPECT_EQ(mir::utilities::isCenterVertex(mir::Shape::Hexahedron, 21), false);
}

//----------------------------------------------------------------------

TEST(mir_shape_tests, determine_central_vertex)
{
  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Triangle), -1);
  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Quad), -1);

  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Tetrahedron), 10);
  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Pyramid), 13);
  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Triangular_Prism), 15);
  EXPECT_EQ(mir::utilities::getCenterVertex(mir::Shape::Hexahedron), 20);
}

//----------------------------------------------------------------------

TEST(mir_compute_averages, float_value)
{
  std::vector<axom::float64> values = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};

  axom::float64 average = mir::utilities::computeAverageFloat(values);

  EXPECT_DOUBLE_EQ(average, 2.5);
}

//----------------------------------------------------------------------

TEST(mir_compute_averages, point_value)
{
  std::vector<mir::Point2> points = {mir::Point2::make_point(0.0, 0.0, 0.0),
                                     mir::Point2::make_point(1.0, 0.0, 0.0),
                                     mir::Point2::make_point(0.0, 1.0, 0.0),
                                     mir::Point2::make_point(0.0, 0.0, 1.0),
                                     mir::Point2::make_point(1.0, 1.0, 0.0),
                                     mir::Point2::make_point(1.0, 0.0, 1.0),
                                     mir::Point2::make_point(0.0, 1.0, 1.0),
                                     mir::Point2::make_point(1.0, 1.0, 1.0)};

  mir::Point2 centroid = mir::utilities::computeAveragePoint(points);

  EXPECT_DOUBLE_EQ(centroid[0], 0.5);
  EXPECT_DOUBLE_EQ(centroid[1], 0.5);
  EXPECT_DOUBLE_EQ(centroid[2], 0.5);
}

//----------------------------------------------------------------------

template <typename ExecSpace, typename IndexT = int>
void test_node_to_zone_relation_builder(const conduit::Node &hostMesh)
{
  namespace bputils = axom::mir::utilities::blueprint;

  // host -> device
  conduit::Node deviceMesh;
  axom::mir::utilities::blueprint::copy<ExecSpace>(deviceMesh, hostMesh);
  conduit::Node &deviceTopo = deviceMesh["topologies/mesh"];

  // Run the algorithm on the device
  conduit::Node deviceRelation;
  axom::mir::NodeToZoneRelationBuilder<ExecSpace> n2z;
  n2z.execute(deviceTopo, deviceRelation);

  // device -> host
  conduit::Node hostRelation;
  axom::mir::utilities::blueprint::copy<seq_exec>(hostRelation, deviceRelation);
#if 0
  // Print the results.
  conduit::Node options;
  options["num_children_threshold"] = 1000;
  options["num_elements_threshold"] = 1000;
  hostRelation.to_summary_string_stream(std::cout, options);
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
  const auto zonesView = bputils::make_array_view<IndexT>(hostRelation["zones"]);
  const auto sizesView = bputils::make_array_view<IndexT>(hostRelation["sizes"]);
  const auto offsetsView = bputils::make_array_view<IndexT>(hostRelation["offsets"]);
  EXPECT_EQ(sizesView.size(), sizeof(sizes)/sizeof(int));
  EXPECT_EQ(offsetsView.size(), sizeof(offsets)/sizeof(int));
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

TEST(mir_utilities, node_to_zone_relation_builder_unstructured)
{
   /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims{{4, 3}};
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

TEST(mir_utilities, node_to_zone_relation_builder_rectilinear)
{
   /*
    8---9--10--11
    |   |   |   |
    4---5---6---7
    |   |   |   |
    0---1---2---3
    */
  conduit::Node mesh;
  axom::StackArray<int, 2> dims{{4, 3}};
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

int main(int argc, char* argv[])
{
  int result = 0;
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;  // create & initialize test logger,

  result = RUN_ALL_TESTS();
  return result;
}

#endif  //  MIR_UTILITIES_TEST_H_
