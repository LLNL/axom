// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"

#include "axom/primal/operators/is_convex.hpp"

#include "axom/core/execution/execution_space.hpp"

// Get permutations for hexahedron vertex order
#include <algorithm>

TEST(primal_is_convex, convex_hexahedron_indices)
{
  const int NUM_VERTS = 8;

  using PointType = axom::primal::Point<double, 3>;
  using HexahedronType = axom::primal::Hexahedron<double, 3>;

  PointType data[NUM_VERTS];
  data[0] = PointType {0, 0, 0};
  data[1] = PointType {1, 0, 0};
  data[2] = PointType {1, 1, 0};
  data[3] = PointType {0, 1, 0};
  data[4] = PointType {0, 0, 1};
  data[5] = PointType {1, 0, 1};
  data[6] = PointType {1, 1, 1};
  data[7] = PointType {0, 1, 1};

  // Indices to permute
  std::vector<int> indices = {0, 1, 2, 3, 4, 5, 6, 7};

  // Verify there are 48 index orderings where a unit hexahedron
  // is convex (4 for each face, 6 faces in a hexahedron, 2
  // clockwise/counter-clockwise winding orientations)
  int num_convex_orderings = 0;
  do
  {
    HexahedronType hex(data[indices[0]],
                       data[indices[1]],
                       data[indices[2]],
                       data[indices[3]],
                       data[indices[4]],
                       data[indices[5]],
                       data[indices[6]],
                       data[indices[7]]);

    num_convex_orderings += axom::primal::is_convex(hex);
  } while(std::next_permutation(indices.begin(), indices.end()));

  EXPECT_EQ(num_convex_orderings, 48);
}

template <typename ExecPolicy>
void check_is_convex_hex_policy()
{
  using PointType = axom::primal::Point<double, 3>;
  using HexahedronType = axom::primal::Hexahedron<double, 3>;

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  axom::Array<bool> convex_val_d(1, 1, kernel_allocator);
  auto convex_val_v = convex_val_d.view();

  axom::Array<bool> non_convex_val_d(1, 1, kernel_allocator);
  auto non_convex_val_v = non_convex_val_d.view();

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      // Use literal constant for Windows
      PointType data[8];
      data[0] = PointType {0, 0, 0};
      data[1] = PointType {1, 0, 0};
      data[2] = PointType {1, 1, 0};
      data[3] = PointType {0, 1, 0};
      data[4] = PointType {0, 0, 1};
      data[5] = PointType {1, 0, 1};
      data[6] = PointType {1, 1, 1};
      data[7] = PointType {0, 1, 1};

      PointType non_planar_datum {0, -1, 0};

      HexahedronType hex_convex(data[0],
                                data[1],
                                data[2],
                                data[3],
                                data[4],
                                data[5],
                                data[6],
                                data[7]);
      convex_val_v[i] = axom::primal::is_convex(hex_convex);

      HexahedronType hex_non_convex(non_planar_datum,
                                    data[1],
                                    data[2],
                                    data[3],
                                    data[4],
                                    data[5],
                                    data[6],
                                    data[7]);
      non_convex_val_v[i] = axom::primal::is_convex(hex_non_convex);
    });

  // Copy convexity checks back to host
  axom::Array<bool> convex_val_h =
    axom::Array<bool>(convex_val_d, host_allocator);

  axom::Array<bool> non_convex_val_h =
    axom::Array<bool>(non_convex_val_d, host_allocator);

  // Verify values
  EXPECT_EQ(convex_val_h[0], true);
  EXPECT_EQ(non_convex_val_h[0], false);
}

//------------------------------------------------------------------------------
TEST(primal_is_convex, convex_hexahedron_seq)
{
  check_is_convex_hex_policy<axom::SEQ_EXEC>();
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #ifdef AXOM_USE_OPENMP
TEST(primal_is_convex, convex_hexahedron_omp)
{
  check_is_convex_hex_policy<axom::OMP_EXEC>();
}
  #endif /* AXOM_USE_OPENMP */

  #ifdef AXOM_USE_CUDA
TEST(primal_is_convex, convex_hexahedron_cuda)
{
  check_is_convex_hex_policy<axom::CUDA_EXEC<256>>();
}
  #endif /* AXOM_USE_CUDA */

  #ifdef AXOM_USE_HIP
TEST(primal_is_convex, convex_hexahedron_hip)
{
  check_is_convex_hex_policy<axom::HIP_EXEC<256>>();
}
  #endif /* AXOM_USE_HIP */

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  int result = RUN_ALL_TESTS();

  return result;
}