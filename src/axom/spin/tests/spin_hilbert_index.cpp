// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file spin_hilbert_index.cpp
 *  \brief Tests the Hilbert curve convention and reversible integer conversions.
 */

#include "axom/spin/HilbertIndex.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/memory_management.hpp"
#include "axom/slic.hpp"
#include "gtest/gtest.h"

#include <cstdint>
#include <cstdlib>
#include <limits>
#include <random>
#include <type_traits>
#include <vector>

namespace
{
template <int DIM>
void exhaustive(int bits)
{
  using H = axom::spin::Hilbertizer<std::int64_t, std::uint64_t, DIM>;
  const std::uint64_t side = std::uint64_t {1} << bits;
  const std::uint64_t count = std::uint64_t {1} << (bits * DIM);
  std::vector<bool> seen(count, false);

  // Enumerate the Cartesian grid independently of the Hilbert decoder.
  for(std::uint64_t ordinal = 0; ordinal < count; ++ordinal)
  {
    SCOPED_TRACE(::testing::Message() << "dimension " << DIM << ", point " << ordinal);
    typename H::PointType point;
    auto remaining = ordinal;
    for(int axis = 0; axis < DIM; ++axis)
    {
      point[axis] = static_cast<std::int64_t>(remaining % side);
      remaining /= side;
    }
    const auto key = H::toIndex(point, bits);
    ASSERT_LT(key, count);
    ASSERT_FALSE(seen[key]);
    seen[key] = true;
    EXPECT_EQ(H::toPoint(key, bits), point);
  }

  typename H::PointType last;
  for(std::uint64_t key = 0; key < count; ++key)
  {
    SCOPED_TRACE(::testing::Message() << "dimension " << DIM << ", key " << key);
    EXPECT_TRUE(seen[key]);
    const auto point = H::toPoint(key, bits);
    ASSERT_TRUE(H::isValidPoint(point, bits));
    EXPECT_EQ(H::toIndex(point, bits), key);
    if(key > 0)
    {
      std::int64_t distance = 0;
      for(int axis = 0; axis < DIM; ++axis)
      {
        distance += std::abs(point[axis] - last[axis]);
      }
      EXPECT_EQ(distance, 1);
    }
    if(bits > 1)
    {
      auto parent = point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        parent[axis] >>= 1;
      }
      EXPECT_EQ(H::toIndex(parent, bits - 1), key >> DIM);
    }
    last = point;
  }
}

template <typename Coord, typename Index, int DIM>
void allDepths()
{
  using H = axom::spin::Hilbertizer<Coord, Index, DIM>;
  constexpr int coord_bits = std::numeric_limits<Coord>::digits;
  constexpr int index_bits = std::numeric_limits<Index>::digits;
  constexpr int capacity = coord_bits < index_bits / DIM ? coord_bits : index_bits / DIM;
  EXPECT_EQ(H::COORD_BITS, coord_bits);
  EXPECT_EQ(H::INDEX_BITS, index_bits);
  EXPECT_EQ(H::MAX_UNIQUE_BITS, capacity);
  EXPECT_EQ(H::maxBitsPerCoord(), capacity);
  std::mt19937_64 random(0x534b494c4c494e47ULL + DIM);
  for(int bits = 1; bits <= capacity; ++bits)
  {
    SCOPED_TRACE(::testing::Message()
                 << "dimension " << DIM << ", depth " << bits << ", coordinate digits "
                 << coord_bits << ", index digits " << index_bits);
    const std::uint64_t max_coord = (std::uint64_t {1} << bits) - 1;
    const auto max_key = std::numeric_limits<std::uint64_t>::max() >> (64 - bits * DIM);
    EXPECT_TRUE(H::isValidDepth(bits));
    EXPECT_TRUE(H::isValidIndex(static_cast<Index>(max_key), bits));
    for(int iteration = 0; iteration < 100; ++iteration)
    {
      typename H::PointType point;
      for(int axis = 0; axis < DIM; ++axis)
      {
        // Check every corner, then sample points with a fixed seed.
        point[axis] =
          static_cast<Coord>(iteration < (1 << DIM) ? (((iteration >> axis) & 1) ? max_coord : 0)
                                                    : (random() & max_coord));
      }
      ASSERT_TRUE(H::isValidPoint(point, bits));
      const auto encoded = H::toIndex(point, bits);
      ASSERT_TRUE(H::isValidIndex(encoded, bits));
      EXPECT_EQ(H::toPoint(encoded, bits), point);
      EXPECT_EQ(axom::spin::convertPointToHilbert<Index>(point, bits), encoded);
      EXPECT_EQ((axom::spin::convertHilbertToPoint<Coord, DIM>(encoded, bits)), point);

      // Independent keys: zero, maximum, every single bit, then random keys.
      const auto key =
        static_cast<Index>(iteration == 0                  ? 0
                             : iteration == 1              ? max_key
                             : iteration <= bits * DIM + 1 ? (std::uint64_t {1} << (iteration - 2))
                                                           : random() & max_key);
      const auto decoded = H::toPoint(key, bits);
      ASSERT_TRUE(H::isValidPoint(decoded, bits));
      EXPECT_EQ(H::toIndex(decoded, bits), key);
    }

    // Validate the first out-of-range value whenever the input type can hold it.
    if(bits < coord_bits)
    {
      for(int axis = 0; axis < DIM; ++axis)
      {
        typename H::PointType outside;
        outside[axis] = static_cast<Coord>(max_coord + 1);
        EXPECT_FALSE(H::isValidPoint(outside, bits));
      }
    }
    if(bits * DIM < index_bits)
    {
      EXPECT_FALSE(H::isValidIndex(static_cast<Index>(max_key + 1), bits));
    }
    if(bits == capacity)
    {
      const typename H::PointType point(static_cast<Coord>(max_coord));
      const auto key = H::toIndex(point);
      EXPECT_TRUE(H::isValidPoint(point));
      EXPECT_TRUE(H::isValidIndex(key));
      EXPECT_EQ(key, H::toIndex(point, bits));
      EXPECT_EQ(H::toPoint(key), point);
      EXPECT_EQ(axom::spin::convertPointToHilbert<Index>(point), key);
      EXPECT_EQ((axom::spin::convertHilbertToPoint<Coord, DIM>(key)), point);
    }
  }
  for(int bits :
      {std::numeric_limits<int>::min(), -1, 0, capacity + 1, std::numeric_limits<int>::max()})
  {
    EXPECT_FALSE(H::isValidDepth(bits));
    EXPECT_FALSE(H::isValidPoint(typename H::PointType(), bits));
    EXPECT_FALSE(H::isValidIndex(0, bits));
  }
  if constexpr(std::is_signed<Coord>::value)
  {
    for(int axis = 0; axis < DIM; ++axis)
    {
      typename H::PointType negative;
      negative[axis] = -1;
      EXPECT_FALSE(H::isValidPoint(negative));
    }
  }
  if constexpr(std::is_signed<Index>::value)
  {
    EXPECT_FALSE(H::isValidIndex(-1));
  }
}

template <typename Index, int DIM>
void everyCoordinateWidth()
{
  allDepths<std::int8_t, Index, DIM>();
  allDepths<std::uint8_t, Index, DIM>();
  allDepths<std::int16_t, Index, DIM>();
  allDepths<std::uint16_t, Index, DIM>();
  allDepths<std::int32_t, Index, DIM>();
  allDepths<std::uint32_t, Index, DIM>();
  allDepths<std::int64_t, Index, DIM>();
  allDepths<std::uint64_t, Index, DIM>();
}

template <int DIM>
void everyIndexWidth()
{
  everyCoordinateWidth<std::int8_t, DIM>();
  everyCoordinateWidth<std::uint8_t, DIM>();
  everyCoordinateWidth<std::int16_t, DIM>();
  everyCoordinateWidth<std::uint16_t, DIM>();
  everyCoordinateWidth<std::int32_t, DIM>();
  everyCoordinateWidth<std::uint32_t, DIM>();
  everyCoordinateWidth<std::int64_t, DIM>();
  everyCoordinateWidth<std::uint64_t, DIM>();
}

template <typename ExecSpace, int DIM>
void checkExecutionSpace()
{
  using H = axom::spin::Hilbertizer<std::uint32_t, std::uint64_t, DIM>;
  constexpr int count = 128;
  auto* results =
    axom::allocate<std::uint64_t>(count * 3, axom::execution_space<ExecSpace>::allocatorID());
  axom::for_all<ExecSpace>(count, [=] AXOM_HOST_DEVICE(axom::IndexType i) {
    const int bits = 1 + static_cast<int>(i % H::MAX_UNIQUE_BITS);
    const auto mask = axom::numeric_limits<std::uint64_t>::max() >> (64 - bits * DIM);
    const auto key = (0x9e3779b97f4a7c15ULL * (i + 1)) & mask;
    const auto point = H::toPoint(key, bits);
    results[3 * i] = H::toIndex(point, bits);
    results[3 * i + 1] = axom::spin::convertPointToHilbert<std::uint64_t>(
      axom::spin::convertHilbertToPoint<std::uint32_t, DIM>(key, bits),
      bits);
    results[3 * i + 2] = H::isValidPoint(point, bits) && H::isValidIndex(key, bits);
  });
  std::uint64_t host[count * 3];
  axom::copy(host, results, sizeof(host));
  axom::deallocate(results);
  for(int i = 0; i < count; ++i)
  {
    const int bits = 1 + i % H::MAX_UNIQUE_BITS;
    const auto mask = std::numeric_limits<std::uint64_t>::max() >> (64 - bits * DIM);
    const auto key = (0x9e3779b97f4a7c15ULL * (i + 1)) & mask;
    EXPECT_EQ(host[3 * i], key);
    EXPECT_EQ(host[3 * i + 1], key);
    EXPECT_EQ(host[3 * i + 2], 1);
  }
}
}  // namespace

TEST(spin_hilbert_index, known_orientation)
{
  using H2 = axom::spin::Hilbertizer<std::int64_t, std::uint64_t, 2>;
  const H2::PointType square[] = {{0, 0}, {0, 1}, {1, 1}, {1, 0}};
  for(std::uint64_t key = 0; key < 4; ++key)
  {
    EXPECT_EQ(H2::toIndex(square[key], 1), key);
    EXPECT_EQ(H2::toPoint(key, 1), square[key]);
  }
  const H2::PointType curve[] = {{0, 0},
                                 {1, 0},
                                 {1, 1},
                                 {0, 1},
                                 {0, 2},
                                 {0, 3},
                                 {1, 3},
                                 {1, 2},
                                 {2, 2},
                                 {2, 3},
                                 {3, 3},
                                 {3, 2},
                                 {3, 1},
                                 {2, 1},
                                 {2, 0},
                                 {3, 0}};
  for(std::uint64_t key = 0; key < 16; ++key)
  {
    EXPECT_EQ(H2::toIndex(curve[key], 2), key);
    EXPECT_EQ(H2::toPoint(key, 2), curve[key]);
  }
  using H3 = axom::spin::Hilbertizer<std::int64_t, std::int64_t, 3>;
  EXPECT_EQ(H3::toIndex({5, 10, 20}, 5), 7865);
  EXPECT_EQ(H3::toPoint(7865, 5), (H3::PointType {5, 10, 20}));
  const H3::PointType cube[] =
    {{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}, {1, 0, 0}};
  for(std::int64_t key = 0; key < 8; ++key)
  {
    EXPECT_EQ(H3::toIndex(cube[key], 1), key);
    EXPECT_EQ(H3::toPoint(key, 1), cube[key]);
  }
}

TEST(spin_hilbert_index, four_dimensional_orientation)
{
  using H = axom::spin::Hilbertizer<std::uint32_t, std::uint64_t, 4>;
  // At depth one, use binary-reflected Gray codes with axis zero first.
  for(std::uint64_t key = 0; key < 16; ++key)
  {
    const auto gray = key ^ (key >> 1);
    H::PointType point;
    for(int axis = 0; axis < 4; ++axis)
    {
      point[axis] = static_cast<std::uint32_t>((gray >> (3 - axis)) & 1);
    }
    EXPECT_EQ(H::toIndex(point, 1), key);
    EXPECT_EQ(H::toPoint(key, 1), point);
  }
}

TEST(spin_hilbert_index, exhaustive_bijection_adjacency_and_prefixes)
{
  for(int bits = 1; bits <= 5; ++bits)
  {
    SCOPED_TRACE(bits);
    exhaustive<2>(bits);
    exhaustive<3>(bits);
  }
  for(int bits = 1; bits <= 4; ++bits)
  {
    SCOPED_TRACE(bits);
    exhaustive<4>(bits);
  }
}

TEST(spin_hilbert_index, every_depth_and_integer_width)
{
  everyIndexWidth<2>();
  everyIndexWidth<3>();
  everyIndexWidth<4>();
}

TEST(spin_hilbert_index, wide_coordinates_are_validated_before_narrowing)
{
  using H = axom::spin::Hilbertizer<std::uint64_t, std::uint8_t, 2>;
  EXPECT_TRUE(H::isValidPoint({15, 15}));
  EXPECT_FALSE(H::isValidPoint({std::uint64_t {1} << 40, 0}));
  EXPECT_FALSE(H::isValidPoint({0, (std::uint64_t {1} << 40) | 1}));
}

#if !defined(NDEBUG) && GTEST_HAS_DEATH_TEST
TEST(spin_hilbert_index, invalid_preconditions)
{
  using H = axom::spin::Hilbertizer<int, std::int32_t, 2>;
  EXPECT_DEATH_IF_SUPPORTED(H::toIndex({-1, 0}, 2), "isValidPoint");
  EXPECT_DEATH_IF_SUPPORTED(H::toIndex({0, 4}, 2), "isValidPoint");
  EXPECT_DEATH_IF_SUPPORTED(H::toPoint(-1, 2), "isValidIndex");
  EXPECT_DEATH_IF_SUPPORTED(H::toPoint(16, 2), "isValidIndex");
  EXPECT_DEATH_IF_SUPPORTED(H::toIndex({0, 0}, 0), "isValidPoint");
  EXPECT_DEATH_IF_SUPPORTED(H::toPoint(0, H::MAX_UNIQUE_BITS + 1), "isValidIndex");
}
#endif

AXOM_CUDA_TEST(spin_hilbert_index, execution_spaces)
{
  checkExecutionSpace<axom::SEQ_EXEC, 2>();
  checkExecutionSpace<axom::SEQ_EXEC, 3>();
  checkExecutionSpace<axom::SEQ_EXEC, 4>();
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(RAJA_ENABLE_CUDA) && \
  defined(AXOM_USE_UMPIRE)
  checkExecutionSpace<axom::CUDA_EXEC<256>, 2>();
  checkExecutionSpace<axom::CUDA_EXEC<256>, 3>();
  checkExecutionSpace<axom::CUDA_EXEC<256>, 4>();
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_HIP) && defined(RAJA_ENABLE_HIP) && \
  defined(AXOM_USE_UMPIRE)
  checkExecutionSpace<axom::HIP_EXEC<256>, 2>();
  checkExecutionSpace<axom::HIP_EXEC<256>, 3>();
  checkExecutionSpace<axom::HIP_EXEC<256>, 4>();
#endif
}

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger;
  return RUN_ALL_TESTS();
}
