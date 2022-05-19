// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/utilities/BitUtilities.hpp"

#include <limits>

namespace
{
constexpr int NRAND_SAMPLES = 1000;

// Returns a random integer value between 0 and the max value
template <typename T>
T random_int()
{
  static_assert(std::is_integral<T>::value, "T must be an integral type");

  constexpr T max_int = std::numeric_limits<T>::max();
  constexpr double max_d = static_cast<double>(max_int);

  const auto val = axom::utilities::random_real(0., max_d);
  return static_cast<T>(val);
}

axom::uint64 shifted(int bit) { return static_cast<axom::uint64>(1) << bit; }

}  // namespace

TEST(core_bit_utilities, trailingZeroes)
{
  constexpr axom::uint64 ZERO = axom::uint64(0);
  constexpr int BITS = axom::utilities::BitTraits<axom::uint64>::BITS_PER_WORD;
  ASSERT_EQ(64, BITS);

  // Axom's trailingZeros will return 64 when given 0
  {
    EXPECT_EQ(BITS, axom::utilities::trailingZeros(ZERO));
  }

  // Test with a known trailing bit
  for(int i = 0; i < BITS; ++i)
  {
    axom::uint64 val = ::shifted(i);
    EXPECT_EQ(i, axom::utilities::trailingZeros(val));

    // Value doesn't change when you set bits to left of trailing bit
    for(int j = i + 1; j < BITS; ++j)
    {
      axom::uint64 val2 = ::shifted(i) + ::shifted(j);
      EXPECT_EQ(axom::utilities::trailingZeros(val),
                axom::utilities::trailingZeros(val2));
    }
  }

  // Test trailing zeros on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    axom::uint64 rand_val = ::random_int<axom::uint64>();

    // Explicitly find first bit and check that it matches
    int bit = 0;
    for(; bit < BITS; ++bit)
    {
      if(rand_val & shifted(bit)) break;
    }
    EXPECT_EQ(bit, axom::utilities::trailingZeros(rand_val));
  }
}

TEST(core_bit_utilities, popCount)
{
  constexpr axom::uint64 ZERO = axom::uint64(0);
  constexpr int BITS = axom::utilities::BitTraits<axom::uint64>::BITS_PER_WORD;
  ASSERT_EQ(64, BITS);

  // Test pop count when zero bits are set
  {
    EXPECT_EQ(0, axom::utilities::popCount(ZERO));
    EXPECT_EQ(BITS, axom::utilities::popCount(~ZERO));
  }

  // Test pop count when one bit is set
  for(int i = 0; i < BITS; ++i)
  {
    axom::uint64 val = ::shifted(i);
    EXPECT_EQ(1, axom::utilities::popCount(val));
    EXPECT_EQ(BITS - 1, axom::utilities::popCount(~val));
  }

  // Test pop count when two bits are set
  for(int i = 0; i < BITS; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      axom::uint64 val = shifted(i) + shifted(j);
      EXPECT_EQ(2, axom::utilities::popCount(val));
      EXPECT_EQ(BITS - 2, axom::utilities::popCount(~val));
    }
  }

  // Test pop count when three bits are set
  for(int i = 0; i < BITS; ++i)
  {
    for(int j = 0; j < i; ++j)
    {
      for(int k = 0; k < j; ++k)
      {
        axom::uint64 val = shifted(i) + shifted(j) + shifted(k);
        EXPECT_EQ(3, axom::utilities::popCount(val));
        EXPECT_EQ(BITS - 3, axom::utilities::popCount(~val));
      }
    }
  }

  // Test pop count on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    auto val = ::random_int<axom::uint64>();

    int bits = 0;
    for(int i = 0; i < BITS; ++i)
    {
      if(val & shifted(i)) ++bits;
    }

    EXPECT_EQ(bits, axom::utilities::popCount(val));
  }
}

TEST(core_bit_utilities, leadingZeros)
{
  constexpr axom::int32 ZERO = axom::int32(0);
  constexpr int BITS = axom::utilities::BitTraits<axom::uint32>::BITS_PER_WORD;
  ASSERT_EQ(32, BITS);

  // Axom's leadingZeros will return 32 when given 0
  EXPECT_EQ(BITS, axom::utilities::leadingZeros(ZERO));

  for(int i = 0; i < BITS; ++i)
  {
    axom::int32 val = ::shifted(i);
    EXPECT_EQ(BITS - i - 1, axom::utilities::leadingZeros(val));

    // Value doesn't change if you set bits to right of leading zero
    for(int j = 0; j < i; ++j)
    {
      axom::int32 val2 = ::shifted(i) + ::shifted(j);
      EXPECT_EQ(axom::utilities::leadingZeros(val),
                axom::utilities::leadingZeros(val2));
    }
  }

  // Test leading zeros on some random numbers
  for(int n = 0; n < NRAND_SAMPLES; ++n)
  {
    axom::int32 rand_val = ::random_int<axom::int32>();

    // Explicitly find first bit and check that it matches
    int bit = 0;
    for(; bit < BITS; ++bit)
    {
      if(rand_val & shifted(BITS - bit - 1)) break;
    }
    EXPECT_EQ(bit, axom::utilities::leadingZeros(rand_val));
  }
}
