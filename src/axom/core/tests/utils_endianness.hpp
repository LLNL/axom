// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <limits>
#include <iostream>

TEST(utils_endianness, print_endianness)
{
  std::cout << "Checking system endianness..." << std::endl;

  const bool isLittleEndian = axom::utilities::isLittleEndian();

  std::cout << "System is " << (isLittleEndian ? "little" : "big") << " endian."
            << std::endl;
}

TEST(utils_endianness, endianness_16)
{
  union SixteenBit
  {
    std::uint8_t raw[2];
    std::int16_t i_val;
    std::uint16_t ui_val;
    short short_val;
    unsigned short ushort_val;
  };

  SixteenBit valOrig = {{0, 1}};
  SixteenBit valSwap = {{1, 0}};

  // Test short
  {
    short v1 = valOrig.short_val;
    short v2 = axom::utilities::swapEndian(v1);
    short v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.short_val, v1);
    EXPECT_EQ(valSwap.short_val, v2);
    EXPECT_EQ(valOrig.short_val, v3);
  }

  // Test unsigned short
  {
    unsigned short v1 = valOrig.ushort_val;
    unsigned short v2 = axom::utilities::swapEndian(v1);
    unsigned short v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.ushort_val, v1);
    EXPECT_EQ(valSwap.ushort_val, v2);
    EXPECT_EQ(valOrig.ushort_val, v3);
  }

  // Test int16
  {
    std::int16_t v1 = valOrig.i_val;
    std::int16_t v2 = axom::utilities::swapEndian(v1);
    std::int16_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.i_val, v1);
    EXPECT_EQ(valSwap.i_val, v2);
    EXPECT_EQ(valOrig.i_val, v3);
  }

  // Test uint16
  {
    std::uint16_t v1 = valOrig.i_val;
    std::uint16_t v2 = axom::utilities::swapEndian(v1);
    std::uint16_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.ui_val, v1);
    EXPECT_EQ(valSwap.ui_val, v2);
    EXPECT_EQ(valOrig.ui_val, v3);
  }
}

TEST(utils_endianness, endianness_32)
{
  union ThirtyTwoBit
  {
    std::uint8_t raw[4];
    std::int32_t i_val;
    std::uint32_t ui_val;
    int int_val;
    unsigned int uint_val;
    float f_val;
  };

  ThirtyTwoBit valOrig = {{0, 1, 2, 3}};
  ThirtyTwoBit valSwap = {{3, 2, 1, 0}};

  // Test int
  {
    int v1 = valOrig.int_val;
    int v2 = axom::utilities::swapEndian(v1);
    int v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.int_val, v1);
    EXPECT_EQ(valSwap.int_val, v2);
    EXPECT_EQ(valOrig.int_val, v3);
  }

  // Test unsigned int
  {
    unsigned int v1 = valOrig.uint_val;
    unsigned int v2 = axom::utilities::swapEndian(v1);
    unsigned int v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.uint_val, v1);
    EXPECT_EQ(valSwap.uint_val, v2);
    EXPECT_EQ(valOrig.uint_val, v3);
  }

  // Test int32
  {
    std::int32_t v1 = valOrig.i_val;
    std::int32_t v2 = axom::utilities::swapEndian(v1);
    std::int32_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.i_val, v1);
    EXPECT_EQ(valSwap.i_val, v2);
    EXPECT_EQ(valOrig.i_val, v3);
  }

  // Test uint32
  {
    std::uint32_t v1 = valOrig.i_val;
    std::uint32_t v2 = axom::utilities::swapEndian(v1);
    std::uint32_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.ui_val, v1);
    EXPECT_EQ(valSwap.ui_val, v2);
    EXPECT_EQ(valOrig.ui_val, v3);
  }

  // Test float
  {
    float v1 = valOrig.f_val;
    float v2 = axom::utilities::swapEndian(v1);
    float v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.f_val, v1);
    EXPECT_EQ(valSwap.f_val, v2);
    EXPECT_EQ(valOrig.f_val, v3);
  }
}

TEST(utils_endianness, endianness_64)
{
  union SixtyFourBit
  {
    std::uint8_t raw[8];
#ifndef AXOM_NO_INT64_T
    std::int64_t i_val;
    std::uint64_t ui_val;
#endif
    double d_val;
  };

  SixtyFourBit valOrig = {{0, 1, 2, 3, 4, 5, 6, 7}};
  SixtyFourBit valSwap = {{7, 6, 5, 4, 3, 2, 1, 0}};

#ifndef AXOM_NO_INT64_T
  // Test int64
  {
    std::int64_t v1 = valOrig.i_val;
    std::int64_t v2 = axom::utilities::swapEndian(v1);
    std::int64_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.i_val, v1);
    EXPECT_EQ(valSwap.i_val, v2);
    EXPECT_EQ(valOrig.i_val, v3);
  }

  // Test uint64
  {
    std::uint64_t v1 = valOrig.i_val;
    std::uint64_t v2 = axom::utilities::swapEndian(v1);
    std::uint64_t v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.ui_val, v1);
    EXPECT_EQ(valSwap.ui_val, v2);
    EXPECT_EQ(valOrig.ui_val, v3);
  }
#endif

  // Test double
  {
    double v1 = valOrig.d_val;
    double v2 = axom::utilities::swapEndian(v1);
    double v3 = axom::utilities::swapEndian(v2);

    EXPECT_EQ(valOrig.d_val, v1);
    EXPECT_EQ(valSwap.d_val, v2);
    EXPECT_EQ(valOrig.d_val, v3);
  }
}

#if 0
// Note: Checks that swapEndian static_asserts for types within invalid sizes.
//       Commented out since it should lead to a compile error
TEST(utils_endianness,invalid_byte_width)
{
  {
    std::int8_t v1 = 5;
    std::int8_t v2 = axom::utilities::swapEndian(v1);

    EXPECT_EQ(v1, v2);
  }

  {
    std::uint8_t v1 = 5;
    std::uint8_t v2 = axom::utilities::swapEndian(v1);

    EXPECT_EQ(v1, v2);
  }
}
#endif

#if 0
// Note: Checks that swapEndian static_asserts for non-native types.
//       Commented out since it should lead to a compile error.
TEST(utils_endianness,invalid_non_native_types)
{
  struct AxomUtilsTestsNonNative
  {
    std::uint16_t a;
    std::uint16_t b;
  };

  AxomUtilsTestsNonNative v1;
  v1.a = 5;
  v1.b = 12;

  AxomUtilsTestsNonNative v2 = axom::utilities::swapEndian(v1);

  EXPECT_NE(v1.a, v2.a);
  EXPECT_NE(v1.b, v2.b);
}
#endif
