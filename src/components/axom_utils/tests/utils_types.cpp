/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "gtest/gtest.h"
#include <limits>

#include "axom/Types.hpp"

TEST(axom_utils_types,check_types_8)
{
  typedef axom::common::int8 SignedType;
  typedef axom::common::uint8 UnsigneType;
  static const std::size_t EXP_BYTES = 1;

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_integer);
  EXPECT_TRUE( std::numeric_limits<UnsigneType>::is_integer);

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_signed);
  EXPECT_FALSE( std::numeric_limits<UnsigneType>::is_signed);

  EXPECT_EQ( 7, std::numeric_limits<SignedType>::digits);
  EXPECT_EQ( 8, std::numeric_limits<UnsigneType>::digits);

  EXPECT_EQ(EXP_BYTES, sizeof(SignedType) );
  EXPECT_EQ(EXP_BYTES, sizeof(UnsigneType) );

}



TEST(axom_utils_types,check_types_16)
{
  typedef axom::common::int16 SignedType;
  typedef axom::common::uint16 UnsigneType;
  static const std::size_t EXP_BYTES = 2;

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_integer);
  EXPECT_TRUE( std::numeric_limits<UnsigneType>::is_integer);

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_signed);
  EXPECT_FALSE( std::numeric_limits<UnsigneType>::is_signed);

  EXPECT_EQ( 15, std::numeric_limits<SignedType>::digits);
  EXPECT_EQ( 16, std::numeric_limits<UnsigneType>::digits);

  EXPECT_EQ(EXP_BYTES, sizeof(SignedType) );
  EXPECT_EQ(EXP_BYTES, sizeof(UnsigneType) );
}


TEST(axom_utils_types,check_types_32)
{
  typedef axom::common::int32 SignedType;
  typedef axom::common::uint32 UnsigneType;
  static const std::size_t EXP_BYTES = 4;

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_integer);
  EXPECT_TRUE( std::numeric_limits<UnsigneType>::is_integer);

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_signed);
  EXPECT_FALSE( std::numeric_limits<UnsigneType>::is_signed);

  EXPECT_EQ( 31, std::numeric_limits<SignedType>::digits);
  EXPECT_EQ( 32, std::numeric_limits<UnsigneType>::digits);

  EXPECT_EQ(EXP_BYTES, sizeof(SignedType) );
  EXPECT_EQ(EXP_BYTES, sizeof(UnsigneType) );
}


TEST(axom_utils_types,check_types_64)
{
#ifndef AXOM_NO_INT64_T
  typedef axom::common::int64 SignedType;
  typedef axom::common::uint64 UnsigneType;
  static const std::size_t EXP_BYTES = 8;

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_integer);
  EXPECT_TRUE( std::numeric_limits<UnsigneType>::is_integer);

  EXPECT_TRUE( std::numeric_limits<SignedType>::is_signed);
  EXPECT_FALSE( std::numeric_limits<UnsigneType>::is_signed);

  EXPECT_EQ( 63, std::numeric_limits<SignedType>::digits);
  EXPECT_EQ( 64, std::numeric_limits<UnsigneType>::digits);

  EXPECT_EQ(EXP_BYTES, sizeof(SignedType) );
  EXPECT_EQ(EXP_BYTES, sizeof(UnsigneType) );
#else
  std::cout<<" Skipping 64-bit tests --"
           <<" 64-bit integer typedefs not defined in this configuration."
           << std::endl;

  EXPECT_TRUE(true);
#endif
}
