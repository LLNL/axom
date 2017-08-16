/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "slic/slic.hpp"


TEST(gtest_slic_interface,initialize_finalize)
{
  EXPECT_FALSE(axom::slic::isInitialized());

  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

TEST(gtest_slic_interface,initialize_finalize_twice)
{
  EXPECT_FALSE(axom::slic::isInitialized());

  // Call initialize and finalize
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());


  // Call initialize and finalize again
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

