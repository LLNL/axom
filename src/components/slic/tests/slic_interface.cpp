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


TEST(slic_interface,initialize_finalize)
{
  EXPECT_FALSE(axom::slic::isInitialized());

  // Initialize and finalize slic
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());

  // Check that we can Initialize and finalize slic a second time
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

