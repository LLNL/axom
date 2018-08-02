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

#include "axom/slic/interface/slic.hpp"


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
