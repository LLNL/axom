// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/nvtx/interface.hpp"

#include "gtest/gtest.h"

//------------------------------------------------------------------------------
// HELPER METHODS
//------------------------------------------------------------------------------
namespace
{
//------------------------------------------------------------------------------
void check_defaults()
{
  EXPECT_EQ(axom::nvtx::DEFAULT_COLOR, axom::nvtx::get_color());
  EXPECT_EQ(axom::nvtx::DEFAULT_CATEGORY, axom::nvtx::get_category());
}

//------------------------------------------------------------------------------
void check_set_color(axom::nvtx::Color newColor)
{
  axom::nvtx::set_color(newColor);
  EXPECT_EQ(newColor, axom::nvtx::get_color());
}

//------------------------------------------------------------------------------
void check_set_category(uint32_t newCategory)
{
  axom::nvtx::set_category(newCategory);
  EXPECT_EQ(newCategory, axom::nvtx::get_category());
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------
TEST(utils_nvtx_settings, default_settings) { check_defaults(); }

//------------------------------------------------------------------------------
TEST(utils_nvtx_settings, set_color)
{
  check_set_color(axom::nvtx::Color::BLACK);
  check_set_color(axom::nvtx::Color::BLUE);
  check_set_color(axom::nvtx::Color::YELLOW);
}

//------------------------------------------------------------------------------
TEST(utils_nvtx_settings, set_category)
{
  constexpr uint32_t MAX_CATEGORY = 42;
  for(uint32_t icategory = 1; icategory < MAX_CATEGORY; ++icategory)
  {
    check_set_category(icategory);
  }
}

//------------------------------------------------------------------------------
TEST(utils_nvtx_settings, set_and_reset)
{
  check_set_color(axom::nvtx::Color::BLUE);

  constexpr uint32_t MY_CATEGORY = 42;
  check_set_category(MY_CATEGORY);

  axom::nvtx::reset();

  check_defaults();
}
