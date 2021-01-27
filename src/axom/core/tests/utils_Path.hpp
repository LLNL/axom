// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Path.hpp"

TEST(utils_Path, trivial_construct)
{
  const std::string trivial_path_name = "foo";
  axom::utilities::Path path(trivial_path_name);
  EXPECT_EQ(static_cast<std::string>(path), trivial_path_name);
}

TEST(utils_Path, two_component_construct)
{
  const std::string two_component_path_name = "foo/bar";
  axom::utilities::Path path(two_component_path_name);
  EXPECT_EQ(static_cast<std::string>(path), two_component_path_name);
}

TEST(utils_Path, two_component_construct_alt_delim)
{
  const std::string two_component_path_name = "foo:bar";
  axom::utilities::Path path(two_component_path_name, ':');
  EXPECT_EQ(static_cast<std::string>(path), two_component_path_name);
}
