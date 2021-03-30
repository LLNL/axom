// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/utilities/Path.hpp"

TEST(utils_Path, trivial_construct)
{
  const std::string path_str = "foo";
  axom::utilities::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(utils_Path, two_component_construct)
{
  const std::string path_str = "foo/bar";
  axom::utilities::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(utils_Path, two_component_construct_alt_delim)
{
  const std::string path_str = "foo:bar";
  axom::utilities::Path path(path_str, ':');
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(utils_Path, join_trivial)
{
  axom::utilities::Path first("foo");
  axom::utilities::Path second("bar");
  auto joined = axom::utilities::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar");
}

TEST(utils_Path, join_trivial_alt_delim)
{
  axom::utilities::Path first("foo");
  axom::utilities::Path second("bar");
  auto joined = axom::utilities::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar");
}

TEST(utils_Path, join_multi_component)
{
  axom::utilities::Path first("foo/bar");
  axom::utilities::Path second("baz/quux");
  auto joined = axom::utilities::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar/baz/quux");
}

TEST(utils_Path, join_multi_component_alt_delim)
{
  axom::utilities::Path first("foo:bar", ':');
  axom::utilities::Path second("baz:quux", ':');
  auto joined = axom::utilities::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar:baz:quux");
}

TEST(utils_Path, parent_basic)
{
  axom::utilities::Path path("foo/bar");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(utils_Path, parent_basic_alt_delim)
{
  axom::utilities::Path path("foo;bar", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(utils_Path, parent_multi_component)
{
  axom::utilities::Path path("foo/bar/baz");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo/bar");
}

TEST(utils_Path, parent_multi_component_alt_delim)
{
  axom::utilities::Path path("foo;bar;baz", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo;bar");
}

TEST(utils_Path, dirname_basename_basic)
{
  axom::utilities::Path path("foo/bar");
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(utils_Path, dirname_basename_basic_alt_delim)
{
  axom::utilities::Path path("foo;bar", ';');
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(utils_Path, dirname_basename_multi_component)
{
  axom::utilities::Path path("foo/bar/baz");
  EXPECT_EQ(path.dirName(), "foo/bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(utils_Path, dirname_basename_multi_component_alt_delim)
{
  axom::utilities::Path path("foo;bar;baz", ';');
  EXPECT_EQ(path.dirName(), "foo;bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(utils_Path, basic_equals)
{
  axom::utilities::Path first("foo/bar/baz");
  axom::utilities::Path second("foo/bar/baz");
  axom::utilities::Path third("foo:bar:baz", ':');
  EXPECT_TRUE(first == second);
  EXPECT_FALSE(first != second);
  EXPECT_FALSE(first == third);
  EXPECT_TRUE(first != third);
}

TEST(utils_Path, parts_empty)
{
  axom::utilities::Path path("");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 0);
}

TEST(utils_Path, parts_single)
{
  axom::utilities::Path path("foo");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 1);
}

TEST(utils_Path, parts_multipart)
{
  axom::utilities::Path path("foo/bar/baz");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 3);
  EXPECT_EQ(parts[0], "foo");
  EXPECT_EQ(parts[1], "bar");
  EXPECT_EQ(parts[2], "baz");
}
