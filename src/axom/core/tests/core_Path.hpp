// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/core/Path.hpp"

TEST(core_Path, trivial_construct)
{
  const std::string path_str = "foo";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct)
{
  const std::string path_str = "foo/bar";
  axom::Path path(path_str);
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, two_component_construct_alt_delim)
{
  const std::string path_str = "foo:bar";
  axom::Path path(path_str, ':');
  EXPECT_EQ(static_cast<std::string>(path), path_str);
}

TEST(core_Path, join_trivial)
{
  axom::Path first("foo");
  axom::Path second("bar");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar");
}

TEST(core_Path, join_trivial_alt_delim)
{
  axom::Path first("foo");
  axom::Path second("bar");
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar");
}

TEST(core_Path, join_multi_component)
{
  axom::Path first("foo/bar");
  axom::Path second("baz/quux");
  auto joined = axom::Path::join({first, second});
  EXPECT_EQ(static_cast<std::string>(joined), "foo/bar/baz/quux");
}

TEST(core_Path, join_multi_component_alt_delim)
{
  axom::Path first("foo:bar", ':');
  axom::Path second("baz:quux", ':');
  auto joined = axom::Path::join({first, second}, ':');
  EXPECT_EQ(static_cast<std::string>(joined), "foo:bar:baz:quux");
}

TEST(core_Path, parent_basic)
{
  axom::Path path("foo/bar");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(core_Path, parent_basic_alt_delim)
{
  axom::Path path("foo;bar", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo");
}

TEST(core_Path, parent_multi_component)
{
  axom::Path path("foo/bar/baz");
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo/bar");
}

TEST(core_Path, parent_multi_component_alt_delim)
{
  axom::Path path("foo;bar;baz", ';');
  auto parent = path.parent();
  EXPECT_EQ(static_cast<std::string>(parent), "foo;bar");
}

TEST(core_Path, dirname_basename_basic)
{
  axom::Path path("foo/bar");
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_basic_alt_delim)
{
  axom::Path path("foo;bar", ';');
  EXPECT_EQ(path.dirName(), "foo");
  EXPECT_EQ(path.baseName(), "bar");
}

TEST(core_Path, dirname_basename_multi_component)
{
  axom::Path path("foo/bar/baz");
  EXPECT_EQ(path.dirName(), "foo/bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, dirname_basename_multi_component_alt_delim)
{
  axom::Path path("foo;bar;baz", ';');
  EXPECT_EQ(path.dirName(), "foo;bar");
  EXPECT_EQ(path.baseName(), "baz");
}

TEST(core_Path, basic_equals)
{
  axom::Path first("foo/bar/baz");
  axom::Path second("foo/bar/baz");
  axom::Path third("foo:bar:baz", ':');
  EXPECT_TRUE(first == second);
  EXPECT_FALSE(first != second);
  EXPECT_FALSE(first == third);
  EXPECT_TRUE(first != third);
}

TEST(core_Path, parts_empty)
{
  axom::Path first("");
  axom::Path second("/");
  axom::Path third(":", ':');

  EXPECT_EQ(first.parts().size(), 0);
  EXPECT_EQ(second.parts().size(), 0);
  EXPECT_EQ(third.parts().size(), 0);
}

TEST(core_Path, parts_single)
{
  axom::Path first("foo");
  axom::Path second("/foo");
  axom::Path third("foo:", ':');

  EXPECT_EQ(first.parts().size(), 1);
  EXPECT_EQ(first.parts()[0], "foo");
  EXPECT_EQ(second.parts().size(), 1);
  EXPECT_EQ(second.parts()[0], "foo");
  EXPECT_EQ(third.parts().size(), 1);
  EXPECT_EQ(third.parts()[0], "foo");
}

TEST(core_Path, parts_multipart)
{
  axom::Path path("foo/bar/baz");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 3);
  EXPECT_EQ(parts[0], "foo");
  EXPECT_EQ(parts[1], "bar");
  EXPECT_EQ(parts[2], "baz");
}

TEST(core_Path, parts_multipart_empty)
{
  axom::Path path("foo//bar/baz");
  auto parts = path.parts();

  EXPECT_EQ(parts.size(), 3);
  EXPECT_EQ(parts[0], "foo");
  EXPECT_EQ(parts[1], "bar");
  EXPECT_EQ(parts[2], "baz");
}
