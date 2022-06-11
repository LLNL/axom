// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/// Tests sidre's ItemCollection hierachy,
/// including MapCollection, ListCollection and IndexexCollection

#include "gtest/gtest.h"

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/sidre.hpp"

namespace sidre = axom::sidre;

template <typename CollectionType>
class ItemCollectionTest : public ::testing::Test
{
protected:
  void SetUp() override { m_coll = new CollectionType; }

  void TearDown() override { delete m_coll; }

protected:
  sidre::ItemCollection<double>* m_coll {nullptr};
};

using MyTypes =
  ::testing::Types<sidre::IndexedCollection<double>, sidre::ListCollection<double>>;
TYPED_TEST_SUITE(ItemCollectionTest, MyTypes);

TYPED_TEST(ItemCollectionTest, TestEmpty)
{
  auto* coll = this->m_coll;

  EXPECT_EQ(0, coll->getNumItems());
}

