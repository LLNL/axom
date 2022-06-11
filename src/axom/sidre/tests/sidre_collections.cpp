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

struct NamedItem
{
  NamedItem(const std::string& name) : m_name(name) { }
  const std::string& getName() const { return m_name; }

  std::string m_name;
};

template <typename CollectionType>
class ItemCollectionTest : public ::testing::Test
{
public:
  using ValueType = typename CollectionType::value_type;

protected:
  void SetUp() override { m_coll = new CollectionType; }

  void TearDown() override { delete m_coll; }

protected:
  template <typename T>
  typename std::enable_if<std::is_same<double, T>::value, double>::type
  create_item(const std::string& str)
  {
    return static_cast<double>(str.size() * 1.111);
  }

  template <typename T>
  typename std::enable_if<std::is_same<NamedItem, T>::value, NamedItem>::type
  create_item(const std::string& str)
  {
    return NamedItem(str);
  }

protected:
  sidre::ItemCollection<ValueType>* m_coll {nullptr};
};

using MyTypes =
  ::testing::Types<sidre::IndexedCollection<double>,
                   sidre::ListCollection<double>,
                   //sidre::MapCollection<double>, // note: invalid since double doesn't have required getName()
                   sidre::IndexedCollection<NamedItem>,
                   sidre::ListCollection<NamedItem>,
                   sidre::MapCollection<NamedItem>>;
TYPED_TEST_SUITE(ItemCollectionTest, MyTypes);

TYPED_TEST(ItemCollectionTest, TestEmpty)
{
  auto* coll = this->m_coll;

  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertItem)
{
  auto* coll = this->m_coll;

  auto val = this->template create_item<typename TestFixture::ValueType>("foo");

  coll->insertItem(&val, "foo");
  EXPECT_EQ(1, coll->getNumItems());
}
