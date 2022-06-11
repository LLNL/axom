// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/// Tests sidre's ItemCollection hierachy,
/// including MapCollection, ListCollection and IndexedCollection

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

  static constexpr bool IsNameBased =
    std::is_same<sidre::MapCollection<ValueType>, CollectionType>::value;

protected:
  void SetUp() override { m_coll = new CollectionType; }

  void TearDown() override
  {
    // Remove and deallocate items
    for(auto idx = m_coll->getFirstValidIndex(); idx != sidre::InvalidIndex;
        idx = m_coll->getNextValidIndex(idx))
    {
      auto* val = m_coll->removeItem(idx);
      delete val;
    }

    delete m_coll;
  }

  // Adds several initial items to the collection
  std::map<std::string, axom::IndexType> addItems(const std::vector<std::string>& strs)
  {
    std::map<std::string, axom::IndexType> insertionMap;

    for(auto& str : strs)
    {
      auto* val = this->template create_item<ValueType>(str);
      // elide warning about non-empty string in ListCollection
      std::string key = IsNameBased ? str : "";

      insertionMap[str] = m_coll->insertItem(val, key);
    }

    return insertionMap;
  }

protected:
  /// SFINAE function to create a double from a string
  /// Note: caller is responsible for deallocating the associated memory
  template <typename T>
  typename std::enable_if<std::is_same<double, T>::value, double*>::type
  create_item(const std::string& str) const
  {
    // creates arbitrary double from the input string's length
    return new double(str.size() * 1.111);
  }

  /// SFINAE function to create a NewItem from a string
  /// Note: caller is responsible for deallocating the associated memory
  template <typename T>
  typename std::enable_if<std::is_same<NamedItem, T>::value, NamedItem*>::type
  create_item(const std::string& str) const
  {
    return new NamedItem(str);
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

TYPED_TEST(ItemCollectionTest, testEmpty)
{
  auto* coll = this->m_coll;

  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertItem)
{
  using ValueType = typename TestFixture::ValueType;

  auto* coll = this->m_coll;

  auto* val = this->template create_item<ValueType>("foo");
  std::string key = TestFixture::IsNameBased ? "foo" : "";

  coll->insertItem(val, key);
  EXPECT_EQ(1, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, insertSeveralItems)
{
  using ValueType = typename TestFixture::ValueType;

  auto* coll = this->m_coll;

  for(auto& str : {"a", "b", "c", "aa", "bb"})
  {
    auto* val = this->template create_item<ValueType>(str);
    std::string key = TestFixture::IsNameBased ? str : "";

    coll->insertItem(val, key);
  }
  EXPECT_EQ(5, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, hasItems)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    EXPECT_TRUE(coll->hasItem(idx));
  }
}

TYPED_TEST(ItemCollectionTest, removeItems)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    EXPECT_TRUE(coll->hasItem(idx));
    auto* val = coll->removeItem(idx);
    delete val;
    EXPECT_FALSE(coll->hasItem(idx));
  }
  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, getItemsAndRemove)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // Tests gettting items -- then remove it
  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    auto* val = coll->getItem(idx);
    coll->removeItem(idx);
    delete val;
  }
  EXPECT_EQ(0, coll->getNumItems());
}

TYPED_TEST(ItemCollectionTest, getItemsAndRemoveAll)
{
  auto* coll = this->m_coll;

  std::vector<std::string> names {"a", "b", "c", "aa", "bb", "cc"};

  auto map = this->addItems(names);
  EXPECT_EQ(names.size(), coll->getNumItems());

  // Tests gettting items -- then remove it
  for(auto kv : map)
  {
    axom::IndexType idx = kv.second;
    auto* val = coll->getItem(idx);
    delete val;
  }
  EXPECT_EQ(names.size(), coll->getNumItems());

  coll->removeAllItems();
  EXPECT_EQ(0, coll->getNumItems());
}
