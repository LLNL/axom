// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/FlatMap.hpp"
#include "axom/core/FlatMapView.hpp"

// gtest includes
#include "gtest/gtest.h"

template <typename Key, typename Value, typename ExecSpaceType>
struct FlatMapTestParams
{
  using KeyType = Key;
  using ValueType = Value;
  using ExecSpace = ExecSpaceType;
};

template <typename FlatMapTestParams>
class core_flatmap_forall : public ::testing::Test
{
public:
  using KeyType = typename FlatMapTestParams::KeyType;
  using ValueType = typename FlatMapTestParams::ValueType;
  using ExecSpace = typename FlatMapTestParams::ExecSpace;

  using MapType = axom::FlatMap<KeyType, ValueType>;
  using MapViewType = typename MapType::View;

  template <typename T>
  KeyType getKey(T input)
  {
    return (KeyType)input;
  }

  template <typename T>
  ValueType getValue(T input)
  {
    return (ValueType)input;
  }

  ValueType getDefaultValue() { return ValueType(); }
};

using ViewTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  FlatMapTestParams<int, double, axom::OMP_EXEC>,
#endif
  FlatMapTestParams<int, double, axom::SEQ_EXEC>>;

TYPED_TEST_SUITE(core_flatmap_forall, ViewTypes);

AXOM_TYPED_TEST(core_flatmap_forall, insert_and_find)
{
  using MapType = typename TestFixture::MapType;
  using MapViewType = typename TestFixture::MapViewType;
  using ExecSpace = typename TestFixture::ExecSpace;

  MapType test_map;

  const int NUM_ELEMS = 100;
  const int EXTRA_THREADS = 100;

  // First do insertions of elements.
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    test_map.insert({key, value});
  }

  MapViewType test_map_view(test_map);

  axom::Array<int> valid_vec(NUM_ELEMS + EXTRA_THREADS, 0);
  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  const auto valid_out = valid_vec.view();
  const auto keys_out = keys_vec.view();
  const auto values_out = values_vec.view();

  // Read values out in a captured lambda.
  axom::for_all<ExecSpace>(
    NUM_ELEMS + EXTRA_THREADS,
    AXOM_LAMBDA(axom::IndexType idx) {
      auto it = test_map_view.find(idx);
      if(it != test_map_view.end())
      {
        keys_out[idx] = it->first;
        values_out[idx] = it->second;
        valid_out[idx] = true;
      }
    });

  // Check contents on the host
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    EXPECT_EQ(valid_out[i], true);
    EXPECT_EQ(keys_out[i], this->getKey(i));
    EXPECT_EQ(values_out[i], this->getValue(i * 10.0 + 5.0));
  }
  for(int i = NUM_ELEMS; i < NUM_ELEMS + EXTRA_THREADS; i++)
  {
    EXPECT_EQ(valid_out[i], false);
  }
}

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  MapType test_map;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Construct a flat map with the key-value pairs.
  test_map = MapType::template create<ExecSpace>(keys_vec, values_vec);

  // Check contents on the host
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val = this->getValue(i * 10.0 + 5.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    EXPECT_EQ(expected_val, test_map.at(expected_key));
  }
}

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched_with_dups)
{
  using MapType = typename TestFixture::MapType;
  using ExecSpace = typename TestFixture::ExecSpace;

  MapType test_map;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS * 2);
  axom::Array<double> values_vec(NUM_ELEMS * 2);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Add some duplicate key values
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 7.0);

    keys_vec[i + NUM_ELEMS] = key;
    values_vec[i + NUM_ELEMS] = value;
  }

  // Construct a flat map with the key-value pairs.
  test_map = MapType::template create<ExecSpace>(keys_vec, values_vec);

  // Check contents on the host. Only one of the duplicate keys should remain.
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val1 = this->getValue(i * 10.0 + 5.0);
    auto expected_val2 = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    EXPECT_TRUE((test_map.at(expected_key) == expected_val1) ||
                (test_map.at(expected_key) == expected_val2));
  }

  // Check that we only have one instance of every key in the map
  axom::Array<std::pair<int, double>> kv_out(NUM_ELEMS);
  int index = 0;
  for(auto &pair : test_map)
  {
    EXPECT_LT(index, NUM_ELEMS);
    kv_out[index++] = {pair.first, pair.second};
  }

  std::sort(kv_out.begin(),
            kv_out.end(),
            [](const std::pair<int, double> &first,
               const std::pair<int, double> &second) -> bool {
              return first.first < second.first;
            });

  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val1 = this->getValue(i * 10.0 + 5.0);
    auto expected_val2 = this->getValue(i * 10.0 + 7.0);
    EXPECT_EQ(kv_out[i].first, expected_key);
    EXPECT_TRUE((kv_out[i].second == expected_val1) ||
                (kv_out[i].second == expected_val2));
  }
}

template <typename KeyType>
struct ConstantHash
{
  using argument_type = KeyType;
  using result_type = axom::IndexType;

  AXOM_HOST_DEVICE axom::IndexType operator()(KeyType) const { return 0; }
};

AXOM_TYPED_TEST(core_flatmap_forall, insert_batched_constant_hash)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KeyType = typename TestFixture::KeyType;
  using ValueType = typename TestFixture::ValueType;

  using MapType = axom::FlatMap<KeyType, ValueType, ConstantHash<KeyType>>;

  const int NUM_ELEMS = 100;

  axom::Array<int> keys_vec(NUM_ELEMS);
  axom::Array<double> values_vec(NUM_ELEMS);
  // Create batch of array elements
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto key = this->getKey(i);
    auto value = this->getValue(i * 10.0 + 5.0);

    keys_vec[i] = key;
    values_vec[i] = value;
  }

  // Construct a flat map with the key-value pairs.
  MapType test_map = MapType::template create<ExecSpace>(keys_vec, values_vec);

  // Check contents on the host
  EXPECT_EQ(NUM_ELEMS, test_map.size());

  // Check that every element we inserted is in the map
  for(int i = 0; i < NUM_ELEMS; i++)
  {
    auto expected_key = this->getKey(i);
    auto expected_val = this->getValue(i * 10.0 + 5.0);
    EXPECT_EQ(1, test_map.count(expected_key));
    EXPECT_EQ(expected_val, test_map.at(expected_key));
  }
}
