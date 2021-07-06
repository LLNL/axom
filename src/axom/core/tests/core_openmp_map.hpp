// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/Map.hpp"
#include "gtest/gtest.h" /* for TEST and EXPECT_* macros */
#include <string>

namespace axom
{
namespace internal
{
template <typename Key, typename T>
experimental::Map<Key, T> init(int N, int len)
{
  experimental::Map<Key, T> test(N, len);
  EXPECT_EQ(N * len, test.max_size());
  EXPECT_EQ(0, test.size());
  EXPECT_EQ(true, test.empty());
  EXPECT_EQ(len, test.bucket_size());
  EXPECT_EQ(N, test.bucket_count());
  return test;
}

template <typename Key, typename T>
void test_storage(experimental::Map<Key, T> &test)
{
 // for(int i = 0; i < test.max_size(); i++)
 // {
 //   Key key = i;
 //   T value = key * 27;
 //   auto ret_test = test.insert(key, value);
 //   EXPECT_EQ(true, ret_test.second);
 // }
  experimental::Map<Key, T> *test2 = &test;
  axom::for_all< axom::OMP_EXEC >(0, test.max_size(), AXOM_LAMBDA(IndexType idx){
    Key key = idx;
    T value = key * 27;
    auto ret_test = test2->insert(key, value);
    EXPECT_EQ(true, ret_test.second);
  } );
  EXPECT_EQ(false, test.empty());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }
  //This should fail, since we're at capacity.
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

template <typename Key, typename T>
void test_brackets(experimental::Map<Key, T> &test)
{
  for(int i = 0; i < test.size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test[key]);
  }
}

template <typename Key, typename T>
void test_insert_assign(experimental::Map<Key, T> &test)
{
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 27;
    auto ret_test = test.insert_or_assign(key, value);
    EXPECT_EQ(true, ret_test.second);
  }

  EXPECT_EQ(false, test.empty());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }

  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 28;
    auto ret_test = test.insert_or_assign(key, value);
    EXPECT_EQ(false, ret_test.second);
    EXPECT_EQ(ret_test.first->key, key);
  }
  EXPECT_EQ(test.size(), test.max_size());
  for(int i = 0; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 28, test.find(key).value);
  }
}

template <typename Key, typename T>
void test_remove(experimental::Map<Key, T> &test)
{
  std::size_t to_erase = test.size();
  for(std::size_t i = 0; i < to_erase; i++)
  {
    Key key = (Key)i;
    bool erased = test.erase(key);
    EXPECT_EQ(erased, true);
    EXPECT_EQ(test.find(key), test.end());
  }
  EXPECT_EQ(test.size(), 0);
  test.insert(0, 900);
  auto ret = test.find(0);
  EXPECT_EQ(900, ret.value);
  EXPECT_EQ(test.size(), 1);
}

template <typename Key, typename T>
void test_rehash(experimental::Map<Key, T> &test, int num, int fact)
{
  auto original_size = test.size();
  test.rehash(num, fact);

  for(int i = 0; i < original_size; i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }

  for(int i = original_size; i < test.max_size(); i++)
  {
    Key key = i;
    T value = key * 27;
    auto ret_test = test.insert(key, value);
    EXPECT_EQ(true, ret_test.second);
    if(ret_test.second == false){
      std::cout << "Failed on step " << i << std::endl;
    }
  }

  for(int i = original_size; i < test.max_size(); i++)
  {
    Key key = i;
    EXPECT_EQ(key * 27, test.find(key).value);
  }
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

}  // namespace internal

TEST(core_map, initialization)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int> test = internal::init<int, int>(i, j);
    }
  }
}

TEST(core_map, insertion)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int> test = internal::init<int, int>(i, j);
      internal::test_storage<int, int>(test);
    }
  }
}

TEST(core_map, insert_or_assign)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int> test = internal::init<int, int>(i, j);
      internal::test_insert_assign<int, int>(test);
    }
  }
}

TEST(core_map, brackets)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int> test = internal::init<int, int>(i, j);
      internal::test_storage<int, int>(test);
      internal::test_brackets<int, int>(test);
    }
  }
}

TEST(core_map, removal)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int> test = internal::init<int, int>(i, j);
      internal::test_storage<int, int>(test);
      internal::test_remove<int, int>(test);
    }
  }
}

TEST(core_map, rehash)
{
  for(int i : {5})
  {
    for(int j : {5})
    {
      for(int k : {2})
      {
        experimental::Map<int, int> test = internal::init<int, int>(i, j);
        internal::test_storage<int, int>(test);
        internal::test_rehash<int, int>(test, -1, k);
        EXPECT_EQ(test.max_size(), k * i * j);
        internal::test_remove<int, int>(test);
      }
    }
  }
//  for(int i : {1, 2, 5, 10, 20, 100})
//  {
//    for(int j : {1, 2, 5, 10})
//    {
//      for(int k = 1; k < 4; k++)
//      {
//        experimental::Map<int, int> test = internal::init<int, int>(i, j);
//        internal::test_storage<int, int>(test);
//        internal::test_rehash<int, int>(test, test.max_size() + 20 * k, -1);
//        EXPECT_EQ((i * j + 20 * k) * j, test.max_size());
//        internal::test_remove<int, int>(test);
//      }
//    }
//  }
}
} /* namespace axom */
