// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/Map.hpp"
#include "gtest/gtest.h" /* for TEST and EXPECT_* macros */
#include <string>
#include "axom/core/utilities/Utilities.hpp"
#include "axom/core/Types.hpp"
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <utility>

namespace axom
{
namespace internal
{
template <typename Key, typename T, typename Hash, typename Policy>
experimental::Map<Key, T, Hash, Policy> init(int N, int len)
{
  experimental::Map<Key, T, Hash, Policy> test(N, len);
  EXPECT_EQ(N * len, test.max_size());
  EXPECT_EQ(0, test.size());
  EXPECT_EQ(true, test.empty());
  EXPECT_EQ(len, test.bucket_size());
  EXPECT_EQ(N, test.bucket_count());
  return test;
}

template <typename Key, typename T, typename Hash, typename Policy>
void test_storage(experimental::Map<Key, T, Hash, Policy> &test)
{
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  axom::for_all<Policy>(
    0,
    test.max_size(),
    AXOM_LAMBDA(IndexType idx) {
      Key key = idx;
      T value = key * 27;
      auto ret_test = test2->insert(key, value);
      EXPECT_EQ(true, ret_test.second);
    });
  EXPECT_EQ(false, test.empty());
  axom::for_all<Policy>(
    0,
    test.max_size(),
    AXOM_LAMBDA(IndexType idx) {
      Key key = idx;
      EXPECT_EQ(key * 27, test2->find(key).value);
    });
  //This should fail, since we're at capacity.
  auto ret = test.insert(test.max_size(), 900);
  EXPECT_EQ(false, ret.second);
}

template <typename Key, typename T, typename Hash, typename Policy>
void insert_perf(experimental::Map<Key, T, Hash, Policy> &test,
                 Key *keys,
                 T *vals,
                 unsigned int count)
{
  experimental::Map<Key, T, Hash, Policy> *test2 = &test;
  axom::for_all<Policy>(
    0,
    count,
    AXOM_LAMBDA(IndexType idx) {
      Key key = keys[idx];
      T value = vals[idx];
      //Key key = idx;
      //T value = key+2;
      auto ret_test = test2->insert(key, value);
//      EXPECT_EQ(true, ret_test.second);
    });
  EXPECT_EQ(false, test.empty());
 // axom::for_all<Policy>(
  //  0,
 //   count,
 //   AXOM_LAMBDA(IndexType idx) {
 //     Key key = keys[idx];
  //    EXPECT_EQ(vals[idx], test2->find(key).value);
 //   });
  //This should fail, since we're at capacity.
  //auto ret = test.insert(test.max_size(), 900);
  //EXPECT_EQ(false, ret.second);
}

template <typename Key, typename T, typename Hash, typename Policy>
void insert_perf_baseline(Key *keys,
                 T *vals,
                 unsigned int count)
{
  std::unordered_map<unsigned int, unsigned int> test(50000000);
  for(unsigned int i = 0; i < count; i++){
      Key key = keys[i];
      T value = vals[i];
      //Key key = idx;
      //T value = key+2;
      test[key] =value;
//      EXPECT_EQ(true, ret_test.second);
    }
}

template <typename Key, typename T>
void gen_data(Key *&keys, T *&vals, unsigned int &count)
{
  //Key * data = axom::allocate<Key>(count, axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType::Unified));
  std::vector<unsigned int> data;
  std::mt19937 gen;
  for(unsigned int i = 0; i < count; i++)
  {
    data.push_back(gen());
  }
  std::sort(data.begin(), data.end());
  std::cout << "Size before unique: " << data.size() << std::endl;
  auto to_del = std::unique(data.begin(), data.end()); 
  std::cout << "Size after unique: " << data.size() << std::endl;
  data.erase(to_del, data.end());
  std::random_shuffle(data.begin(), data.end());
  count = data.size();
  keys = axom::allocate<Key>(count);
  axom::copy(keys, data.data(), count*sizeof(unsigned int));
  vals = axom::allocate<T>(count);
  for(unsigned int i = 0; i < count; i++)
  {
    vals[i] = gen();
  }
}

template <typename Key>
void terrible_test(Key * keys, unsigned int count)
{
  for(unsigned int i = 0; i < count; i++){
    for(unsigned int j = i+1; j < count; j++){
      if(keys[i] == keys[j]){
        std::cout << "LITERALLY THE WORST" << std::endl;
      }
    }
  }
}

}  // namespace internal

TEST(core_map, initialization)
{
  for(int i : {1, 2, 5, 10, 20, 100})
  {
    for(int j : {1, 2, 5, 10})
    {
      experimental::Map<int, int, std::hash<int>, axom::OMP_EXEC> test =
        internal::init<int, int, std::hash<int>, axom::OMP_EXEC>(i, j);
    }
  }
}

TEST(core_map, baseline)
{
  for(unsigned int i : {50000000})
  {
    for(unsigned int j : {100})
    {
      for(float load : {1.03})
      {
        unsigned int *keys;
        unsigned int *vals;
        unsigned int count = i;
        std::cout << "Entering generator" << std::endl;
        internal::gen_data<unsigned int, unsigned int>(keys, vals, count);
        std::cout << "Exiting generator" << std::endl;
        //internal::terrible_test<unsigned int>(keys, count);
        internal::insert_perf_baseline<unsigned int, unsigned int, std::hash<unsigned int>, axom::OMP_EXEC>(keys, vals, count);
      }
    }
  }

}

TEST(core_map, omp_insertion_perf)
{
  for(unsigned int i : {50000000})
  {
    for(unsigned int j : {100})
    {
      for(float load : {1.03})
      {
        unsigned int *keys;
        unsigned int *vals;
        unsigned int count = i;
        std::cout << "Entering generator" << std::endl;
        internal::gen_data<unsigned int, unsigned int>(keys, vals, count);
        std::cout << "Exiting generator" << std::endl;
        //internal::terrible_test<unsigned int>(keys, count);
        experimental::Map<unsigned int, unsigned int, std::hash<unsigned int>, axom::OMP_EXEC> test =
          internal::init<unsigned int, unsigned int, std::hash<unsigned int>, axom::OMP_EXEC>((unsigned int) ceil(count*load), j);
        internal::insert_perf<unsigned int, unsigned int, std::hash<unsigned int>, axom::OMP_EXEC>(test, keys, vals, count);
      }
    }
  }
}

} /* namespace axom */
