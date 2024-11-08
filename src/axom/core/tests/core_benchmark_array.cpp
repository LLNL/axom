// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "benchmark/benchmark.h"
#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/fmt.hpp"

#include <tuple>
#include <vector>
#include <type_traits>
#include <cstdint>

//-----------------------------------------------------------------------------
// Arguments and types to test
//-----------------------------------------------------------------------------

void CustomArgs(benchmark::internal::Benchmark* b)
{
  b->Arg(1 << 3);   // small
  b->Arg(1 << 16);  // larger than  32K L1 cache
  b->Arg(1 << 19);  // larger than 256K L2 cache
  //b->Arg(1 << 25);  // larger than  25M L3 cache
}

template <typename T>
struct Wrapper
{
  T j;
};

using Types = std::tuple<int, std::pair<int, int>, Wrapper<int>, std::string>;

//-----------------------------------------------------------------------------
// Helper functions for generating type names and values
//-----------------------------------------------------------------------------

template <typename T>
struct get_value_impl
{
  static T get(int i) { return static_cast<T>(i); }
};

template <typename T>
struct get_value_impl<std::pair<T, T>>
{
  static std::pair<T, T> get(int i)
  {
    return std::make_pair(static_cast<T>(i), static_cast<T>(i + 1));
  }
};

template <typename T>
struct get_value_impl<Wrapper<T>>
{
  static Wrapper<T> get(int i) { return Wrapper<T> {static_cast<T>(i)}; }
};

template <>
struct get_value_impl<std::string>
{
  static std::string get(int i)
  {
    return axom::fmt::format("This is a somewhat long string -- {}", i);
  }
};

template <typename T>
T get_value(int i)
{
  return get_value_impl<T>::get(i);
}

template <typename T>
std::string get_type_name();

template <>
std::string get_type_name<int>()
{
  return "int";
}

template <>
std::string get_type_name<std::pair<int, int>>()
{
  return "std::pair<int, int>";
}

template <>
std::string get_type_name<Wrapper<int>>()
{
  return "Wrapper<int>";
}

template <>
std::string get_type_name<std::string>()
{
  return "std::string";
}

//-----------------------------------------------------------------------------
// Benchmarks for array and vector construction
//-----------------------------------------------------------------------------
template <typename T>
static void Array_ctor(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    axom::Array<T> arr(size);
    assert(arr.size() == size);
    assert(arr.capacity() >= size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
static void Vector_ctor(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    std::vector<T> arr(size);
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) >= size);
    benchmark::DoNotOptimize(arr);
  }
}

//-----------------------------------------------------------------------------
// Benchmarks for array and vector push_back and emplace_back
//-----------------------------------------------------------------------------
template <typename T>
void Array_push_back_initialSize(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    axom::Array<T> arr(0, size);
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == size);

    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) == size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Array_emplace_back_initialSize(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    axom::Array<T> arr(0, size);
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == size);

    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) == size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Array_push_back_startEmpty(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    axom::Array<T> arr;
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == 0);

    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) >= size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Array_emplace_back_startEmpty(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    axom::Array<T> arr;
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == 0);
    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) >= size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Vector_push_back_initialSize(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    std::vector<T> arr;
    arr.reserve(size);
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == size);
    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) == size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Vector_emplace_back_initialSize(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    std::vector<T> arr;
    arr.reserve(size);
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == size);
    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) == size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Vector_push_back_startEmpty(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    std::vector<T> arr;
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == 0);
    for(int i = 0; i < size; ++i)
    {
      arr.push_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) >= size);
    benchmark::DoNotOptimize(arr);
  }
}

template <typename T>
void Vector_emplace_back_startEmpty(benchmark::State& state)
{
  const int size = state.range(0);
  for(auto _ : state)
  {
    std::vector<T> arr;
    assert(static_cast<int>(arr.size()) == 0);
    assert(static_cast<int>(arr.capacity()) == 0);
    for(int i = 0; i < size; ++i)
    {
      arr.emplace_back(get_value<T>(i));
    }
    assert(static_cast<int>(arr.size()) == size);
    assert(static_cast<int>(arr.capacity()) >= size);
    benchmark::DoNotOptimize(arr);
  }
}
//-----------------------------------------------------------------------------
// Register all the tests
//-----------------------------------------------------------------------------

template <typename T>
void RegisterBenchmark()
{
  auto tname = [](const std::string& n) {
    return axom::fmt::format("{}<{}>", n, get_type_name<T>());
  };

  // clang-format off
  benchmark::RegisterBenchmark(tname("Array::ctor"), &Array_ctor<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("vector::ctor"), &Vector_ctor<T>)->Apply(CustomArgs);


  benchmark::RegisterBenchmark(tname("Array::push_back_startEmpty"), &Array_push_back_startEmpty<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("Array::emplace_back_startEmpty"), &Array_emplace_back_startEmpty<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("Array::push_back_initialSize"), &Array_push_back_initialSize<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("Array::emplace_back_initialSize"), &Array_emplace_back_initialSize<T>)->Apply(CustomArgs);

  benchmark::RegisterBenchmark(tname("vector::push_back_startEmpty"), &Array_push_back_startEmpty<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("vector::emplace_back_startEmpty"), &Array_emplace_back_startEmpty<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("vector::push_back_initialSize"), &Array_push_back_initialSize<T>)->Apply(CustomArgs);
  benchmark::RegisterBenchmark(tname("vector::emplace_back_initialSize"), &Array_emplace_back_initialSize<T>)->Apply(CustomArgs);
  // clang-format on
}

//-----------------------------------------------------------------------------
// Main and helper functions to consistently register templated test types
//
// Note: This can likely be simplified/streamlined a bit when we move to C++17,
// e.g. using std::apply and paraeter pack fold expressions
//-----------------------------------------------------------------------------

template <typename Tuple, std::size_t... Is>
void RegisterBenchmarksImpl(std::index_sequence<Is...>)
{
  using expander = int[];
  (void)expander {
    0,
    (RegisterBenchmark<typename std::tuple_element<Is, Tuple>::type>(), 0)...};
}

template <typename Tuple>
void RegisterBenchmarks()
{
  RegisterBenchmarksImpl<Tuple>(
    std::make_index_sequence<std::tuple_size<Tuple>::value> {});
}

int main(int argc, char* argv[])
{
  ::benchmark::Initialize(&argc, argv);
  axom::slic::SimpleLogger logger;
  RegisterBenchmarks<Types>();
  ::benchmark::RunSpecifiedBenchmarks();
  return 0;
}
