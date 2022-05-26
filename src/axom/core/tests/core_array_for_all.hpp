// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/Array.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/synchronize.hpp"
#include "axom/core/execution/for_all.hpp"

// gtest includes
#include "gtest/gtest.h"

namespace testing
{
//------------------------------------------------------------------------------
//  This test harness defines some types that are useful for the tests below
//------------------------------------------------------------------------------
template <typename TheExecSpace>
class core_array_for_all : public ::testing::Test
{
public:
  using ExecSpace = TheExecSpace;

  // Define some memory spaces
#ifdef AXOM_USE_UMPIRE
  static constexpr axom::MemorySpace host_memory = axom::MemorySpace::Host;
#else
  static constexpr axom::MemorySpace host_memory = axom::MemorySpace::Dynamic;
#endif

  static constexpr axom::MemorySpace exec_space_memory =
    axom::execution_space<ExecSpace>::memory_space;

  // Define some Array type aliases
  template <typename T>
  using HostTArray = axom::Array<T, 1, host_memory>;
  template <typename T>
  using DynamicTArray = axom::Array<T, 1, axom::MemorySpace::Dynamic>;
  using HostArray = axom::Array<int, 1, host_memory>;
  using DynamicArray = axom::Array<int, 1, axom::MemorySpace::Dynamic>;
  using KernelArray = axom::Array<int, 1, exec_space_memory>;
  using KernelArrayView = axom::ArrayView<int, 1, exec_space_memory>;
};

// Generate a list of available execution types
using MyTypes = ::testing::Types<
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_OPENMP)
  axom::OMP_EXEC,
#endif
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  axom::CUDA_EXEC<100>,
  axom::CUDA_EXEC<256>,
  axom::CUDA_EXEC<256, axom::ASYNC>,
#endif
  axom::SEQ_EXEC>;

TYPED_TEST_SUITE(core_array_for_all, MyTypes);

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, explicit_ArrayView)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KernelArray = typename TestFixture::KernelArray;
  using KernelArrayView = typename TestFixture::KernelArrayView;
  using HostArray = typename TestFixture::HostArray;

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr int N = 374;
  KernelArray arr(N);

  // Modify array using lambda and ArrayView
  KernelArrayView arr_view(arr);
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_view[idx] = N - idx; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr = arr;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, auto_ArrayView)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KernelArray = typename TestFixture::KernelArray;
  using HostArray = typename TestFixture::HostArray;

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr int N = 374;
  KernelArray arr(N);

  // Modify array using lambda and ArrayView
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_view[idx] = N - idx; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr = arr;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, auto_ArrayView_const)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using KernelArray = typename TestFixture::KernelArray;
  using HostArray = typename TestFixture::HostArray;

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr int N = 374;
  KernelArray arr(N);

  // Populate an array on the host...
  HostArray source(N);
  for(int i = 0; i < N; ++i)
  {
    source[i] = i;
  }

  // Then copy it over to the device
  KernelArray kernelSource = source;
  auto kernelSourceView = kernelSource.view();

  // First, modify array using lambda and KernelArray::ArrayView operator[] const
  auto arrData = arr.data();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) {
      arrData[idx] = N - kernelSourceView[idx];
    });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr = arr;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }

  // Then modify array using lambda and KernelArray::ConstArrayView operator[] const
  KernelArray arrConst(N);
  auto arrConstData = arrConst.data();

  const KernelArray& kernelSourceCref = kernelSource;
  auto kernelSourceConstView = kernelSourceCref.view();

  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) {
      arrConstData[idx] = N - kernelSourceConstView[idx];
    });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArrConst = arrConst;
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, dynamic_array)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  // Test that array has been initialized to default values (0 for int)
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N; ++i)
    {
      int default_value {0};
      EXPECT_EQ(localArr[i], default_value);
    }
  }

  // Modify array using lambda and ArrayView
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_view[idx] = N - idx; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(localArr[i], N - i);
    }
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, dynamic_array_insert)
{
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;
  using ExecSpace = typename TestFixture::ExecSpace;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif

  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);
  auto arr_v = arr.view();

  // Set some elements
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_v[idx] = idx - 5 * idx + 7; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  const axom::IndexType N_insert = 10;
  // Let's push back some elements
  for(axom::IndexType i = 0; i < N_insert; i++)
  {
    arr.push_back(100 + i);
  }
  EXPECT_EQ(arr.size(), N + N_insert);
  {
    // Check elements
    HostArray host_arr = arr;
    for(axom::IndexType i = 0; i < N; i++)
    {
      EXPECT_EQ(host_arr[i], i - 5 * i + 7);
    }
    for(axom::IndexType i = 0; i < N_insert; i++)
    {
      EXPECT_EQ(host_arr[i + N], 100 + i);
    }
  }

  // Do the same, but with an iterator at the beginning
  for(axom::IndexType i = 0; i < N_insert; i++)
  {
    arr.insert(arr.begin(), 200 + i);
  }
  EXPECT_EQ(arr.size(), N + N_insert * 2);
  {
    // Check elements
    HostArray host_arr = arr;
    for(axom::IndexType i = 0; i < N; i++)
    {
      EXPECT_EQ(host_arr[i + N_insert], i - 5 * i + 7);
    }
    for(axom::IndexType i = 0; i < N_insert; i++)
    {
      EXPECT_EQ(host_arr[i + N + N_insert], 100 + i);
      // elements pushed at beginning are pushed in reverse order
      EXPECT_EQ(host_arr[i], 200 + N_insert - i - 1);
    }
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, dynamic_array_range_insert)
{
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;
  using ExecSpace = typename TestFixture::ExecSpace;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif

  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);
  auto arr_v = arr.view();

  // Set some elements
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_v[idx] = idx - 5 * idx + 7; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Create a host range to set
  const axom::IndexType N_range = 100;
  HostArray range_vals(N_range);

  for(axom::IndexType i = 0; i < N_range; i++)
  {
    range_vals[i] = 100 + i;
  }
  // Insert range at the end
  arr.insert(arr.end(), range_vals.size(), range_vals.data());
  EXPECT_EQ(arr.size(), N + N_range);
  {
    // Check elements
    HostArray host_arr = arr;
    for(axom::IndexType i = 0; i < N; i++)
    {
      EXPECT_EQ(host_arr[i], i - 5 * i + 7);
    }
    for(axom::IndexType i = 0; i < N_range; i++)
    {
      EXPECT_EQ(host_arr[i + N], 100 + i);
    }
  }

  // Insert range at beginning
  arr.insert(arr.begin(), range_vals.size(), range_vals.data());
  EXPECT_EQ(arr.size(), N + N_range * 2);
  {
    // Check elements
    HostArray host_arr = arr;
    for(axom::IndexType i = 0; i < N; i++)
    {
      EXPECT_EQ(host_arr[i + N_range], i - 5 * i + 7);
    }
    for(axom::IndexType i = 0; i < N_range; i++)
    {
      EXPECT_EQ(host_arr[i + N + N_range], 100 + i);
      // elements pushed at beginning are pushed in reverse order
      EXPECT_EQ(host_arr[i], 100 + i);
    }
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, dynamic_array_range_set)
{
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;
  using ExecSpace = typename TestFixture::ExecSpace;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif

  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);
  auto arr_v = arr.view();

  // Set some elements
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_v[idx] = idx - 5 * idx + 7; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Create a host range to set
  const axom::IndexType N_range = 100;
  HostArray range_vals(N_range);

  for(axom::IndexType i = 0; i < N_range; i++)
  {
    range_vals[i] = 100 + i;
  }
  // Overwrite first N_range elements in beginning of range
  arr.set(range_vals.data(), N_range, 0);
  EXPECT_EQ(arr.size(), N);
  {
    // Check elements
    HostArray host_arr = arr;
    for(axom::IndexType i = N_range; i < N; i++)
    {
      EXPECT_EQ(host_arr[i], i - 5 * i + 7);
    }
    for(axom::IndexType i = 0; i < N_range; i++)
    {
      EXPECT_EQ(host_arr[i], 100 + i);
    }
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, dynamic_array_initializer_list)
{
  using DynamicArray = typename TestFixture::DynamicArray;
  using HostArray = typename TestFixture::HostArray;
  using ExecSpace = typename TestFixture::ExecSpace;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif

  // Construct array with an initializer list
  {
    DynamicArray arr({1, 2, 3, 4, 5}, kernelAllocID);
    EXPECT_EQ(arr.size(), 5);
    EXPECT_EQ(arr.capacity(), 5);
    EXPECT_EQ(arr.getAllocatorID(), kernelAllocID);

    HostArray arr_host = arr;
    for(axom::IndexType i = 0; i < 5; i++)
    {
      EXPECT_EQ(arr_host[i], i + 1);
    }
  }

  // Assign an initializer list to an array
  {
    constexpr axom::IndexType N = 10;
    DynamicArray arr(N, N, kernelAllocID);
    arr.fill(6);

    EXPECT_EQ(arr.size(), 10);

    arr = {1, 2, 3, 4, 5};
    // after assignment, array should contain just the initializer list values
    EXPECT_EQ(arr.size(), 5);
    EXPECT_EQ(arr.capacity(), 10);
    EXPECT_EQ(arr.getAllocatorID(), kernelAllocID);
    HostArray arr_host = arr;
    for(axom::IndexType i = 0; i < 5; i++)
    {
      EXPECT_EQ(arr_host[i], i + 1);
    }
  }
}

//------------------------------------------------------------------------------
constexpr int MAGIC_DEFAULT_CTOR = 222;

struct NonTrivialDefaultCtor
{
  NonTrivialDefaultCtor() = default;

  int m_val {MAGIC_DEFAULT_CTOR};
};

AXOM_TYPED_TEST(core_array_for_all, nontrivial_default_ctor_obj)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialDefaultCtor>;
  using HostArray =
    typename TestFixture::template HostTArray<NonTrivialDefaultCtor>;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  // Check array contents on host
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(localArr[i].m_val, MAGIC_DEFAULT_CTOR);
    }
  }

  const int MAGIC_PREFILL = 111;
  // Fill with placeholder value
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_view[idx].m_val = MAGIC_PREFILL; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Fill with instance of copy-constructed type
  arr.fill(NonTrivialDefaultCtor {});

  // Check array contents on host
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(localArr[i].m_val, MAGIC_DEFAULT_CTOR);
    }
  }

  // Resize array, adding new default-constructed elements
  constexpr axom::IndexType N2 = 500;
  arr.resize(N2);

  // Check array contents on host
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N2; ++i)
    {
      EXPECT_EQ(localArr[i].m_val, MAGIC_DEFAULT_CTOR);
    }
  }
}

//------------------------------------------------------------------------------
struct NonTrivialCtor
{
  NonTrivialCtor(int val) : m_val(val) { }

  int m_val;
};

AXOM_TYPED_TEST(core_array_for_all, nontrivial_ctor_obj)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialCtor>;
  using HostArray = typename TestFixture::template HostTArray<NonTrivialCtor>;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  const int MAGIC_PREFILL = 111;
  const int MAGIC_FILL = 555;
  // Fill with placeholder value
  auto arr_view = arr.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType idx) { arr_view[idx].m_val = MAGIC_PREFILL; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Fill with instance of copy-constructed type
  arr.fill(NonTrivialCtor {MAGIC_FILL});

  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr(arr, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i].m_val, MAGIC_FILL);
  }
}

//------------------------------------------------------------------------------
struct NonTrivialDtor
{
  constexpr static int MAGIC_DTOR {555};

  ~NonTrivialDtor()
  {
    m_val = MAGIC_DTOR;
    NonTrivialDtor::dtor_calls++;
  }

  int m_val {0};

  static int dtor_calls;
};

constexpr int NonTrivialDtor::MAGIC_DTOR;
int NonTrivialDtor::dtor_calls {0};

AXOM_TYPED_TEST(core_array_for_all, nontrivial_dtor_obj)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialDtor>;
  using HostArray = typename TestFixture::template HostTArray<NonTrivialDtor>;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  NonTrivialDtor::dtor_calls = 0;
  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  // Initialization should not invoke the destructor
  EXPECT_EQ(NonTrivialDtor::dtor_calls, 0);

  NonTrivialDtor::dtor_calls = 0;
  // Array::resize(N-100) should invoke the destructor on 100 elements
  arr.resize(N - 100);
  EXPECT_EQ(NonTrivialDtor::dtor_calls, 100);

  // Resize to original size - this should leave destructed memory uninitialized
  arr.resize(axom::ArrayOptions::Uninitialized {}, N);

  // Check array contents on host
  {
    HostArray localArr(arr, hostAllocID);
    // Non-destructed elements should be in original state
    for(int i = 0; i < N - 100; i++)
    {
      EXPECT_EQ(localArr[i].m_val, 0);
    }
    // Destroyed elements should be set to magic flag value
    for(int i = N - 100; i < N; i++)
    {
      EXPECT_EQ(localArr[i].m_val, NonTrivialDtor::MAGIC_DTOR);
    }
  }

  // Reset objects in array
  arr.fill(NonTrivialDtor {});

  // Array::clear should invoke the destructor N times
  NonTrivialDtor::dtor_calls = 0;
  arr.clear();
  EXPECT_EQ(NonTrivialDtor::dtor_calls, N);

  // Resize to original size - this should leave destructed memory uninitialized
  arr.resize(axom::ArrayOptions::Uninitialized {}, N);

  // All elements should be in the destructed state
  {
    HostArray localArr(arr, hostAllocID);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(localArr[i].m_val, NonTrivialDtor::MAGIC_DTOR);
    }
  }
}

//------------------------------------------------------------------------------
constexpr static int MAGIC_COPY_CTOR {333};
struct NonTrivialCopyCtor
{
  NonTrivialCopyCtor() { }
  NonTrivialCopyCtor(const NonTrivialCopyCtor&) { m_val *= MAGIC_COPY_CTOR; }
  NonTrivialCopyCtor& operator=(const NonTrivialCopyCtor&)
  {
    m_val *= MAGIC_COPY_CTOR;
    return *this;
  }
  NonTrivialCopyCtor(NonTrivialCopyCtor&& other) = default;
  NonTrivialCopyCtor& operator=(NonTrivialCopyCtor&& other) = default;
  ~NonTrivialCopyCtor() = default;

  int m_val {1};
};

AXOM_TYPED_TEST(core_array_for_all, nontrivial_copy_ctor_obj)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialCopyCtor>;
  using IntArray = typename TestFixture::DynamicArray;
  using IntHostArray = typename TestFixture::HostArray;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif
  int hostAllocID = axom::execution_space<axom::SEQ_EXEC>::allocatorID();

  // Helper function to check all values in the array for consistency
  auto check_array_values = [=](const DynamicArray& arr, int expected) -> bool {
    // Copy device side values into int array
    IntArray values(arr.size(), arr.size(), kernelAllocID);
    const auto values_v = values.view();
    const auto arr_v = arr.view();
    axom::for_all<ExecSpace>(
      arr.size(),
      AXOM_LAMBDA(axom::IndexType i) { values_v[i] = arr_v[i].m_val; });

    // handles synchronization, if necessary
    if(axom::execution_space<ExecSpace>::async())
    {
      axom::synchronize<ExecSpace>();
    }

    // Check array contents on host
    IntHostArray values_host(values, hostAllocID);
    bool allValuesEq = true;
    for(int i = 0; i < arr.size(); ++i)
    {
      allValuesEq = allValuesEq && (values_host[i] == expected);
    }
    return allValuesEq;
  };

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 374;
  DynamicArray arr(N, N, kernelAllocID);

  // Check some cases where the copy constructor shouldn't be invoked
  {
    // Create an array of N items using default MemorySpace for ExecSpace
    constexpr axom::IndexType N = 374;
    DynamicArray arr(N, N, kernelAllocID);

    // Array default-construction should not invoke copy constructor
    EXPECT_TRUE(check_array_values(arr, 1));

    // Array resize should not invoke copy constructor
    arr.resize(2 * N);
    EXPECT_TRUE(check_array_values(arr, 1));

    // Array move-construction should not invoke copy constructor
    DynamicArray arr2(std::move(arr));
    EXPECT_TRUE(check_array_values(arr2, 1));

    // Array move-assignment should not invoke copy constructor
    DynamicArray arr3 = std::move(arr2);
    EXPECT_TRUE(check_array_values(arr3, 1));
  }

  // Check some cases where the copy constructor should be invoked
  {
    // Create an array of N items using default MemorySpace for ExecSpace
    constexpr axom::IndexType N = 374;
    DynamicArray arr(N, N, kernelAllocID);

    // Array copy-construction should invoke each element's copy constructor
    DynamicArray arr2(arr);
    EXPECT_TRUE(check_array_values(arr2, MAGIC_COPY_CTOR));

    // Array copy-assignment should invoke each element's copy constructor
    DynamicArray arr3 = arr;
    EXPECT_TRUE(check_array_values(arr3, MAGIC_COPY_CTOR));

    // Transfers between memory spaces should invoke each element's copy constructor
    DynamicArray arr4(arr, hostAllocID);
    EXPECT_TRUE(check_array_values(arr4, MAGIC_COPY_CTOR));

    // Fill with instance of copy-constructed type - each element should be
    // copy-constructed from the argument
    arr.fill(NonTrivialCopyCtor {});
    EXPECT_TRUE(check_array_values(arr, MAGIC_COPY_CTOR));

    // Second fill should be idempotent - i.e. isn't affected by the data
    // already present in the array
    arr.fill(NonTrivialCopyCtor {});
    EXPECT_TRUE(check_array_values(arr, MAGIC_COPY_CTOR));
  }
}

//------------------------------------------------------------------------------
AXOM_TYPED_TEST(core_array_for_all, nontrivial_emplace)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialCopyCtor>;
  using IntArray = typename TestFixture::DynamicArray;
  using HostIntArray = typename TestFixture::HostArray;

  int kernelAllocID = axom::execution_space<ExecSpace>::allocatorID();
#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
  if(axom::execution_space<ExecSpace>::onDevice())
  {
    kernelAllocID = axom::getUmpireResourceAllocatorID(
      umpire::resource::MemoryResourceType::Device);
  }
#endif

  // Helper function to copy device values to a host array
  auto convert_to_host_array = [=](const DynamicArray& arr) -> HostIntArray {
    IntArray arr_ints(arr.size(), arr.size(), kernelAllocID);
    const auto arr_v = arr.view();
    const auto arr_ints_v = arr_ints.view();
    axom::for_all<ExecSpace>(
      arr.size(),
      AXOM_LAMBDA(axom::IndexType i) { arr_ints_v[i] = arr_v[i].m_val; });

    // handles synchronization, if necessary
    if(axom::execution_space<ExecSpace>::async())
    {
      axom::synchronize<ExecSpace>();
    }

    return arr_ints;
  };

  // Create an array of N items using default MemorySpace for ExecSpace
  constexpr axom::IndexType N = 10;
  DynamicArray arr(N, N, kernelAllocID);

  // Check default-constructed array contents on host
  {
    HostIntArray localArr = convert_to_host_array(arr);
    for(int i = 0; i < N; ++i)
    {
      EXPECT_EQ(localArr[i], 1);
    }
  }

  // Emplace some elements
  // emplace of xvalue - should bind to rvalue emplace
  arr.emplace_back(NonTrivialCopyCtor {});
  constexpr int MAGIC_MOVE_EXPLICIT = 222;
  // emplace of explicitly-moved object - should bind to rvalue emplace
  {
    NonTrivialCopyCtor explicitMove;
    explicitMove.m_val = MAGIC_MOVE_EXPLICIT;
    arr.emplace_back(std::move(explicitMove));
  }
  {
    NonTrivialCopyCtor lvalue;
    // emplace of object as lvalue - should bind to copying emplace
    arr.emplace_back(lvalue);
  }
  EXPECT_EQ(arr.size(), N + 3);
  {
    HostIntArray localArr = convert_to_host_array(arr);
    for(int i = 0; i < N + 1; ++i)
    {
      EXPECT_EQ(localArr[i], 1);
    }
    EXPECT_EQ(localArr[N + 1], MAGIC_MOVE_EXPLICIT);
    // our one copied element should be set to the copy ctor value
    EXPECT_EQ(localArr[N + 2], MAGIC_COPY_CTOR);
  }

  // Emplace some elements in the front
  // emplace of xvalue - should bind to rvalue emplace
  arr.emplace(arr.begin(), NonTrivialCopyCtor {});
  // emplace of explicitly-moved object - should bind to rvalue emplace
  {
    NonTrivialCopyCtor explicitMove;
    explicitMove.m_val = MAGIC_MOVE_EXPLICIT;
    arr.emplace(arr.begin(), std::move(explicitMove));
  }
  {
    NonTrivialCopyCtor lvalue;
    // emplace of object as lvalue - should bind to copying emplace
    arr.emplace(arr.begin(), lvalue);
  }

  EXPECT_EQ(arr.size(), N + 6);
  {
    HostIntArray localArr = convert_to_host_array(arr);
    for(int i = 2; i < N + 4; ++i)
    {
      EXPECT_EQ(localArr[i], 1);
    }
    // check our explicitly-moved values
    EXPECT_EQ(localArr[1], MAGIC_MOVE_EXPLICIT);
    EXPECT_EQ(localArr[N + 4], MAGIC_MOVE_EXPLICIT);
    // copied objects should be at the beginning and end
    EXPECT_EQ(localArr[0], MAGIC_COPY_CTOR);
    EXPECT_EQ(localArr[N + 5], MAGIC_COPY_CTOR);
  }
}

}  // end namespace testing
