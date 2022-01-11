// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
  HostArray localArr(arr, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i], N - i);
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

  // construct a view of the array before we invoke the destructor
  auto arr_v = arr.view();

  // Array::clear should invoke the destructor N times
  NonTrivialDtor::dtor_calls = 0;
  arr.clear();
  EXPECT_EQ(NonTrivialDtor::dtor_calls, N);

  // Copy post-destructor array
  DynamicArray arr_after_dtor(N, N, kernelAllocID);
  auto arr_after_dtor_v = arr_after_dtor.view();
  axom::for_all<ExecSpace>(
    N,
    AXOM_LAMBDA(axom::IndexType i) { arr_after_dtor_v[i] = arr_v[i]; });

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr(arr_after_dtor_v, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i].m_val, NonTrivialDtor::MAGIC_DTOR);
  }
}

//------------------------------------------------------------------------------
struct NonTrivialCopyCtor
{
  constexpr static int MAGIC_COPY_CTOR {333};

  NonTrivialCopyCtor() { }
  NonTrivialCopyCtor(const NonTrivialCopyCtor&) { m_val *= MAGIC_COPY_CTOR; }

  int m_val {1};
};

constexpr int NonTrivialCopyCtor::MAGIC_COPY_CTOR;

AXOM_TYPED_TEST(core_array_for_all, nontrivial_copy_ctor_obj)
{
  using ExecSpace = typename TestFixture::ExecSpace;
  using DynamicArray =
    typename TestFixture::template DynamicTArray<NonTrivialCopyCtor>;
  using HostArray = typename TestFixture::template HostTArray<NonTrivialCopyCtor>;

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

  // Fill with instance of copy-constructed type
  arr.fill(NonTrivialCopyCtor {});

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  HostArray localArr(arr, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i].m_val, NonTrivialCopyCtor::MAGIC_COPY_CTOR);
  }

  // Second fill should be idempotent - i.e. isn't affected by the data already
  // present in the array
  arr.fill(NonTrivialCopyCtor {});

  // handles synchronization, if necessary
  if(axom::execution_space<ExecSpace>::async())
  {
    axom::synchronize<ExecSpace>();
  }

  // Check array contents on host
  localArr = HostArray(arr, hostAllocID);
  for(int i = 0; i < N; ++i)
  {
    EXPECT_EQ(localArr[i].m_val, NonTrivialCopyCtor::MAGIC_COPY_CTOR);
  }
}

}  // end namespace testing
