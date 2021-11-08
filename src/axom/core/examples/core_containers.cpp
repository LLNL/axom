// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file core_containers.cpp
 *  \brief This example code is a demonstration of the Axom Core containers.
 */

/* This example code contains snippets used in the Core Sphinx documentation.
  * They begin and end with comments such as
  *
  * timer_start
  * timer_end
  *
  * each prepended with an underscore.
  */

// Axom includes
#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"
#include "axom/core/Macros.hpp"
#include "axom/core/memory_management.hpp"

#ifdef AXOM_USE_RAJA
  #include "axom/core/execution/execution_space.hpp"
  #include "axom/core/execution/for_all.hpp"
#endif

// C/C++ includes
#include <iostream>
#include <vector>

void showArray(axom::Array<int>& a, const char* name)
{
  std::cout << "Array " << name << " = " << a << std::endl;
}

void showTupleArray(axom::MCArray<int>& a, const char* name)
{
  const auto numComponents = a.shape()[1];
  std::cout << "MCArray " << name << " with " << a.shape()[0] << " "
            << numComponents << "-tuples = [" << std::endl;
  for(int i = 0; i < a.shape()[0]; ++i)
  {
    std::cout << "  "
              << axom::ArrayView<int>(a.data() + (i * numComponents),
                                      numComponents)
              << std::endl;
  }
  std::cout << "]" << std::endl;
}

void showTupleArrayView(axom::MCArrayView<int>& a, const char* name)
{
  const auto numComponents = a.shape()[1];
  std::cout << "MCArrayView " << name << " with " << a.shape()[0] << " "
            << numComponents << "-tuples = [" << std::endl;
  for(int i = 0; i < a.shape()[0]; ++i)
  {
    std::cout << "  "
              << axom::ArrayView<int>(a.data() + (i * numComponents),
                                      numComponents)
              << std::endl;
  }
  std::cout << "]" << std::endl;
}

void demoArrayBasic()
{
  // _arraybasic_start
  // Here is an Array of ints with length three.
  axom::Array<int> a(3);
  std::cout << "Length of a = " << a.size() << std::endl;
  a[0] = 2;
  a[1] = 5;
  a[2] = 11;

  // An Array increases in size if a value is pushed back.
  a.push_back(4);
  std::cout << "After pushing back a value, a's length = " << a.size()
            << std::endl;

  // You can also insert a value in the middle of the Array.
  // Here we insert value 6 at position 2 and value 1 at position 4.
  showArray(a, "a");
  a.insert(2, 6);
  a.insert(4, 1);
  std::cout << "After inserting two values, ";
  showArray(a, "a");
  // _arraybasic_end

  // _arraytuple_start
  // Here is an MCArray of ints, containing two triples.
  const int numTuples = 2;
  const int numComponents = 3;
  axom::MCArray<int> b(numTuples, numComponents);
  // Set tuple 0 to (1, 4, 2).
  b(0, 0) = 1;
  b(0, 1) = 4;
  b(0, 2) = 2;
  // Set tuple 1 to one tuple, (8, 0, -1).
  // The first argument to set() is the buffer to copy into the MCArray, the
  // second is the number of tuples in the buffer, and the third argument
  // is the first tuple to fill from the buffer.
  int ival[3] = {8, 0, -1};
  b.set(ival, 3, 3);

  showTupleArray(b, "b");

  // Now, insert two tuples, (0, -1, 1), (1, -1, 0), into the MCArray, directly
  // after tuple 0.
  int jval[6] = {0, -1, 1, 1, -1, 0};
  b.insert(1, numTuples * numComponents, jval);

  showTupleArray(b, "b");
  // _arraytuple_end

  // _extbuffer_start
  // The internal buffer maintained by an MCArray is accessible.
  int* pa = a.data();
  // An MCArray can be constructed with a pointer to an external buffer.
  // Here's an Array interpreting the memory pointed to by pa as three 2-tuples.
  axom::MCArrayView<int> c(pa, 3, 2);

  showArray(a, "a");
  showTupleArrayView(c, "c");

  // Since c is an alias to a's internal memory, changes affect both Arrays.
  a[0] = 1;
  c(1, 1) = 9;

  std::cout
    << "Array a and MCArrayView c use the same memory, a's internal buffer."
    << std::endl;
  showArray(a, "a");
  showTupleArrayView(c, "c");
  // _extbuffer_end

  // _iteration_start
  // Iteration over multidimensional arrays uses the shape() method
  // to retrieve the extents in each dimension.
  for(int i = 0; i < c.shape()[0]; i++)
  {
    for(int j = 0; j < c.shape()[1]; j++)
    {
      // Note that c's operator() accepts two arguments because it is two-dimensional
      std::cout << "In ArrayView c, index (" << i << ", " << j << ") yields "
                << c(i, j) << std::endl;
    }
  }

  // To iterate over the "flat" data in an Array, regardless of dimension,
  // use a range-based for loop.
  std::cout << "Range-based for loop over ArrayView c yields: ";
  for(const int value : c)
  {
    std::cout << value << " ";
  }
  std::cout << std::endl;

  // Alternatively, the "flat" data can be iterated over with operator[]
  // from 0 -> size().
  std::cout << "Standard for loop over ArrayView c yields: ";
  for(int i = 0; i < c.size(); i++)
  {
    std::cout << c[i] << " ";
  }
  std::cout << std::endl;
  // _iteration_end
}

// The following example requires CUDA + Umpire + unified memory
#if defined(AXOM_USE_UMPIRE) && defined(AXOM_USE_CUDA) && \
  defined(__CUDACC__) && defined(UMPIRE_ENABLE_UM)
  #define AXOM_CONTAINERS_EXAMPLE_ON_DEVICE
#endif

#ifdef AXOM_CONTAINERS_EXAMPLE_ON_DEVICE

// _cuda_kernel_start
// Aliases used for convenience
using UnifiedIntArrayView = axom::ArrayView<int, 1, axom::MemorySpace::Unified>;
using DeviceIntArrayView = axom::ArrayView<int, 1, axom::MemorySpace::Device>;

__global__ void add(const UnifiedIntArrayView A,
                    const UnifiedIntArrayView B,
                    DeviceIntArrayView C)
{
  for(int i = 0; i < A.size(); i++)
  {
    C[i] = A[i] + B[i];
  }
}
// _cuda_kernel_end

#endif

void demoArrayDevice()
{
#ifdef AXOM_CONTAINERS_EXAMPLE_ON_DEVICE
  // _cuda_array_create_start
  constexpr int N = 10;
  const int allocator_id = axom::getUmpireResourceAllocatorID(
    umpire::resource::MemoryResourceType::Unified);

  // The last template parameter specifies a memory space.
  // Its default value is Dynamic, which lets the user specify the
  // memory space at runtime with a memory allocator ID.  The
  // third constructor parameter specifies the allocator.
  // If this argument is not provided host memory will be allocated.
  axom::Array<int> A_dynamic(N, N, allocator_id);

  // We also have the option to "lock down" the memory space to allow for
  // compile-time guarantees against dereferencing pointers in the wrong memory space.
  axom::Array<int, 1, axom::MemorySpace::Unified> B_unified(N);

  // Despite having different types, both of these arrays are in unified memory.
  for(int i = 0; i < N; i++)
  {
    A_dynamic[i] = i * 5;
    B_unified[i] = i * 2;
  }

  // Since our kernel requires that its arguments to be in unified memory at compile time,
  // we lock down the dynamic array.  In the general case this will result in a transfer,
  // but because the dynamic array happens to be in unified memory a direct copy is possible.
  axom::Array<int, 1, axom::MemorySpace::Unified> A_unified = A_dynamic;

  // The result array is allocated in device memory
  axom::Array<int, 1, axom::MemorySpace::Device> C_device(N);

  // _cuda_array_create_end
  // _cuda_array_call_start

  // Passing by reference is not possible for CUDA kernels, so the three arrays
  // are converted to corresponding ArrayViews that are "shallow copies" of the
  // original Array.
  add<<<1, 1>>>(A_unified, B_unified, C_device);

  // Since our result array is in device memory, we copy it to host memory so we can view it.
  axom::Array<int, 1, axom::MemorySpace::Host> C_host = C_device;
  std::cout << "Array C_host = " << C_host << std::endl;

  // Since by default allocations happen in host memory, we could have also used a dynamic array (the default)
  axom::Array<int> C_dynamic = C_device;
  std::cout << "Array C_dynamic = " << C_dynamic << std::endl;
  // _cuda_array_call_end

  #ifdef AXOM_USE_RAJA
  // _array_w_raja_start
  // To use a lambda as a kernel, we create the ArrayViews explicitly.
  const UnifiedIntArrayView A_view = A_unified;
  const UnifiedIntArrayView B_view = B_unified;
  // Create a new array for our RAJA result
  axom::Array<int, 1, axom::MemorySpace::Device> C_device_raja(N);
  DeviceIntArrayView C_view = C_device_raja;

  // Declare the lambda mutable so our copy of C_view (captured by value) is mutable
  axom::for_all<axom::CUDA_EXEC<1>>(
    0,
    N,
    [=] AXOM_HOST_DEVICE(axom::IndexType i) mutable {
      C_view[i] = A_view[i] + B_view[i] + 1;
    });

  // Finally, copy things over to host memory so we can display the data
  axom::Array<int, 1, axom::MemorySpace::Host> C_host_raja = C_view;
  std::cout << "Array C_host_raja = " << C_host_raja << std::endl;
  // _array_w_raja_end
  #endif
#endif
}

int main(int AXOM_UNUSED_PARAM(argc), char** AXOM_UNUSED_PARAM(argv))
{
  demoArrayBasic();
  demoArrayDevice();
  return 0;
}
