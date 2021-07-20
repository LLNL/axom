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
#include "axom/core/StackArray.hpp"
#include "axom/core/execution/execution_space.hpp"
#include "axom/core/execution/for_all.hpp"
#include "axom/core/execution/synchronize.hpp"
#include "axom/core/memory_management.hpp"

#ifdef WIN32
  #include "windows.h"
void sleep(int numSeconds)
{
  int numMilliSecs = numSeconds * 1000;
  Sleep(numMilliSecs);
}
#else
  #include <unistd.h>  // for sleep()
#endif

// C/C++ includes
#include <iostream>
#include <vector>

void demoMemoryManageBasic()
{
  // _membasic_start
  int *dynamic_memory_array;
  int *dyn_array_dst;
  int len = 20;

  //Allocation looks similar to use of malloc() in C -- just template
  //return type instead of casting.
  dynamic_memory_array = axom::allocate<int>(len);

  for(int i = 0; i < len; i++)
  {
    dynamic_memory_array[i] = i;
  }

  for(int i = 0; i < len; i++)
  {
    std::cout << i << " Current value: " << dynamic_memory_array[i] << std::endl;
  }

  //Now, a copy operation.
  dyn_array_dst = axom::allocate<int>(len);

  //It's used exactly like memcpy -- destination, soruce, number of bytes.
  axom::copy(dyn_array_dst, dynamic_memory_array, sizeof(int) * len);

  for(int i = 0; i < len; i++)
  {
    std::cout << i << " Current value: " << dyn_array_dst[i] << std::endl;
    std::cout << "Matches old value? " << std::boolalpha
              << (dynamic_memory_array[i] == dyn_array_dst[i]) << std::endl;
  }

  //Deallocate is exactly like free. Of course, we won't try to access the now-illegal memory after this:
  axom::deallocate(dyn_array_dst);

  //Reallocate is like realloc -- copies existing contents into a larger memory space.
  //Slight deviation from realloc() in that it asks for item count, rather than bytes.
  dynamic_memory_array = axom::reallocate(dynamic_memory_array, len * 2);
  for(int i = 20; i < len * 2; i++)
  {
    dynamic_memory_array[i] = i;
  }
  for(int i = 0; i < len * 2; i++)
  {
    std::cout << i << " Current value: " << dynamic_memory_array[i] << std::endl;
  }
  // _membasic_end
}

void demoAxomExecution()
{
  // _exebasic_start
  //This part of the code works regardless of Umpire's presence, allowing for generic use of axom::allocate across C and C++
  //codes.
  int *A = axom::allocate<int>(1000);
  int *B = axom::allocate<int>(1000);
  int *C = axom::allocate<int>(1000);

  for(int i = 0; i < 1000; i++)
  {
    A[i] = i * 5;
    B[i] = i * 2;
    C[i] = 0;
  }

  //Axom provides an API for the most basic usage of RAJA, the for_all loop.
  //This'll just be a simple reduce.
  axom::for_all<axom::SEQ_EXEC>(
    0,
    1000,
    AXOM_LAMBDA(axom::IndexType i) { C[i] = A[i] + B[i]; });

  std::cout << "Sums: " << std::endl;
  for(int i = 0; i < 1000; i++)
  {
    std::cout << C[i] << " ";
    C[i] = 0;
  }

// _exebasic_end

//Now, let's say we want to try out use of CUDA. We just change that execution space.
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE) && \
  defined(AXOM_USE_CUDA) && defined(__CUDACC__)
  // _cudaexebasic_start
  //This example requires Umpire to be in use, and Unified memory available.
  A = axom::allocate<int>(1000,
                          axom::getUmpireResourceAllocatorID(
                            umpire::resource::MemoryResourceType::Unified));
  B = axom::allocate<int>(1000,
                          axom::getUmpireResourceAllocatorID(
                            umpire::resource::MemoryResourceType::Unified));
  C = axom::allocate<int>(1000,
                          axom::getUmpireResourceAllocatorID(
                            umpire::resource::MemoryResourceType::Unified));

  for(int i = 0; i < 1000; i++)
  {
    A[i] = i * 5;
    B[i] = i * 2;
    C[i] = 0;
  }

  axom::for_all<axom::CUDA_EXEC<256>>(
    0,
    1000,
    AXOM_LAMBDA(axom::IndexType i) { C[i] = A[i] + B[i]; });

  std::cout << "Sums: " << std::endl;
  for(int i = 0; i < 1000; i++)
  {
    std::cout << C[i] << " ";
  }
  std::cout << std::endl;
// _cudaexebasic_end
#endif
}

int main(int AXOM_NOT_USED(argc), char **AXOM_NOT_USED(argv))
{
  demoMemoryManageBasic();
  demoAxomExecution();
  return 0;
}
