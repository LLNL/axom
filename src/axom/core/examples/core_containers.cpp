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
#include "axom/core/Macros.hpp"
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
    // FIXME: Replace with ArrayView
    axom::ArrayView<int> temp(a.data() + (i * numComponents), numComponents);
    std::cout << "  " << temp << std::endl;
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
    // FIXME: Replace with ArrayView
    axom::ArrayView<int> temp(a.data() + (i * numComponents), numComponents);
    std::cout << "  " << temp << std::endl;
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
  b.set(ival, 1, 1);

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

  std::cout << "MCArrays a and c use the same memory, a's internal buffer."
            << std::endl;
  showArray(a, "a");
  showTupleArrayView(c, "c");
  // _extbuffer_end
}

int main(int AXOM_NOT_USED(argc), char** AXOM_NOT_USED(argv))
{
  demoArrayBasic();
  return 0;
}
