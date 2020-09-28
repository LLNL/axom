// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file
 * This file tests the ability to disable warnings about unused local typedefs
 * on all supported compilers using the AXOM_DISABLE_UNUSED_LOCAL_TYPEDEF
 * build variable.
 */

#include <iostream>

int main()
{
  typedef int IntT;

  std::cout << "I have defined type IntT, but am not using it." << std::endl;
  return 0;
}
