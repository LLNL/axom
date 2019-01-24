/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file
 * This file tests the ability to disable warnings about strict aliasing
 * on all supported compilers using the AXOM_DISABLE_ALIASING_WARNINGS
 * build variable.
 *
 * The example is modeled after a known example in python's headers
 * (see http://legacy.python.org/dev/peps/pep-3123) and might need
 * to be built with optimization enabled to trigger the warning.
 */

#include <iostream>

struct Foo {int i;  Foo* ob; };
struct Bar {int i;  Foo* ob; };


int main()
{
  Foo foo = { 1, NULL};
  ((Bar*)(&foo))->i++;        // violates strict aliasing

  std::cout << " foo.i: " << foo.i << std::endl;
  return 0;
}
