/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
