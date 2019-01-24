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
 * This file tests the ability to disable warnings about uninitialized variables
 * on all supported compilers using the AXOM_DISABLE_UNINITIALIZED_WARNINGS
 * build variable.
 */

#include <iostream>
#include <cstdlib>


int main()
{
  int* result;          // Note: variable not allocated or initialized

  if( rand()%2 == 0 )
    *result = 5;

  return 0;
}
