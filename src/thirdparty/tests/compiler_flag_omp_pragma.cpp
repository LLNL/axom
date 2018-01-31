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
 * This file tests the ability to disable warnings about unknown [openmp]
 * pragmas on all supported compilers using the
 * AXOM_DISABLE_OMP_PRAGMA_WARNINGS build variable.
 * It should be compiled without the openmp flag.
 */

#include <iostream>

int main()
{
  const int SZ = 100;
  int arr[SZ];

  #pragma omp parallel for
  for(int i=0 ; i<SZ ; ++i)
    arr[i] = i;

  std::cout <<"Value of array element 0 is " << arr[0] << std::endl;

  return 0;
}
