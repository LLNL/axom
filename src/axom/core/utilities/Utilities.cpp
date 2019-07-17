// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file
 *
 * \brief   Implementation file for utility functions.
 *
 */

#include "axom/config.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include <cstdlib>  // for exit, EXIT_SUCCESS, EXIT_FAILURE

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#endif

namespace axom
{
namespace utilities
{
int binomial_coefficient(int n, int k)
{
  if(k > n - k)
  {
    k = n - k;
  }
  int val = 1;
  for(int i = 1; i <= k; ++i)
  {
    val *= (n - k + i);
    val /= i;
  }
  return val;
}

void processAbort()
{
#ifndef AXOM_USE_MPI
  abort();
#else
  int mpi = 0;
  MPI_Initialized(&mpi);
  if(mpi)
  {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
  abort();
#endif
}

}  // end namespace utilities
}  // end namespace axom
