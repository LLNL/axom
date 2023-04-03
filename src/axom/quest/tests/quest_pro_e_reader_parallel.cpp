// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#include "axom/quest/readers/PProEReader.hpp"

#include "gtest/gtest.h"

#include "mpi.h"


//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);

  // finalized when exiting main scope
  result = RUN_ALL_TESTS();

  MPI_Finalize();

  return result;
}