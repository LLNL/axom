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

#include "gtest/gtest.h"

#include "mpi.h"

#include "axom/lumberjack/RootCommunicator.hpp"

TEST(lumberjack_RootCommunicator, case01)
{
  const int ranksLimit = 5;
  axom::lumberjack::RootCommunicator c;
  c.initialize(MPI_COMM_WORLD, ranksLimit);

  int commRank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
  EXPECT_EQ(c.rank(), commRank);
  
  c.finalize();
}

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
