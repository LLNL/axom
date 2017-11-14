/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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

#include "quest/ANNQuery.hpp"

#include "slic/slic.hpp"

template < typename T >
void verify_array(T* standard, T* expt, int n)
{
  int mismatches = 0;

  for (int i = 0; i < n; ++i) {
    if (standard[i] != expt[i]) {
      ++mismatches;
      SLIC_INFO("i " << i << " standard " << standard[i] << " expt " << expt[i]);
    }
  }

  if (mismatches > 0) {
    ADD_FAILURE() << " with " << mismatches << " mismatches.";
  } else {
    SUCCEED();
  }
}

//----------------------------------------------------------------------
TEST(quest_ann, simple_2D_query)
{
  SLIC_INFO("*** This test verifies a simple 2D all-nearest-neighbors query.");

  double x[] = {-1.2, -1.0, -0.8, -1.0, 0.8,  1.0, 1.2, 1.0};
  double y[] = { 0.0, -0.2,  0.0,  0.2, 0.0, -0.2, 0.0, 0.2};
  double z[] = { 0.0,  0.0,  0.0,  0.0, 0.0,  0.0, 0.0, 0.0};
  int region[] = {0, 0, 0, 0, 1, 1, 1, 1};
  int n = 8;
  double limit = 1.9;
  int neighbor[] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int expneighbor[] = {-1, 4, 4, 4, 2, 2, -1, 2};

  {
    SCOPED_TRACE("brute force limit 1.9");
    axom::quest::all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, neighbor);
    verify_array(expneighbor, neighbor, n);
  }
  {
    SCOPED_TRACE("indexed limit 1.9");
    axom::quest::all_nearest_neighbors_index1(x, y, z, region, n, limit, neighbor);
    verify_array(expneighbor, neighbor, n);
  }
}

//----------------------------------------------------------------------
TEST(quest_ann, simple_3D_query)
{
  SLIC_INFO("*** This test verifies a simple 3D all-nearest-neighbors query.");

  double x[] = {-1.2, -1.0, -0.8, -1.0, 0.8, 1.0, 1.2, 1.0};
  double y[] = { 0.0, -0.2,  0.0, -0.1, 0.0, 0.2, 0.0, 0.1};
  double z[] = { 0.0,  0.0,  0.0,  0.2, 0.0, 0.0, 0.0, 0.2};
  int region[] = {0, 0, 0, 0, 1, 1, 1, 1};
  int n = 8;
  double limit = 1.9;
  int neighbor[] = {-1, -1, -1, -1, -1, -1, -1, -1};
  int expneighbor[] = {-1, 4, 4, 4, 2, 2, -1, 2};

  {
    SCOPED_TRACE("brute force limit 1.9");
    axom::quest::all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, neighbor);
    verify_array(expneighbor, neighbor, n);
  }
  {
    SCOPED_TRACE("indexed limit 1.9");
    axom::quest::all_nearest_neighbors_index1(x, y, z, region, n, limit, neighbor);
    verify_array(expneighbor, neighbor, n);
  }
}

// //----------------------------------------------------------------------
// TEST(quest_ann, cplx_13region_query)
// {
//   SLIC_INFO("*** 13-region closely-packed query.");

// }


//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
