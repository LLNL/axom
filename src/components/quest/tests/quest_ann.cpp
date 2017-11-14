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

//----------------------------------------------------------------------
TEST(quest_ann, cplx_13region_query)
{
  SLIC_INFO("*** 13-region closely-packed query.");

  // Pre-formatted lines, not to be disturbed by uncrustify: *INDENT-OFF*
  double x[] = { -2.7, -2.3, -1.5, -1.2, -0.8, -0.9, -1.8,
                 -0.8, -0.3,  0.4,  1.4,  1.5,  0.9,
                 -2.6, -2.5, -2.0, -1.7, -1.4, -1.7, -2.0, -1.3,
                 -1.6, -1.3, -0.9, -0.8, -0.9, -1.1, -1.3, -1.4,
                 -0.9, -0.3,  0.2,  0.9,  0.9,  0.5, -0.5, -0.7,
                  1.0,  1.1,  1.3,  1.6,  2.0,  2.0,  2.3,  1.6,
                 -2.5, -1.9, -1.3, -1.2, -0.9, -1.1, -1.5, -2.1, -2.3,
                 -1.0, -0.9, -0.4,  0.0, -0.1, -0.6, -1.0,
                  0.1,  0.5,  1.1,  1.3,  0.9,  0.4,
                  1.3,  1.4,  2.1,  2.4,  2.3,  1.9,
                 -1.0, -0.3,  0.0, -0.3, -0.8,
                 -0.1,  0.0,  0.3,  0.5,  0.4,  0.3,  0.1,
                  0.7,  1.1,  1.8,  2.0,  1.8,  1.4,  1.0,  0.7,
                  3.5,  3.7,  4.0,  3.6 };

  double y[] = {  1.3,  1.2,  1.2,  1.2,  1.7,  2.3,  1.8,
                  2.0,  1.3,  0.9,  1.3,  2.1,  2.9,
                  0.8,  0.2,  0.4,  0.5,  0.8,  1.1,  1.0,  1.1,
                  0.5,  0.2,  0.3,  0.9,  1.4,  1.2,  1.0,  0.6,
                  0.0, -0.2, -0.4, -0.1,  0.7,  0.9,  1.3,  0.6,
                  0.7,  0.3, -0.1,  0.4,  0.0,  0.6,  0.8,  1.4,
                 -0.6, -0.9, -1.7, -0.8, -0.2,  0.1,  0.2,  0.3,  0.0,
                 -1.3, -1.9, -1.8, -1.0, -0.7, -0.2, -0.7,
                 -0.8, -1.3, -1.3, -0.5, -0.4, -0.5,
                 -0.3, -1.1, -1.5, -1.0,  0.4, -0.3,
                 -2.5, -2.8, -2.4, -2.0, -2.0,
                 -1.5, -2.2, -2.6, -1.9, -1.6, -1.3, -1.1,
                 -2.2, -2.6, -2.3, -1.7, -1.4, -1.3, -1.5,  0.7,
                  1.0,  0.8,  0.9,  1.5 };

  double z[97];

  int region[] = {  1,  1,  1,  1,  1,  1,  1,
                    2,  2,  2,  2,  2,  2,
                    3,  3,  3,  3,  3,  3,  3,  3,
                    4,  4,  4,  4,  4,  4,  4,  4,
                    5,  5,  5,  5,  5,  5,  5,  5,
                    6,  6,  6,  6,  6,  6,  6,  6,
                    7,  7,  7,  7,  7,  7,  7,  7,  7,
                    8,  8,  8,  8,  8,  8,  8,
                    9,  9,  9,  9,  9,  9,
                   10, 10, 10, 10, 10, 10,
                   11, 11, 11, 11, 11,
                   12, 12, 12, 12, 12, 12, 12,
                   13, 13, 13, 13, 13, 13, 13, 13,
                   14, 14, 14, 14 };
  // *INDENT-ON*

  int n = 97;
  double limit = 1.4;
  int bfneighbor[97];
  int idxneighbor[97];

  for (int i = 0; i < n; ++i) {
    z[i] = 0.;
    bfneighbor[i] = -1;
    idxneighbor[i] = -1;
  }

  {
    SCOPED_TRACE("Comparing brute force with indexed, limit 1.4");
    axom::quest::all_nearest_neighbors_bruteforce(x, y, z, region, n, limit, bfneighbor);
    axom::quest::all_nearest_neighbors_index1(x, y, z, region, n, limit, idxneighbor);
    verify_array(bfneighbor, idxneighbor, n);
  }
}


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
