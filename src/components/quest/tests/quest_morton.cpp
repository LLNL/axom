/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include "gtest/gtest.h"

#include "quest/Point.hpp"
#include "quest/BitTwiddle.hpp"

#include <cstdlib>

// Uncomment the line below for true randomized points
#ifndef MORTON_TESTER_SHOULD_SEED
//  #define MORTON_TESTER_SHOULD_SEED
#endif

#ifdef MORTON_TESTER_SHOULD_SEED
  #include <ctime>      // for time() used by srand()
#endif

namespace {

    static const int MAX_ITER = 1000000;

    // Generate a random integer in the range [beg, end)
    int randomInt(int beg, int end)
    {
        int range = end-beg;
        return (std::rand() % range) + beg;
    }

    template<int DIM>
    quest::Point<int, DIM> randomPoint(int beg, int end)
    {
        quest::Point<int,DIM> pt;
        for(int i=0; i< DIM; ++i)
            pt[i] = randomInt(beg,end);

        return pt;
    }
}

//----------------------------------------------------------------------
TEST( quest_point, test_morton_2D)
{
  using namespace quest;

  static const int DIM = 2;
  typedef int CoordType;
  typedef Point<CoordType, DIM> GridPoint;

  for(int i=0; i< MAX_ITER; ++i)
  {
      GridPoint origPt = randomPoint<DIM>(0, 1<<16);

      MortonIndex mortonIdx = convertPointToMorton( origPt);
      GridPoint convertedPt = convertMortonToPoint2D<CoordType>( mortonIdx);

      EXPECT_EQ( origPt, convertedPt );
  }

}


TEST( quest_point, test_morton_3D)
{
  using namespace quest;

  static const int DIM = 3;
  typedef int CoordType;
  typedef Point<CoordType, DIM> GridPoint;


  GridPoint origPt = randomPoint<DIM>(0, 1<<10);

  MortonIndex mortonIdx = convertPointToMorton( origPt);
  GridPoint convertedPt = convertMortonToPoint3D<CoordType>( mortonIdx);

  EXPECT_EQ( origPt, convertedPt );

}



//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

#ifdef MORTON_TESTER_SHOULD_SEED
  std::srand( std::time(0) );
#else
  std::srand( 105);
#endif

  result = RUN_ALL_TESTS();

  return result;
}


