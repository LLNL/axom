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
#include "quest/MortonIndex.hpp"

#include <cstdlib>
#include <limits>

// Uncomment the line below for true randomized points
#ifndef MORTON_TESTER_SHOULD_SEED
//  #define MORTON_TESTER_SHOULD_SEED
#endif

#ifdef MORTON_TESTER_SHOULD_SEED
  #include <ctime>      // for time() used by srand()
#endif

namespace {

    static const int MAX_ITER = 10000;

    // Generate a random integer in the range [beg, end)
    template<typename CoordType>
    CoordType randomInt(CoordType beg, CoordType end)
    {
        CoordType range = end-beg;

        if(range == 0)
            range = std::numeric_limits<CoordType>::max();

        return (std::rand() % range) + beg;
    }

    template<typename CoordType, int DIM>
    quest::Point<CoordType, DIM> randomPoint(CoordType beg, CoordType end)
    {
        quest::Point<CoordType,DIM> pt;
        for(int i=0; i< DIM; ++i)
            pt[i] = randomInt(beg,end);

        return pt;
    }
}


TEST( quest_point, test_max_set_bit)
{
    std::cout <<" This test checks that MortonBase's maxSetBit function works properly"
            << std::endl;


    typedef int CoordType;

    quest::Mortonizer<CoordType,2> morton2;
    EXPECT_EQ( morton2.maxSetBit( 0), 0);


    int maxBit = std::numeric_limits<CoordType>::digits;
    for(int i=0; i<= maxBit; ++i)
    {
        int val = 1 << i;
        EXPECT_EQ( morton2.maxSetBit( val), i);

        for(int j=0; j< MAX_ITER; ++j)
        {
          int randVal = randomInt<int>(0, val);
          EXPECT_EQ( morton2.maxSetBit( val + randVal), i);
        }
    }

}


TEST( quest_point, test_mortonizer)
{
  using namespace quest;

  Point<int,2> pt2(2);  // (0b10, 0b10)
  Mortonizer<int,2> morton2;
  MortonIndex mIdx2 = morton2.mortonize(pt2);
  EXPECT_EQ( mIdx2, static_cast<MortonIndex>(0b1100) ); // interleaved bits

  MortonIndex mIdx2Alt = morton2.mortonize(pt2[0],pt2[1]);
  EXPECT_EQ( mIdx2, mIdx2Alt );
  EXPECT_EQ( morton2.demortonize(mIdx2Alt), pt2 );

  // Example from doxygen documentation Point(6,3) has MortonIndex 29
  Point<int,2> ptDoxygenExample;
  ptDoxygenExample[0] = 6;
  EXPECT_EQ(ptDoxygenExample[0], 0b0110 );

  ptDoxygenExample[1] = 3;
  EXPECT_EQ(ptDoxygenExample[1], 0b0011 );

  MortonIndex mEx = morton2.mortonize(ptDoxygenExample);
  EXPECT_EQ( mEx, static_cast<MortonIndex>(30) ); // interleaved bits
  EXPECT_EQ(30, 0b00011110 );


  Point<int,3> pt3(2);  // (0b10, 0b10, 0b10)
  Mortonizer<int,3> morton3;
  MortonIndex mIdx3 = morton3.mortonize(pt3);
  EXPECT_EQ( mIdx3, static_cast<MortonIndex>(0b111000) ); // interleaved bits

  MortonIndex mIdx3Alt = morton3.mortonize(pt3[0],pt3[1],pt3[2]);
  EXPECT_EQ( mIdx3, mIdx3Alt );
  EXPECT_EQ( morton3.demortonize(mIdx3Alt), pt3 );


  // The following will not compile -- static_assert CoordType must be integral
  //Mortonizer<double,2> mortonD;
}


template<typename CoordType>
void testMortonizer2D()
{
    using namespace quest;

    static const int DIM = 2;
    typedef Point<CoordType, DIM> GridPoint;

    int maxBits = quest::Mortonizer<CoordType,DIM>::uniqueCoordBits();
    std::cout<<"\tMax bits is: " << maxBits << std::endl;


    for(int i=0; i< MAX_ITER; ++i)
    {
        GridPoint origPt = randomPoint<CoordType, DIM>(0, 1 << maxBits);

        MortonIndex mortonIdx = convertPointToMorton2D( origPt);
        GridPoint convertedPt = convertMortonToPoint2D<CoordType>( mortonIdx);
        EXPECT_EQ( origPt, convertedPt );

        MortonIndex convertedMortonIdx = convertPointToMorton2D( convertedPt );
        EXPECT_EQ( mortonIdx, convertedMortonIdx );
    }
}


TEST( quest_point, test_integral_types_2D)
{
    std::cout <<"*** Testing morton indexing in 2D with different coord types" << std::endl;

    std::cout <<"Testing char -- ";
    testMortonizer2D<char>();

    std::cout <<"Testing uchar -- ";
    testMortonizer2D<unsigned char>();

    // --
    std::cout <<"Testing short -- ";
    testMortonizer2D<short>();

    std::cout <<"Testing ushort -- ";
    testMortonizer2D<unsigned short>();

    // --
    std::cout <<"Testing int-- ";
    testMortonizer2D<int>();

    std::cout <<"Testing uint-- ";
    testMortonizer2D<unsigned int>();

    // --
    std::cout <<"Testing long-- ";
    testMortonizer2D<long long int>();

    std::cout <<"Testing ull-- ";
    testMortonizer2D<unsigned long long int>();

}


template<typename CoordType>
void testMortonizer3D()
{
    using namespace quest;

    static const int DIM = 3;
    typedef Point<CoordType, DIM> GridPoint;

    int maxBits = quest::Mortonizer<CoordType,DIM>::uniqueCoordBits();
    std::cout<<"\tMax bits is: " << maxBits << std::endl;


    for(int i=0; i< MAX_ITER; ++i)
    {
        GridPoint origPt = randomPoint<CoordType, DIM>(0, 1 << maxBits);

        MortonIndex mortonIdx = convertPointToMorton3D( origPt);
        GridPoint convertedPt = convertMortonToPoint3D<CoordType>( mortonIdx);
        EXPECT_EQ( origPt, convertedPt );

        MortonIndex convertedMortonIdx = convertPointToMorton3D( convertedPt );
        EXPECT_EQ( mortonIdx, convertedMortonIdx );
    }
}


TEST( quest_point, test_integral_types_3D)
{
    std::cout <<"*** Testing morton indexing in 3D with different coord types" << std::endl;

    std::cout <<"Testing char -- ";
    testMortonizer3D<char>();

    std::cout <<"Testing uchar -- ";
    testMortonizer3D<unsigned char>();

    // --
    std::cout <<"Testing short -- ";
    testMortonizer3D<short>();

    std::cout <<"Testing ushort -- ";
    testMortonizer3D<unsigned short>();

    // --
    std::cout <<"Testing int-- ";
    testMortonizer3D<int>();

    std::cout <<"Testing uint-- ";
    testMortonizer3D<unsigned int>();

    // --
    std::cout <<"Testing long-- ";
    testMortonizer3D<long long int>();

    std::cout <<"Testing ull-- ";
    testMortonizer3D<unsigned long long int>();

}



TEST( quest_point, test_point_hasher)
{
    using namespace quest;

    std::cout<<"** Here we test the point hasher which can be used e.g. in an unordered_map" << std::endl;

    typedef int CoordType;
    PointHash<CoordType> ptHash;
    Point<CoordType,1> p1(2);
    std::size_t exp = 0b10;
    EXPECT_EQ( ptHash(p1), exp );

    Point<CoordType,2> p2(2);
    exp = 0b1100;
    EXPECT_EQ( ptHash(p2), exp);
    p2[0] = 0b111000;
    p2[1] = 0b001100;
    exp = 0b010111100000;
    EXPECT_EQ( ptHash(p2), exp);

    Point<CoordType,3> p3(2);
    exp = 0b111000;
    EXPECT_EQ( ptHash(p3), exp);
    p3[0] = 0b011111;
    p3[1] = 0b000100;
    p3[2] = 0b010000;
    exp = 051311;                   // in octal (read bits bottom up, left to right)
    EXPECT_EQ( ptHash(p3), exp);

    Point<CoordType,4> p4(2);
    exp = 0b11110000;
    EXPECT_EQ( ptHash(p4), exp);
    p4[0] = 0b111111;
    p4[1] = 0b001000;
    p4[2] = 0b010000;
    p4[3] = 0b100000;
    exp = 0x953111;                 // in hex (read bits bottom up, left to right)
    EXPECT_EQ( ptHash(p4), exp);

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


