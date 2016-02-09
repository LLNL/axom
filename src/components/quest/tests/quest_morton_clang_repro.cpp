/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/// KW: This file demonstrates a bug in clang's optimizer (at level -O2 and above).
///     Specifically, the loop in expandBits() erroneously get optimized away in the 'DIM=3' example,
//      but is fine in the 'DIM=2' example.  The only difference is the data in the B and S arrays
//      I tested this against the clang 3.7.0 version on LC.
//      I have also tried to manually unroll the loop, and get the same problem.
//      Note that the gcc version works as expected.
//      I have a working fix for my code (but not the compiler bug), but I am not thrilled with it.
//      My fix is to set the loop variable to 'volatile'
//		This file has no dependencies on the rest of the CS toolkit.

// Running the code:
// 2D example (everything is fine).
// >  clang++ -DDIM=2 -O2 quest_morton_clang_repro.cpp -o morton2D
// Buggy 3D example
// > clang++ -DDIM=3 -O2 quest_morton_clang_repro.cpp -o morton3D_buggy
// Buggy 3D example, manually unrolled loop
// > clang++ -DDIM=3 -O2 -DATTEMPT_HAND_UNROLL quest_morton_clang_repro.cpp -o morton3D_stillBuggy
// 3D example with naive bugfix
// > clang++ -DDIM=3 -O2 -DAPPLY_FIX quest_morton_clang_repro.cpp -o morton3D_bugFix

#include <cassert>
#include <cstdlib>
#include <limits>
#include <iostream>
#include <sstream>

// Try to set dimension if not already set
// Otherwise code will assume D==2 version
#if not defined(DIM) || DIM == 3
#define USE_3D
#endif

#if defined(__clang__) && defined(HAS_OPTS) && defined(APPLY_FIX)
  #define CIRCUMVENT_CLANG_OPTIMIZATION_BUG
#endif

namespace
{

// Simple macro for logging to cout
#define LOG( msg )                                                                \
    do {                                                                          \
        std::ostringstream oss;                                                   \
        oss << msg;                                                               \
        std::cout << oss.str() << std::endl;                                      \
    } while ( 0 )

#define NOLOG( ignore_EXP ) ( (void)0 )

typedef std::size_t IndexType;

#ifdef USE_3D
    const IndexType B[] = {
            0x9249249249249249,     // 0010'0100'1001'0010'0100'1001
            0x30C30C30C30C30C3,     // 0000'1100'0011'0000'1100'0011
            0xF00F00F00F00F00F,     // 0000'0000'1111'0000'0000'1111
            0x00FF0000FF0000FF,     // 0000'0000'0000'0000'1111'1111
            0xFFFF00000000FFFF,     // x16
            0x00000000FFFFFFFF };     // x32

    const int S[] = { 2, 4, 8, 16, 32, 64 };

#else   // 2D version
    const IndexType B[] = {
            0x5555555555555555,     // 0101'0101
            0x3333333333333333,     // 0011'0011
            0x0F0F0F0F0F0F0F0F,     // 0000'1111
            0x00FF00FF00FF00FF,     // 0x8  1x8
            0x0000FFFF0000FFFF,     // 0x16 1x16
            0x00000000FFFFFFFF};    //  0x32 1x32

    const int S[] = { 1, 2, 4, 8, 16, 32};
#endif

    /**
     * This function expands the input parameter x with zeros between its bits
     * E.g. if x is defined by the following four bits: b0.b1.b2.b3
     *      the returned value will be b0.0.b1.0.b2.0.b3 in 2D mode (a zero b/w each bit),
     *      and b0.0.0.b1.0.0.b2.0.0.b3 in 3D mode (two zeros between each bit)
     */
    IndexType expandBits (IndexType x)
    {
        LOG("In function expandBits");
        LOG("\t input value: " << x);

#ifndef ATTEMPT_HAND_UNROLL

  #if defined(CIRCUMVENT_CLANG_OPTIMIZATION_BUG) && defined(USE_3D)
        for(volatile int i=5; i >= 0; --i)
  #else
        for (int i = 5; i >= 0; --i)
#endif
        {
            x = (x | (x << S[i])) & B[i];
            //LOG("\t -- step " << i << " -- value " << x);
        }

#else   // Hand unrolled version exhibits same bug
        x = (x | (x << S[5])) & B[5];
        x = (x | (x << S[4])) & B[4];
        x = (x | (x << S[3])) & B[3];
        x = (x | (x << S[2])) & B[2];
        x = (x | (x << S[1])) & B[1];
        x = (x | (x << S[0])) & B[0];
#endif

        LOG("\t output value: " << x);

        return x;
    }



#ifdef USE_3D
    /**
     * Bit interleave a 3D point
     */
    IndexType mortonize3D (const int & x, const int & y, const int & z)
    {
        return (expandBits (x) | (expandBits (y) << 1) | (expandBits (z) << 2));
    }
#else
    /**
     * Bit interleave a 2D point
     */
    IndexType mortonize2D (const int & x, const int & y)
    {
        return (expandBits (x) | (expandBits (y) << 1));
    }
#endif

}

//-----------------------------------------------------------------------------
int main ()
{
    static const int RESULT_AS_EXPECTED = 0;
    static const int UNEXPECTED_RESULT = 1;

    int result = 0;



    // Marked as volatile so compiler will not
    // overoptimize the calls to expandBits() for the case where x==2
    volatile int x = 2;  // 0b10 in binary

    // Print out a message indicating the defined values
#ifdef USE_3D
    LOG("Using 3D data arrays");
#ifdef CIRCUMVENT_CLANG_OPTIMIZATION_BUG
    LOG("\t applying bugfix -- loop variable is marked volatile");
#elif defined(ATTEMPT_HAND_UNROLL)
    LOG("\t attempting to fix by manually unrolling the loop");
#else
    LOG("\t w/o bugfix");
#endif
#else
    LOG("Using 2D data arrays");
#endif


#ifdef USE_3D
    // 3D example -- find the bit interleaved Morton code for a 3D point
    int pt3D[] =
        { x, x, x };
    IndexType idx3D = mortonize3D (pt3D[0], pt3D[1], pt3D[2]);

    LOG("Orig pt 3D: (" << pt3D[0] <<", " << pt3D[1] <<", " << pt3D[2] <<")");
    LOG("Mortenized 3D: " << idx3D);
    LOG("Expand bit("<<x<<") = " << expandBits(x));

    if (static_cast<IndexType> (56) == idx3D)
    {
        LOG("Got the correct result -- 56");

      #if defined(__clang__)
        #if defined(CIRCUMVENT_CLANG_OPTIMIZATION_BUG) or not defined(HAS_OPTS)
            result = RESULT_AS_EXPECTED;
        #else
          LOG("... but we expected the wrong result due to a clang bug.");
          LOG("Please look into why the test is passing now. Perhaps the bug is fixed.");
          result = UNEXPECTED_RESULT;
        #endif
      #else
        result = RESULT_AS_EXPECTED;
      #endif
    }
    else
    {
        // Note: 10540996613548315209 == B[0] == 0x9249249249249249
        LOG("Got the incorrect result. Expected 56, got " << idx3D);

        #if defined(__clang__)
          #if defined(CIRCUMVENT_CLANG_OPTIMIZATION_BUG)
            LOG("Workaround is no longer working.");
            result = UNEXPECTED_RESULT;
          #elif not defined(HAS_OPTS)
            result = UNEXPECTED_RESULT;
          #else
            LOG("... but we expected this due to a known clang bug.");
            result = RESULT_AS_EXPECTED;
          #endif
        #else
          result = UNEXPECTED_RESULT;
        #endif
    }


#else
    // 2D example -- find the bit interleaved Morton code for a 2D point

    int pt2D[] =
    {   x,x};
    IndexType idx2D = mortonize2D(pt2D[0], pt2D[1]);

    LOG("Orig pt 2D: (" << pt2D[0] <<", " << pt2D[1] <<")");
    LOG("Mortenized 2D: " << idx2D);
    LOG("Expand bit("<<x<<") = " << expandBits(x));

    if(static_cast<IndexType>(12) == idx2D)
    {
        LOG("Got the correct result -- 12");
        result = RESULT_AS_EXPECTED;
    }
    else
    {
        LOG("Got the incorrect result. Expected 12, got " << idx2D);
        result = UNEXPECTED_RESULT;
    }

#endif

    LOG("Test " << (result==0 ? "PASSED" : "FAILED") );

    return result;
}
