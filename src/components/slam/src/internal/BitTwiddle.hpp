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
* \file BitTwiddle
*
* \brief Some bit twiddling operations in suppport of slam::BitSet
*
*/


#ifndef SLAM_BIT_TWIDDLE_H_
#define SLAM_BIT_TWIDDLE_H_

#include "axom/Types.hpp"

namespace axom
{
namespace slam
{
namespace internal
{

/** Counts the number of trailing zeros in a word */
int trailingZeros(axom::common::uint64 word)
{
    // adapted from bit twiddling hacks, modified for 64 bits:
    // https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel

    typedef axom::common::uint64 Word;

    int cnt = 64;   // cnt tracks the number of zero bits on the right
    word &= -static_cast<axom::common::int64>(word);
    if (word) cnt--;
    if (word & 0x00000000FFFFFFFF) cnt -= 32;
    if (word & 0x0000FFFF0000FFFF) cnt -= 16;
    if (word & 0x00FF00FF00FF00FF) cnt -= 8;
    if (word & 0x0F0F0F0F0F0F0F0F) cnt -= 4;
    if (word & 0x3333333333333333) cnt -= 2;
    if (word & 0x5555555555555555) cnt -= 1;

    return cnt;
}

/** Counts the number of set bits in a word */
int popCount(axom::common::uint64 word)
{
    // 64 bit popcount implementation from: 
    // http://chessprogramming.wikispaces.com/Population+Count#SWARPopcount

    typedef axom::common::uint64 Word;

    const Word masks[] = {
        0x5555555555555555,              //  0 --  -1/3
        0x3333333333333333,              //  1 --  -1/5
        0x0F0F0F0F0F0F0F0F,              //  2 --  -1/17
        0x0101010101010101,              //  3 --  -1/255
    };

    // aggregate counts of 2-,4-,8- and 16-bits
    word = word - ((word >> 1)  & masks[0]);
    word = (word & masks[1]) + ((word >> 2)  & masks[1]);
    word = (word + (word >> 4)) & masks[2];
    return static_cast<int>((word * masks[3]) >> 56);
}

} // end namespace internal
} // end namespace slam
} // end namespace axom

#endif //  SLAM_BIT_TWIDDLE_H_
