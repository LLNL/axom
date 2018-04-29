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
    typedef axom::common::uint64 Word;

    for (int i = 0; i < 64; ++i)
    {
        if ( (word & (Word(1) << i)) != 0)
            return i;
    }
    return 0;
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
