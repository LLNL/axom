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
    typedef axom::common::uint64 Word;
    int cnt = 0;
    for (int i = 0; i < 64; ++i)
    {
        if ( (word & (Word(1) << i)) != 0)
            ++cnt;
    }
    return cnt;
}

} // end namespace internal
} // end namespace slam
} // end namespace axom

#endif //  SLAM_BIT_TWIDDLE_H_
