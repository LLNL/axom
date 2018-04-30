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
* \note Using builtin compiler intrinsics when possible, 
* with fallback to efficient bit-twiddling operations.
*/


#ifndef SLAM_BIT_TWIDDLE_H_
#define SLAM_BIT_TWIDDLE_H_

#include "axom/Types.hpp"

#undef _SLAM_USE_INTRINSICS
#undef _SLAM_USE_INTRINSICS_WIN
#undef _SLAM_USE_INTRINSICS_GCC
#undef _SLAM_USE_INTRINSICS_PPC

// Setup intrinsics -- adapted from https://stackoverflow.com/a/22291538
#ifdef WIN32
    #if(_MSC_VER >= 1600)
    #define _SLAM_USE_INTRINSICS
    #define _SLAM_USE_INTRINSICS_WIN
    #include <intrin.h>
    #endif
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    #define _SLAM_USE_INTRINSICS
    #define _SLAM_USE_INTRINSICS_GCC
    #include <x86intrin.h>
//#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
//    #define _SLAM_USE_INTRINSICS
//    #define _SLAM_USE_INTRINSICS_PPC
//    #include <altivec.h>
#endif


namespace axom
{
namespace slam
{
namespace internal
{

#ifdef _SLAM_USE_INTRINSICS_WIN

/** Utility function to use intrinsic trailing zeros function on windows */
inline int intrinsicTrailingZeros(axom::common::uint64 word)
{
    unsigned long cnt = 0;

#ifdef _M_X64
    return _BitScanForward64(&cnt, word) ? cnt : 64;
#else
    typedef union 
    { 
        axom::common::uint64 ull_type; 
        axom::common::uint32 i_type[2]; 
    } UnsignedUnionType;

    UnsignedUnionType val = { word };

    return _BitScanForward(&cnt, val.i_type[0])
        ? cnt
        : _BitScanForward(&cnt, val.i_type[1])
        ? 32 + cnt
        : 64;
#endif // _M_X64

}

/** Utility function to use intrinsic trailing zeros function on windows */
inline int intrinsicPopCount(axom::common::uint64 word)
{
#ifdef _M_X64
    return static_cast<int>(__popcnt64(word));
#else
    typedef union
    {
        axom::common::uint64 ull_type;
        axom::common::uint32 i_type[2];
    } UnsignedUnionType;

    UnsignedUnionType val = { word };

    return static_cast<int>(__popcnt(val.i_type[0]) + __popcnt(val.i_type[1]));
#endif
}

#endif // _SLAM_USE_INTRINSICS_WIN


/** Counts the number of trailing zeros in a word */
inline int trailingZeros(axom::common::uint64 word)
{
#ifdef _SLAM_USE_INTRINSICS
  #if defined(_SLAM_USE_INTRINSICS_WIN)
    return intrinsicTrailingZeros(word);
  #elif defined(_SLAM_USE_INTRINSICS_GCC)
    return __builtin_ctzl(word);
  #endif       
#else
    // Explicit implementation adapted from bit twiddling hacks
    // https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel
    // and modified for 64 bits:

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
#endif // _SLAM_USE_INTRINSICS
}

/** Counts the number of set bits in a word */
inline int popCount(axom::common::uint64 word)
{
#ifdef _SLAM_USE_INTRINSICS
  #if defined(_SLAM_USE_INTRINSICS_WIN)
    return intrinsicPopCount(word);
  #elif defined(_SLAM_USE_INTRINSICS_GCC)
    return __builtin_popcountl(word);
  #endif       
#else
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
#endif // _SLAM_USE_INTRINSICS
}

#undef _SLAM_USE_INTRINSICS
#undef _SLAM_USE_INTRINSICS_WIN
#undef _SLAM_USE_INTRINSICS_GCC
#undef _SLAM_USE_INTRINSICS_PPC

} // end namespace internal
} // end namespace slam
} // end namespace axom

#endif //  SLAM_BIT_TWIDDLE_H_
