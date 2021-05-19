// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file BitTwiddle
 *
 * \brief Some bit twiddling operations in support of slam::BitSet
 */

#ifndef SLAM_BIT_TWIDDLE_H_
#define SLAM_BIT_TWIDDLE_H_

#include "axom/core/Types.hpp"

namespace axom
{
namespace slam
{
namespace internal
{
/**
 * Helper traits template for determining the number of bits and bytes
 * required to represent a given type \a T.
 */
template <typename T>
struct BitTraits;

template <>
struct BitTraits<axom::uint64>
{
  enum
  {
    NUM_BYTES = 8,
    BITS_PER_WORD = NUM_BYTES << 3,
    LG_BITS_PER_WORD = 6
  };
};

template <>
struct BitTraits<axom::uint32>
{
  enum
  {
    NUM_BYTES = 4,
    BITS_PER_WORD = NUM_BYTES << 3,
    LG_BITS_PER_WORD = 5
  };
};

template <>
struct BitTraits<axom::uint16>
{
  enum
  {
    NUM_BYTES = 2,
    BITS_PER_WORD = NUM_BYTES << 3,
    LG_BITS_PER_WORD = 4
  };
};

template <>
struct BitTraits<axom::uint8>
{
  enum
  {
    NUM_BYTES = 1,
    BITS_PER_WORD = NUM_BYTES << 3,
    LG_BITS_PER_WORD = 3
  };
};

/**
 * Counts the number of trailing zeros in \a word
 *
 * \return The number of zeros to the right of the first set bit in \word,
 * starting with the least significant bit, or 64 if \a word == 0.
 */
/* clang-format off */
inline int trailingZeros(axom::uint64 word)
{
  // Explicit implementation adapted from bit twiddling hacks
  // https://graphics.stanford.edu/~seander/bithacks.html#ZerosOnRightParallel
  // and modified for 64 bits:

  // cnt tracks the number of zero bits on the right of the first set bit
  int cnt = BitTraits<axom::uint64>::BITS_PER_WORD;

  word &= -static_cast<axom::int64>(word);
  if (word)                      { cnt--;     }
  if (word & 0x00000000FFFFFFFF) { cnt -= 32; }
  if (word & 0x0000FFFF0000FFFF) { cnt -= 16; }
  if (word & 0x00FF00FF00FF00FF) { cnt -=  8; }
  if (word & 0x0F0F0F0F0F0F0F0F) { cnt -=  4; }
  if (word & 0x3333333333333333) { cnt -=  2; }
  if (word & 0x5555555555555555) { cnt -=  1; }

  return cnt;
}
/* clang-format on */

/** Counts the number of set bits in \a word */
inline int popCount(axom::uint64 word)
{
  // 64 bit popcount implementation from:
  // http://chessprogramming.wikispaces.com/Population+Count#SWARPopcount

  using Word = axom::uint64;

  const Word masks[] = {
    0x5555555555555555,  //  0 --  -1/3
    0x3333333333333333,  //  1 --  -1/5
    0x0F0F0F0F0F0F0F0F,  //  2 --  -1/17
    0x0101010101010101,  //  3 --  -1/255
  };

  // aggregate counts of 2-,4-,8- and 16-bits
  word = word - ((word >> 1) & masks[0]);
  word = (word & masks[1]) + ((word >> 2) & masks[1]);
  word = (word + (word >> 4)) & masks[2];
  return static_cast<int>((word * masks[3]) >> 56);
}

}  // end namespace internal
}  // end namespace slam
}  // end namespace axom

#endif  //  SLAM_BIT_TWIDDLE_H_
