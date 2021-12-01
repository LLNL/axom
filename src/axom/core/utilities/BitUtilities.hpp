// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *
 * \file BitUtilities.hpp
 *
 * \brief Header file containing bitwise utility functions.
 *
 */

#ifndef AXOM_BIT_UTILITIES_HPP
#define AXOM_BIT_UTILITIES_HPP

#include "axom/core/Macros.hpp"  // for AXOM_DEVICE_CODE
#include "axom/core/Types.hpp"

namespace axom
{
namespace utilities
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
  constexpr static int NUM_BYTES = 8;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 6;
};

template <>
struct BitTraits<axom::uint32>
{
  constexpr static int NUM_BYTES = 4;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 5;
};

template <>
struct BitTraits<axom::uint16>
{
  constexpr static int NUM_BYTES = 2;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 4;
};

template <>
struct BitTraits<axom::uint8>
{
  constexpr static int NUM_BYTES = 1;
  constexpr static int BITS_PER_WORD = NUM_BYTES << 3;
  constexpr static int LG_BITS_PER_WORD = 3;
};

/**
 * Counts the number of trailing zeros in \a word
 *
 * \return The number of zeros to the right of the first set bit in \word,
 * starting with the least significant bit, or 64 if \a word == 0.
 */
/* clang-format off */
AXOM_HOST_DEVICE inline int trailingZeros(axom::uint64 word)
{
#ifdef AXOM_DEVICE_CODE
  return word != axom::uint64(0) ? __ffsll(word) - 1 : 64;
#else
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
#endif
}
/* clang-format on */

/** Counts the number of set bits in \a word */
AXOM_HOST_DEVICE inline int popCount(axom::uint64 word)
{
#ifdef AXOM_DEVICE_CODE
  // Use CUDA intrinsic for popcount
  return __popcll(word);
#else
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
#endif
}

//------------------------------------------------------------------------------
//
// count leading zeros
//
AXOM_HOST_DEVICE inline axom::int32 leadingZeros(axom::int32 x)
{
#ifdef AXOM_DEVICE_CODE
  // Use CUDA intrinsic for count leading zeros
  return __clz(x);
#else
  axom::int32 y;
  axom::int32 n = 32;
  y = x >> 16;
  if(y != 0)
  {
    n = n - 16;
    x = y;
  }
  y = x >> 8;
  if(y != 0)
  {
    n = n - 8;
    x = y;
  }
  y = x >> 4;
  if(y != 0)
  {
    n = n - 4;
    x = y;
  }
  y = x >> 2;
  if(y != 0)
  {
    n = n - 2;
    x = y;
  }
  y = x >> 1;
  if(y != 0) return axom::int32(n - 2);
  return axom::int32(n - x);
#endif
}

}  // namespace utilities
}  // namespace axom

#endif  // AXOM_BIT_UTILITIES_HPP
