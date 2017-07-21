/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */


/*!
 ******************************************************************************
 *
 * \file Utilities.hpp
 *
 * \brief Header file containing utility functions.
 *
 ******************************************************************************
 */

#ifndef AXOM_UTILITIES_HPP_
#define AXOM_UTILITIES_HPP_

#include "axom/config.hpp"
#include "axom/Types.hpp"
#include "axom/Macros.hpp"

#include <cmath>
#include <algorithm>

#ifdef AXOM_USE_CXX11
  #include <type_traits>
#endif


namespace axom {
namespace utilities {

/*!
 * \brief Gracefully aborts the application
 */
  void processAbort();

  /*!
   * \brief Returns the absolute value of x
   * \param x value whose absolute value is computed
   * \return abs(x) the absolute value of x
   */
  template < typename T >
  inline T abs( const T& x) { return ( (x < 0)? -x : x ); }

  /*!
   * \brief Swaps the values of a, b
   * \param [in,out] a 1st object to swap
   * \param [in,out] b 2nd object to swap
   */
  template < typename T >
  inline void swap( T& a, T& b )
  {
     T tmp = a;
     a = b; b = tmp;
  }

  /*!
   * \brief Clamps an input val to a given range.
   * \param [in] val  The value to clamp
   * \param [in] lower The lower range
   * \param [in] upper The upper range
   * \return The clamped value.
   * \post lower <= returned value <= upper.
   */
  template < typename T >
  inline T clampVal( T val, T lower, T upper )
  {
    return (val < lower) ? lower
        : (val > upper) ? upper : val;
  }

  /*!
   * Tests the endianness of the system
   *
   * \return True, if the system is little endian, false otherwise
   */
  inline bool isLittleEndian()
  {
    // Note: Endianness test adapted from: https://stackoverflow.com/a/2103095

    enum {
      O32_LITTLE_ENDIAN = 0x03020100ul,
      O32_BIG_ENDIAN    = 0x00010203ul,
      O32_PDP_ENDIAN    = 0x01000302ul };

    const union {
      axom::common::uint8 raw[4];
      axom::common::uint32 value;
    }  host_order = { { 0, 1, 2, 3 } };

    return host_order.value == O32_LITTLE_ENDIAN;
  }

  /*!
   * \brief Swaps the endianness of the input value
   *
   * \param val The input value
   * \return The value with endianness swapped
   * \note Assumes endianness is either little or big (not PDP)
   * \pre T is a native arithmetic type (i.e. integral or floating point)
   * \pre sizeof(T) must be 2, 4, or 8 bytes
   */
  template < typename T >
  T swapEndian(T val)
  {
    const int NBYTES = sizeof(T);

    AXOM_STATIC_ASSERT_MSG( NBYTES == 2 || NBYTES == 4 || NBYTES == 8,
        "swapEndian only valid for types of size 2, 4 or 8 bytes.");

  #ifdef AXOM_USE_CXX11
    AXOM_STATIC_ASSERT_MSG( std::is_arithmetic<T>::value,
       "swapEndian only valid for native arithmetic types");
  #endif

    union
    {
      axom::common::uint8 raw[ NBYTES ];
      T val;
    } swp;

    axom::common::uint8* src = reinterpret_cast<axom::common::uint8*>(&val);

    // Reverse the bytes
    for(int i=0; i < NBYTES; ++i)
    {
      swp.raw[i] = src[(NBYTES-1)-i];
    }

    return swp.val;
  }

/*!
 * \brief Fuzzy comparison of two real valued quantities
 *
 * \param a The first real valued quantities we are comparing
 * \param b The second real valued quantities we are comparing
 * \param thresh The threshold of the fuzzy comparison.  Default is 1.0e-8
 * \return \c true if the absolute value of the difference is less than thresh and false otherwise
 */
  template<typename RealType>
  bool isNearlyEqual(RealType a, RealType b, RealType thresh = 1.0e-8)
  {
    return std::fabs(a-b) <= thresh;
  }

/*!
 * \brief Fuzzy comparison of two real valued quantities
 *
 * \param a The first real valued quantities we are comparing
 * \param b The second real valued quantities we are comparing
 * \param relThresh The relative threshold of the fuzzy comparison.  Default is 1.0e-6
 * \param absThresh The absolute threshold of the fuzzy comparison.  Default is 1.0e-8
 * \return \c true if the absolute value of the difference is less than the sum of absThresh
 *         and the relative difference (relThresh times the absolute max of a and b)
 */
  template<typename RealType>
  bool isNearlyEqualRelative(RealType a, RealType b, RealType relThresh = 1.0e-6,
                            RealType absThresh = 1.0e-8)
  {
    RealType maxFabs = std::max(std::fabs(a), std::fabs(b) );
    return std::fabs(a-b) <= ( maxFabs * relThresh + absThresh);

    // Equation from Real-Time Collsion Detection book -- http://realtimecollisiondetection.net/pubs/Tolerances/
    // Note: If we use this, we must update the doxygen
    // return std::fabs(a-b) <=  std::max(absThresh, relThresh * maxFabs );
  }


}  // namespace utilities
}  // namespace axom

#endif  // AXOM_UTILITIES_HPP_
