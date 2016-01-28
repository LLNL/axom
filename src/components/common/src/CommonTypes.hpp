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
 *  \file CommonTypes.hpp
 *
 *  \brief File exposing some common types used by toolkit components.
 *
 */

#ifndef COMMONTYPES_HPP_
#define COMMONTYPES_HPP_


#ifndef USE_CXX11
  #include <cstddef>            // brings in NULL
#endif

#include "common/cstdint_wrapper.hpp"  // for fixed bitwidth integer types

namespace asctoolkit {

namespace common {

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
#define ATK_NULLPTR NULL
#endif


  typedef detail::int8_t int8;      /** Eight bit signed integer type */
  typedef detail::uint8_t uint8;    /** Eight bit unsigned integer type */

  typedef detail::int16_t int16;    /** Sixteen bit signed integer type */
  typedef detail::uint16_t uint16;  /** Sixteen bit unsigned integer type */

  typedef detail::int32_t int32;    /** Thirty-two bit signed integer type */
  typedef detail::uint32_t uint32;  /** Thirty-two bit unsigned integer type */

  // Note: KW -- We assume that ATK_NO_INT64_T will be defined
  // on systems/compilers that do not support 64 bit integer types
  #ifndef ATK_NO_INT64_T
  typedef detail::int64_t int64;      /** Sixty-four bit signed integer type */
  typedef detail::uint64_t uint64;    /** Sixty-four bit unsigned integer type */
  #endif


} // end namespace common
} // end namespace asctoolkit


#endif // COMMONTYPES_HPP_
