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
 *  \brief File containing common types used in toolkit components.
 *
 */

#ifndef COMMONTYPES_HPP_
#define COMMONTYPES_HPP_

#ifdef USE_CXX11
  #include <cstdint>            // for fixed width types in c++11
#else
  #include <boost/cstdint.hpp>  // for fixed width types -- fallback to boost when c++11 unavailable
  #include <cstddef>            // brings in NULL
#endif

namespace asctoolkit
{

namespace common
{

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
#define ATK_NULLPTR NULL
#endif

// Adds support to common for fixed width integer types in the toolkit
// The code below assumes that ATK_NO_INT64_T will be defined
// (e.g. by the build system or host_config file)
// if 64 bit integer types are not defined by the compiler

#ifdef USE_CXX11
typedef std::int8_t int8;
typedef std::uint8_t uint8;

typedef std::int16_t int16;
typedef std::uint16_t uint16;

typedef std::int32_t int32;
typedef std::uint32_t uint32;

  #ifndef  ATK_NO_INT64_T
typedef std::int64_t int64;
typedef std::uint64_t uint64;
  #endif

#else
typedef boost::int8_t int8;
typedef boost::uint8_t uint8;

typedef boost::int16_t int16;
typedef boost::uint16_t uint16;

typedef boost::int32_t int32;
typedef boost::uint32_t uint32;

  #if defined(BOOST_NO_INT64_T) && !defined(ATK_NO_INT64_T)
    #define ATK_NO_INT64_T 1
  #endif

  #ifndef ATK_NO_INT64_T
typedef boost::int64_t int64;
typedef boost::uint64_t uint64;
  #endif

#endif



} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* COMMONTYPES_HPP_ */
