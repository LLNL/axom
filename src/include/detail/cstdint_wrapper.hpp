
#ifndef COMMON_CSTDINT_WRAPPER_HPP_
#define COMMON_CSTDINT_WRAPPER_HPP_

#include "axom/config.hpp"    // defines AXOM_USE_CXX11

#ifdef AXOM_USE_CXX11
  #include <cstdint>            // for fixed width types in c++11
#endif


/**
 * \file
 * \brief This file typedefs the fixed-width types that we want to expose in
 * COMMON
 *
 * These types are defined in the C++11 standard (actually in C99),
 * When C++11 is not available, we attempt to define the types ourselves
 * using a code snippet adapted from http://snipplr.com/view/18199/stdinth
 * (see copyright info for that snippet below).
 *
 * We place these types in the detail namespace of common
 */



namespace axom
{
namespace common
{
namespace detail
{


#ifdef AXOM_USE_CXX11
typedef std::int8_t int8_t;
typedef std::uint8_t uint8_t;

typedef std::int16_t int16_t;
typedef std::uint16_t uint16_t;

typedef std::int32_t int32_t;
typedef std::uint32_t uint32_t;

  #ifndef  AXOM_NO_INT64_T
typedef std::int64_t int64_t;
typedef std::uint64_t uint64_t;
  #endif

#else

/* Note (KW 1/2016): code adapted from: http://snipplr.com/view/18199/stdinth
 *
 * ISO C9x  7.18  Integer types <stdint.h>
 * Based on ISO/IEC SC22/WG14 9899 Committee draft (SC22 N2794)
 *
 *  THIS SOFTWARE IS NOT COPYRIGHTED
 *
 *  Contributor: Danny Smith <danny_r_smith_2001@yahoo.co.nz>
 *
 *  This source code is offered for use in the public domain. You may
 *  use, modify or distribute it freely.
 *
 *  This code is distributed in the hope that it will be useful but
 *  WITHOUT ANY WARRANTY. ALL WARRANTIES, EXPRESS OR IMPLIED ARE HEREBY
 *  DISCLAIMED. This includes but is not limited to warranties of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 *  Date: 2000-12-02
 *
 * mwb: This was modified in the following ways:
 *
 *      - make it compatible with Visual C++ 6 (which uses
 *          non-standard keywords and suffixes for 64-bit types)
 *      - some environments need stddef.h included (for wchar stuff?)
 *      - handle the fact that Microsoft's limits.h header defines
 *          SIZE_MAX
 *      - make corrections for SIZE_MAX, INTPTR_MIN, INTPTR_MAX, UINTPTR_MAX,
 *          PTRDIFF_MIN, PTRDIFF_MAX, SIG_ATOMIC_MIN, and SIG_ATOMIC_MAX
 *          to be 64-bit aware.
 *
 * KW for axom project:
 *      - Only took parts of file that we needed for fixed-bitwidth types in the
 *        toolkit
 *      - Appended 'AXOM_COMMON_' to macros to avoid possible collisions with
 *        stdint
 */


#if _MSC_VER && (_MSC_VER < 1300)
/* using MSVC 6 or earlier - no "long long" type, but might have _int64 type */
  #define __AXOM_COMMON_STDINT_LONGLONG           __int64
  #define __AXOM_COMMON_STDINT_LONGLONG_SUFFIX    i64
#else
  #define __AXOM_COMMON_STDINT_LONGLONG           long long
  #define __AXOM_COMMON_STDINT_LONGLONG_SUFFIX    LL
#endif


/* 7.18.1.1  Exact-width integer types */
typedef signed char int8_t;
typedef unsigned char uint8_t;

typedef short int16_t;
typedef unsigned short uint16_t;

typedef int int32_t;
typedef unsigned uint32_t;

typedef __AXOM_COMMON_STDINT_LONGLONG int64_t;
typedef unsigned __AXOM_COMMON_STDINT_LONGLONG uint64_t;


/* 7.18.1.2  Minimum-width integer types */
typedef signed char int_least8_t;
typedef unsigned char uint_least8_t;

typedef short int_least16_t;
typedef unsigned short uint_least16_t;

typedef int int_least32_t;
typedef unsigned uint_least32_t;

typedef __AXOM_COMMON_STDINT_LONGLONG int_least64_t;
typedef unsigned __AXOM_COMMON_STDINT_LONGLONG uint_least64_t;


/*  7.18.1.3  Fastest minimum-width integer types
 *  Not actually guaranteed to be fastest for all purposes
 *  Here we use the exact-width types for 8 and 16-bit ints.
 */
typedef char int_fast8_t;
typedef unsigned char uint_fast8_t;

typedef short int_fast16_t;
typedef unsigned short uint_fast16_t;

typedef int int_fast32_t;
typedef unsigned int uint_fast32_t;

typedef __AXOM_COMMON_STDINT_LONGLONG int_fast64_t;
typedef unsigned __AXOM_COMMON_STDINT_LONGLONG uint_fast64_t;


  #undef __AXOM_COMMON_STDINT_LONGLONG
  #undef __AXOM_COMMON_STDINT_LONGLONG_SUFFIX

#endif  // AXOM_USE_CXX11
}     // end namespace detail
}     // end namespace common
}     // end namespace axom


#endif // COMMON_CSTDINT_WRAPPER_HPP_
