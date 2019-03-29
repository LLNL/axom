/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

/*!
 *  \file Types.hpp
 *
 *  \brief File exposing some common types used by axom components.
 *
 */

#ifndef AXOM_TYPES_HPP_
#define AXOM_TYPES_HPP_

// Axom includes
#include "axom/config.hpp"  // for compile time definitions

// C/C++ includes
#include <cstdint>          // for c++11 fixed with types

namespace axom
{

using int8  = std::int8_t;      /*!< 8-bit signed integer type      */
using uint8 = std::uint8_t;     /*!< 8-bit unsigned integer type    */

using int16  = std::int16_t;    /*!< 16-bit signed integer type     */
using uint16 = std::uint16_t;   /*!< 16-bit unsigned integer type   */

using int32  = std::int32_t;    /*!< 32-bit signed integer type     */
using uint32 = std::uint32_t;   /*!< 32-bit unsigned integer type   */

// Note: KW -- We assume that AXOM_NO_INT64_T will be defined
// on systems/compilers that do not support 64 bit integer types
#ifndef AXOM_NO_INT64_T
using int64  = std::int64_t;    /*!< 64-bit signed integer type     */
using uint64 = std::uint64_t;   /*!< 64-bit unsigned integer type   */
#endif

using float32 = float;
using float64 = double;


#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
using IndexType = int64;
#else
using IndexType = int32;
#endif

} // end namespace axom


#endif // AXOM_TYPES_HPP_
