// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

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
#include <cstdint>  // for c++11 fixed with types

#ifdef AXOM_USE_MPI
  #include <mpi.h>  // for MPI types
#endif

namespace axom
{
/*
  Axom integer types are deprecated.
  Their use can trigger an warning or an error, or be quietly allowed,
  by using cmake -DAXOM_DEPRECATED_TYPES=<WARN|ERROR|ALLOW>
  Eventually, these types will be removed.
*/
#if AXOM_DEPRECATED_TYPES_N == 1 || AXOM_DEPRECATED_TYPES_N == 2
  #if AXOM_DEPRECATED_TYPES_N == 1
    #warning \
      "Using deprecated Axom types.  Please see cmake variable AXOM_DEPRECATED_TYPES"
  #endif
using int8 = std::int8_t;   /*!< 8-bit signed integer type      */
using uint8 = std::uint8_t; /*!< 8-bit unsigned integer type    */

using int16 = std::int16_t;   /*!< 16-bit signed integer type     */
using uint16 = std::uint16_t; /*!< 16-bit unsigned integer type   */

using int32 = std::int32_t;   /*!< 32-bit signed integer type     */
using uint32 = std::uint32_t; /*!< 32-bit unsigned integer type   */

  // Note: KW -- We assume that AXOM_NO_INT64_T will be defined
  // on systems/compilers that do not support 64 bit integer types
  #ifndef AXOM_NO_INT64_T
using int64 = std::int64_t;   /*!< 64-bit signed integer type     */
using uint64 = std::uint64_t; /*!< 64-bit unsigned integer type   */
  #endif
#endif

using float32 = float;
using float64 = double;

#if defined(AXOM_USE_64BIT_INDEXTYPE) && !defined(AXOM_NO_INT64_T)
using IndexType = std::int64_t;
#else
using IndexType = std::int32_t;
#endif

#ifdef AXOM_USE_MPI

/*!
 * \brief Traits class to map Axom types to their corresponding MPI type.
 *
 * \note The mpi_traits are initialized in Types.cpp
 * \see Types.cpp
 */
template <class AxomType>
struct mpi_traits
{
  static const MPI_Datatype type;
};

/// \name Specialization of mpi_traits
/// @{

//------------------------------------------------------------------------------
template <>
struct mpi_traits<float64>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<float32>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::int8_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::uint8_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::int16_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::uint16_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::int32_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::uint32_t>
{
  static const MPI_Datatype type;
};

  //------------------------------------------------------------------------------
  #ifndef AXOM_NO_INT64_T
template <>
struct mpi_traits<std::int64_t>
{
  static const MPI_Datatype type;
};

//------------------------------------------------------------------------------
template <>
struct mpi_traits<std::uint64_t>
{
  static const MPI_Datatype type;
};

  #endif /* end AXOM_NO_INT64_T */

  /// @}

#endif /* end AXOM_USE_MPI */

}  // end namespace axom

#endif  // AXOM_TYPES_HPP_
