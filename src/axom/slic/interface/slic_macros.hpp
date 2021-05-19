// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SLIC_MACROS_HPP_
#define AXOM_SLIC_MACROS_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"  // for AXOM_HOST_DEVICE macros

/*!
 * \file slic_macros.hpp
 */

/// \name ERROR MACROS
/// @{

/*!
 * \def SLIC_ERROR( msg )
 * \brief Logs an error and aborts the application.
 * \param [in] msg user-supplied message
 * \note The SLIC_ERROR macro is always active.
 * \warning This macro calls processAbort().
 *
 * Usage:
 * \code
 *   SLIC_ERROR( "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_ERROR(msg)                                           \
  do                                                              \
  {                                                               \
    std::ostringstream __oss;                                     \
    __oss << msg;                                                 \
    axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__); \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_ERROR_IF( EXP, msg )
 * \brief Logs an error iff EXP is true and aborts the application.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 * \note The SLIC_ERROR_IF macro is always active.
 * \warning This macro calls processAbort() if EXP is true.
 *
 * Usage:
 * \code
 *   SLIC_ERROR_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_ERROR_IF(EXP, msg)                                     \
  do                                                                \
  {                                                                 \
    if(EXP)                                                         \
    {                                                               \
      std::ostringstream __oss;                                     \
      __oss << msg;                                                 \
      axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__); \
    }                                                               \
  } while(axom::slic::detail::false_value)

/// @}

/// \name WARNING MACROS
/// @{

/*!
 * \def SLIC_WARNING( msg )
 * \brief Logs a warning message.
 * \param [in] msg user-supplied message
 * \note The SLIC_WARNING macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_WARNING( "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_WARNING(msg)                                           \
  do                                                                \
  {                                                                 \
    std::ostringstream __oss;                                       \
    __oss << msg;                                                   \
    axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_WARNING_IF( EXP, msg )
 * \brief Logs an error iff EXP is true and aborts the application.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 * \note The SLIC_WARNING_IF macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_WARNING_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_WARNING_IF(EXP, msg)                                     \
  do                                                                  \
  {                                                                   \
    if(EXP)                                                           \
    {                                                                 \
      std::ostringstream __oss;                                       \
      __oss << msg;                                                   \
      axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
    }                                                                 \
  } while(axom::slic::detail::false_value)

/// @}

// Use complete debug macros when not on device
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

  //-----------------------------------------------------------------------------
  /// \name ASSERT MACROS
  /// @{

  /*!
 * \def SLIC_ASSERT( EXP )
 * \brief Asserts that a given expression is true. If the expression is not true
 *  an error will be logged and the application will be aborted.
 * \param [in] EXP user-supplied boolean expression.
 * \note This macro is only active when debugging is turned on.
 * \warning This macro calls processAbort() if EXP is false.
 *
 * Usage:
 * \code
 *   SLIC_ASSERT( my_val >= 0 );
 * \endcode
 *
 */
  #define SLIC_ASSERT(EXP)                                            \
    do                                                                \
    {                                                                 \
      if(!(EXP))                                                      \
      {                                                               \
        std::ostringstream __oss;                                     \
        __oss << "Failed Assert: " << #EXP << std::ends;              \
        axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__); \
      }                                                               \
    } while(axom::slic::detail::false_value)

  /*!
 * \def SLIC_ASSERT_MSG( EXP, msg )
 * \brief Same as SLIC_ASSERT, but with a custom error message.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 * \note This macro is only active when debugging is turned on.
 * \warning This macro calls processAbort() if EXP is false.
 * \see SLIC_ASSERT( EXP )
 *
 * Usage:
 * \code
 *   SLIC_ASSERT_MSG( my_val >= 0, "my_val must always be positive" );
 * \endcode
 *
 */
  #define SLIC_ASSERT_MSG(EXP, msg)                                          \
    do                                                                       \
    {                                                                        \
      if(!(EXP))                                                             \
      {                                                                      \
        std::ostringstream __oss;                                            \
        __oss << "Failed Assert: " << #EXP << std::endl << msg << std::ends; \
        axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__);        \
      }                                                                      \
    } while(axom::slic::detail::false_value)

  /// @}

  //-----------------------------------------------------------------------------
  /// \name DEBUG MACROS
  /// @{

  /*!
 * \def SLIC_CHECK( EXP )
 * \brief Checks that a given expression is true. If the expression is not true
 *  a warning is logged, but, in contrast to the similar SLIC_ASSERT macro the
 *  application is not aborted.
 * \param [in] EXP user-supplied boolean expression.
 * \note This macro is only active when debugging is turned on.
 *
 * Usage:
 * \code
 *   SLIC_CHECK( my_val >= 0 );
 * \endcode
 *
 */
  #define SLIC_CHECK(EXP)                                                 \
    do                                                                    \
    {                                                                     \
      if(!(EXP))                                                          \
      {                                                                   \
        std::ostringstream __oss;                                         \
        __oss << "Failed Check: " << #EXP << std::ends;                   \
        if(axom::slic::debug::checksAreErrors)                            \
        {                                                                 \
          axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__);   \
        }                                                                 \
        else                                                              \
        {                                                                 \
          axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
        }                                                                 \
      }                                                                   \
    } while(axom::slic::detail::false_value)

  /*!
 * \def SLIC_CHECK_MSG( EXP, msg )
 * \brief Same as SLIC_CHECK, but with a custom error message.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 * \note This macro is only active when debugging is turned on.
 * \see SLIC_DEBUG( EXP )
 *
 * Usage:
 * \code
 *   SLIC_CHECK_MSG( my_val >= 0, "my_val must always be positive" );
 * \endcode
 *
 */
  #define SLIC_CHECK_MSG(EXP, msg)                                          \
    do                                                                      \
    {                                                                       \
      if(!(EXP))                                                            \
      {                                                                     \
        std::ostringstream __oss;                                           \
        __oss << "Failed Check: " << #EXP << std::endl << msg << std::ends; \
        if(axom::slic::debug::checksAreErrors)                              \
        {                                                                   \
          axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__);     \
        }                                                                   \
        else                                                                \
        {                                                                   \
          axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__);   \
        }                                                                   \
      }                                                                     \
    } while(axom::slic::detail::false_value)

/// @}

// Use assert when on device
#elif defined(AXOM_DEBUG) && defined(AXOM_DEVICE_CODE)

  #define SLIC_ASSERT(EXP) assert(EXP)
  #define SLIC_ASSERT_MSG(EXP, msg) assert(EXP)
  #define SLIC_CHECK(EXP) assert(EXP)
  #define SLIC_CHECK_MSG(EXP, msg) assert(EXP)

#else  // turn off debug macros and asserts

  #define SLIC_ASSERT(ignore_EXP) ((void)0)
  #define SLIC_ASSERT_MSG(ignore_EXP, ignore_msg) ((void)0)
  #define SLIC_CHECK(ignore_EXP) ((void)0)
  #define SLIC_CHECK_MSG(ignore_EXP, ignore_msg) ((void)0)

#endif /* END ifdef AXOM_DEBUG */

/*!
 * \def SLIC_INFO( msg )
 * \brief Logs an Info message.
 * \param [in] msg user-supplied message
 * \note The SLIC_INFO macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_INFO( "informative text goes here" );
 * \endcode
 *
 */
#define SLIC_INFO(msg)                                \
  do                                                  \
  {                                                   \
    std::ostringstream __oss;                         \
    __oss << msg;                                     \
    axom::slic::logMessage(axom::slic::message::Info, \
                           __oss.str(),               \
                           __FILE__,                  \
                           __LINE__);                 \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_INFO_IF( EXP, msg )
 * \brief Logs an Info message iff EXP is true
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 * \note The SLIC_INFO_IF macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_INFO_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_INFO_IF(EXP, msg)                          \
  do                                                    \
  {                                                     \
    if(EXP)                                             \
    {                                                   \
      std::ostringstream __oss;                         \
      __oss << msg;                                     \
      axom::slic::logMessage(axom::slic::message::Info, \
                             __oss.str(),               \
                             __FILE__,                  \
                             __LINE__);                 \
    }                                                   \
  } while(axom::slic::detail::false_value)

#ifdef AXOM_DEBUG

  /*!
 * \def SLIC_DEBUG( msg )
 * \brief Logs a Debug message.
 * \param [in] msg user-supplied message
 * \note The SLIC_Debug macro is active in debug mode.
 *
 * Usage:
 * \code
 *   SLIC_DEBUG( "debug message goes here" );
 * \endcode
 *
 */
  #define SLIC_DEBUG(msg)                                \
    do                                                   \
    {                                                    \
      std::ostringstream __oss;                          \
      __oss << msg;                                      \
      axom::slic::logMessage(axom::slic::message::Debug, \
                             __oss.str(),                \
                             __FILE__,                   \
                             __LINE__);                  \
    } while(axom::slic::detail::false_value)

  /*!
 * \def SLIC_DEBUG_IF( EXP, msg )
 * \brief Logs an Debug message iff EXP is true
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 * \note The SLIC_DEBUG_IF macro is active in debug mode.
 *
 * Usage:
 * \code
 *   SLIC_DEBUG_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
  #define SLIC_DEBUG_IF(EXP, msg)                          \
    do                                                     \
    {                                                      \
      if(EXP)                                              \
      {                                                    \
        std::ostringstream __oss;                          \
        __oss << msg;                                      \
        axom::slic::logMessage(axom::slic::message::Debug, \
                               __oss.str(),                \
                               __FILE__,                   \
                               __LINE__);                  \
      }                                                    \
    } while(axom::slic::detail::false_value)

#else  // turn off debug macros

  #define SLIC_DEBUG(ignore_EXP) ((void)0)
  #define SLIC_DEBUG_IF(ignore_EXP, ignore_msg) ((void)0)

#endif

namespace axom
{
namespace slic
{
namespace detail
{
/*!
 * \brief Variable of a type that evaluates as false.
 *
 * \note Workaround for warnings about constant expressions in slic macros.
 */
struct FalseType
{
  AXOM_HOST_DEVICE
  FalseType() { }

  AXOM_HOST_DEVICE
  inline operator bool() const { return false; }
};

static const FalseType false_value;

} /* namespace detail */
} /* namespace slic */
} /* namespace axom */

#endif /* AXOM_SLIC_MACROS_HPP_ */
