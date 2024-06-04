// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SLIC_MACROS_HPP_
#define AXOM_SLIC_MACROS_HPP_

#include "axom/config.hpp"
#include "axom/core/Macros.hpp"

/*!
 * \file slic_macros.hpp
 */

///@{
//! \name SLIC_ERROR MACROS
//!
//! \collective
//! \attention These error macros are collective operations.
//! All ranks in the user-supplied communicator must call the macro
//! when used within an MPI distributed environment, and slic::enableAbortOnError()
//! is called for the current active logger (default is enabled
//! for loggers)
//! \sa axom::slic::isAbortOnErrorsEnabled()
//! \sa axom::slic::setAbortOnError(bool status)
//!

/*!
 * \def SLIC_ERROR( msg )
 * \brief Logs an error and aborts the application.
 *
 * \param [in] msg user-supplied message
 *
 * \warning This macro calls processAbort().
 * \note The SLIC_ERROR macro is always active.
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
    if(axom::slic::isAbortOnErrorsEnabled())                      \
    {                                                             \
      axom::slic::abort();                                        \
    }                                                             \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_ERROR_IF( EXP, msg )
 * \brief Logs an error iff EXP is true and aborts the application.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \warning This macro calls processAbort() iff EXP is true.
 * \note The SLIC_ERROR_IF macro is always active.
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
      if(axom::slic::isAbortOnErrorsEnabled())                      \
      {                                                             \
        axom::slic::abort();                                        \
      }                                                             \
    }                                                               \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_ERROR_ROOT( msg )
 * \brief Macro that logs given error message only on root.
 *
 * \param [in] msg user-supplied message.
 *
 * \warning This macro calls processAbort() iff EXP is true.
 * \note The SLIC_ERROR_ROOT macro is always active.
 * \note By default, all ranks are considered to be root.
 *       Must call `axom::slic::initialize(is_root={true|false})`
 *       or set via `axom::slic::setIsRoot({true|false})` to filter based on root.
 *
 * Usage:
 * \code
 *   SLIC_ERROR_ROOT( "An error has occurred!" );
 * \endcode
 *
 */
#define SLIC_ERROR_ROOT(msg)                                        \
  do                                                                \
  {                                                                 \
    if(axom::slic::isRoot())                                        \
    {                                                               \
      std::ostringstream __oss;                                     \
      __oss << msg;                                                 \
      axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__); \
      if(axom::slic::isAbortOnErrorsEnabled())                      \
      {                                                             \
        axom::slic::abort();                                        \
      }                                                             \
    }                                                               \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_ERROR_ROOT_IF( EXP, msg )
 * \brief Macro that logs given error message only on root iff EXP is true.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_ERROR_ROOT_IF macro is always active.
 * \note By default, all ranks are considered to be root.
 *       Must call `axom::slic::initialize(is_root={true|false})`
 *       or set via `axom::slic::setIsRoot({true|false})` to filter based on root.
 *
 * Usage:
 * \code
 *   SLIC_ERROR_ROOT_IF( (my_val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_ERROR_ROOT_IF(EXP, msg)                                  \
  do                                                                  \
  {                                                                   \
    if(EXP)                                                           \
    {                                                                 \
      if(axom::slic::isRoot())                                        \
      {                                                               \
        std::ostringstream __oss;                                     \
        __oss << msg;                                                 \
        axom::slic::logErrorMessage(__oss.str(), __FILE__, __LINE__); \
        if(axom::slic::isAbortOnErrorsEnabled())                      \
        {                                                             \
          axom::slic::abort();                                        \
        }                                                             \
      }                                                               \
    }                                                                 \
  } while(axom::slic::detail::false_value)

///@}

///@{
//! \name SLIC_WARNING MACROS
//!
//! \collective
//! \attention These warning macros can be set as collective operations.
//! These warning macros are collective if slic::enableAbortOnWarning()
//! is called for the current active logger (default is disabled
//! for loggers).
//! These warning macros must then be called by all ranks in the
//! user-supplied communicator when used within an MPI distributed
//! environment.
//!
//! \sa axom::slic::isAbortOnWarningsEnabled()
//! \sa axom::slic::setAbortOnWarning(bool status)
//!

/*!
 * \def SLIC_WARNING( msg )
 * \brief Logs a warning message.
 *
 * \param [in] msg user-supplied message
 *
 * \note The SLIC_WARNING macro is always active.
 * \note Aborts the application when `slic::enableAbortOnWarning()`
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
    if(axom::slic::isAbortOnWarningsEnabled())                      \
    {                                                               \
      axom::slic::abort();                                          \
    }                                                               \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_WARNING_IF( EXP, msg )
 * \brief Logs a warning iff EXP is true
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_WARNING_IF macro is always active.
 * \note Aborts the application when `slic::enableAbortOnWarning()`
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
      if(axom::slic::isAbortOnWarningsEnabled())                      \
      {                                                               \
        axom::slic::abort();                                          \
      }                                                               \
    }                                                                 \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_WARNING_ROOT( msg )
 * \brief Macro that logs given warning message only on root.
 *
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_WARNING_ROOT macro is always active.
 * \note By default, all ranks are considered to be root.
 *       Must call `axom::slic::initialize(is_root={true|false})`
 *       or set via `axom::slic::setIsRoot({true|false})` to filter based on root.
 *
 * Usage:
 * \code
 *   SLIC_WARNING_ROOT( "A warning has occurred!" );
 * \endcode
 *
 */
#define SLIC_WARNING_ROOT(msg)                                        \
  do                                                                  \
  {                                                                   \
    if(axom::slic::isRoot())                                          \
    {                                                                 \
      std::ostringstream __oss;                                       \
      __oss << msg;                                                   \
      axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
      if(axom::slic::isAbortOnWarningsEnabled())                      \
      {                                                               \
        axom::slic::abort();                                          \
      }                                                               \
    }                                                                 \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_WARNING_ROOT_IF( EXP, msg )
 * \brief Macro that logs given warning message only on root iff EXP is true.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_WARNING_ROOT_IF macro is always active.
 * \note By default, all ranks are considered to be root.
 *       Must call `axom::slic::initialize(is_root={true|false})`
 *       or set via `axom::slic::setIsRoot({true|false})` to filter based on root.
 *
 * Usage:
 * \code
 *   SLIC_WARNING_ROOT_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_WARNING_ROOT_IF(EXP, msg)                                  \
  do                                                                    \
  {                                                                     \
    if(EXP)                                                             \
    {                                                                   \
      if(axom::slic::isRoot())                                          \
      {                                                                 \
        std::ostringstream __oss;                                       \
        __oss << msg;                                                   \
        axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
        if(axom::slic::isAbortOnWarningsEnabled())                      \
        {                                                               \
          axom::slic::abort();                                          \
        }                                                               \
      }                                                                 \
    }                                                                   \
  } while(axom::slic::detail::false_value)

///@}

// Use complete debug macros when not on device
#if defined(AXOM_DEBUG) && !defined(AXOM_DEVICE_CODE)

  //-----------------------------------------------------------------------------
  /// @{
  //! \name SLIC_ASSERT MACROS
  //!
  //! \collective
  //! \attention These assert macros are collective operations.
  //! All ranks in the user-supplied communicator must call the macro
  //! when used within an MPI distributed environment, and slic::enableAbortOnError()
  //! is called for the current active logger (default is enabled
  //! for loggers)
  //!
  //! \sa axom::slic::isAbortOnErrorsEnabled()
  //! \sa axom::slic::setAbortOnError(bool status)
  //!

  /*!
 * \def SLIC_ASSERT( EXP )
 * \brief Asserts that a given expression is true. If the expression is not true
 *  an error will be logged and the application will be aborted.
 *
 * \param [in] EXP user-supplied boolean expression.
 *
 * \warning This macro calls processAbort() iff EXP is false.
 * \note This macro is only active when AXOM_DEBUG is defined.
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
        if(axom::slic::isAbortOnErrorsEnabled())                      \
        {                                                             \
          axom::slic::abort();                                        \
        }                                                             \
      }                                                               \
    } while(axom::slic::detail::false_value)

  /*!
 * \def SLIC_ASSERT_MSG( EXP, msg )
 * \brief Same as SLIC_ASSERT, but with a custom error message.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 *
 * \warning This macro calls processAbort() iff EXP is false.
 * \note This macro is only active when AXOM_DEBUG is defined.
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
        if(axom::slic::isAbortOnErrorsEnabled())                             \
        {                                                                    \
          axom::slic::abort();                                               \
        }                                                                    \
      }                                                                      \
    } while(axom::slic::detail::false_value)

  ///@}

  //-----------------------------------------------------------------------------
  /// @{
  //! \name SLIC_CHECK MACROS
  //!
  //! \collective
  //! \attention These check macros can be set as collective operations.
  //! These check macros are collective if either:
  //! - slic::debug::checksAreErrors is set to true (default is false) and
  //! slic::enableAbortOnError() is called for the current active logger (default is
  //! enabled for loggers)
  //! - slic::debug::checksAreErrors is set to false (default is false) and
  //! slic::enableAbortOnWarning() is called for the current active logger (default is
  //! disabled for loggers)
  //!
  //! These check macros must then be called by all ranks in the
  //! user-supplied communicator when used within an MPI distributed
  //! environment.
  //!
  //! \sa axom::slic::isAbortOnErrorsEnabled()
  //! \sa axom::slic::setAbortOnError(bool status)
  //! \sa axom::slic::isAbortOnWarningsEnabled()
  //! \sa axom::slic::setAbortOnWarning(bool status)
  //!

  /*!
 * \def SLIC_CHECK( EXP )
 * \brief Checks that a given expression is true. If the expression is not true
 *  a warning is logged, but, in contrast to the similar SLIC_ASSERT macro the
 *  application is not aborted.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \note This macro is only active when AXOM_DEBUG is defined.
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
          if(axom::slic::isAbortOnErrorsEnabled())                        \
          {                                                               \
            axom::slic::abort();                                          \
          }                                                               \
        }                                                                 \
        else                                                              \
        {                                                                 \
          axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__); \
          if(axom::slic::isAbortOnWarningsEnabled())                      \
          {                                                               \
            axom::slic::abort();                                          \
          }                                                               \
        }                                                                 \
      }                                                                   \
    } while(axom::slic::detail::false_value)

  /*!
 * \def SLIC_CHECK_MSG( EXP, msg )
 * \brief Same as SLIC_CHECK, but with a custom error message.
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 *
 * \note This macro is only active when AXOM_DEBUG is defined.
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
          if(axom::slic::isAbortOnErrorsEnabled())                          \
          {                                                                 \
            axom::slic::abort();                                            \
          }                                                                 \
        }                                                                   \
        else                                                                \
        {                                                                   \
          axom::slic::logWarningMessage(__oss.str(), __FILE__, __LINE__);   \
          if(axom::slic::isAbortOnWarningsEnabled())                        \
          {                                                                 \
            axom::slic::abort();                                            \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } while(axom::slic::detail::false_value)

/// @}

// Use assert when on device (HIP does not yet support assert())
#elif defined(AXOM_DEBUG) && defined(AXOM_DEVICE_CODE)
  #if !defined(__HIP_DEVICE_COMPILE__)
    #define SLIC_ASSERT(EXP) assert(EXP)
    #define SLIC_ASSERT_MSG(EXP, msg) assert(EXP)
    #define SLIC_CHECK(EXP) assert(EXP)
    #define SLIC_CHECK_MSG(EXP, msg) assert(EXP)
  #else
// AXOM_DEBUG is defined so the expression and message cannot be ignored
// as they probably contain parameters that were decorated with AXOM_UNUSED_PARAM
// which will require the expression and message to be used to avoid warnings.
// Here we use them in a way that will not really have a runtime effect.
// The msg code block may contain << operators so we need a stream-like object
// but we enclose it in "if(false)" so it never executes.
namespace axom
{
namespace slic
{
namespace internal
{
class blackhole
{ };
// Write an object to a blackhole.
template <typename T>
AXOM_HOST_DEVICE inline blackhole &operator<<(blackhole &bh, T)
{
  return bh;
}
}  // namespace internal
}  // namespace slic
}  // namespace axom

    #define SLIC_ASSERT(EXP) ((void)(EXP))
    #define SLIC_ASSERT_MSG(EXP, msg)            \
      {                                          \
        if(false)                                \
        {                                        \
          ((void)(EXP));                         \
          axom::slic::internal::blackhole __oss; \
          __oss << msg;                          \
        }                                        \
      }
    #define SLIC_CHECK(EXP) ((void)(EXP))
    #define SLIC_CHECK_MSG(EXP, msg)             \
      {                                          \
        if(false)                                \
        {                                        \
          ((void)(EXP));                         \
          axom::slic::internal::blackhole __oss; \
          __oss << msg;                          \
        }                                        \
      }
  #endif
#else  // turn off debug macros and asserts

  #define SLIC_ASSERT(ignore_EXP) ((void)0)
  #define SLIC_ASSERT_MSG(ignore_EXP, ignore_msg) ((void)0)
  #define SLIC_CHECK(ignore_EXP) ((void)0)
  #define SLIC_CHECK_MSG(ignore_EXP, ignore_msg) ((void)0)

#endif /* END ifdef AXOM_DEBUG */

/*!
 * \def SLIC_INFO( msg )
 * \brief Logs an Info message.
 *
 * \param [in] msg user-supplied message
 *
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
 * \def SLIC_INFO_TAGGED( msg, tag )
 * \brief Logs an Info message to a tagged stream
 *
 * \param [in] msg user-supplied message
 * \param [in] tag user-supplied tag
 *
 * \note The SLIC_INFO_TAGGED macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_INFO_TAGGED("informative text goes here", "tag");
 * \endcode
 *
 */
#define SLIC_INFO_TAGGED(msg, tag)                    \
  do                                                  \
  {                                                   \
    std::ostringstream __oss;                         \
    __oss << msg;                                     \
    axom::slic::logMessage(axom::slic::message::Info, \
                           __oss.str(),               \
                           tag,                       \
                           __FILE__,                  \
                           __LINE__,                  \
                           false,                     \
                           true);                     \
  } while(axom::slic::detail::false_value)

/*!
 * \def SLIC_INFO_IF( EXP, msg )
 * \brief Logs an Info message iff EXP is true
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
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

/*!
 * \def SLIC_INFO_ROOT( msg )
 * \brief Logs an Info message if on root
 *
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_INFO_ROOT macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_INFO_ROOT( "informative text goes here" );
 * \endcode
 *
 */
#define SLIC_INFO_ROOT(msg) SLIC_INFO_IF(axom::slic::isRoot(), msg)

/*!
 * \def SLIC_INFO_ROOT_IF( EXP, msg )
 * \brief Logs an Info message if on root and iff EXP is true
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_INFO_ROOT_IF macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_INFO_ROOT_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
#define SLIC_INFO_ROOT_IF(EXP, msg) \
  SLIC_INFO_IF((EXP) && (axom::slic::isRoot()), msg)

#ifdef AXOM_DEBUG

  /*!
 * \def SLIC_DEBUG( msg )
 * \brief Logs a Debug message.
 *
 * \param [in] msg user-supplied message
 *
 * \note The SLIC_Debug macro is active when AXOM_DEBUG is defined.
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
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_DEBUG_IF macro is active when AXOM_DEBUG is defined.
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

  /*!
 * \def SLIC_DEBUG_ROOT( msg )
 * \brief Logs a Debug message if on root
 *
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_DEBUG_ROOT macro is active when AXOM_DEBUG is defined.
 *
 * Usage:
 * \code
 *   SLIC_DEBUG_ROOT( "informative text goes here" );
 * \endcode
 *
 */
  #define SLIC_DEBUG_ROOT(msg) SLIC_DEBUG_IF(axom::slic::isRoot(), msg)

  /*!
 * \def SLIC_DEBUG_ROOT_IF( EXP, msg )
 * \brief Logs a Debug message if on root and iff EXP is true
 *
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message.
 *
 * \note The SLIC_DEBUG_ROOT_IF macro is active when AXOM_DEBUG is defined.
 *
 * Usage:
 * \code
 *   SLIC_DEBUG_ROOT_IF( (val < 0), "my_val should always be positive" );
 * \endcode
 *
 */
  #define SLIC_DEBUG_ROOT_IF(EXP, msg) \
    SLIC_DEBUG_IF((EXP) && (axom::slic::isRoot()), msg)

#else  // turn off debug macros

  #define SLIC_DEBUG(ignore_EXP) ((void)0)
  #define SLIC_DEBUG_IF(ignore_EXP, ignore_msg) ((void)0)
  #define SLIC_DEBUG_ROOT(ignore_EXP) ((void)0)
  #define SLIC_DEBUG_ROOT_IF(ignore_EXP, ignore_msg) ((void)0)

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
