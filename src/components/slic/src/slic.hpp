/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file slic.hpp
 *
 * \date May 31, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef SLIC_HPP_
#define SLIC_HPP_

#include "slic/Logger.hpp"
#include "common/Utilities.hpp" // for utilities::processAbort()

// C/C++ includes
#include <iostream> // for std::endl, std::ends
#include <sstream>  // for std::ostringstream


/// \name ERROR MACROS
/// @{

/*!
 ******************************************************************************
 * \def SLIC_ERROR( EXP, msg )
 * \brief Logs an error iff EXP is true and aborts the application.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 * \note The SLIC_ERROR macro is always active.
 * \warning This macro calls processAbort() if EXP is false.
 *
 * Usage:
 * \code
 *   SLIC_ERROR( my_val < 0, "my_val should always be positive" );
 * \endcode
 *
 ******************************************************************************
 */
#define SLIC_ERROR( EXP, msg )                                                \
do {                                                                          \
  if ( EXP ) {                                                                \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Fatal,msg,__FILE__, __LINE__ );            \
    asctoolkit::utilities::processAbort();                                    \
  }                                                                           \
} while ( 0 )

/// @}

/// \name WARNING MACROS
/// @{

/*!
 ******************************************************************************
 * \def SLIC_WARNING( EXP, msg )
 * \brief Logs a warning iff EXP is true.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 * \note The SLIC_WARNING macro is always active.
 *
 * Usage:
 * \code
 *   SLIC_WARNING( my_val < 0, "my_val should always be positive" );
 * \endcode
 *
 ******************************************************************************
 */
#define SLIC_WARNING( EXP, msg )                                              \
do {                                                                          \
  if ( EXP ) {                                                                \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Warning,msg,__FILE__, __LINE__ );          \
  }                                                                           \
} while ( 0 )

/// @}

#ifdef ATK_DEBUG

//-----------------------------------------------------------------------------
/// \name ASSERT MACROS
/// @{

/*!
 ******************************************************************************
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
 ******************************************************************************
 */
#define SLIC_ASSERT( EXP )                                                    \
do {                                                                          \
  if ( !(EXP) ) {                                                             \
    std::ostringstream oss;                                                   \
    oss << "Failed Assert: " << # EXP << std::ends;                           \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Fatal,oss.str(),__FILE__,__LINE__ );       \
    asctoolkit::utilities::processAbort();                                    \
  }                                                                           \
} while ( 0 )

/*!
 ******************************************************************************
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
 ******************************************************************************
 */
#define SLIC_ASSERT_MSG( EXP, msg )                                           \
do {                                                                          \
  if ( !(EXP) ) {                                                             \
    std::ostringstream oss;                                                   \
    oss << "Failed Assert: " << # EXP << std::endl << msg << std::ends;       \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Fatal,oss.str(),__FILE__,__LINE__ );       \
    asctoolkit::utilities::processAbort();                                    \
  }                                                                           \
} while ( 0 )

/// @}

//-----------------------------------------------------------------------------
/// \name DEBUG MACROS
/// @{

/*!
 ******************************************************************************
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
 ******************************************************************************
 */
#define SLIC_CHECK( EXP )                                                     \
do {                                                                          \
  if ( !(EXP) ) {                                                             \
    std::ostringstream oss;                                                   \
    oss << "Failed Check: " << # EXP << std::ends;                            \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Warning,oss.str(),__FILE__,__LINE__ );     \
  }                                                                           \
} while ( 0 )

/*!
 ******************************************************************************
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
 ******************************************************************************
 */
#define SLIC_CHECK_MSG( EXP, msg )                                            \
do {                                                                          \
  if ( !(EXP) ) {                                                             \
    std::ostringstream oss;                                                   \
    oss << "Failed Check: " << # EXP << std::endl << msg <<  std::ends;       \
    asctoolkit::slic::Logger::log(                                            \
        asctoolkit::slic::message::Warning,oss.str(),__FILE__,__LINE__ );     \
  }                                                                           \
} while ( 0 )

/// @}

#else // turn off debug macros and asserts

#define SLIC_ASSERT( ignore_EXP ) ( (void)0 )
#define SLIC_ASSERT_MSG( ignore_EXP, ignore_msg ) ( (void)0 )
#define SLIC_CHECK( ignore_EXP ) ( (void)0 )
#define SLIC_CHECK_MSG( ignore_EXP, ignore_msg ) ( (void)0 )

#endif /* END ifdef ATK_DEBUG */

#endif /* SLIC_HPP_ */
