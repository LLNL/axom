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
#include "slic/LogStream.hpp"
#include "slic/MessageLevel.h"

#include "common/CommonTypes.hpp" // for ATK_NULLPTR
#include "common/Utilities.hpp"   // for utilities::processAbort()

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
    asctoolkit::slic::logMessage(                                             \
        asctoolkit::slic::message::Fatal,msg,__FILE__, __LINE__ );            \
    if ( asctoolkit::slic::getAbortOnError() ) {                             \
       asctoolkit::utilities::processAbort();                                 \
    }                                                                         \
  }                                                                           \
} while ( 0 )

/// @}

/*!
 ******************************************************************************
 * \def SLIC_ERROR_MSG( EXP, msg )
 * \brief Same as SLIC_ERROR, but with a custom error message.
 * \param [in] EXP user-supplied boolean expression.
 * \param [in] msg user-supplied message
 * \note The SLIC_ERROR_MSG is always active.
 * \warning This macro calls processAbort() if EXP is false.
 * \see SLIC_ERROR( EXP )
 *
 * Usage:
 * \code
 *   SLIC_ERROR_MSG( my_val >= 0, "my_val must always be positive" );
 * \endcode
 *
 ******************************************************************************
 */
#define SLIC_ERROR_MSG( EXP, msg )                                            \
do {                                                                          \
  if ( !(EXP) ) {                                                             \
    std::ostringstream oss;                                                   \
    oss << "Failed Error: " << # EXP << std::endl << msg << std::ends;        \
    asctoolkit::slic::logMessage(                                             \
        asctoolkit::slic::message::Fatal,oss.str(),__FILE__,__LINE__ );       \
    if ( asctoolkit::slic::getAbortOnError() ) {                             \
       asctoolkit::utilities::processAbort();                                 \
    }                                                                         \
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
    asctoolkit::slic::logMessage(                                             \
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
    asctoolkit::slic::logMessage(                                             \
        asctoolkit::slic::message::Fatal,oss.str(),__FILE__,__LINE__ );       \
    if ( asctoolkit::slic::getAbortOnAssert() ) {                            \
       asctoolkit::utilities::processAbort();                                 \
    }                                                                         \
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
    asctoolkit::slic::logMessage(                                             \
        asctoolkit::slic::message::Fatal,oss.str(),__FILE__,__LINE__ );       \
    if ( asctoolkit::slic::getAbortOnAssert() ) {                            \
       asctoolkit::utilities::processAbort();                                 \
    }                                                                         \
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
    asctoolkit::slic::logMessage(                                             \
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
    asctoolkit::slic::logMessage(                                             \
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

namespace asctoolkit {


namespace slic {

struct RuntimeAbortBehavior
{
   static bool willAbortOnAssert;
   static bool willAbortOnError;
};


/*!
 *******************************************************************************
 * \brief Initializes the SLIC logging environment.
 *******************************************************************************
 */
void initialize();

/*!
 *******************************************************************************
 * \brief Checks if the SLIC logging environment is initialized.
 * \return status true if initialized, else, false.
 *******************************************************************************
 */
bool isInitialized();

/*!
 *******************************************************************************
 * \brief Creates a new logger associated with the given name.
 * \param [in] name the name to associate with the new logger.
 * \param [in] imask inheritance mask, indicates the log level streams that
 *  will be inherited from the "root" logger. By default, nothing is inherited.
 *******************************************************************************
 */
void createLogger( const std::string& name,
                   char imask=inherit::nothing );

/*!
 *******************************************************************************
 * \brief Activates the logger associated with the given name.
 * \param [in] name the name of the logger to activate.
 *******************************************************************************
 */
void activateLogger( const std::string& name );

/*!
 *******************************************************************************
 * \brief Sets desired logging level.
 * \param [in] level user-supplied level to log.
 *******************************************************************************
 */
void setLoggingLevel( message::Level level );

/*!
 *******************************************************************************
 * \brief Sets abort behavior when a SLIC_ASSERT is evaluated to false.  The
 *  default setting is for SLIC_ASSERT to abort the process.
 * \param [in] sets whether a SLIC_ASSERT failed evaluation will abort the
 * process.
 *******************************************************************************
 */
void setAbortOnAssert( bool willAbort );

/*!
 *******************************************************************************
 * \brief Gets abort behavior when a SLIC_ASSERT is evaluated to false.  The
 *  default setting is for SLIC_ASSERT to abort the process.
 * \return if true, a SLIC_ASSERT failed evaluation will abort the process.
 *******************************************************************************
 */
bool getAbortOnAssert();

/*!
 *******************************************************************************
 * \brief Sets abort behavior when a SLIC_ERROR is evaluated to false.  The
 * default behavior is for SLIC_ERROR to abort a process.
 * \param [in] if true, a SLIC_ERROR failed evaluation will abort the process.
 *******************************************************************************
 */
void setAbortOnError( bool willAbort );

/*!
 *******************************************************************************
 * \brief Gets abort behavior when a SLIC_ERROR is evaluated to false.  The
 * default behavior is for SLIC_ERROR to abort a process.
 * \return if true, a SLIC_ERROR failed evaluation will abort the process.
 *******************************************************************************
 */
bool getAbortOnError();

/*!
 *******************************************************************************
 * \brief Adds the given stream to the the given level.
 * \param [in] ls pointer to the log stream.
 * \param [in] level the level to log.
 * \pre ls != ATK_NULLPTR
 *******************************************************************************
 */
void addStreamToLevel( LogStream* ls, message::Level level);

/*!
 *******************************************************************************
 * \brief Adds the given stream to all levels.
 * \param [in] ls pointer to the log stream.
 * \pre ls != ATK_NULLPTR.
 *******************************************************************************
 */
void addStreamToAllLevels( LogStream* ls );

/*!
 *******************************************************************************
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] filter_dulicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 *******************************************************************************
 */
void logMessage( message::Level level,
                 const std::string& message,
                 bool filter_duplicates=false );

/*!
 *******************************************************************************
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied associated with this message.
 * \param [in] filter_dulicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 *******************************************************************************
 */
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag,
                 bool filter_duplicates=false );

/*!
 *******************************************************************************
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_dulicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 *******************************************************************************
 */
void logMessage( message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates=false );

/*!
 *******************************************************************************
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied tag associated with the message.
 * \param [in] fileName the name of the file this message is logged form.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_dulicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 *******************************************************************************
 */
void logMessage( message::Level level,
                 const std::string& message,
                 const std::string& tag,
                 const std::string& fileName,
                 int line,
                 bool filter_duplicates=false );

/*!
 *******************************************************************************
 * \brief Flushes all streams.
 * \see Logger::flushStreams.
 *******************************************************************************
 */
void flushStreams();

/*!
 *******************************************************************************
 * \brief Finalizes the slic logging environment.
 *******************************************************************************
 */
void finalize();

} /* namespace slic */

} /* namespace asctoolkit */

#endif /* SLIC_HPP_ */
