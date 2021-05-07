// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file slic.hpp
 */

#ifndef SLIC_HPP_
#define SLIC_HPP_

#include "axom/config.hpp"
#include "axom/slic/core/Logger.hpp"
#include "axom/slic/core/LogStream.hpp"
#include "axom/slic/core/MessageLevel.hpp"
#include "axom/slic/interface/slic_macros.hpp"

// C/C++ includes
#include <iostream>  // for std::endl, std::ends
#include <sstream>   // for std::ostringstream

namespace axom
{
namespace slic
{
struct debug
{
  AXOM_EXPORT static bool checksAreErrors;
};

/*!
 * \brief Initializes the SLIC logging environment.
 */
void initialize();

/*!
 * \brief Checks if the SLIC logging environment is initialized.
 * \return status true if initialized, else, false.
 */
bool isInitialized();

/*!
 * \brief Creates a new logger associated with the given name.
 * \param [in] name the name to associate with the new logger.
 * \param [in] imask inheritance mask, indicates the log level streams that
 *  will be inherited from the "root" logger. By default, nothing is inherited.
 */
void createLogger(const std::string& name, char imask = inherit::nothing);

/*!
 * \brief Activates the logger associated with the given name.
 * \param [in] name the name of the logger to activate.
 * \return True if the named logger was activated, False otherwise
 */
bool activateLogger(const std::string& name);

/*!
 * \brief Returns the name of the active logger.
 * \return s a string corresponding to the name of the active logger.
 */
std::string getActiveLoggerName();

/*!
 * \brief Sets desired logging level.
 * \param [in] level user-supplied level to log.
 */
void setLoggingMsgLevel(message::Level level);

/*!
 * \brief Gets the current logging level.
 */
message::Level getLoggingMsgLevel();

/*!
 * \brief Toggles the abort behavior for errors on the current active logger.
 * \param [in] status user-supplied flag.
 */
void setAbortOnError(bool status);

/*!
 * \brief Enables aborts on error messages for the current active logger.
 * \note This is equivalent to calling slic::setAbortOnError( true )
 * \post slic::isAbortOnErrorsEnabled() == true.
 * \pre slic::isInitialized() == true
 */
void enableAbortOnError();

/*!
 * \brief Disables aborts on error messages for the current active logger.
 * \note this is equivalent to calling slic::setAbortOnError( false )
 * \post slic::isAbortOnErrorsEnabled() == false.
 * \pre slic::isInitialized() == true
 */
void disableAbortOnError();

/*!
 * \brief Checks whether aborts on errors are enabled for the current logger.
 * \return status true if aborts on errors are enabled, otherwise, false.
 * \pre slic::isInitialized() == true.
 */
bool isAbortOnErrorsEnabled();

/*!
 * \brief Toggles the abort behavior for warnings on the current active logger.
 * \param [in] status user-supplied flag.
 */
void setAbortOnWarning(bool status);

/*!
 * \brief Enables aborts on warnings messages for the current active logger.
 * \note This is equivalent to calling slic::setAbortOnWarning( true )
 * \post slic::isAbortOnWarningsEnabled() == true.
 * \pre slic::isInitialized() == true.
 */
void enableAbortOnWarning();

/*!
 * \brief Disables aborts on warnings messages for the current active logger.
 * \note This is equivalent to calling slic::setAbortOnWarnings( false )
 * \post slic::isAbortOnWarnigsEnabled() == true.
 * \pre slic::isInitialized() == true.
 */
void disableAbortOnWarning();

/*!
 * \brief Checks whether aborts on warnings are enabled for the current logger.
 * \return status true if aborts on warnings are enabled, otherwise, false.
 * \pre slic::isInitialized() == true.
 */
bool isAbortOnWarningsEnabled();

/*!
 * \brief Sets the function to call when program abort is requested
 * \param [in] abort_func The user-specified function to call
 * \pre slic::isInitialized() == true.
 */
void setAbortFunction(AbortFunctionPtr abort_func);

/*!
 * \brief Adds the given stream to the the given level.
 * \param [in] ls pointer to the log stream.
 * \param [in] level the level to log.
 * \pre ls != nullptr
 */
void addStreamToMsgLevel(LogStream* ls, message::Level level);

/*!
 * \brief Adds the given stream to all levels.
 * \param [in] ls pointer to the log stream.
 * \pre ls != nullptr.
 */
void addStreamToAllMsgLevels(LogStream* ls);

/*!
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                bool filter_duplicates = false);

/*!
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied associated with this message.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                bool filter_duplicates = false);

/*!
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& fileName,
                int line,
                bool filter_duplicates = false);

/*!
 * \brief Logs the given message to all registered streams.
 * \param [in] level the level of the message being logged.
 * \param [in] message user-supplied message.
 * \param [in] tag user-supplied tag associated with the message.
 * \param [in] fileName the name of the file this message is logged form.
 * \param [in] line the line number within the file this message is logged.
 * \param [in] filter_duplicates optional parameter that indicates whether
 * duplicate messages resulting from running in parallel will be filtered out.
 * Default is false.
 */
void logMessage(message::Level level,
                const std::string& message,
                const std::string& tag,
                const std::string& fileName,
                int line,
                bool filter_duplicates = false);

/*!
 * \brief Convenience method to log an error message.
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file that the message is logged.
 */
void logErrorMessage(const std::string& message,
                     const std::string& fileName,
                     int line);

/*!
 * \brief Convenience method to log warning messages.
 * \param [in] message user-supplied message.
 * \param [in] fileName the name of the file this message is logged from.
 * \param [in] line the line number within the file that the message is logged.
 */
void logWarningMessage(const std::string& message,
                       const std::string& fileName,
                       int line);

/*!
 * \brief Flushes all streams.
 * \see Logger::flushStreams.
 */
void flushStreams();

/*!
 * \brief Pushes all streams.
 * \see Logger::pushStreams.
 */
void pushStreams();

/*!
 * \brief Finalizes the slic logging environment.
 */
void finalize();

/*!
 * \brief Uses glibc's backtrace() functionality to return a stacktrace
 * \return s a string corresponding to the stacktrace.
 */
std::string stacktrace();

} /* namespace slic */

} /* namespace axom */

#endif /* SLIC_HPP_ */
