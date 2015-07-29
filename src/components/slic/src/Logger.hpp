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
 * \file Logger.hpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include "slic/MessageLevel.h"

// C/C++ includes
#include <string> // for STL string
#include <vector> // for STL vector
#include <map>    // for STL map

#include "common/ATKMacros.hpp"

namespace asctoolkit {

namespace slic {


// Forward declarations
class LogStream;

/*!
 *******************************************************************************
 * \class Logger
 *
 * \brief A singleton class for logging error, debug and warning messages.
 *
 *  Applications using this logging facility must call initialize() and
 *  finalize(), typically in the beginning and end of the application's main()
 *  respectively. Logging messages can then be accomplished through calls to
 *  'logMessage(type,msg,file,line)' or through the convenience macro layer,
 *  which allows for logging messages to be compiled out.
 *
 *  The Logger appends the messages to the set of registered LogStreams.
 *  Messages are logged using the LogStream abstraction layer. This allows an
 *  application to use some of the predefined LogStream mechanisms or implement
 *  a custom one by implementing a derivative of the LogStream class.
 *
 * \see LogStream MessageType
 *******************************************************************************
 */
class Logger
{

public:

  /*!
   *****************************************************************************
   * \brief Sets the logging level to the given level. This controls which
   *  messages are logged based on severity. All messages with equal or higher
   *  severity to the given level will be logged.
   * \param [in] level the logging level.
   *****************************************************************************
   */
  void setLoggingLevel( message::Level level );

  /*!
   *****************************************************************************
   * \brief Binds the given stream to the given level for this Logger instance.
   * \param [in] ls pointer to the user-supplied LogStream object.
   * \param [in] level the level that this stream will be associated with.
   * \param [in] pass_ownership flag that indicates whether the given logger
   *  instance owns the supplied LogStream object. This parameter is optional.
   *  Default is true.
   * \note The Logger takes ownership of the LogStream object.
   * \pre ls != NULL.
   *****************************************************************************
   */
  void addStreamToLevel( LogStream* ls, message::Level level,
                         bool pass_ownership=true );

  /*!
   *****************************************************************************
   * \brief Binds the given stream to all the levels for this Logger instance.
   * \param [in] ls pointer to the user-supplied LogStream object.
   * \note The Logger takes ownership of the LogStream object.
   * \pre ls != NULL.
   *****************************************************************************
   */
  void addStreamToAllLevels( LogStream* ls );

  /*!
   *****************************************************************************
   * \brief Returns the number of streams at the given level.
   * \param [in] level the level in query.
   * \return N the number of streams at the given level.
   * \post N >= 0
   *****************************************************************************
   */
  int getNumStreamsAtLevel( message::Level level );

  /*!
   *****************************************************************************
   * \brief Returns the ith stream at the given level.
   * \param [in] level the level in query.
   * \param [in] i the index of the stream in query.
   * \return stream_ptr pointer to the stream.
   * \pre i >= 0 && i < this->getNumStreamsAtLevel( level )
   * \post stream_ptr != NULL.
   *****************************************************************************
   */
  LogStream* getStream( message::Level level, int i );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] filter_dulicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   * Default is false.
   *****************************************************************************
   */
  void logMessage( message::Level level, const std::string& message,
                   bool filter_duplicates=false );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] tagName user-supplied tag to associated with the given message.
   * \param [in] filter_dulicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   * Default is false.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& tagName,
                   bool filter_duplicates=false );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] fileName name of the file this call is made from.
   * \param [in] line line within the file that this call is made from.
   * \param [in] filter_dulicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   * Default is false.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& fileName,
                   int line,
                   bool filter_duplicates=false );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] tagName user-supplied tag to associated with the given message.
   * \param [in] fileName name of the file this call is made from.
   * \param [in] line line within the file that this call is made from.
   * \param [in] filter_dulicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   * Default is false.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& tagName,
                   const std::string& fileName,
                   int line,
                   bool filter_duplicates=false );

  /*!
   *****************************************************************************
   * \brief Flushes all streams.
   * \note When used within an MPI distributed environment, flushAllStreams is
   *  a collective operation. All ranks in the user-supplied communicator must
   *  call this method.
   *****************************************************************************
   */
  void flushStreams();

  /// \name Static Methods
  ///@{

  /*!
   *****************************************************************************
   * \brief Initializes the logging environment, with the root logger.
   * \post Logger::getActiveLogger() != NULL.
   *****************************************************************************
   */
  static void initialize();

  /*!
   *****************************************************************************
   * \brief Creates a new logger associated with the given name.
   * \param [in] name the name to associate with the new logger.
   * \param [in] imask inheritance mask, indicates the log level(s), which will
   *  be inherited from the "root" logger. By default, nothing is inherited.
   * \return status return status, true if the logger is created, else false.
   * \note False is returned if a logger associated with the given name
   *  already exists.
   *****************************************************************************
   */
  static bool createLogger(const std::string& name,
                           char imask=inherit::nothing );

  /*!
   *****************************************************************************
   * \brief Activates the logger with the associate name.
   * \param [in] name the name of the logger to activate.
   * \return status return status, true if the logger is activated, else false.
   * \note False is returned if the logger with the given name does not exist.
   *****************************************************************************
   */
  static bool activateLogger(const std::string& name);

  /*!
   *****************************************************************************
   * \brief Finalizes the logging environment.
   * \post Logger::getActiveLogger() == NULL.
   *****************************************************************************
   */
  static void finalize();

  /*!
   *****************************************************************************
   * \brief Returns a pointer to the logger instance.
   * \return logger pointer to the logger instance.
   * \pre s_Logger != NULL
   * \post logger != NULL
   *****************************************************************************
   */
  static Logger* getActiveLogger();

  /*!
   *****************************************************************************
   * \brief Returns the root logger
   * \return logger pointer to the root logger instance.
   *****************************************************************************
   */
  static Logger* getRootLogger();

  ///@}

private:

  /*!
   *****************************************************************************
   * \brief Default constructor, made private since this is a singleton.
   *****************************************************************************
   */
  Logger();

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~Logger();

  /// \name Private class members
  ///@{

  bool m_isEnabled[ message::Num_Levels ];
  std::map< LogStream*, LogStream* > m_streamObjectsManager;
  std::vector< LogStream* > m_logStreams[ message::Num_Levels ];

  ///@}

  /// \name Static Members
  ///@{

  static Logger* s_Logger;
  static std::map< std::string, Logger* > s_loggers;

  ///@}

  DISABLE_COPY_AND_ASSIGNMENT(Logger);

};

} /* namespace slic */

} /* namespace asctoolkit */

#endif /* LOGGER_HPP_ */
