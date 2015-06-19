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
   * \note The Logger takes ownership of the LogStream object.
   * \pre ls != NULL.
   *****************************************************************************
   */
  void addStreamToLevel( LogStream* ls, message::Level level );

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
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   *****************************************************************************
   */
  void logMessage( message::Level level, const std::string& message );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] tagName user-supplied tag to associated with the given message.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& tagName );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] fileName name of the file this call is made from.
   * \param [in] line line within the file that this call is made from.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& fileName,
                   int line );

  /*!
   *****************************************************************************
   * \brief Logs the given message to all registered streams.
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] tagName user-supplied tag to associated with the given message.
   * \param [in] fileName name of the file this call is made from.
   * \param [in] line line within the file that this call is made from.
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
                   const std::string& tagName,
                   const std::string& fileName,
                   int line );

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
   * \brief Initializes the logging environment.
   * \post Logger::getInstance() != NULL.
   *****************************************************************************
   */
  static void initialize();

  /*!
   *****************************************************************************
   * \brief Finalizes the logging environment.
   * \post Logger::getInstance() == NULL.
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
  static Logger* getInstance();

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

  ///@}

  /// \name Disabled Methods
  ///@{

  Logger( const Logger& ); // Not implemented
  Logger& operator=( const Logger& ); // Not implemented

  ///@}

};

} /* namespace slic */

} /* namespace asctoolkit */

#endif /* LOGGER_HPP_ */
