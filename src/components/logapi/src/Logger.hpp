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

#include "logapi/MessageLevel.h"

// C/C++ includes
#include <string> // for STL string
#include <vector> // for STL vector

namespace asctoolkit {

namespace logapi {


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
 *  The Logger delegates the messages to the corresponding LogStream according
 *  to the type of the message. Messages can filtered out through calls to
 *  enable()/disable(). Messages are logged using the LogStream abstraction
 *  layer. This allows an application to use some of the predefined LogStream
 *  mechanisms or implement a custom one by deriving from the LogStream class.
 *  Each message type is associated with a concrete LogStream instance set by
 *  the application with an invocation to the setLogStream(type,logStream)
 *  method. Each message type may point to a different LogStream instance or
 *  all can point to the same stream.
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
   * \brief Adds the given log stream to this Logger instance.
   * \param [in] ls pointer to a user-supplied LogStream object.
   * \pre ls != NULL.
   *****************************************************************************
   */
  void addLogStream( LogStream* ls );


  /*!
   *****************************************************************************
   * \brief Logs the given message
   * \param [in] level the level of the given message.
   * \param [in] message the user-supplied message to log.
   * \param [in] fileName the name of the file that calls this log message.
   * \param [in] line the line number within the file that logs this message.
   * \pre type >= FATAL && type < Num_Msg_Types
   * \pre m_Streams[ type ] != NULL
   *****************************************************************************
   */
  void logMessage( message::Level level,
                   const std::string& message,
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
  void flushAllStreams();

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
  std::vector< LogStream* > m_logStreams;

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

} /* namespace logapi */

} /* namespace asctoolkit */

#endif /* LOGGER_HPP_ */
