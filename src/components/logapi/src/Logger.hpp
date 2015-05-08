/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
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

// C/C++ includes
#include <string> // for STL string

namespace asctoolkit {

namespace logapi {

/*!
 *******************************************************************************
 * \enum MessageType
 *
 * \brief Enumerates the different types of messaging.
 *
 * \note The ordering of the enumerators reflects the level of severity of
 *  a message.
 *
 * \see Logger
 *******************************************************************************
 */
enum MessageType {
  Fatal,        //!< FATAL log a non-recoverable event.
  Error,        //!< ERROR log an error that *may* be recoverable.
  Warning,      //!< WARNING log a warning.
  Info,         //!< INFO log information that is useful for users & developers.
  Debug,        //!< DEBUG log information that is useful for developers.

  Num_Msg_Types //!< Num_Msg_Types
};

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
   * \brief Enables log messages of the given type.
   * \param [in] type message type.
   * \pre type >= FATAL && type < Num_Msg_Types
   *****************************************************************************
   */
  void enable( int type );

  /*!
   *****************************************************************************
   * \brief Disables log messages of the given type.
   * \param [in] type message type.
   * \pre type >= FATAL && type < Num_Msg_types
   *****************************************************************************
   */
  void disable( int type );

  /*!
   *****************************************************************************
   * \brief Sets the LogStream for the given message type.
   * \param [in] type the message type/level.
   * \param [in] ls pointer to an application log stream.
   * \pre type >= FATAL && type < Num_Msg_Types
   * \pre ls != NULL
   *****************************************************************************
   */
  void setLogStream( int type, LogStream* ls);

  /*!
   *****************************************************************************
   * \brief Logs the given message
   * \param [in] type
   * \param [in] message
   * \param [in] fileName
   * \param [in] line
   * \pre type >= FATAL && type < Num_Msg_Types
   * \pre m_Streams[ type ] != NULL
   *****************************************************************************
   */
  void logMessage( int type,
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

  bool m_StreamState[ Num_Msg_Types ];
  LogStream* m_Streams[ Num_Msg_Types ];

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
