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
 * \file LogStream.hpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef LOGSTREAM_HPP_
#define LOGSTREAM_HPP_

#include "logapi/MessageLevel.h"

/// \name Wildcards
/// @{

#define MSG_IGNORE_TAG ""
#define MSG_IGNORE_FILE ""
#define MSG_IGNORE_LINE -1

/// @}

// C/C++ includes
#include <string> // For STL string

namespace asctoolkit {

namespace logapi {

/*!
 *******************************************************************************
 * \class LogStream
 *
 * \brief Abstract base class defining the interface of all LogStreams.
 *
 *  Concrete instances implement the append() method accordingly, e.g., to
 *  append the message to a file, write it to stdout, etc.
 *
 * \see Logger Console
 *******************************************************************************
 */
class LogStream
{

public:
  LogStream();
  virtual ~LogStream();

  /*!
   *****************************************************************************
   * \brief Sets the format string.
   * \param fmt a format string
   * \note The following keywords in the format string are replaced:
   *  <ul>
   *    <li> <LEVEL> with the message type, e.g, ERROR, FATAL, etc. </li>
   *    <li> <MESSAGE> with the user-supplied message </li>
   *    <li> <TAG> user-supplied tag </li>
   *    <li> <FILE> with the filename </li>
   *    <li> <LINE> with the line number </li>
   *  </ul>
   *****************************************************************************
   */
  void setFormatString(const std::string& format) { m_formatString = format; }

  /*!
   *****************************************************************************
   * \brief Appends the given message to the stream.
   *
   * \param [in] msgLevel the level of the message.
   * \param [in] message the user-supplied message.
   * \param [in] tagName user-supplied tag to associate with the given message.
   * \param [in] fileName the file where this message is appended
   * \param [in] line the line within the file at which the message is appended.
   *
   * \note The following wildcards may be used to ignore a particular field:
   * <ul>
   *   <li> MSG_IGNORE_TAG  </li>
   *   <li> MSG_IGNORE_FILE </li>
   *   <li> MSG_IGNORE_LINE </li>
   * </ul>
   *****************************************************************************
   */
  virtual void append( message::Level msgLevel,
                       const std::string& message,
                       const std::string& tagName,
                       const std::string& fileName,
                       int line
                       ) = 0;

  /*!
   *****************************************************************************
   * \brief Flushes the log stream. It's a NO-OP by default.
   * \note The intent of this method is to be overridden by concrete
   *  implementations. This is primarily useful for applications running
   *  in a distributed MPI environment, where the flush is a collective
   *  operation intended for a synchronization checkpoint.
   *****************************************************************************
   */
  virtual void flush() { };

protected:

  /*!
   *****************************************************************************
   * \brief Returns the formatted message as a single string.
   * \param [in] msgLevel the level of the given message.
   * \param [in] message the user-supplied message.
   * \param [in] tagName user-supplied tag, may be MSG_IGNORE_TAG
   * \param [in] fileName filename where this message is logged, may be
   *  MSG_IGNORE_FILE to ignore this field.
   * \param [in] line the line number within the file where the message is
   *  logged. Likewise, may be set to MSG_IGNORE_LINE to ignore this field.
   * \return str the formatted message string.
   * \post str != "".
   *****************************************************************************
   */
   std::string getFormatedMessage( const std::string& msgLevel,
                                   const std::string& message,
                                   const std::string& tagName,
                                   const std::string& fileName,
                                   int line );

private:

  std::string m_formatString;

  /*!
   *****************************************************************************
   * \brief Replaces the given key in the message string with the given value.
   * \param [in,out] msg the message string that will be modified.
   * \param [in] key the key in the message that will be replace.
   * \param [in] value the value to replace it with.
   *****************************************************************************
   */
  void replaceKey( std::string& msg,
                   const std::string& key,
                   const std::string& value );

  /// \name Disabled Methods
  ///@{

  LogStream( const LogStream& ); // Not implemented
  LogStream& operator=( const LogStream& ); // Not implemented

  ///@}
};

} /* namespace logapi */

} /* namespace asctoolkit */

#endif /* LOGSTREAM_HPP_ */
