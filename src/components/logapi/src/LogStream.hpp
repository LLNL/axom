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

#include "logapi/MessageType.h"

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
   *    <li> <MTYPE> with the message type, e.g, ERROR, FATAL, etc. </li>
   *    <li> <FILE> with the filename </li>
   *    <li> <LINE> with the line number </li>
   *    <li> <MESSAGE> with the user-supplied message </li>
   *  </ul>
   *****************************************************************************
   */
  void setFormatString(const std::string& fmt) { m_fmtString = fmt; }

  /*!
   *****************************************************************************
   * \brief Appends the given message to the stream.
   * \param [in] msgType the type of the message.
   * \param [in] msgTypeName string representation of the message type.
   * \param [in] message the user-supplied message.
   * \param [in] fileName the file where this message is appended
   * \param [in] line the line within the file at which the message is appended.
   *****************************************************************************
   */
  virtual void append( MessageType msgType,
                       const std::string& msgTypeName,
                       const std::string& message,
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
   * \return
   *****************************************************************************
   */
   std::string getFormatedMessage( const std::string& msgTypeName,
                                   const std::string& message,
                                   const std::string& fileName,
                                   int line );

private:

  std::string m_fmtString;

  /// \name Disabled Methods
  ///@{

  LogStream( const LogStream& ); // Not implemented
  LogStream& operator=( const LogStream& ); // Not implemented

  ///@}
};

} /* namespace logapi */

} /* namespace asctoolkit */

#endif /* LOGSTREAM_HPP_ */
