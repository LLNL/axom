// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file LogStream.hpp
 *
 */

#ifndef LOGSTREAM_HPP_
#define LOGSTREAM_HPP_

#include "axom/slic/core/MessageLevel.hpp"
#include "axom/core/Macros.hpp"

/// \name Wildcards
/// @{

#define MSG_IGNORE_TAG ""
#define MSG_IGNORE_FILE ""
#define MSG_IGNORE_LINE -1

/// @}

// C/C++ includes
#include <string>  // For STL string

namespace axom
{
namespace slic
{
/*!
 * \class LogStream
 *
 * \brief Abstract base class defining the interface of all LogStreams.
 *
 *  Concrete instances implement the append() method accordingly, e.g., to
 *  append the message to a file, write it to stdout, etc.
 *
 * \see Logger Console
 */
class LogStream
{
public:
  LogStream();
  virtual ~LogStream();

  /*!
   * \brief Sets the format string.
   * \param format a format string
   * \note The following keywords in the format string are replaced:
   *  <ul>
   *    <li> <LEVEL> with the message type, e.g, ERROR, FATAL, etc. </li>
   *    <li> <MESSAGE> with the user-supplied message </li>
   *    <li> <TAG> user-supplied tag </li>
   *    <li> <FILE> with the filename </li>
   *    <li> <LINE> with the line number </li>
   *    <li> <LINE> with the MPI rank </li>
   *    <li> <TIMESTAMP> date/time the message is logged </li>
   *  </ul>
   *
   *  An example illustrating how an application can setup the format string
   *  is shown below:
   * \code
   *  std::string format =
   *         std::string( "***********************************\n" )+
   *         std::string( "* <TIMESTAMP>\n\n" ) +
   *         std::string( "* LEVEL=<LEVEL>\n" ) +
   *         std::string( "* MESSAGE=<MESSAGE>\n" ) +
   *         std::string( "* FILE=<FILE>\n" ) +
   *         std::string( "* LINE=<LINE>\n" ) +
   *         std::string( "* RANK=<RANK>\n" ) +
   *         std::string( "***********************************\n" );
   * \endcode
   */
  void setFormatString(const std::string& format) { m_formatString = format; }

  /*!
   * \brief Appends the given message to the stream.
   *
   * \param [in] msgLevel the level of the message.
   * \param [in] message the user-supplied message.
   * \param [in] tagName user-supplied tag to associate with the given message.
   * \param [in] fileName the file where this message is appended
   * \param [in] line the line within the file at which the message is appended.
   * \param [in] filter_duplicates optional parameter that indicates whether
   * duplicate messages resulting from running in parallel will be filtered out.
   *
   * \note The following wildcards may be used to ignore a particular field:
   * <ul>
   *   <li> MSG_IGNORE_TAG  </li>
   *   <li> MSG_IGNORE_FILE </li>
   *   <li> MSG_IGNORE_LINE </li>
   * </ul>
   */
  virtual void append(message::Level msgLevel,
                      const std::string& message,
                      const std::string& tagName,
                      const std::string& fileName,
                      int line,
                      bool filter_duplicates) = 0;

  /*!
   * \brief Flushes the log stream. It's a NO-OP by default.
   * \note The intent of this method is to be overridden by concrete
   *  implementations. This is primarily useful for applications running
   *  in a distributed MPI environment, where the flush is a collective
   *  operation intended for a synchronization checkpoint.
   */
  virtual void flush() {};

  /*!
   * \brief Pushes messages incrementally up the log stream. NO-OP by default.
   * \note The intent of this method is to be overridden by concrete
   *  implementations that need to be incrementally advanced. This is primarily
   *  useful for applications running in a distributed MPI environment, where
   *  the push is a collective operation intended for a incrementally advancing
   *  messages through the log stream.
   */
  virtual void push() {};

protected:
  /*!
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
   */
  std::string getFormatedMessage(const std::string& msgLevel,
                                 const std::string& message,
                                 const std::string& tagName,
                                 const std::string& rank,
                                 const std::string& fileName,
                                 int line);

  /*!
   * \brief Returns a time-stamp.
   * \return str a textual representation of the current time.
   */
  std::string getTimeStamp();

private:
  std::string m_formatString;

  /*!
   * \brief Replaces the given key in the message string with the given value.
   * \param [in,out] msg the message string that will be modified.
   * \param [in] key the key in the message that will be replace.
   * \param [in] value the value to replace it with.
   */
  void replaceKey(std::string& msg,
                  const std::string& key,
                  const std::string& value,
                  std::size_t pos = std::string::npos);

  DISABLE_COPY_AND_ASSIGNMENT(LogStream);
  DISABLE_MOVE_AND_ASSIGNMENT(LogStream);
};

} /* namespace slic */

} /* namespace axom */

#endif /* LOGSTREAM_HPP_ */
