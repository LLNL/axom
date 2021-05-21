// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file SynchronizedStream.hpp
 *
 */

#ifndef SYNCHRONIZEDSTREAM_HPP_
#define SYNCHRONIZEDSTREAM_HPP_

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/Macros.hpp"

// C/C++ includes
#include <iostream>  // for std::ostream

// MPI
#include <mpi.h>  // For MPI

namespace axom
{
namespace slic
{
/*!
 * \class SynchronizedStream
 *
 * \brief A concrete instance of LogStream that dumps messages to a C++
 *  std::ostream object.
 *
 * \note The intent of this class is to illustrate how to using the Logging
 *  facility within an MPI distributed environment and provide a utility that
 *  could be useful for debugging problems at small scales.
 *
 * \warning Do not use this for large-scale production runs.
 * \warning The intent of this class is to be used primarily with std::cout,
 *  std::cerr, etc. It is suggested that applications do not use this class
 *  with an std::ofstream object.
 */
class SynchronizedStream : public LogStream
{
public:
  SynchronizedStream(std::ostream* stream, MPI_Comm comm);
  SynchronizedStream(std::ostream* stream,
                     MPI_Comm comm,
                     const std::string& format);

  virtual ~SynchronizedStream();

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
   * \note This method doesn't put anything to the console. Instead the
   *  messages are cached locally to each ranks and are dumped to the console
   *  in rank order when flush is called.
   */
  virtual void append(message::Level msgLevel,
                      const std::string& message,
                      const std::string& tagName,
                      const std::string& fileName,
                      int line,
                      bool filter_duplicates);

  /*!
   * \brief Dumps the messages to the console in rank-order.
   */
  virtual void flush();

private:
  /// Forward declarations
  struct MessageCache;

  /// \name Private Members
  /// @{

  MPI_Comm m_comm;
  MessageCache* m_cache;
  std::ostream* m_stream;
  /// @}

  /*!
   * \brief Default constructor. Made private to prevent applications from
   *  using it. Instead the constructor that passes the underlying MPI comm
   *  should be used.
   */
  SynchronizedStream()
    : m_comm(MPI_COMM_NULL)
    , m_cache(static_cast<MessageCache*>(nullptr))
    , m_stream(static_cast<std::ostream*>(nullptr)) {};

  DISABLE_COPY_AND_ASSIGNMENT(SynchronizedStream);
  DISABLE_MOVE_AND_ASSIGNMENT(SynchronizedStream);
};

} /* namespace slic */
} /* namespace axom */

#endif /* SYNCHRONIZEDSTREAM_HPP_ */
