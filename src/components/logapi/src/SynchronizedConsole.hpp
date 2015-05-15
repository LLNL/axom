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
 * \file SynchronizedConsole.hpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef SYNCHRONIZEDCONSOLE_HPP_
#define SYNCHRONIZEDCONSOLE_HPP_

#include "logapi/LogStream.hpp"

// MPI
#include <mpi.h> // For MPI


namespace asctoolkit {
namespace logapi {

/*!
 *******************************************************************************
 * \class SynchronizedConsole
 *
 * \brief A concrete instance of LogStream that dumps messages to the console
 *  in rank order.
 *
 * \note The intent of this class is to illustrate how to using the Logging
 *  facility within an MPI distributed environment and provide a utility that
 *  could be useful for debugging problems at small scales.
 *
 * \warning Do not use this for large-scale production runs.
 *******************************************************************************
 */
class SynchronizedConsole : public LogStream
{
public:
  SynchronizedConsole();
  virtual ~SynchronizedConsole();

  /*!
   *****************************************************************************
   * \brief Sets the MPI communicator.
   * \param [in] comm the MPI communicator to use by this instance.
   * \pre comm != MPI_COMM_NULL
   *****************************************************************************
   */
  void setCommunicator(MPI_Comm comm);

  /*!
   *****************************************************************************
   * \brief Appends the given message to the stream.
   *
   * \param [in] msgType the type of the message.
   * \param [in] msgTypeName string representation of the message type.
   * \param [in] message the user-supplied message.
   * \param [in] fileName the file where this message is appended
   * \param [in] line the line within the file at which the message is appended.
   *
   * \note This method doesn't put anything to the console. Instead the
   *  messages are cached locally to each ranks and are dumped to the console
   *  in rank order when flush is called.
   *****************************************************************************
   */
  virtual void append( MessageType msgType,
                       const std::string& msgTypeName,
                       const std::string& message,
                       const std::string& fileName,
                       int line );

  /*!
   *****************************************************************************
   * \brief Dumps the messages to the console in rank-order.
   *****************************************************************************
   */
  virtual void flush();

private:

  /// Forward declarations
  struct MessageCache;

  /// \name Private Members
  /// @{

  MPI_Comm m_comm;
  MessageCache* m_cache;

  /// @}

  /// \name Disabled Methods
  /// @{

  SynchronizedConsole( const SynchronizedConsole& ); // Not implemented
  SynchronizedConsole& operator=( const SynchronizedConsole& ); // Not implemented

  /// @}
};

} /* namespace logapi */
} /* namespace asctoolkit */

#endif /* SYNCHRONIZEDCONSOLE_HPP_ */
