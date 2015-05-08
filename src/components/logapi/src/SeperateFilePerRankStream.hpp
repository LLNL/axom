/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file SeperateFilePerRankStream.hpp
 *
 * \date May 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef SEPERATEFILEPERRANKSTREAM_HPP_
#define SEPERATEFILEPERRANKSTREAM_HPP_

#include "logapi/LogStream.hpp"

// MPI
#include <mpi.h>

namespace asctoolkit {

namespace logapi {

class SeperateFilePerRankStream : public LogStream
{
public:
  SeperateFilePerRankStream( const std::string& prefix,
                           MPI_Comm c );

  virtual ~SeperateFilePerRankStream();

  /// \brief See LogStream::append
  virtual void append( MessageType msgType,
                       const std::string& msgTypeName,
                       const std::string& message,
                       const std::string& fileName,
                       int line );

private:

  /*!
   *****************************************************************************
   * \brief Default constructor. Made private to prevent users from using it.
   *****************************************************************************
   */
  SeperateFilePerRankStream() { m_comm=MPI_COMM_NULL; m_fstream=0; };

  // Forward Declarations
  struct FileStream;

  /// \name Private Members
  /// @{

  MPI_Comm m_comm;
  FileStream* m_fstream;
  /// @}

  /// \name Disabled Methods
  /// @{

  SeperateFilePerRankStream( const SeperateFilePerRankStream& ); // Not implemented
  SeperateFilePerRankStream& operator=( const SeperateFilePerRankStream& ); // Not implemented

  /// @}

};

} /* namespace logapi */

} /* namespace asctoolkit */

#endif /* SEPERATEFILEPERRANKSTREAM_HPP_ */
