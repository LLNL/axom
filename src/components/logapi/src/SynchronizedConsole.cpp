/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file SynchronizedConsole.cpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "SynchronizedConsole.hpp"

// C/C++ includes
#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

namespace asctoolkit {
namespace logapi {

struct SynchronizedConsole::MessageCache
{

  std::vector< std::string > messages;

  void printMessages( const int rank )
  {
    const unsigned N = messages.size();

    if( N==0 ) {
      /* short-circuit */
      return;
    }

    std::cout << "======\n";
    std::cout << "RANK[" << rank << "] NMSGS=" << N << std::endl;
    std::cout << "======\n";

    for ( unsigned i=0; i < N; ++i ) {
      std::cout << messages[ i ] << std::endl;
    } // END for all messages

    std::cout.flush();

    messages.clear();
  }

};

//------------------------------------------------------------------------------
SynchronizedConsole::SynchronizedConsole():
    m_comm(MPI_COMM_NULL),
    m_cache( new MessageCache() )
{

}

//------------------------------------------------------------------------------
SynchronizedConsole::~SynchronizedConsole()
{

  if ( m_cache != NULL ) {
    delete m_cache;
    m_cache = NULL;
  }

}

//------------------------------------------------------------------------------
void SynchronizedConsole::setCommunicator(MPI_Comm comm)
{
  assert( "pre: NULL MPI Communicator" && (comm != MPI_COMM_NULL)  );
  m_comm = comm;
}

//------------------------------------------------------------------------------
void SynchronizedConsole::append( int msgType,
                                  const std::string& msgTypeName,
                                  const std::string& message,
                                  const std::string& fileName,
                                  int line )
{
  assert( "pre: null message cache!" && (m_cache != NULL) );

  // STEP 0: form the message
  std::ostringstream oss;
  oss << "[" << msgTypeName << "]: " << message << std::endl;
  oss << "FILE: " << fileName << std::endl;
  oss << "LINE: " << line << std::endl;

  // STEP 1: cache it
  m_cache->messages.push_back( oss.str() );
}

//------------------------------------------------------------------------------
void SynchronizedConsole::flush()
{
  assert( "pre: null message cache!" && (m_cache != NULL) );
  assert( "pre: null MPI communicator!" && (m_comm != MPI_COMM_NULL) );

  int rank   = -1;
  int nranks =  0;
  MPI_Comm_rank( m_comm, &rank   );
  MPI_Comm_size( m_comm, &nranks );

  const int prevrank = rank-1;
  const int nextrank = rank+1;

  MPI_Request null_request = MPI_REQUEST_NULL;

  if ( rank == 0 ) {

    /* rank 0 */

    // print messages at this rank
    m_cache->printMessages( rank );

    if ( nranks > 1 ) {

      /* signal next rank */
      MPI_Isend(NULL,0,MPI_INT,1,0,m_comm,&null_request);

    } // END if

  } else if( rank == nranks-1 ) {

    /* last rank */

    // Wait for signal from previous rank
    MPI_Recv(NULL,0,MPI_INT,prevrank,MPI_ANY_TAG,m_comm,MPI_STATUSES_IGNORE);

    // print messages at this rank
    m_cache->printMessages( rank );

  } else {

    // Wait for signal from previous rank
    MPI_Recv(NULL,0,MPI_INT,prevrank,MPI_ANY_TAG,m_comm,MPI_STATUSES_IGNORE);

    // print messages at this rank
    m_cache->printMessages( rank );

    // signal next rank
    MPI_Isend(NULL,0,MPI_INT,nextrank,0,m_comm,&null_request);

  } // END else

}

} /* namespace logapi */
} /* namespace asctoolkit */
