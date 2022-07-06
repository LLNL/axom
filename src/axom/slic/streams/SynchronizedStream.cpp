// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/SynchronizedStream.hpp"

#include <vector>

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slic
{
struct SynchronizedStream::MessageCache
{
  std::vector<std::string> messages;

  void printMessages(std::ostream* stream)
  {
    if(stream == nullptr)
    {
      std::cerr << "ERROR: cannot write to NULL stream!\n";
      return;
    }

    const unsigned N = messages.size();

    if(N == 0)
    {
      /* short-circuit */
      return;
    }

    for(unsigned i = 0; i < N; ++i)
    {
      (*stream) << messages[i];
    }  // END for all messages

    messages.clear();
  }
};

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(std::ostream* stream, MPI_Comm comm)
  : m_comm(comm)
  , m_cache(new MessageCache())
  , m_stream(stream)
{ }

//------------------------------------------------------------------------------
SynchronizedStream::SynchronizedStream(std::ostream* stream,
                                       MPI_Comm comm,
                                       const std::string& format)
  : m_comm(comm)
  , m_cache(new MessageCache)
  , m_stream(stream)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
SynchronizedStream::~SynchronizedStream()
{
  delete m_cache;
  m_cache = static_cast<MessageCache*>(nullptr);
}

//------------------------------------------------------------------------------
void SynchronizedStream::append(message::Level msgLevel,
                                const std::string& message,
                                const std::string& tagName,
                                const std::string& fileName,
                                int line,
                                bool AXOM_UNUSED_PARAM(filter_duplicates))
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return;
  }

  int rank = -1;
  MPI_Comm_rank(m_comm, &rank);

  // STEP 1: cache formatted message
  m_cache->messages.push_back(
    this->getFormatedMessage(message::getLevelAsString(msgLevel),
                             message,
                             tagName,
                             std::to_string(rank),
                             fileName,
                             line));
}

//------------------------------------------------------------------------------
void SynchronizedStream::flush()
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return;
  }

  if(m_comm == MPI_COMM_NULL)
  {
    std::cerr << "ERROR: NULL communicator!\n";
    return;
  }

  int rank = -1;
  int nranks = 0;
  MPI_Comm_rank(m_comm, &rank);
  MPI_Comm_size(m_comm, &nranks);

  const int prevrank = rank - 1;
  const int nextrank = rank + 1;

  if(rank > 0)
  {
    // wait for signal from previous rank
    MPI_Recv(nullptr, 0, MPI_INT, prevrank, MPI_ANY_TAG, m_comm, MPI_STATUSES_IGNORE);
  }

  // print messages for this rank
  m_cache->printMessages(m_stream);

  // signal next rank
  if(nranks > 1 && nextrank < nranks)
  {
    MPI_Request null_request = MPI_REQUEST_NULL;
    MPI_Isend(nullptr, 0, MPI_INT, nextrank, 0, m_comm, &null_request);
    MPI_Request_free(&null_request);
  }
}

//------------------------------------------------------------------------------
// void SynchronizedStream::determineAbortState()
bool SynchronizedStream::confirmAbort()
{
  if(m_cache == nullptr)
  {
    std::cerr << "ERROR: NULL cache!\n";
    return false;
  }

  if(m_comm == MPI_COMM_NULL)
  {
    std::cerr << "ERROR: NULL communicator!\n";
    return false;
  }

  int rank = -1;
  int nranks = 0;
  MPI_Comm_rank(m_comm, &rank);
  MPI_Comm_size(m_comm, &nranks);

  bool rank_abort = getAbortFlag();
  bool all_abort = false;
  MPI_Allreduce(&rank_abort, &all_abort, 1, MPI_C_BOOL, MPI_LOR, m_comm);

  return all_abort;
}

} /* namespace slic */
} /* namespace axom */
