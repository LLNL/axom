// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file NonCollectiveRootCommunicator.cpp
 *
 * \brief Implementation of the NonCollectiveRootCommunicator class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/NonCollectiveRootCommunicator.hpp"
#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace lumberjack
{
void NonCollectiveRootCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  static int mpiTag = 32767;
  m_mpiComm = comm;
  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;
  m_mpiTag = mpiTag;
  ++mpiTag;
}

void NonCollectiveRootCommunicator::finalize() { }

int NonCollectiveRootCommunicator::rank() { return m_mpiCommRank; }

void NonCollectiveRootCommunicator::ranksLimit(int value)
{
  m_ranksLimit = value;
}

int NonCollectiveRootCommunicator::ranksLimit() { return m_ranksLimit; }

int NonCollectiveRootCommunicator::numPushesToFlush() { return 1; }

void NonCollectiveRootCommunicator::push(
  const char* packedMessagesToBeSent,
  std::vector<const char*>& receivedPackedMessages)
{
  if(m_mpiCommRank == 0)
  {
    const char* currPackedMessages = nullptr;
    bool receive_messages = true;
    while(receive_messages)
    {
      currPackedMessages = mpiNonBlockingReceiveMessages(m_mpiComm, m_mpiTag);

      if(isPackedMessagesEmpty(currPackedMessages))
      {
        if(currPackedMessages == nullptr)
        {
          receive_messages = false;
        }
        else
        {
          delete[] currPackedMessages;
        }
      }
      else
      {
        receivedPackedMessages.push_back(currPackedMessages);
      }

      currPackedMessages = nullptr;
    }
  }
  else
  {
    if(isPackedMessagesEmpty(packedMessagesToBeSent) == false)
    {
      mpiNonBlockingSendMessages(m_mpiComm, 0, packedMessagesToBeSent, m_mpiTag);
    }
  }
}

bool NonCollectiveRootCommunicator::isOutputNode()
{
  if(m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

int NonCollectiveRootCommunicator::mpiTag() const { return m_mpiTag; }

}  // end namespace lumberjack
}  // end namespace axom
