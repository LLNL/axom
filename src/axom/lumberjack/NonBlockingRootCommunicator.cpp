// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file NonBlockingRootCommunicator.cpp
 *
 * \brief Implementation of the NonBlockingRootCommunicator class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/NonBlockingRootCommunicator.hpp"
#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace lumberjack
{
void NonBlockingRootCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  m_mpiComm = comm;
  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;
}

void NonBlockingRootCommunicator::finalize() { }

int NonBlockingRootCommunicator::rank() { return m_mpiCommRank; }

void NonBlockingRootCommunicator::ranksLimit(int value) { m_ranksLimit = value; }

int NonBlockingRootCommunicator::ranksLimit() { return m_ranksLimit; }

int NonBlockingRootCommunicator::numPushesToFlush() { return 1; }

void NonBlockingRootCommunicator::push(const char* packedMessagesToBeSent,
                            std::vector<const char*>& receivedPackedMessages)
{
  constexpr int mpiTag = 32767; 
  if(m_mpiCommRank == 0)
  {
    const char* currPackedMessages = nullptr;
    bool receive_messages = true;
    while(receive_messages)
    {
      currPackedMessages = mpiNonBlockingReceiveMessages(m_mpiComm, mpiTag);

      if(isPackedMessagesEmpty(currPackedMessages))
      {
        if (currPackedMessages == nullptr )
        {
          receive_messages = false;
        } else {
          delete [] currPackedMessages;
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
      mpiNonBlockingSendMessages(m_mpiComm, 0, packedMessagesToBeSent, mpiTag);
    }
  }
}

bool NonBlockingRootCommunicator::isOutputNode()
{
  if(m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

}  // end namespace lumberjack
}  // end namespace axom