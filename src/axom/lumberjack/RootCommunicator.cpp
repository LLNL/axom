/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*!
 ******************************************************************************
 *
 * \file RootCommunicator.cpp
 *
 * \brief Implementation of the RootCommunicator class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/RootCommunicator.hpp"

#include <cstdlib>

#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace lumberjack
{

void RootCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  m_mpiComm = comm;
  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;
}

void RootCommunicator::finalize()
{}

int RootCommunicator::rank()
{
  return m_mpiCommRank;
}

void RootCommunicator::ranksLimit(int value)
{
  m_ranksLimit = value;
}

int RootCommunicator::ranksLimit()
{
  return m_ranksLimit;
}

int RootCommunicator::numPushesToFlush()
{
  return 1;
}

void RootCommunicator::push(const char* packedMessagesToBeSent,
                            std::vector<const char*>& receivedPackedMessages)
{
  MPI_Barrier(m_mpiComm);
  if (m_mpiCommRank == 0)
  {
    const char* currPackedMessages;
    int ranksDoneCount = 0;
    while(ranksDoneCount < (m_mpiCommSize-1))
    {
      currPackedMessages = mpiBlockingReceiveMessages(m_mpiComm);
      if (isPackedMessagesEmpty(currPackedMessages))
      {
        if ((currPackedMessages != nullptr) || (currPackedMessages[0] == '\0'))
        {
          delete [] currPackedMessages;
        }
      }
      else
      {
        receivedPackedMessages.push_back(currPackedMessages);
      }
      ++ranksDoneCount;
    }
  }
  else
  {
    if (isPackedMessagesEmpty(packedMessagesToBeSent))
    {
      mpiNonBlockingSendMessages(m_mpiComm, 0, zeroMessage);
    }
    else
    {
      mpiNonBlockingSendMessages(m_mpiComm, 0, packedMessagesToBeSent);
    }
  }
  MPI_Barrier(m_mpiComm);
}

bool RootCommunicator::isOutputNode()
{
  if (m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

} // end namespace lumberjack
} // end namespace axom
