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
 * \file BinaryTreeCommunicator.cpp
 *
 * \brief Implementation of the BinaryTreeCommunicator class.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/BinaryTreeCommunicator.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include <cstdlib>
#include <cmath>

#include "axom/lumberjack/MPIUtility.hpp"

namespace axom
{
namespace lumberjack
{

void BinaryTreeCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
  m_mpiComm = comm;
  MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
  MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
  m_ranksLimit = ranksLimit;

  // Calculate tree information about this rank
  m_treeHeight = int(axom::utilities::log2(m_mpiCommSize) + 1);
  m_parentRank = (m_mpiCommRank - 1) >> 1;
  m_leftChildRank = (m_mpiCommRank * 2) + 1;
  m_rightChildRank = (m_mpiCommRank * 2) + 2;
  m_childCount = 0;
  if (m_leftChildRank < m_mpiCommSize)
  {
    ++m_childCount;
  }
  else
  {
    m_leftChildRank = -1;
  }
  if (m_rightChildRank < m_mpiCommSize)
  {
    ++m_childCount;
  }
  else
  {
    m_rightChildRank = -1;
  }

}

void BinaryTreeCommunicator::finalize()
{}

int BinaryTreeCommunicator::rank()
{
  return m_mpiCommRank;
}

void BinaryTreeCommunicator::ranksLimit(int value)
{
  m_ranksLimit = value;
}

int BinaryTreeCommunicator::ranksLimit()
{
  return m_ranksLimit;
}

int BinaryTreeCommunicator::numPushesToFlush()
{
  return m_treeHeight-1;
}

void BinaryTreeCommunicator::push(const char* packedMessagesToBeSent,
                                  std::vector<const char*>& receivedPackedMessages)
{
  MPI_Barrier(m_mpiComm);
  if (m_mpiCommRank != 0)
  {
    mpiNonBlockingSendMessages(m_mpiComm, m_parentRank, packedMessagesToBeSent);
  }

  int childrenDoneCount = 0;
  const char* currPackedMessages;
  while(childrenDoneCount < m_childCount)
  {
    currPackedMessages = mpiBlockingRecieveMessages(m_mpiComm);
    if (currPackedMessages != nullptr)
    {
      receivedPackedMessages.push_back(currPackedMessages);
    }
    ++childrenDoneCount;
  }

  MPI_Barrier(m_mpiComm);
}

bool BinaryTreeCommunicator::isOutputNode()
{
  if (m_mpiCommRank == 0)
  {
    return true;
  }
  return false;
}

} // end namespace lumberjack
} // end namespace axom
