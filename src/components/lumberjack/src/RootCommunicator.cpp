/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file RootCommunicator.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the class implementation of the RootCommunicator.
 *******************************************************************************
 */

#include "lumberjack/RootCommunicator.hpp"

#include "common/CommonTypes.hpp"

#include <cstdlib>

#include "lumberjack/MPIUtility.hpp"
#include "lumberjack/Utility.hpp"

namespace asctoolkit {
namespace lumberjack {

void RootCommunicator::initialize(MPI_Comm comm, int ranksLimit)
{
    m_mpiComm = comm;
    MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
    MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
    m_ranksLimit = ranksLimit;
}

void RootCommunicator::finalize()
{

}

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
    if (m_mpiCommRank == 0){
        const char* currPackedMessages;
        int ranksDoneCount = 0;
        while(ranksDoneCount < (m_mpiCommSize-1)){
            currPackedMessages = mpiBlockingRecieveMessages(m_mpiComm);
            if (currPackedMessages != ATK_NULLPTR) {
                receivedPackedMessages.push_back(currPackedMessages);
            }
            ++ranksDoneCount;
        }
    }
    else {
        mpiNonBlockingSendMessages(m_mpiComm, 0, packedMessagesToBeSent);
    }
    MPI_Barrier(m_mpiComm);
}

bool RootCommunicator::shouldMessagesBeOutputted()
{
    if (m_mpiCommRank == 0){
        return true;
    }
    return false;
}

} // end namespace lumberjack
} // end namespace asctoolkit
