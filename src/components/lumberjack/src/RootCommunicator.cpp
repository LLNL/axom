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

bool RootCommunicator::shouldMessagesBeOutputted()
{
    if (m_mpiCommRank == 0){
        return true;
    }
    return false;
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

void RootCommunicator::pushMessagesOnce(std::vector<Message*>& messages,
                                        std::vector<Combiner*>& combiners)
{
    MPI_Barrier(m_mpiComm);
    if (m_mpiCommRank == 0){
        Message* message;
        int ranksDoneCount = 0;
        while(ranksDoneCount < (m_mpiCommSize-1)){
            message = mpiBlockingRecieveAnyMessage(m_mpiComm, m_ranksLimit);
            if (message == ATK_NULLPTR) {
                ++ranksDoneCount;
            }
            else {
                messages.push_back(message);
            }
        }
        combineMessages(messages, combiners, m_ranksLimit);
    }
    else {
        mpiNonBlockingSendMessages(m_mpiComm, 0, messages);
    }
    MPI_Barrier(m_mpiComm);
}

void RootCommunicator::pushMessagesFully(std::vector<Message*>& messages,
                                         std::vector<Combiner*>& combiners)
{
    pushMessagesOnce(messages, combiners);
}

} // end namespace lumberjack
} // end namespace asctoolkit
