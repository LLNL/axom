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

void RootCommunicator::pushMessagesOnce(std::vector<Message*>& messages)
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
    }
    else {
        char outOfMessagesChar = '0';
        MPI_Request mpiRequest;
        for(int i=0; i<(int)messages.size(); ++i){
            std::string packedMessage = messages[i]->pack();
            MPI_Isend(const_cast<char*>(packedMessage.c_str()),
                     packedMessage.size(), MPI_CHAR, 0, 0, m_mpiComm, &mpiRequest);
            delete messages[i];
        }
        // Send that we are done sending messages from this rank
        MPI_Isend(&outOfMessagesChar, 1, MPI_CHAR, 0, 0, m_mpiComm, &mpiRequest);
        messages.clear();
    }
    MPI_Barrier(m_mpiComm);
}

void RootCommunicator::pushMessagesFully(std::vector<Message*>& messages)
{
    pushMessagesOnce(messages);
}

} // end namespace lumberjack
} // end namespace asctoolkit
