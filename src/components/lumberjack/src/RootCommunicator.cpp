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

#include <cstdlib>
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

void RootCommunicator::pushMessageInfosOnce(std::vector<MessageInfo*>& messageInfos)
{
    MPI_Request mpiRequests[m_mpiCommSize];
    MPI_Status mpiStatuses[m_mpiCommSize];

    MPI_Barrier(m_mpiComm);
    if (m_mpiCommRank == 0){
        char* charArray = ATK_NULLPTR;
        int messageSize;
        MPI_Status mpiStatus;
        for(int i=1; i<m_mpiCommSize; ++i){
            messageSize = -1;
            MPI_Probe(i, 0, m_mpiComm, &mpiStatus);
            MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);
            charArray = (char*)malloc((messageSize+1)*sizeof(char));
            MPI_Irecv(charArray, messageSize, MPI_CHAR, i, 0, m_mpiComm, &mpiRequests[i]);
            charArray[messageSize] = '\0';
            std::string messageString(charArray);
            free(charArray);
            MessageInfo* mi = new MessageInfo();
            mi->unpack(messageString, m_ranksLimit);
            messageInfos.push_back(mi);
        }
    }
    else {
        for(int i=0; i<(int)messageInfos.size(); ++i){
            std::string packedMessageInfo = messageInfos[i]->pack();
            MPI_Isend(const_cast<char*>(packedMessageInfo.c_str()),
                     packedMessageInfo.size(), MPI_CHAR, 0, 0, m_mpiComm, &mpiRequests[i]);
            delete messageInfos[i];
        }
        messageInfos.clear();
    }
    MPI_Waitall(m_mpiCommSize, mpiRequests, mpiStatuses);
}

void RootCommunicator::pushMessageInfosFully(std::vector<MessageInfo*>& messageInfos)
{
    pushMessageInfosOnce(messageInfos);
}

} // end namespace lumberjack
} // end namespace asctoolkit
