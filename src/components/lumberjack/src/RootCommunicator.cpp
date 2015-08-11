#include "lumberjack/RootCommunicator.hpp"

#include <cstdlib>

#include <iostream>
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

void RootCommunicator::pushMessagesOnce(std::vector<MessageInfo*>& messages)
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
            messages.push_back(mi);
        }
    }
    else {
        for(int i=0; i<(int)messages.size(); ++i){
            std::string message = messages[i]->pack();
            MPI_Isend(const_cast<char*>(message.c_str()),
                     message.size(), MPI_CHAR, 0, 0, m_mpiComm, &mpiRequests[i]);
            delete messages[i];
        }
        messages.clear();
    }
    MPI_Waitall(m_mpiCommSize, mpiRequests, mpiStatuses);
}

void RootCommunicator::pushMessagesFully(std::vector<MessageInfo*>& messages)
{
    pushMessagesOnce(messages);
}

}
}