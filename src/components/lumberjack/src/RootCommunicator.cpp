#include "lumberjack/RootCommunicator.hpp"

#include <cstdlib>

namespace asctoolkit {
namespace lumberjack {

void RootCommunicator::initialize(MPI_Comm comm)
{
    m_mpiComm = comm;
    MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
    MPI_Comm_size(m_mpiComm, &m_mpiCommSize);
}

void RootCommunicator::finalize()
{

}

void RootCommunicator::pushMessagesOnce()
{
    if (m_mpiCommRank == 0){
        char* charArray = ATK_NULLPTR;
        int messageSize;
        MPI_Status mpiStatus;
        for(int i=1; i<m_mpiCommSize; ++i){
            messageSize = -1;
            MPI_Probe(i, 0, m_mpiComm, &mpiStatus);
            MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);
            charArray = (char*)malloc((messageSize+1)*sizeof(char));
            MPI_Recv(charArray, messageSize, MPI_CHAR, i, 0, m_mpiComm, &mpiStatus);
            charArray[messageSize] = '\0';
            std::string messageString(charArray);
            free(charArray);
            MessageInfo mi;
            mi.unpack(messageString, 5);
            m_messages.push_back(mi);
        }
    }
    else {
        for(int i=0; i<(int)m_messages.size(); ++i){
            std::string message = m_messages[i].pack();
            MPI_Send(const_cast<char*>(message.c_str()),
                     message.size(), MPI_CHAR, 0, 0, m_mpiComm);
        }
        m_messages.clear();
    }
}

void RootCommunicator::pushMessagesFully()
{
    pushMessagesOnce();
}

std::vector<MessageInfo>* RootCommunicator::getMessages()
{
    if (m_mpiCommRank == 0){
        std::vector<MessageInfo>* returnedVector = new std::vector<MessageInfo>;
        returnedVector->swap(m_messages);
        return returnedVector;
    }
    return ATK_NULLPTR;
}

void RootCommunicator::queueMessage(const std::string& message, const std::string& fileName, const int lineNumber)
{
    MessageInfo messageInfo(message, m_mpiCommRank, fileName, lineNumber);
    m_messages.push_back(messageInfo);
}

}
}