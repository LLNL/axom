#include "lumberjack/BinaryTreeCommunicator.hpp"

namespace asctoolkit {
namespace lumberjack {

void BinaryTreeCommunicator::setCommunicator(MPI_Comm comm)
{

}

void BinaryTreeCommunicator::pushMessagesOnce()
{

}

void BinaryTreeCommunicator::pushMessagesFully()
{

}

std::vector<MessageInfos> BinaryTreeCommunicator::getMessages()
{

}

void BinaryTreeCommunicator::queueMessage(const std::string& message,
                                    const std::string& fileName,
                                    const int lineNumber)
{
    MessageInfo mi(message, m_commRank, fileName, lineNumber);
    m_messages.append(mi);
}

void BinaryTreeCommunicator::queueMessage(const std::string& message)
{
    queueMessage(message, "", -1);
}

}
}
