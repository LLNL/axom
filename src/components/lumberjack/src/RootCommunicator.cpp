#include "lumberjack/RootCommunicator.hpp"

namespace asctoolkit {
namespace lumberjack {

void RootCommunicator::setCommunicator(MPI_Comm* comm)
{
	m_comm = comm;
	MPI_Comm_rank(&m_comm, &m_commRank)
}

void RootCommunicator::pushMessagesOnce()
{
	if (m_commRank == 0){
		receiveMessages();
	}
	else {
		
	}
}

void RootCommunicator::receiveMessages()
{

}

void RootCommunicator::pushMessagesFully()
{
	pushMessagesOnce();
}

std::vector<MessageInfos> RootCommunicator::getMessages()
{
	if (m_commRank == 0){
		std::vector<MessageInfos>* returnedVector = new std::vector<MessageInfos>;
		returnedVector.swap(m_messages);
		return returnedVector;
	}
	return nullptr;
}

}
}