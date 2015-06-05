#include "lumberjack/RootCommunicator.hpp"

#include "common/CommonTypes.hpp"

namespace asctoolkit {
namespace lumberjack {

void RootCommunicator::setCommunicator(MPI_Comm comm)
{
	m_comm = comm;
	MPI_Comm_rank(m_comm, &m_commRank);
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

std::vector<MessageInfo>* RootCommunicator::getMessages()
{
	if (m_commRank == 0){
		std::vector<MessageInfo>* returnedVector = new std::vector<MessageInfo>;
		returnedVector->swap(m_messages);
		return returnedVector;
	}
	return ATK_NULLPTR;
}

}
}