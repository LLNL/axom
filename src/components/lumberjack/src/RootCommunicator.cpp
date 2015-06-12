#include "lumberjack/RootCommunicator.hpp"

#include "common/CommonTypes.hpp"

#include <sstream>

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
		int incomingMpiCommRank = -1;
		std::stringstream ss;
		MPI_Status recieveStatus;
		for(int i=1; i<m_mpiCommSize; ++i){
			MPI_Recv(&incomingMpiCommRank, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &recieveStatus);
			//ToDo: Recieve entire message
			ss << "This message came from rank " << incomingMpiCommRank;
			MessageInfo mi(ss.str(), incomingMpiCommRank, "", -1);
			ss.str(std::string());
			m_messages.push_back(mi);
		}
	}
	else {
		for(int i=0; i<(int)m_messages.size(); ++i){
			//ToDO: Send entire message
			MPI_Send(&m_mpiCommRank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
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

void RootCommunicator::queueMessage(MessageInfo messageInfo)
{
	m_messages.push_back(messageInfo);
}

}
}