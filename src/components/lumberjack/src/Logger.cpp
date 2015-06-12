#include "lumberjack/Logger.hpp"

#include "common/CommonTypes.hpp"

namespace asctoolkit {
namespace lumberjack {

void Logger::initialize(MPI_Comm comm, Communicator* communicator)
{
	m_mpiComm = comm;
	MPI_Comm_rank(m_mpiComm, &m_mpiCommRank);
	m_communicator = communicator;
	m_communicator->initialize(m_mpiComm);
}

void Logger::finalize()
{
	m_communicator->finalize();
	m_communicator = ATK_NULLPTR;
}

void Logger::setOutputStream()
{

}

void Logger::setErrorStream()
{

}

void Logger::flushOutputStream()
{

}

void Logger::flushErrorStream()
{

}

std::vector<MessageInfo>* Logger::getMessages()
{
	return m_communicator->getMessages();
}

void Logger::queueMessage(const std::string& message)
{
	queueMessage(message, "", -1);
}

void Logger::queueMessage(const std::string& message, const std::string& fileName, const int lineNumber)
{
	MessageInfo messageInfo(message, m_mpiCommRank, fileName, lineNumber);
	m_communicator->queueMessage(messageInfo);
}

void Logger::pushMessagesOnce()
{
	m_communicator->pushMessagesOnce();
}

void Logger::pushMessagesFully()
{
	m_communicator->pushMessagesFully();
}

}
}
