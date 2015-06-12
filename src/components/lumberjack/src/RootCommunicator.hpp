#ifndef ROOTCOMMUNICATOR_HPP
#define ROOTCOMMUNICATOR_HPP

#include <string>

#include "mpi.h"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class RootCommunicator: public Communicator {
	public:
		void initialize(MPI_Comm comm);
		void finalize();
		void pushMessagesOnce();
		void pushMessagesFully();
		std::vector<MessageInfo>* getMessages();
		void queueMessage(MessageInfo messageInfo);
	private:
		MPI_Comm m_mpiComm;
		int m_mpiCommRank;
		int m_mpiCommSize;
		std::vector<MessageInfo> m_messages;
};

}
}

#endif