#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

#include "mpi.h"
#include "lumberjack/Communicator.hpp"

namespace asctoolkit {
namespace lumberjack {

class RootCommunicator: public Communicator {
	public:
		void setCommunicator(MPI_Comm* comm);
		void pushMessagesOnceUpTree();
		void pushMessagesFullyUpTree();
		std::vector<MessageInfos> getMessages();
	private:
		void receiveMessages();
		
		MPI_Comm* m_comm;
		int m_commRank
		std::vector<int> m_ranks;
		std::vector<MessageInfo> m_messages;
};

}
}

#endif