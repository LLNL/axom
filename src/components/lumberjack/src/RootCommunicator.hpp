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
		void setCommunicator(MPI_Comm comm);
		void pushMessagesOnce();
		void pushMessagesFully();
		std::vector<MessageInfo>* getMessages();
	private:
		void receiveMessages();
		
		MPI_Comm m_comm;
		int m_commRank;
		std::vector<int> m_ranks;
		std::vector<MessageInfo> m_messages;
};

}
}

#endif