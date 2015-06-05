#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <vector>

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Communicator {
	public:
		virtual void setCommunicator(MPI_Comm comm) = 0;
		virtual void pushMessagesOnce() = 0;
		virtual void pushMessagesFully() = 0;
		virtual std::vector<MessageInfo>* getMessages() = 0;
};

}
}

#endif