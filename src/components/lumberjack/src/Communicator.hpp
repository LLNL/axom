#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <vector>

namespace asctoolkit {
namespace lumberjack {

class Communicator {
	public:
		virtual void setCommunicator(MPI_Comm comm) = 0;
		virtual void pushMessagesOnce() = 0;
		virtual void pushMessagesFully() = 0;
		virtual std::vector<MessageInfos> getMessages() = 0;
};

}
}

#endif