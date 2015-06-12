#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

#include "mpi.h"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Logger {
	public:
		void initialize(MPI_Comm comm, Communicator* communicator);
		void finalize();

		void setOutputStream();
		void setErrorStream();

		void flushOutputStream();
		void flushErrorStream();

		std::vector<MessageInfo>* getMessages();

		void queueMessage(const std::string& message);
		void queueMessage(const std::string& message, const std::string& fileName, const int lineNumber);

		void pushMessagesOnce();
		void pushMessagesFully();
	private:
		MPI_Comm m_mpiComm;
		int m_mpiCommRank;
		Communicator* m_communicator;
};

}
}

#endif