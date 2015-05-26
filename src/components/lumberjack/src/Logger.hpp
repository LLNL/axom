#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

namespace asctoolkit {
namespace lumberjack {

class Logger {
	public:
		void initialize(/*MPI_COMM*/);
		void finalize();

		void setOutputStream();

		void addMessage(const std::string& message);
		void addMessage(const std::string& message, const std::string& fileName, int lineNumber);

		void pushMessagesOnceUpTree();
		void pushMessagesFullyUpTree();

		void outputMessages();
};

}
}

#endif