#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <vector>

#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Communicator {
    public:
        virtual void initialize(MPI_Comm comm) = 0;
        virtual void finalize() = 0;
        virtual void pushMessagesOnce() = 0;
        virtual void pushMessagesFully() = 0;
        virtual std::vector<MessageInfo>* getMessages() = 0;
        virtual void queueMessage(const std::string& message, const std::string& fileName, const int lineNumber) = 0;
};

}
}

#endif