#ifndef LOGGER_HPP
#define LOGGER_HPP

#include <string>

#include "mpi.h"
#include "lumberjack/Combiner.hpp"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageEqualityCombiner.hpp"
#include "lumberjack/MessageInfo.hpp"

namespace asctoolkit {
namespace lumberjack {

class Logger {
    public:
        void initialize(Communicator* communicator, int ranksLimit);
        void finalize();

        void setOutputStream();
        void setErrorStream();

        void flushOutputStream();
        void flushErrorStream();

        void addMessageCombiner(Combiner* combiner);
        void removeMessageCombiner(const std::string& combinerIdentifier);
        void clearMessageCombiners();

        std::vector<MessageInfo*>* getMessages();

        void ranksLimit(int ranksLimit);
        int ranksLimit();

        void queueMessage(const std::string& message);
        void queueMessage(const std::string& message, const std::string& fileName, const int lineNumber);

        void pushMessagesOnce();
        void pushMessagesFully();
    private:
        void combineMessages();

        Communicator* m_communicator;
        int m_ranksLimit;
        MessageEqualityCombiner* m_messageEqualityCombiner;
        std::vector<Combiner*> m_combiners;
        std::vector<MessageInfo*> m_messages;
};

}
}

#endif