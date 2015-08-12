/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

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

        void addMessageCombiner(Combiner* combiner);
        void removeMessageCombiner(const std::string& combinerIdentifier);
        void clearMessageCombiners();

        void getMessageInfos(std::vector<MessageInfo*>& filledVector);

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