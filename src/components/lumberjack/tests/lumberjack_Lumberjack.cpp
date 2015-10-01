/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "gtest/gtest.h"

#include "lumberjack/Lumberjack.hpp"

#include "common/CommonTypes.hpp"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/Message.hpp"
#include "lumberjack/Utility.hpp"

#include <stdlib.h>
#include <time.h>

class TestCommunicator: public asctoolkit::lumberjack::Communicator {
    public:
        void initialize(MPI_Comm comm, int ranksLimit)
        {
            m_mpiComm = comm;
            m_ranksLimit = ranksLimit;
            srand(time(ATK_NULLPTR));
        }

        void finalize()
        {

        }

        void ranksLimit(int value)
        {
            m_ranksLimit = value;
        }

        int ranksLimit()
        {
            return m_ranksLimit;
        }

        void pushMessagesOnce(std::vector<asctoolkit::lumberjack::Message*>& messages,
                              std::vector<asctoolkit::lumberjack::Combiner*>& combiners)
        {
            asctoolkit::lumberjack::combineMessages(messages, combiners, m_ranksLimit);
        }

        void pushMessagesFully(std::vector<asctoolkit::lumberjack::Message*>& messages,
                               std::vector<asctoolkit::lumberjack::Combiner*>& combiners)
        {
            pushMessagesOnce(messages, combiners);
        }

        bool shouldMessagesBeOutputted()
        {
            return true;
        }

        int rank()
        {
            return rand() % (m_ranksLimit*4);
        }

    private:
        MPI_Comm m_mpiComm;
        int m_ranksLimit;
};


TEST(lumberjack_Lumberjack, combineMessages01)
{
    int ranksLimit = 5;
    TestCommunicator communicator;
    communicator.initialize(MPI_COMM_NULL, ranksLimit);
    asctoolkit::lumberjack::Lumberjack lumberjack;
    lumberjack.initialize(&communicator, ranksLimit);

    lumberjack.queueMessage("Should be combined.");
    lumberjack.queueMessage("Should be combined.");
    lumberjack.queueMessage("Should be combined.");
    lumberjack.queueMessage("Should be combined.");
    lumberjack.queueMessage("Should be combined.");
    lumberjack.queueMessage("Should be combined.");

    lumberjack.pushMessagesFully();

    std::vector<asctoolkit::lumberjack::Message*> messages;
    lumberjack.getMessages(messages);

    EXPECT_EQ((int)messages.size(), 1);
    EXPECT_EQ(messages[0]->text(), "Should be combined.");
    EXPECT_EQ(messages[0]->ranksCount(), 6);

    lumberjack.finalize();
    communicator.finalize();
}
