#include "gtest/gtest.h"

#include "lumberjack/Logger.hpp"

#include "common/CommonTypes.hpp"
#include "lumberjack/Communicator.hpp"
#include "lumberjack/MessageInfo.hpp"

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

        void pushMessagesOnce(std::vector<asctoolkit::lumberjack::MessageInfo*>& messages)
        {
        	messages = messages;
        }

        void pushMessagesFully(std::vector<asctoolkit::lumberjack::MessageInfo*>& messages)
        {
        	messages = messages;
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


TEST(lumberjack_Logger, combineMessages01)
{
	int ranksLimit = 5;
	TestCommunicator communicator;
	communicator.initialize(ATK_NULLPTR, ranksLimit);
	asctoolkit::lumberjack::Logger logger;
	logger.initialize(&communicator, ranksLimit);

	logger.queueMessage("Should be combined.");
	logger.queueMessage("Should be combined.");
	logger.queueMessage("Should be combined.");
	logger.queueMessage("Should be combined.");
	logger.queueMessage("Should be combined.");
	logger.queueMessage("Should be combined.");

	logger.pushMessagesFully();

	std::vector<asctoolkit::lumberjack::MessageInfo*>* messages = logger.getMessages();
	asctoolkit::lumberjack::MessageInfo* mi = messages->at(0);
	EXPECT_EQ((int)messages->size(), 1);
	EXPECT_EQ(mi->message(), "Should be combined.");
	EXPECT_EQ(mi->rankCount(), 6);

	logger.finalize();
	communicator.finalize();
}
