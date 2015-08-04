#include "lumberjack/Logger.hpp"
#include "lumberjack/RootCommunicator.hpp"
#include "lumberjack/MessageInfo.hpp"

#include "common/CommonTypes.hpp"

#include <mpi.h>
#include <iostream>

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int commRank = -1;
    int commSize = -1;

    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);
    int ranksLimit = commSize/2;

    asctoolkit::lumberjack::RootCommunicator communicator;
    communicator.initialize(MPI_COMM_WORLD, ranksLimit);
    asctoolkit::lumberjack::Logger logger;
    logger.initialize(&communicator, ranksLimit);

    if (commRank == 0){
        logger.queueMessage("This message is not important and will not combined");
    }
    else {
        logger.queueMessage("This message is not important and will be combined");
    }
    logger.pushMessagesOnce();
    std::vector<asctoolkit::lumberjack::MessageInfo*>* messageInfos = logger.getMessages();

    if (commRank == 0){
        std::cout << "Rank 0: Printing Messages Recieved!!" << std::endl;
        for(int i=0; i<(int)(messageInfos->size()); ++i){
            asctoolkit::lumberjack::MessageInfo* currMessageInfo = messageInfos->at(i);
            std::cout << "(" << currMessageInfo->stringOfRanks() << ") " << currMessageInfo->rankCount() <<
                         " '" << currMessageInfo->message() << "'" << std::endl;
            delete currMessageInfo;
        }
    }
    else {
        if (messageInfos == ATK_NULLPTR){
            std::cout << "Rank " << commRank << ": Success! No messages left in queue!" << std::endl;
        }
        else {
            std::cout << "Rank " << commRank << ": Failure! Messages left in queue!" << std::endl;
        }
    }

    delete messageInfos;
    MPI_Barrier(MPI_COMM_WORLD);
    logger.finalize();
    communicator.finalize();
    MPI_Finalize();

    return 0;
}
