/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

#include "lumberjack/Logger.hpp"
#include "lumberjack/RootCommunicator.hpp"
#include "lumberjack/MessageInfo.hpp"

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
        logger.queueMessage("This message will not combined");
    }
    else {
        logger.queueMessage("This message will be combined");
    }
    logger.pushMessagesOnce();

    std::vector<asctoolkit::lumberjack::MessageInfo*> messageInfos;
    logger.getMessageInfos(messageInfos);
    for(int i=0; i<(int)(messageInfos.size()); ++i){
        std::cout << "(" << messageInfos[i]->stringOfRanks() << ") " << messageInfos[i]->rankCount() <<
                     " '" << messageInfos[i]->message() << "'" << std::endl;
        delete messageInfos[i];
    }

    logger.finalize();
    communicator.finalize();
    MPI_Finalize();

    return 0;
}
