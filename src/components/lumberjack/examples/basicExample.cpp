/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file basicExample.cpp
 * \author Chris White (white238@llnl.gov)
 *******************************************************************************
 */

#include "lumberjack/Logger.hpp"
#include "lumberjack/BinaryTreeCommunicator.hpp"
#include "lumberjack/Message.hpp"

#include <mpi.h>
#include <iostream>

//------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    // Initialize MPI and get rank and comm size
    MPI_Init(&argc, &argv);

    int commRank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);

    // Determine how many ranks we want to individually track per message
    int ranksLimit = commSize/2;

    // Initialize which lumberjack communicator we want
    asctoolkit::lumberjack::BinaryTreeCommunicator communicator;
    communicator.initialize(MPI_COMM_WORLD, ranksLimit);

    // Initialize lumberjack logger
    asctoolkit::lumberjack::Logger logger;
    logger.initialize(&communicator, ranksLimit);

    // Queue messages into lumberjack
    if (commRank == 0){
        logger.queueMessage("This message will not combined");
    }
    else {
        logger.queueMessage("This message will be combined");
        logger.queueMessage("This message will be combined");
        logger.queueMessage("This message will be combined");
    }
    // Push messages fully through lumberjack's communicator
    logger.pushMessagesFully();

    // Get messages back out of lumberjack since they have been pushed.
    std::vector<asctoolkit::lumberjack::Message*> messages;
    logger.getMessages(messages);
    for(int i=0; i<(int)(messages.size()); ++i){
        std::cout << "(" << messages[i]->stringOfRanks() << ") " << messages[i]->ranksCount() <<
                     " '" << messages[i]->text() << "'" << std::endl;
        delete messages[i];
    }

    // Finalize the lumberjack logger
    logger.finalize();
    // Finalize the lumberjack communicator
    communicator.finalize();
    // Finalize MPI
    MPI_Finalize();

    return 0;
}
