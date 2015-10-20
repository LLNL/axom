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
 * \file speedTest.cpp
 * \author Chris White (white238@llnl.gov)
 *******************************************************************************
 */

#include "lumberjack/Lumberjack.hpp"
#include "lumberjack/BinaryTreeCommunicator.hpp"
#include "lumberjack/RootCommunicator.hpp"
#include "lumberjack/Message.hpp"
#include "lumberjack/Utility.hpp"

#include <mpi.h>

#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

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
    asctoolkit::lumberjack::Communicator* communicator;
    if (std::string(argv[1]) == "b") {
        communicator = new asctoolkit::lumberjack::BinaryTreeCommunicator;
    } else if (std::string(argv[1]) == "r") {
        communicator = new asctoolkit::lumberjack::RootCommunicator;
    }
    communicator->initialize(MPI_COMM_WORLD, ranksLimit);

    // Initialize lumberjack logger
    asctoolkit::lumberjack::Lumberjack lj;
    lj.initialize(communicator, ranksLimit);

    // Read lines from file
    std::string currMessage;
    std::vector<std::string> lines;
    std::ifstream file(argv[3]);
    while(std::getline(file, currMessage)){
        currMessage += '\n';
        lines.push_back(currMessage);
    }
    file.close();

    // Start clock
    std::clock_t begin = clock();

    // Queue messages into lumberjack
    int cycleCount = 0;
    int cycleLimit = asctoolkit::lumberjack::stringToInt(argv[2]);
    int linesSize = (int)lines.size();
    for (int i = 0; i < linesSize; ++i){
        lj.queueMessage(lines[i]);
        ++cycleCount;
        if (cycleCount > cycleLimit) {
            lj.pushMessagesOnce();
            cycleCount = 0;
        }
    }

    // Push messages fully through lumberjack's communicator
    lj.pushMessagesFully();

    // End clock
    std::clock_t end = clock();

    // Get messages back out of lumberjack since they have been pushed.
    std::vector<asctoolkit::lumberjack::Message*> messages;
    lj.getMessages(messages);

    if (commRank == 0) {
       std::ofstream outFile;
       outFile.open("speedTestOutput");
       for(int i=0; i<(int)(messages.size()); ++i){
           outFile << messages[i]->text();
           delete messages[i];
       }
       outFile.close();
   }

    // Finalize the lumberjack 
    lj.finalize();
    // Finalize the lumberjack communicator
    communicator->finalize();
    delete communicator;

    // Output elapsed time
    if (commRank == 0) {
        std::cout << "Elapsed time: " << ((double)(end - begin)*1000)/CLOCKS_PER_SEC << std::endl;
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
