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
 * \file MPIUtility.cpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the implementation of utility functions.
 *******************************************************************************
 */

#include "lumberjack/MPIUtility.hpp"

#include "common/CommonTypes.hpp"

#include <cstdlib>

namespace asctoolkit {
namespace lumberjack {

Message* mpiBlockingRecieveAnyMessage(MPI_Comm comm, int ranksLimit)
{
    char* charArray = ATK_NULLPTR;
    int messageSize = -1;
    MPI_Status mpiStatus;

    // Get size and source of MPI message
    MPI_Probe(MPI_ANY_SOURCE, 0, comm, &mpiStatus);
    MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);

    // Recieve packed Message
    charArray = (char*)malloc((messageSize+1)*sizeof(char));
    MPI_Recv(charArray, messageSize, MPI_CHAR, mpiStatus.MPI_SOURCE, 0, comm, &mpiStatus);

    // If count is 1, this means the sending node is out of messages to send
    if (messageSize == 1) {
        return ATK_NULLPTR;
    }

    // Unpack Message
    charArray[messageSize] = '\0';
    std::string packedMessage(charArray);
    free(charArray);
    Message* message = new Message();
    message->unpack(packedMessage, ranksLimit);

    return message;
}

} // end namespace lumberjack
} // end namespace asctoolkit
