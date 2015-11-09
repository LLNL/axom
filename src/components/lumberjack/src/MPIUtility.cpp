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
#include <cstring>
#include <iostream>

namespace asctoolkit {
namespace lumberjack {

const char* mpiBlockingRecieveMessages(MPI_Comm comm)
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

    if (messageSize == 1) {
        delete charArray;
        return ATK_NULLPTR;
    }

    return charArray;
}

void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank, const char* packedMessagesToBeSent)
{
    MPI_Request mpiRequest;
    MPI_Isend(const_cast<char*>(packedMessagesToBeSent),
             strlen(packedMessagesToBeSent), MPI_CHAR, destinationRank, 0, comm, &mpiRequest);
}

} // end namespace lumberjack
} // end namespace asctoolkit
