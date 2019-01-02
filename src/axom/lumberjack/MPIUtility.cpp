/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*!
 ******************************************************************************
 *
 * \file MPIUtility.cpp
 *
 * \brief Implementation of the MPI utility methods.
 *
 ******************************************************************************
 */

#include "axom/lumberjack/MPIUtility.hpp"

#include <cstdlib>
#include <cstring>

namespace axom
{
namespace lumberjack
{

const char* mpiBlockingRecieveMessages(MPI_Comm comm)
{
  char* charArray = nullptr;
  int messageSize = -1;
  MPI_Status mpiStatus;

  // Get size and source of MPI message
  MPI_Probe(MPI_ANY_SOURCE, 0, comm, &mpiStatus);
  MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);

  // Setup where to recieve the char array
  charArray = new char[messageSize+1];
  charArray[messageSize] = '\0';

  // Receive packed Message
  MPI_Recv(charArray, messageSize, MPI_CHAR, mpiStatus.MPI_SOURCE, 0, comm,
           &mpiStatus);

  return charArray;
}

void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank,
                                const char* packedMessagesToBeSent)
{
  MPI_Request mpiRequest;
  MPI_Isend(const_cast<char*>(packedMessagesToBeSent),
            strlen(packedMessagesToBeSent), MPI_CHAR,
            destinationRank, 0, comm, &mpiRequest);
}

} // end namespace lumberjack
} // end namespace axom
