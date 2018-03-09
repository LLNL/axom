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

#include "lumberjack/MPIUtility.hpp"

#include "axom/Types.hpp"

#include <cstdlib>
#include <cstring>

namespace axom
{
namespace lumberjack
{

const char* mpiBlockingRecieveMessages(MPI_Comm comm)
{
  char* charArray = AXOM_NULLPTR;
  int messageSize = -1;
  MPI_Status mpiStatus;

  // Get size and source of MPI message
  MPI_Probe(MPI_ANY_SOURCE, 0, comm, &mpiStatus);
  MPI_Get_count(&mpiStatus, MPI_CHAR, &messageSize);

  // Receive packed Message
  charArray = new char[messageSize+1];
  MPI_Recv(charArray, messageSize, MPI_CHAR, mpiStatus.MPI_SOURCE, 0, comm,
           &mpiStatus);
  charArray[messageSize] = '\0';

  if (messageSize == 1)
  {
    delete [] charArray;
    return AXOM_NULLPTR;
  }

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
