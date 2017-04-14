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
 * \file MPIUtility.hpp
 * \author Chris White (white238@llnl.gov)
 *
 * \brief This file contains the definitions of utility functions.
 *******************************************************************************
 */

#ifndef MPIUTILITY_HPP
#define MPIUTILITY_HPP

#include "mpi.h"

#include "lumberjack/Message.hpp"

namespace axom {
namespace lumberjack {

/*!
*****************************************************************************
* \brief Recieves any Message sent to this rank. Returns null if terminating
* message is sent.
*
* \param [in] comm The MPI Communicator.
*****************************************************************************
*/
const char* mpiBlockingRecieveMessages(MPI_Comm comm);

/*!
*****************************************************************************
* \brief Sends all Message sent to the given rank.
*
* Clears and deletes all Message classes when done.
*
* \param [in] comm The MPI Communicator.
* \param [in] destinationRank Where the Message classes is being sent.
* \param [in,out] packedMessagesToBeSent All of the Message classes to be sent packed together.
*****************************************************************************
*/
void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank, const char* packedMessagesToBeSent);
} // end namespace lumberjack
} // end namespace axom

#endif
