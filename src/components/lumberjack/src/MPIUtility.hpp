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

namespace asctoolkit {
namespace lumberjack {

/*!
*****************************************************************************
* \brief Recieves any Message sent to this rank. Returns null if terminating
* message is sent.
*
* \param [in] comm The MPI Communicator.
* \param [in] ranksLimit Limit on how many ranks are individually tracked per Message.
*****************************************************************************
*/
Message* mpiBlockingRecieveAnyMessage(MPI_Comm comm, int ranksLimit);

/*!
*****************************************************************************
* \brief Sends all Message sent to the given rank.
*
* Clears and deletes all Message classes when done.
*
* \param [in] comm The MPI Communicator.
* \param [in] destinationRank Where the Message classes is being sent.
* \param [in,out] messages All of the Message classes to be sent.
*****************************************************************************
*/
void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank, std::vector<Message*>& messages);

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
