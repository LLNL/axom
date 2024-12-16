// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file MPIUtility.hpp
 *
 * \brief This file contains the definitions of MPI utility functions.
 *******************************************************************************
 */

#ifndef MPIUTILITY_HPP
#define MPIUTILITY_HPP

#include "mpi.h"

namespace axom
{
namespace lumberjack
{
/*!
 *****************************************************************************
 * \brief Receives any Message sent to this rank. Returns null if terminating
 *  message is sent.
 *
 * \param [in] comm The MPI Communicator.
 *****************************************************************************
 */
const char* mpiBlockingReceiveMessages(MPI_Comm comm);

/*!
 *****************************************************************************
 * \brief Receives any Message sent to this rank, if there are any messages 
 *  that are sent. Returns null if no messages are sent.
 *
 * \param [in] comm The MPI Communicator.
 * \param [in] tag The MPI tag to use for communication.  When set to zero, 
 *  MPI communication uses default LJ_Tag. 
 *****************************************************************************
 */
const char* mpiNonBlockingReceiveMessages(MPI_Comm comm, int tag = 0);

/*!
 *****************************************************************************
 * \brief Sends all Message sent to the given rank.
 *
 * Clears and deletes all Message classes when done.
 *
 * \param [in] comm The MPI Communicator.
 * \param [in] destinationRank Where the Message classes is being sent.
 * \param [in,out] packedMessagesToBeSent All of the Message classes to be sent
 *  packed together.
 * \param [in] tag The MPI tag to use for communication.  When set to zero, 
 *  MPI communication uses default LJ_Tag.
 *****************************************************************************
 */
void mpiNonBlockingSendMessages(MPI_Comm comm,
                                int destinationRank,
                                const char* packedMessagesToBeSent,
                                int tag = 0);
}  // end namespace lumberjack
}  // end namespace axom

#endif
