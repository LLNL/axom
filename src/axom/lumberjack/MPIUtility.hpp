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
 * \brief Sends all Message sent to the given rank.
 *
 * Clears and deletes all Message classes when done.
 *
 * \param [in] comm The MPI Communicator.
 * \param [in] destinationRank Where the Message classes is being sent.
 * \param [in,out] packedMessagesToBeSent All of the Message classes to be sent
 *  packed together.
 *****************************************************************************
 */
void mpiNonBlockingSendMessages(MPI_Comm comm, int destinationRank,
                                const char* packedMessagesToBeSent);
} // end namespace lumberjack
} // end namespace axom

#endif
