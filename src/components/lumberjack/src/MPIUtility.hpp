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
*****************************************************************************
*/
Message* mpiBlockingRecieveAnyMessage(MPI_Comm comm, int ranksLimit);

} // end namespace lumberjack
} // end namespace asctoolkit

#endif
