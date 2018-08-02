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

#ifndef QUEST_MPICOMM_WRAPPER_HPP_
#define QUEST_MPICOMM_WRAPPER_HPP_

#include "axom/config.hpp" // for Axom compile-time definitions

/*!
 * \file
 *
 * \brief Provides aliases for MPI_Comm and MPI_COMM_SELF, which are used in
 *  the interface, when Axom is compiled without MPI. This avoids complicated
 *  ifdef logic at the user-facing API, which was otherwise, necessary.
 *
 * \note When Axom is compiled with MPI, this file just includes the MPI header.
 */

#ifdef AXOM_USE_MPI
  #include <mpi.h>
#else
using MPI_Comm = int;
constexpr int MPI_COMM_SELF = -1;
#endif

#endif /* QUEST_MPICOMM_WRAPPER_HPP_ */
