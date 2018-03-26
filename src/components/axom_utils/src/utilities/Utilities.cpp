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
 *
 * \file
 *
 * \brief   Implementation file for utility functions.
 *
 */

#include "axom/config.hpp"
#include "axom_utils/Utilities.hpp"

#include <cstdlib> // for exit, EXIT_SUCCESS, EXIT_FAILURE

#ifdef AXOM_USE_MPI
#include <mpi.h>
#endif

namespace axom
{
namespace utilities
{

void processAbort()
{
#ifndef AXOM_USE_MPI
  abort();
#else
  int mpi = 0;
  MPI_Initialized( &mpi );
  if ( mpi )
  {
    MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
  }
  else
  {
    abort();
  }
#endif
}

}   // end namespace utilities
}   // end namespace axom
