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
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for utility functions.
 *
 ******************************************************************************
 */

#include "axom_common/Utilities.hpp"

#include <cstdlib> // for exit, EXIT_SUCCESS, EXIT_FAILURE

#ifdef AXOM_USE_MPI
#include <mpi.h>
#endif

namespace axom {
namespace utilities {

  void processAbort()
  {
#ifndef AXOM_USE_MPI
    exit( EXIT_FAILURE );
#else
    int mpi = 0;
    MPI_Initialized( &mpi );
    if ( mpi )
    {
      MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }
    else
    {
      exit( EXIT_FAILURE );
    }
#endif
  }

}   // end namespace utilities
}   // end namespace axom
