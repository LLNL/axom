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
 * \brief   Header file containing utility functions.
 *
 ******************************************************************************
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

#ifdef USE_MPI
#include <mpi.h>
#endif

namespace asctoolkit
{
namespace utilities
{

   /*!
    * \brief Gracefully aborts the application
    */
   inline void processAbort()
   {
#ifndef USE_MPI
     exit( -1 );
#else
     int mpi = 0;
     MPI_Initialized( &mpi );
     if ( mpi ) {

       MPI_Abort( MPI_COMM_WORLD, -1 );

     } else {

       exit( -1 );

     }
#endif
   }

}  // ending brace for utilities namespace
}  // ending brace for asctoolkit namespace

#endif  // header include guard
