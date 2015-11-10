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


#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cmath>

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

   /*!
    * \brief Fuzzy comparison of two real valued quantities
    *
    * \param a,b The real valued quantities we are comparing
    * \param thresh The threshold of the fuzzy comparison.  Default is 1.0e-8
    * \return \c true if the absolute value of the difference is less than \param thresh and false otherwise
    */
   template<typename RealType>
   bool compareReals(RealType a, RealType b, RealType thresh = 1.0e-8)
   {
       return std::fabs(a-b) <= thresh;
   }

   /*!
    * \brief Fuzzy comparison of two real valued quantities
    *
    * \param a,b The real valued quantities we are comparing
    * \param relThresh The relative threshold of the fuzzy comparison.  Default is 1.0e-6
    * \param absThresh The absolute threshold of the fuzzy comparison.  Default is 1.0e-8
    * \return \c true if the absolute value of the difference is less than the sum of \param absThresh
    *         and the relative difference (\param relThresh times the absolute max of a and b)
    */
   template<typename RealType>
   bool compareRealsRelative(RealType a, RealType b, RealType relThresh = 1.0e-6, RealType absThresh = 1.0e-8)
   {
       RealType maxFabs = std::max(std::fabs(a), std::fabs(b) );
       return std::fabs(a-b) <= ( maxFabs * relThresh + absThresh);
   }


}  // ending brace for utilities namespace
}  // ending brace for asctoolkit namespace

#endif  // header include guard
