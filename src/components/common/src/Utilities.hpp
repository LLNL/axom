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

#include <cmath>
#include <algorithm>

namespace asctoolkit {
namespace utilities {

/*!
 * \brief Gracefully aborts the application
 */
  void processAbort();

/*!
 * \brief Fuzzy comparison of two real valued quantities
 *
 * \param a The first real valued quantities we are comparing
 * \param b The second real valued quantities we are comparing
 * \param thresh The threshold of the fuzzy comparison.  Default is 1.0e-8
 * \return \c true if the absolute value of the difference is less than thresh and false otherwise
 */
  template<typename RealType>
  bool isNearlyEqual(RealType a, RealType b, RealType thresh = 1.0e-8)
  {
    return std::fabs(a-b) <= thresh;
  }

/*!
 * \brief Fuzzy comparison of two real valued quantities
 *
 * \param a The first real valued quantities we are comparing
 * \param b The second real valued quantities we are comparing
 * \param relThresh The relative threshold of the fuzzy comparison.  Default is 1.0e-6
 * \param absThresh The absolute threshold of the fuzzy comparison.  Default is 1.0e-8
 * \return \c true if the absolute value of the difference is less than the sum of absThresh
 *         and the relative difference (relThresh times the absolute max of a and b)
 */
  template<typename RealType>
  bool isNearlyEqualRelative(RealType a, RealType b, RealType relThresh = 1.0e-6,
                            RealType absThresh = 1.0e-8)
  {
    RealType maxFabs = std::max(std::fabs(a), std::fabs(b) );
    return std::fabs(a-b) <= ( maxFabs * relThresh + absThresh);

    // Equation from Real-Time Collsion Detection book -- http://realtimecollisiondetection.net/pubs/Tolerances/
    // Note: If we use this, we must update the doxygen
    // return std::fabs(a-b) <=  std::max(absThresh, relThresh * maxFabs );
  }


}  // ending brace for utilities namespace
}  // ending brace for asctoolkit namespace

#endif  // header include guard
