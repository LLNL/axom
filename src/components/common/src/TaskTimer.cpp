/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file TaskTimer.cpp
 *
 * \date Feb 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "TaskTimer.hpp"

// C/C++ includes
#ifdef USE_CXX11
#include <chrono>
#else
#include <cstddef>
#include <sys/time.h>
#endif

namespace asctoolkit {

namespace utilities {


TaskTimer::TaskTimer() :
        m_startTime(0.0),
        m_endTime(0.0)
{

}

//------------------------------------------------------------------------------
TaskTimer::~TaskTimer()
{

}

//------------------------------------------------------------------------------
double TaskTimer::getCurrentTime()
{
  double t = 0.0;

#ifdef USE_CXX11

  /* aliases */
  using sec   = std::chrono::duration< double >;
  using clock = std::chrono::high_resolution_clock;

  sec ct;
  ct = std::chrono::duration_cast< sec >( clock::now().time_since_epoch() );
  t = ct.count();

#else

  struct timeval ct;
  gettimeofday( &ct, NULL );
  t = ct.tv_sec + (ct.tv_usec / 1000000.);

#endif

  return ( t );
}

} /* namespace utilities */

} /* namespace asctoolkit */
