// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TimeOfDayTimer.hpp
 *
 * \brief A glibc based timer implementation for axom's Timer class
 *
 * \note TimeOfDayTimer is an internal helper class, not meant for external use.
 * It is intended to be used by axom's Timer class in unix-based configurations.
 *******************************************************************************
 */

#ifndef TIMEOFDAY_TIMER_HPP_
#define TIMEOFDAY_TIMER_HPP_

#include <sys/time.h>       // for gettimeofday() and timeval
                            // Note: located in <time> on some systems

namespace axom
{
namespace utilities
{
namespace detail
{

/*!
 * \class
 * \brief A simple timer utility based on the glibc gettimeofday() function
 *
 * \note This is a simple class without any checks to ensure proper usage of the
 *  timer. It is meant to be used as a base class for the Timer class in
 *  axom/core/Timer.hpp Specifically, we do not check that start() was called
 *  before stop(), or if stop() was called before attempting to find the elapsed
 *  time.
 */
class TimeofdayTimer
{
private:
  typedef timeval TimeStruct;
  typedef long int TimeDiff;

  enum { TIMER_ONE      = 1,
         TIMER_THOUSAND = 1000,
         TIMER_MILLION  = 1000000 };

public:
  /*! \brief Constructor for TimeOfDayTimer instance */
  TimeofdayTimer() { reset(); }

  /*! \brief Sets the start time of the timer */
  void start() { gettimeofday(&m_startTime, nullptr); }

  /*! \brief Sets the stop time of the timer */
  void stop()  { gettimeofday(&m_stopTime, nullptr); }

  /*!  \brief Resets the timer */
  void reset()
  {
    m_startTime = (struct timeval){0,0};
    m_stopTime =  (struct timeval){0,0};
  }

  /*! \brief Returns the number of seconds between start() and stop() */
  double elapsedTimeInSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_MILLION);
  }

  /*! \brief Returns the number of milliseconds between start() and stop() */
  double elapsedTimeInMilliSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_THOUSAND);
  }

  /*! \brief Returns the number of microseconds between start() and stop() */
  double elapsedTimeInMicroSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_ONE);
  }

private:
  /*! \brief Computes the time difference between start() and stop() */
  TimeDiff clockDiff() const
  {
    const TimeDiff sDiff = m_stopTime.tv_sec - m_startTime.tv_sec;
    const TimeDiff uDiff = m_stopTime.tv_usec - m_startTime.tv_usec;
    return sDiff * TIMER_MILLION + uDiff;
  }

private:
  TimeStruct m_startTime;
  TimeStruct m_stopTime;
};


} /* namespace detail */
} /* namespace utilities */
} /* namespace axom */

#endif // TIMEOFDAY_TIMER_HPP_
