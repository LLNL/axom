// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file TickCountTimer.hpp
 *
 * \brief A Windows-based timer implementation for axom's Timer class
 *
 * \note TickCountTimer is an internal helper class, not meant for external use.
 * It is intended to be used by axom's Timer class in Windows configurations.
 *******************************************************************************
 */

#ifndef TICK_COUNT_TIMER_HPP_
#define TICK_COUNT_TIMER_HPP_

#ifdef WIN32

#define NOMINMAX
#include <Windows.h>

namespace axom
{
namespace utilities
{
namespace detail
{

/**
 * \class
 * \brief A simple timer utility based on the Windows GetTickCount64 function
 *
 * \note This is a simple class without any checks to ensure proper usage of the
 *  timer. It is meant to be used as a base class for the Timer class in
 *  axom/core/Timer.hpp Specifically, we do not check that start() was called
 *  before stop(), or if stop() was called before attempting to find the elapsed
 *  time.
 */
class TickCountTimer
{
private:
  typedef ULONGLONG TimeStruct;
  typedef long int TimeDiff;

  enum { TIMER_ONE      = 1,
         TIMER_THOUSAND = 1000,
         TIMER_MILLION  = 1000000 };

public:
  TickCountTimer() { reset(); }

  /** \brief Sets the start time of the timer */
  void start() {m_startTime = GetTickCount64(); }

  /** \brief Sets the stop time of the timer */
  void stop() {m_stopTime = GetTickCount64(); }

  /**  \brief Resets the timer */
  void reset()
  {
    m_startTime=0;
    m_stopTime=0;
  }

  /** \brief Returns the number of seconds between start() and stop() */
  double elapsedTimeInSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_MILLION);
  }

  /** \brief Returns the number of milliseconds between start() and stop() */
  double elapsedTimeInMilliSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_THOUSAND);
  }

  /** \brief Returns the number of microseconds between start() and stop() */
  double elapsedTimeInMicroSec() const
  {
    return clockDiff() / static_cast<double>(TIMER_ONE);
  }

private:
  /** \brief Computes the time difference between start() and stop() */
  TimeDiff clockDiff() const
  {
    const TimeDiff sDiff = static_cast<TimeDiff>(m_stopTime - m_startTime);
    return sDiff * TIMER_THOUSAND;
  }

private:
  TimeStruct m_startTime;
  TimeStruct m_stopTime;
};


} /* namespace detail */
} /* namespace utilities */
} /* namespace axom */

#endif // WIN32

#endif // TICK_COUNT_TIMER_HPP_
