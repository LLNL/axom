// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ChronoTimer.hpp
 *
 * \brief A C++11 chrono-based timer implementation for axom's Timer class
 *
 * \note ChronoTimer is an internal helper class, not meant for external usage.
 * It is intended to be used by axom's Timer class in C++11 configurations.
 *******************************************************************************
 */

#ifndef CHRONO_TIMER_HPP_
#define CHRONO_TIMER_HPP_

#include <chrono>

namespace axom
{
namespace utilities
{
namespace detail
{

/*!
 * \class
 *
 * \brief A simple timer utility based on the C++11 chrono library
 *
 * \note This is a simple class without any checks to ensure proper usage of the
 *  timer. It is meant to be used as a base class for the Timer class in
 *  axom/core/Timer.hpp. Specifically, we do not check that start() was called
 *  before stop(), or if stop() was called before attempting to find the elapsed
 *  time.
 */
class ChronoTimer
{
private:
  typedef std::chrono::high_resolution_clock ClockType;
  typedef std::chrono::time_point<ClockType>        TimeStruct;
  typedef std::chrono::duration<double>             TimeDiff;
  typedef std::chrono::duration<double, std::milli> MilliTimeDiff;
  typedef std::chrono::duration<double, std::micro> MicroTimeDiff;

public:
  /*!  \brief Constructor for a ChronoTimer instance */
  ChronoTimer() : m_startTime( TimeStruct() ), m_stopTime( TimeStruct() ) {}

  /*!  \brief Sets the start time of the timer */
  void start() { m_startTime = ClockType::now(); }

  /*!  \brief Sets the stop time of the timer */
  void stop()  { m_stopTime = ClockType::now(); }

  /*!  \brief Resets the timer */
  void reset() { m_startTime = m_stopTime = TimeStruct(); }

  /*!  \brief Returns the number of seconds between start() and stop() */
  double elapsedTimeInSec() const
  {
    return clockDiff().count();
  }

  /*!  \brief Returns the number of milliseconds between start() and stop() */
  double elapsedTimeInMilliSec() const
  {
    return std::chrono::duration_cast< MilliTimeDiff >( clockDiff() ).count();
  }

  /*! \brief Returns the number of microseconds between start() and stop() */
  double elapsedTimeInMicroSec() const
  {
    return std::chrono::duration_cast< MicroTimeDiff >( clockDiff() ).count();
  }

private:

  /*! \brief Computes the difference between start() and stop() */
  TimeDiff clockDiff() const { return m_stopTime - m_startTime; }

private:
  TimeStruct m_startTime;
  TimeStruct m_stopTime;
};

} /* namespace detail */
} /* namespace utilities */
} /* namespace axom */

#endif // CHRONO_TIMER_HPP_
