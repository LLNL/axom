// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 * \file Timer.hpp
 *
 * \brief Defines a simple Timer class to measure execution time.
 *
 ******************************************************************************
 */

#ifndef TIMER_HPP_
#define TIMER_HPP_

#include "axom/config.hpp"

#include <chrono>

namespace axom
{
namespace utilities
{
/*!
 * \brief A simple Timer class to measure execution time.
 *
 *  \note We might want to extend the functionality of the timer class
 *   by making HighPrecisionTimer a template parameter.
 *        API requirements for HighPrecisionTimer class
 *        -- must be default constructible
 *        -- timing functions:
 *              void start()
 *              void stop()
 *              void reset()
 *        -- elapsed time functions:
 *              double elapsedTimeInSec()
 *              double elapsedTimeInMilliSec()
 *              double elapsedTimeInMicroSec()
 *
 *  \note We might want to add support for pausing and resuming the timer while
 *    accumulating the time differences
 *
 *  Example Usage:
 *  \code
 *
 *     utilities::Timer t;
 *     t.start();
 *
 *     // code to measure its execution time
 *
 *     t.stop();
 *     std::cout << "Elapsed Time: << t.elapsed() << std::endl;
 *
 *     t.reset();
 *
 *  \endcode
 */
class Timer
{
private:
  using ClockType = std::chrono::high_resolution_clock;
  using TimeStruct = std::chrono::time_point<ClockType>;
  using TimeDiff = std::chrono::duration<double>;
  using MilliTimeDiff = std::chrono::duration<double, std::milli>;
  using MicroTimeDiff = std::chrono::duration<double, std::micro>;

public:
  /*!
   * \brief Default constructor.
   * \param startRunning Indicates whether to start the timer
   *        during construction (default is false)
   */
  Timer(bool startRunning = false) : m_running(startRunning)
  {
    if(m_running) start();
  }

  /*!
   * \brief Starts the timer.Sets the start time of this Timer instance.
   */
  void start()
  {
    m_running = true;
    m_startTime = ClockType::now();
  }

  /*!
   * \brief Stops the timer. Sets the end time of this Timer instance.
   */
  void stop()
  {
    m_stopTime = ClockType::now();
    m_running = false;
  }

  /*!
   * \brief Returns the elapsed time in seconds.
   * \return t the elapsed time in seconds.
   */
  double elapsed() { return elapsedTimeInSec(); };

  /*!
   * \brief Returns the elapsed time in seconds.
   * \return t the elapsed time in seconds.
   */
  double elapsedTimeInSec()
  {
    if(m_running) stop();
    return clockDiff().count();
  }

  /*!
   * \brief Returns the elapsed time in milliseconds.
   * \return t the elapsed time in milliseconds.
   */
  double elapsedTimeInMilliSec()
  {
    if(m_running) stop();
    return std::chrono::duration_cast<MilliTimeDiff>(clockDiff()).count();
  }

  /*!
   * \brief Returns the elapsed time in microseconds.
   * \return t the elapsed time in microseconds.
   */
  double elapsedTimeInMicroSec()
  {
    if(m_running) stop();
    return std::chrono::duration_cast<MicroTimeDiff>(clockDiff()).count();
  }

  /*!
   * \brief Resets the timer.
   * \post this->elapsed()==0.0
   */
  void reset()
  {
    m_running = false;
    m_startTime = m_stopTime = TimeStruct();
  }

private:
  /*! \brief Computes the difference between start() and stop() */
  TimeDiff clockDiff() const { return m_stopTime - m_startTime; }

  TimeStruct m_startTime;
  TimeStruct m_stopTime;
  bool m_running;
};

}  // namespace utilities
}  // namespace axom

#endif  // TIMER_HPP_
