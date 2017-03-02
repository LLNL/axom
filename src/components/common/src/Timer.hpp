/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef TIMER_HPP_
#define TIMER_HPP_

#include "common/config.hpp"
#ifdef ATK_USE_CXX11
  #include "common/ChronoTimer.hpp"
#else
  #include "common/TimeofdayTimer.hpp"
#endif

namespace {
#ifdef ATK_USE_CXX11
  typedef axom::utilities::detail::ChronoTimer HighPrecisionTimer;
#else
  typedef axom::utilities::detail::TimeofdayTimer HighPrecisionTimer;
#endif
}



namespace axom {
namespace utilities {

/*!
 *******************************************************************************
 * \brief A simple Timer class to measure execution time.
 *
 * \note The actual timing functionality is implemented using a HighPrecisionTimer
 *  instance.  These are located in the detail namespace using the chrono library in C++11
 *  and glibc gettimeofday() otherwise.
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
 *******************************************************************************
 */

class Timer
{
public:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   * \param startRunning Indicates whether to start the timer
   *        during construction (default is false)
   *****************************************************************************
   */
  Timer(bool startRunning = false): m_running(startRunning)
  {
      if(m_running)
          m_hpTimer.start();
  }

  /*!
   *****************************************************************************
   * \brief Starts the timer.Sets the start time of this Timer instance.
   *****************************************************************************
   */
  void start() { m_running = true; m_hpTimer.start(); }

  /*!
   *****************************************************************************
   * \brief Stops the timer. Sets the end time of this Timer instance.
   *****************************************************************************
   */
  void stop() { m_hpTimer.stop(); m_running = false; }

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in seconds.
   * \return t the elapsed time in seconds.
   *****************************************************************************
   */
  double elapsed() { return elapsedTimeInSec();};

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in seconds.
   * \return t the elapsed time in seconds.
   *****************************************************************************
   */
  double elapsedTimeInSec()
  {
      if(m_running)
          m_hpTimer.stop();
      return m_hpTimer.elapsedTimeInSec();
  }

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in milliseconds.
   * \return t the elapsed time in milliseconds.
   *****************************************************************************
   */
  double elapsedTimeInMilliSec()
  {
      if(m_running)
          m_hpTimer.stop();
      return m_hpTimer.elapsedTimeInMilliSec();
  }

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in microseconds.
   * \return t the elapsed time in microseconds.
   *****************************************************************************
   */
  double elapsedTimeInMicroSec()
  {
      if(m_running)
          m_hpTimer.stop();
      return m_hpTimer.elapsedTimeInMicroSec();
  }

  /*!
   *****************************************************************************
   * \brief Resets the timer.
   * \post this->elapsed()==0.0
   *****************************************************************************
   */
  void reset() { m_running = false; m_hpTimer.reset(); }

private:
    HighPrecisionTimer  m_hpTimer;
    bool                m_running;
};

} // namespace utilities 
} // namespace axom 

#endif // TIMER_HPP_ 
