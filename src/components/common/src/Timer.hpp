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
 * \file TaskTimer.hpp
 *
 * \date Feb 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */


#ifdef USE_CXX11
  #include <chrono>
#else
  #include <sys/time.h>
#endif



#ifndef TIMER_HPP_
#define TIMER_HPP_

namespace asctoolkit {

namespace utilities {


/*!
 *******************************************************************************
 * \brief A simple Timer class used to measure execution time.
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
 *  \endcode
 *******************************************************************************
 */
class Timer
{
private:
  // Some typedefs for internal implementation
  #ifdef USE_CXX11
    typedef std::chrono::high_resolution_clock  ClockType;
    typedef std::chrono::time_point<ClockType>  TimeStruct;
    typedef std::chrono::duration<double>       TimeDiff;

    static constexpr TimeStruct s_defaultTime = TimeStruct();
  #else
    typedef timeval                             TimeStruct;
    typedef long int                            TimeDiff;

    static const TimeStruct s_defaultTime;

    enum { TIMER_ONE      = 1
         , TIMER_THOUSAND = 1000
         , TIMER_MILLION  = 1000000
    };

  #endif

public:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   * \param startRunning Indicates whether to start the timer
   *        during construction (default is false)
   *****************************************************************************
   */
  Timer(bool startRunning = false);

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~Timer();

  /*!
   *****************************************************************************
   * \brief Starts the timer.Sets the start time of this Timer instance.
   *****************************************************************************
   */
  void start();
  /*!
   *****************************************************************************
   * \brief Stops the timer. Sets the end time of this Timer instance.
   *****************************************************************************
   */
  void stop();

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
  double elapsedTimeInSec();

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in milliseconds.
   * \return t the elapsed time in milliseconds.
   *****************************************************************************
   */
  double elapsedTimeInMilliSec();

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in microseconds.
   * \return t the elapsed time in microseconds.
   *****************************************************************************
   */
  double elapsedTimeInMicroSec();

  /*!
   *****************************************************************************
   * \brief Resets the timer.
   * \post this->elapsed()==0.0
   *****************************************************************************
   */
  void reset();

private:
    TimeDiff clockDiff();

private:

    TimeStruct m_startTime;
    TimeStruct m_stopTime;
    bool       m_running;

};

} /* namespace utilities */

} /* namespace asctoolkit */

#endif /* TIMER_HPP_ */
