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

#ifndef TASKTIMER_HPP_
#define TASKTIMER_HPP_

namespace asctoolkit {

namespace utilities {


/*!
 *******************************************************************************
 * \brief A simple TaskTimer class used to measure execution time.
 *
 *  Example Usage:
 *  \code
 *
 *     utilities::TaskTimer t;
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
class TaskTimer
{
public:

  /*!
   *****************************************************************************
   * \brief Default constructor.
   *****************************************************************************
   */
  TaskTimer();

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~TaskTimer();

  /*!
   *****************************************************************************
   * \brief Starts the timer.Sets the start time of this TaskTimer instance.
   *****************************************************************************
   */
  void start() { m_startTime = this->getCurrentTime(); };

  /*!
   *****************************************************************************
   * \brief Stops the timer. Sets the end time of this TaskTimer instance.
   *****************************************************************************
   */
  void stop() { m_endTime = this->getCurrentTime(); };

  /*!
   *****************************************************************************
   * \brief Returns the elapsed time in seconds.
   * \return t the elapsed time in seconds.
   *****************************************************************************
   */
  double elapsed() { return (m_endTime-m_startTime);};

  /*!
   *****************************************************************************
   * \brief Resets the timer.
   * \post this->elapsed()==0.0
   *****************************************************************************
   */
  void reset() { m_startTime = m_endTime = 0.0; };

private:
    double m_startTime;
    double m_endTime;

    double getCurrentTime();
};

} /* namespace utilities */

} /* namespace asctoolkit */

#endif /* TASKTIMER_HPP_ */
