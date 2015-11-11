#include <sys/time.h>
#include "common/CommonTypes.hpp"

namespace tinyHydro {

/**
 * \class
 * \brief A simple timer class with functionality to start and stop the timer, and to return the time difference in seconds, milliseconds and microseconds
 * Support
 */
  class Timer {
  public:
    Timer(bool startRunning = false) : m_running(startRunning)
    {
      if(m_running)
        start();
    }

    /**
     * Starts the timer (overwriting all previous start times)
     */
    void start()
    {
      m_running = true;
      gettimeofday(&m_startTime, ATK_NULLPTR);
    }

    /**
     * Stops the timer (overwriting all previous stop times)
     */
    void stop()
    {
      gettimeofday(&m_endTime, ATK_NULLPTR);
      m_running = false;
    }

    /**
     * Gets the time difference, in fractions of a second, between the start time and stop time.
     */
    double  getElapsedTime()              { return getElapsedTimeInSec(); }

    /**
     * Gets the time difference, in fractions of a second, between the start time and stop time.
     */
    double  getElapsedTimeInSec()         { return clockDiff() / static_cast<double>(MILLION); }

    /**
     * Gets the time difference, in fractions of a millisecond, between the start time and stop time.
     */
    double  getElapsedTimeInMilliSec()    { return clockDiff() / static_cast<double>(THOUSAND); }

    /**
     * Gets the time difference, in microseconds, between the start time and stop time.
     */
    double  getElapsedTimeInMicroSec()    { return clockDiff() / static_cast<double>(ONE); }
  private:

    static const long int ONE = 1;
    static const long int THOUSAND = 1000;
    static const long int MILLION = 1000000;

    /**
     * Computes the number of microseconds between the start and stop times.
     * Stops the timer if it is running at invocation time.
     */
    long int clockDiff()
    {
      if(m_running)
        stop();

      const long int sDiff = m_endTime.tv_sec - m_startTime.tv_sec;
      const long int uDiff = m_endTime.tv_usec - m_startTime.tv_usec;
      return sDiff * MILLION + uDiff;
    }
  private:
    timeval m_startTime, m_endTime;
    bool m_running;
  };

} // end namespace tinyHydro
