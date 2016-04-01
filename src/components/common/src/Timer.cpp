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

#include "Timer.hpp"
#include "common/CommonTypes.hpp"   // For ATK_NULLPTR



namespace asctoolkit {

namespace utilities {


#ifdef USE_CXX11
    // Define static constexpr member variable s_defaultTime
    constexpr Timer::TimeStruct Timer::s_defaultTime;
#else
    // Initialize static const member variable s_defaultTime
    const Timer::TimeStruct Timer::s_defaultTime = (struct timeval){0,0};
#endif




//------------------------------------------------------------------------------
Timer::Timer(bool startRunning) :
          m_startTime(s_defaultTime)
        , m_stopTime(s_defaultTime)
        , m_running(startRunning)
{
    if(m_running)
        start();
}

//------------------------------------------------------------------------------
Timer::~Timer()
{

}

//------------------------------------------------------------------------------
void Timer::start()
{
    m_running = true;

#ifdef USE_CXX11
    m_startTime = ClockType::now();
#else
    gettimeofday(&m_startTime, ATK_NULLPTR);
#endif
}

//------------------------------------------------------------------------------
void Timer::stop()
{
#ifdef USE_CXX11
    m_stopTime = ClockType::now();
#else
    gettimeofday(&m_stopTime, ATK_NULLPTR);
#endif

    m_running = false;
}


//------------------------------------------------------------------------------
void Timer::reset()
{
    m_running = false;
    m_startTime = s_defaultTime;
    m_stopTime = s_defaultTime;
}


//------------------------------------------------------------------------------
Timer::TimeDiff Timer::clockDiff()
{
    if(m_running)
        stop();

#ifdef USE_CXX11
    return m_stopTime - m_startTime;
#else
    const long int sDiff = m_stopTime.tv_sec - m_startTime.tv_sec;
    const long int uDiff = m_stopTime.tv_usec - m_startTime.tv_usec;
    return sDiff * TIMER_MILLION + uDiff;
#endif
}


//------------------------------------------------------------------------------
double Timer::elapsedTimeInSec()
{
#ifdef USE_CXX11
    return clockDiff().count();
#else
    return clockDiff() / static_cast<double>(TIMER_MILLION);
#endif
}

double Timer::elapsedTimeInMilliSec()
{
#ifdef USE_CXX11
    typedef std::chrono::duration<double, std::milli> MilliTimeDiff;
    return std::chrono::duration_cast< MilliTimeDiff >( clockDiff() ).count();
#else
    return clockDiff() / static_cast<double>(TIMER_THOUSAND);
#endif
}

double Timer::elapsedTimeInMicroSec()
{
#ifdef USE_CXX11
    typedef std::chrono::duration<double, std::micro> MicroTimeDiff;
    return std::chrono::duration_cast< MicroTimeDiff >( clockDiff() ).count();
#else
    return clockDiff() / static_cast<double>(TIMER_ONE);
#endif
}


} /* namespace utilities */

} /* namespace asctoolkit */
