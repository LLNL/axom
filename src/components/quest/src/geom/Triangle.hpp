/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file Triangle.hpp
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_

#include "quest/Point.hpp"

namespace quest
{

template < typename T, int DIM >
class Triangle
{
public:
    typedef Point<T,DIM> PointType;

public:

  /*!
   *****************************************************************************
   * \brief Constructor -- creates a triangle from the 3 points A,B,C.
   *****************************************************************************
   */
  Triangle( const PointType& A,
            const PointType& B,
            const PointType& C )
        : m_A (A), m_B(B), m_C(C)    {}


  /*!
   *****************************************************************************
   * \brief Destructor
   *****************************************************************************
   */
   ~Triangle() { }

  /*!
   *****************************************************************************
   * \brief
   * \return A
   *****************************************************************************
   */
  const PointType& A() const { return m_A; };

  /*!
   *****************************************************************************
   * \brief
   * \return B
   *****************************************************************************
   */
  const PointType& B() const { return m_B; };

  /*!
   *****************************************************************************
   * \brief
   * \return C
   *****************************************************************************
   */
  const PointType& C() const { return m_C; };

private:

  /*!
   *****************************************************************************
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use.
   *****************************************************************************
   */
  Triangle() { }

  PointType m_A;
  PointType m_B;
  PointType m_C;
};

} /* namespace quest */

#endif /* TRIANGLE_HPP_ */
