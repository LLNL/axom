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

template < typename T, int ndims >
class Triangle
{
public:

  /*!
   *****************************************************************************
   * \brief Constructor -- creates a triangle from the 3 points A,B,C.
   *****************************************************************************
   */
  Triangle( const Point<T,ndims>& A,
            const Point<T,ndims>& B,
            const Point<T,ndims>& C )
    {
      m_A = A;
      m_B = B;
      m_C = C;
    }


  /*!
   *****************************************************************************
   * \brief Destructor
   *****************************************************************************
   */
  virtual ~Triangle() { }

  /*!
   *****************************************************************************
   * \brief
   * \return A
   *****************************************************************************
   */
  const Point< T,ndims >& A() const { return m_A; };

  /*!
   *****************************************************************************
   * \brief
   * \return B
   *****************************************************************************
   */
  const Point< T,ndims >& B() const { return m_B; };

  /*!
   *****************************************************************************
   * \brief
   * \return C
   *****************************************************************************
   */
  const Point< T,ndims >& C() const { return m_C; };

private:

  /*!
   *****************************************************************************
   * \brief Default constructor. Does nothing.
   * \note Made private to prevent its use.
   *****************************************************************************
   */
  Triangle() { }

  Point< T, ndims > m_A;
  Point< T, ndims > m_B;
  Point< T, ndims > m_C;
};

} /* namespace quest */

#endif /* TRIANGLE_HPP_ */
