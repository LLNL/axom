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
#include "quest/Vector.hpp"

namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Triangle;

/*!
 * \brief Overloaded output operator for triangles
 */
template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Triangle<T,DIM> & tri);



template < typename T, int DIM >
class Triangle
{
public:
    typedef Point<T,DIM>  PointType;
    typedef Vector<T,DIM> VectorType;

public:

  /*!
   *****************************************************************************
   * \brief Default constructor. Creates a degenerate triangle.
   *****************************************************************************
   */
  Triangle() { }

  /*!
   *****************************************************************************
   * \brief Custom Constructor. Creates a triangle from the 3 points A,B,C.
   * \param [in] A point instance corresponding to vertex A of the triangle.
   * \param [in] B point instance corresponding to vertex B of the triangle.
   * \param [in] C point instance corresponding to vertex C of the triangle.
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
   * \brief Returns a reference to vertex A of the triangle.
   * \return A reference to vertex A of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& A( ) { return m_A; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex A of the triangle.
   * \return A const reference to vertex A of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& A() const { return m_A; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex B of the triangle.
   * \return B reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& B( ) { return m_B; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex B of the triangle.
   * \return B const reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& B() const { return m_B; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex C of the triangle.
   * \return C reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& C( ) { return m_C; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex C of the triangle.
   * \return C const reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& C() const { return m_C; };


  /*!
   * \brief Returns the normal of the triangle (not normalized)
   */
  VectorType normal() const
  {
      return VectorType::cross_product( VectorType(m_A,m_B), VectorType(m_A,m_C));
  }


  /*!
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
      os <<"{"
         << m_A <<" "
         << m_B <<" "
         << m_C <<"}";

      return os;
  }


private:

  PointType m_A;
  PointType m_B;
  PointType m_C;
};

//------------------------------------------------------------------------------
/// Free functions implementing Triangle's operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Triangle<T,DIM> & tri)
{
    tri.print(os);
    return os;
}


} /* namespace quest */

#endif /* TRIANGLE_HPP_ */
