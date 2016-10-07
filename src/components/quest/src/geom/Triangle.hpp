/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_

#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include <cmath> // for acos()

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

    enum {
        NUM_TRI_VERTS = 3
    };

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
            const PointType& C );


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
  PointType& A( ) { return m_points[0]; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex A of the triangle.
   * \return A const reference to vertex A of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& A() const { return m_points[0]; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex B of the triangle.
   * \return B reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& B( ) { return m_points[1]; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex B of the triangle.
   * \return B const reference to vertex B of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& B() const { return m_points[1]; };

  /*!
   *****************************************************************************
   * \brief Returns a reference to vertex C of the triangle.
   * \return C reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  PointType& C( ) { return m_points[2]; };

  /*!
   *****************************************************************************
   * \brief Returns a const reference to vertex C of the triangle.
   * \return C const reference to vertex C of the triangle.
   * \see quest::Point
   *****************************************************************************
   */
  const PointType& C() const { return m_points[2]; };

  /**
   *****************************************************************************
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   *****************************************************************************
   */
  PointType& operator[](int idx)
  {
      SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);
      return m_points[ idx ];
  }

  /**
   *****************************************************************************
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   *****************************************************************************
   */
  const PointType& operator[](int idx) const
  {
      SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);
      return m_points[ idx ];
  }

  /*!
   *****************************************************************************
   * \brief Returns the normal of the triangle (not normalized)
   * \pre This function is only valid when DIM = 3
   * \return The normal vector to the triangle, zero vector otherwise
   *****************************************************************************
   */
  VectorType normal() const
  {
      SLIC_CHECK_MSG(DIM==3, "Triangle::normal() is only valid in 3D.");

      return (DIM==3)
              ? VectorType::cross_product( VectorType(m_points[0],m_points[1]),
                                           VectorType(m_points[0],m_points[2]))
              : VectorType();
  }

  /**
   *****************************************************************************
   * \brief Returns the area of the triangle
   * \pre Only defined when dimension DIM is 2 or 3
   *****************************************************************************
   */
  double area() const
  {
      SLIC_CHECK_MSG( DIM == 2 || DIM == 3,
            "Triangle::area() is only valid in 2D or 3D");

      VectorType v(m_points[0], m_points[1]);
      VectorType w(m_points[0], m_points[2]);

      return (DIM==2)
                  ? 0.5 * std::fabs(v[0]*w[1] - v[1]*w[0])
                  : 0.5 * VectorType::cross_product( v,w).norm();
  }

  /*!
   *****************************************************************************
   * \brief Computes the request angle corresponding to the given vertex ID.
   * \param [in] idx the index of the corresponding vertex
   * \return alpha the incidence angle.
   * \pre idx >= 0 && idx < NUM_TRI_VERTS
   *****************************************************************************
   */
  double angle( int idx ) const;

  /*!
   *****************************************************************************
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const
  {
      os <<"{"
         << m_points[0] <<" "
         << m_points[1] <<" "
         << m_points[2] <<"}";

      return os;
  }


private:

  PointType m_points[3];
};

} /* namespace quest */

//------------------------------------------------------------------------------
//  Triangle implementation
//------------------------------------------------------------------------------
namespace quest {

template <typename T, int DIM>
Triangle< T,DIM >::Triangle( const PointType& A,
                             const PointType& B,
                             const PointType& C  )
{
  m_points[0] = A;
  m_points[1] = B;
  m_points[2] = C;
}

//------------------------------------------------------------------------------
template <typename T, int DIM>
inline double Triangle< T,DIM >::angle( int idx ) const
{
  SLIC_ASSERT( idx >= 0 && idx < NUM_TRI_VERTS );

  const int idx1 = (idx+1)%NUM_TRI_VERTS;
  const int idx2 = (idx+2)%NUM_TRI_VERTS;

  PointType pt = m_points[idx];
  VectorType V1( pt, m_points[idx1] );
  VectorType V2( pt, m_points[idx2] );
  V1 /= V1.norm();
  V2 /= V2.norm();

  double dotprod = VectorType::dot_product( V1, V2 );
  return ( acos(dotprod) );
}

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
