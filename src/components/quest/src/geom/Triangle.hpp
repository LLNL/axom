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
#include "common/Utilities.hpp"

#include "slic/slic.hpp"

#include <cmath> // for acos()
#include <ostream>

namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Triangle;

/**
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

  /**
   *****************************************************************************
   * \brief Default constructor. Creates a degenerate triangle.
   *****************************************************************************
   */
  Triangle() { }

  /**
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


  /**
   *****************************************************************************
   * \brief Destructor
   *****************************************************************************
   */
   ~Triangle() { }

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

  /**
   *****************************************************************************
   * \brief Returns the normal of the triangle (not normalized)
   * \pre This function is only valid when DIM = 3
   * \return n triangle normal when DIM=3, zero vector otherwise
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

      return 0.5 * VectorType::cross_product(v, w).norm();
  }

  /**
   * \brief Returns the barycentric coordinates of a point within a triangle
   * \return The barycentric coordinates of the triangle inside a Point<T,3>
   *
   * \post The barycentric coordinates sum to 1.
   *
   * Adapted from Real Time Collision Detection by Christer Ericson.
   */
  Point<T,3> computeBarycenterCoords(const PointType& p) const
  {
    Point<T,3> bary;

    VectorType u= VectorType::cross_product(VectorType(m_points[0],m_points[1]),
                                            VectorType(m_points[0],m_points[2]));
    const T x= std::abs(u[0]);
    const T y= std::abs(u[1]);
    const T z= std::abs(u[2]);

    T ood = 1.0 / u[2];      // compute in xy plane by default 
    int c0 = 0;
    int c1 = 1;

    if (x>=y && x>= z)       // compute in yz plane
    {
        c0 = 1;
        c1 = 2;
        ood=1.0/u[0];

    }
    else if (y>=x && y>=z)  // compute in xz plane
    {
        c0 = 0;
        c1 = 2;
        ood=-1.0/u[1];

    }

    // References to triangle vertices for convenience
    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    // Compute ood * area of each sub-triangle
    bary[0] = ood * math::determinant(p[c0] - B[c0], p[c1] - B[c1], 
                                      B[c0] - C[c0], B[c1] - C[c1]);
    bary[1] = ood * math::determinant(p[c0] - C[c0], p[c1] - C[c1],
                                      C[c0] - A[c0], C[c1] - A[c1]);
    bary[2] = 1. - bary[0] - bary[1];

    return bary;
  }

  /**
   *****************************************************************************
   * \brief Returns whether the triangle is degenerate
   * \return true iff the triangle is degenerate (0 area)
   * \see quest::Point
   *****************************************************************************
   */
  bool degenerate() const
  {
    return asctoolkit::utilities::isNearlyEqual(area(),  0.0, 1.0e-12);
  }

  /**
   *****************************************************************************
   * \brief Returns whether Point P is in the triangle for some 3d Triangle
   * \return true iff P is in the triangle
   * \see quest::Point
   *****************************************************************************
   */
  bool checkInTriangle(const Point<double, DIM>& P) const{
    Point<T,3> bC= barycenterCoords(P);
    return ((bC[0]>=0.0) && (bC[1] >= 0.0) && (bC[2]>=0.0) &&
	    (bC[0]<=1.0) && (bC[1]<=1.0) && (bC[2]<=1.0) &&
	    (bC[0]+bC[1]+bC[2]<=1.0));
  }


  /**
   *****************************************************************************
   * \brief Computes the request angle corresponding to the given vertex ID.
   * \param [in] idx the index of the corresponding vertex
   * \return alpha the incidence angle.
   * \pre idx >= 0 && idx < NUM_TRI_VERTS
   *****************************************************************************
   */
  double angle( int idx ) const;

  /**
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

  const PointType& pt = m_points[idx];
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
