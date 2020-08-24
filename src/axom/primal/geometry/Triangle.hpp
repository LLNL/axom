// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef TRIANGLE_HPP_
#define TRIANGLE_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/numerics/Determinants.hpp"
#include "axom/core/utilities/Utilities.hpp"

#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"


#include <cmath>   // for acos()
#include <ostream> // for std::ostream

namespace axom
{
namespace primal
{

// Forward declare the templated classes and operator functions
template < typename T, int NDIMS >
class Triangle;

/**
 * \brief Overloaded output operator for triangles
 */
template < typename T,int NDIMS >
std::ostream& operator<<(std::ostream & os, const Triangle< T,NDIMS > & tri);

/*!
 * \class Triangle
 *
 * \brief Represents a triangular geometric shape defined by three points.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template < typename T,int NDIMS >
class Triangle
{
public:
  typedef Point< T,NDIMS >  PointType;
  typedef Vector< T,NDIMS > VectorType;

  enum
  {
    NUM_TRI_VERTS = 3
  };

public:

  /*!
   * \brief Default constructor. Creates a degenerate triangle.
   */
  AXOM_HOST_DEVICE
  Triangle();

  /*!
   * \brief Custom Constructor. Creates a triangle from the 3 points A,B,C.
   * \param [in] A point instance corresponding to vertex A of the triangle.
   * \param [in] B point instance corresponding to vertex B of the triangle.
   * \param [in] C point instance corresponding to vertex C of the triangle.
   */
  AXOM_HOST_DEVICE
  Triangle( const PointType& A,
            const PointType& B,
            const PointType& C );

  /*!
   * \brief Destructor
   */
  AXOM_HOST_DEVICE
  ~Triangle() { }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  AXOM_HOST_DEVICE
  PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);
    return m_points[ idx ];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  AXOM_HOST_DEVICE
  const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >=0 && idx < NUM_TRI_VERTS);
    return m_points[ idx ];
  }

  /*!
   * \brief Returns the normal of the triangle (not normalized)
   * \pre This function is only valid when NDIMS = 3
   * \return n triangle normal when NDIMS=3, zero vector otherwise
   */
  AXOM_HOST_DEVICE
  VectorType normal() const
  {
    SLIC_CHECK_MSG(NDIMS==3, "Triangle::normal() is only valid in 3D.");

    return (NDIMS==3)
           ? VectorType::cross_product( VectorType(m_points[0],m_points[1]),
                                        VectorType(m_points[0],m_points[2]))
           : VectorType();
  }

  /*!
   * \brief Returns the area of the triangle
   * \pre Only defined when dimension NDIMS is 2 or 3
   */
  AXOM_HOST_DEVICE
  double area() const
  {
    SLIC_CHECK_MSG( NDIMS == 2 || NDIMS == 3,
                    "Triangle::area() is only valid in 2D or 3D");

    VectorType v(m_points[0], m_points[1]);
    VectorType w(m_points[0], m_points[2]);

    // While this code is correct and clear, in the 2D case it may be doing
    // too much work by taking the square root (norm()) of a just-computed
    // square (cross_product()).  This should be revisited if this turns out
    // to be a bottleneck.
    return 0.5 * VectorType::cross_product(v, w).norm();
  }

private:
  /*!
   * \brief Return the volume of the parallelepiped defined by this triangle
   *  instance and the given point.
   * \param [in] p user-supplied point.
   * \note If the volume is approx. zero, all four points are (nearly) coplanar.
   * \return vol the volume or 0.0 iff NDIMS < 3
   */
  double ppedVolume(const PointType& p) const
  {
    /* This method returns double (instead of T) and explicitly specializes
       determinant() on type double to avoid confusion of deduced template
       types. */
    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    if (NDIMS < 3)
    {
      return 0.;
    }
    else
    {
      return numerics::determinant< double > ( A[0], A[1], A[2], 1.,
                                               B[0], B[1], B[2], 1.,
                                               C[0], C[1], C[2], 1.,
                                               p[0], p[1], p[2], 1.  );
    }
  }

public:

  /*!
   * \brief Returns the barycentric coordinates of a point within a triangle
   * \return The barycentric coordinates of the triangle inside a Point<T,3>
   * \pre The point lies in this triangle's plane.
   * \post The barycentric coordinates sum to 1.
   * Adapted from Real Time Collision Detection by Christer Ericson.
   */
  Point< double, 3 > physToBarycentric(const PointType& p) const
  {
    SLIC_CHECK(axom::utilities::isNearlyEqual(ppedVolume(p), 0.));

    Point< double, 3 > bary;

    Vector< double, 3 > u =
      VectorType::cross_product( VectorType(m_points[0],m_points[1]),
                                 VectorType(m_points[0],m_points[2]) );
    const double x = std::abs(u[0]);
    const double y = std::abs(u[1]);
    const double z = std::abs(u[2]);

    double ood = 1.0 / u[2];      // compute in xy plane by default
    int c0 = 0;
    int c1 = 1;

    if (x>=y && x>= z)
    {
      // compute in yz plane
      c0 = 1;
      c1 = 2;
      ood=1.0/u[0];
    }
    else if (y>=x && y>=z)
    {
      // compute in xz plane
      c0 = 0;
      c1 = 2;
      ood=-1.0/u[1];
    }

    // References to triangle vertices for convenience
    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    // Compute ood * area of each sub-triangle
    bary[0] = ood * numerics::determinant( p[c0] - B[c0], p[c1] - B[c1],
                                           B[c0] - C[c0], B[c1] - C[c1]);
    bary[1] = ood * numerics::determinant( p[c0] - C[c0], p[c1] - C[c1],
                                           C[c0] - A[c0], C[c1] - A[c1]);
    bary[2] = 1. - bary[0] - bary[1];

    return bary;
  }

  /*!
   * \brief Returns the physical coordinates of a barycentric point
   * \param [in] bary Barycentric coordinates relative to this triangle
   * \return Physical point represented by bary
   */
  PointType baryToPhysical(const Point<double, 3> & bary) const
  {
    SLIC_CHECK_MSG( axom::utilities::isNearlyEqual(1., bary[0]+bary[1]+bary[2]),
                    "Barycentric coordinates must sum to (near) one." );

    PointType res;
    for (int i = 0 ; i < NDIMS ; ++i)
    {
      res[i] = bary[0] * m_points[0][i] +
               bary[1] * m_points[1][i] +
               bary[2] * m_points[2][i];
    }

    return res;
  }

  /*!
   * \brief Returns whether the triangle is degenerate
   * \return true iff the triangle is degenerate (0 area)
   * \see primal::Point
   */
  AXOM_HOST_DEVICE
  bool degenerate(double eps = 1.0e-12) const
  {
    return axom::utilities::isNearlyEqual(area(),  0.0, eps);
  }

  /*!
   * \brief Returns whether Point P is in the triangle for some 3d Triangle
   * \return true iff P is in the triangle
   * \see primal::Point
   */
  bool checkInTriangle(const PointType& p, double eps = 1.0e-8) const
  {
    if (!axom::utilities::isNearlyEqual(ppedVolume(p), 0., eps))
    {
      return false;
    }

    Point< double,3 > bC= physToBarycentric(p);
    return ((bC[0]>=(0.0-eps)) && (bC[1] >= (0.0-eps)) && (bC[2]>=(0.0-eps)) &&
            (bC[0]<=(1.0+eps)) && (bC[1] <= (1.0+eps)) && (bC[2]<=(1.0+eps)));
  }

  /*!
   * \brief Computes the request angle corresponding to the given vertex ID.
   * \param [in] idx the index of the corresponding vertex
   * \return alpha the incidence angle in the range [0, pi].
   * \pre idx >= 0 && idx < NUM_TRI_VERTS
   */
  double angle( int idx ) const;

  /*!
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
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
  PointType m_points[ 3 ];

};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
//  Triangle implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{

template < typename T, int NDIMS >
Triangle< T,NDIMS >::Triangle()
   : m_points{ PointType(), PointType(), PointType() }
{}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Triangle< T,NDIMS >::Triangle( const PointType& A,
                               const PointType& B,
                               const PointType& C)
   : m_points{A,B,C} 
{}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline double Triangle< T,NDIMS >::angle( int idx ) const
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
template < typename T, int NDIMS >
std::ostream& operator<<(std::ostream & os, const Triangle< T,NDIMS > & tri)
{
  tri.print(os);
  return os;
}

} /* namespace primal */

} /* namespace axom */

#endif /* TRIANGLE_HPP_ */
