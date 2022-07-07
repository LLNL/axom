// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_TRIANGLE_HPP_
#define AXOM_PRIMAL_TRIANGLE_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Sphere.hpp"

#include <cmath>
#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Triangle;

/**
 * \brief Overloaded output operator for triangles
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Triangle<T, NDIMS>& tri);

/*!
 * \accelerated
 * \class Triangle
 *
 * \brief Represents a triangular geometric shape defined by three points.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Triangle
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using SphereType = Sphere<T, NDIMS>;

  static constexpr int DIM = NDIMS;
  static constexpr int NUM_TRI_VERTS = 3;

public:
  /// \brief Default constructor. Creates a degenerate triangle.
  AXOM_HOST_DEVICE
  Triangle() : m_points {PointType(), PointType(), PointType()} { }

  /*!
   * \brief Custom Constructor. Creates a triangle from the 3 points A,B,C.
   * \param [in] A point instance corresponding to vertex A of the triangle.
   * \param [in] B point instance corresponding to vertex B of the triangle.
   * \param [in] C point instance corresponding to vertex C of the triangle.
   */
  AXOM_HOST_DEVICE
  Triangle(const PointType& A, const PointType& B, const PointType& C)
    : m_points {A, B, C}
  { }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  AXOM_HOST_DEVICE
  PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_TRI_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  AXOM_HOST_DEVICE
  const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_TRI_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Returns the normal of the triangle (not normalized)
   * \pre This function is only valid when NDIMS = 3
   * \return n triangle normal when NDIMS=3, zero vector otherwise
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, VectorType>::type normal() const
  {
    return VectorType::cross_product(m_points[1] - m_points[0],
                                     m_points[2] - m_points[0]);
  }

  /// \brief Returns the area of the triangle (3D specialization)
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, double>::type area() const
  {
    return 0.5 * normal().norm();
  }

  /// \brief Returns the area of the triangle (2D specialization)
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type area() const
  {
    return axom::utilities::abs(signedArea());
  }

  /**
   * \brief Returns the signed area of a 2D triangle
   *
   * The area is positive when the vertices are oriented counter-clockwise.
   * \note This function is only available for triangles in 2D since signed
   * areas don't make sense for 3D triangles
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedArea() const
  {
    using axom::numerics::determinant;
    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    // clang-format off
    return 0.5 * determinant(B[0]-A[0], C[0]-A[0],
                             B[1]-A[1], C[1]-A[1]);
    // clang-format on
  }

  /*!
   * \brief Returns the volume of the triangle (synonym for area())
   * \sa area()
   */
  AXOM_HOST_DEVICE double volume() const { return area(); }

  /// \brief Returns the signed volume of a 2D triangle
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedVolume() const
  {
    return signedArea();
  }

  /**
   * \brief Returns the circumsphere (circumscribing circle) of the triangle
   *
   * Implements formula from https://mathworld.wolfram.com/Circumcircle.html
   * \note This function is only available for triangles in 2D
   */
  template <int TDIM = NDIMS>
  typename std::enable_if<TDIM == 2, SphereType>::type circumsphere() const
  {
    using axom::numerics::determinant;
    using axom::numerics::dot_product;
    using axom::utilities::abs;
    using NumericArrayType = primal::NumericArray<T, NDIMS>;

    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    // clang-format off
    const Point<T, 4> sq { dot_product(A.data(), A.data(), NDIMS),
                           dot_product(B.data(), B.data(), NDIMS),
                           dot_product(C.data(), C.data(), NDIMS)};

    const double  a =  determinant(A[0], A[1], 1.,
                                   B[0], B[1], 1.,
                                   C[0], C[1], 1.);

    const double bx = -determinant(sq[0], A[1], 1.,
                                   sq[1], B[1], 1.,
                                   sq[2], C[1], 1.);

    const double by =  determinant(sq[0], A[0], 1.,
                                   sq[1], B[0], 1.,
                                   sq[2], C[0], 1.);

    const double  c = -determinant(sq[0], A[0], A[1],
                                   sq[1], B[0], B[1],
                                   sq[2], C[0], C[1]);
    // clang-format on

    const auto center = NumericArrayType {-bx, -by} / (2 * a);
    const T radius = sqrt(bx * bx + by * by - 4 * a * c) / (2 * abs(a));
    return SphereType(center.data(), radius);
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
       determinant() on type double to avoid confusion of deduced template types. */
    if(NDIMS < 3)
    {
      return 0.;
    }
    else
    {
      const VectorType pA = m_points[0] - p;
      const VectorType pB = m_points[1] - p;
      const VectorType pC = m_points[2] - p;

      // clang-format off
      return numerics::determinant< double > (pA[0], pA[1], pA[2],
                                              pB[0], pB[1], pB[2],
                                              pC[0], pC[1], pC[2]);
      // clang-format on
    }
  }

public:
  /*!
   * \brief Returns the barycentric coordinates of a point within a triangle
   *
   * \return The barycentric coordinates of the triangle
   * \pre The point lies in this triangle's plane.
   * \post The barycentric coordinates sum to 1.
   * Adapted from Real Time Collision Detection by Christer Ericson.
   */
  Point<double, 3> physToBarycentric(const PointType& p) const
  {
    // Query point needs to be in triangle's plane
    SLIC_CHECK(axom::utilities::isNearlyEqual(ppedVolume(p), 0.));

    Point<double, 3> bary;

    auto triArea2D =
      [](double x1, double y1, double x2, double y2, double x3, double y3)
      -> double { return (x1 - x2) * (y2 - y3) - (x2 - x3) * (y1 - y2); };

    // References to triangle vertices for convenience
    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];

    // unnormalized triangle normal
    const auto u = VectorType::cross_product(B - A, C - A);

    double ood;     // one over denomenator
    double nu, nv;  // numerators for 2D projection

    // Find best projection plane for computing weights: xy, yz, xz
    const int projectionDim = DIM == 2 ? 2 : primal::abs(u.array()).argMax();
    switch(projectionDim)
    {
    case 0:  // compute in yz plane
      nu = triArea2D(p[1], p[2], B[1], B[2], C[1], C[2]);
      nv = triArea2D(p[1], p[2], C[1], C[2], A[1], A[2]);
      ood = 1.0 / u[0];
      break;
    case 1:  // compute in xz plane
      nu = triArea2D(p[0], p[2], B[0], B[2], C[0], C[2]);
      nv = triArea2D(p[0], p[2], C[0], C[2], A[0], A[2]);
      ood = -1.0 / u[1];
      break;
    case 2:
    default:  // compute in xy plane
      nu = triArea2D(p[0], p[1], B[0], B[1], C[0], C[1]);
      nv = triArea2D(p[0], p[1], C[0], C[1], A[0], A[1]);
      ood = 1.0 / u[2];
      break;
    }

    // Return barycentric coordinates: ood * area of each sub-triangle
    bary[0] = ood * nu;
    bary[1] = ood * nv;
    bary[2] = 1. - bary[0] - bary[1];

    return bary;
  }

  /*!
   * \brief Returns the physical coordinates of a barycentric point
   * \param [in] bary Barycentric coordinates relative to this triangle
   * \return Physical point represented by bary
   */
  PointType baryToPhysical(const Point<double, 3>& bary) const
  {
    SLIC_CHECK_MSG(
      axom::utilities::isNearlyEqual(1., bary[0] + bary[1] + bary[2]),
      "Barycentric coordinates must sum to (near) one.");

    PointType res;
    for(int i = 0; i < NUM_TRI_VERTS; ++i)
    {
      res.array() += bary[i] * m_points[i].array();
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
    return axom::utilities::isNearlyEqual(area(), 0.0, eps);
  }

  /*!
   * \brief Returns whether Point P is in the triangle for some 3d Triangle
   * \return true iff P is in the triangle
   * \see primal::Point
   */
  bool checkInTriangle(const PointType& p, double eps = 1.0e-8) const
  {
    if(!axom::utilities::isNearlyEqual(ppedVolume(p), 0., eps))
    {
      return false;
    }

    Point<double, 3> bC = physToBarycentric(p);
    return ((bC[0] >= (0.0 - eps)) && (bC[1] >= (0.0 - eps)) &&
            (bC[2] >= (0.0 - eps)) && (bC[0] <= (1.0 + eps)) &&
            (bC[1] <= (1.0 + eps)) && (bC[2] <= (1.0 + eps)));
  }

  /*!
   * \brief Computes the request angle corresponding to the given vertex ID.
   * \param [in] idx the index of the corresponding vertex
   * \return alpha the incidence angle in the range [0, pi].
   * \pre idx >= 0 && idx < NUM_TRI_VERTS
   */
  AXOM_HOST_DEVICE
  double angle(int idx) const;

  /*!
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << "}";

    return os;
  }

private:
  PointType m_points[3];
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
template <typename T, int NDIMS>
AXOM_HOST_DEVICE inline double Triangle<T, NDIMS>::angle(int idx) const
{
  SLIC_ASSERT(idx >= 0 && idx < NUM_TRI_VERTS);

  const int idx1 = (idx + 1) % NUM_TRI_VERTS;
  const int idx2 = (idx + 2) % NUM_TRI_VERTS;

  const PointType& pt = m_points[idx];
  VectorType V1 = (m_points[idx1] - pt).unitVector();
  VectorType V2 = (m_points[idx2] - pt).unitVector();

  double dotprod = VectorType::dot_product(V1, V2);

  // Account for floating point error in (rare) degenerate cases
  return acos(axom::utilities::clampVal(dotprod, -1.0, 1.0));
}

//------------------------------------------------------------------------------
/// Free functions implementing Triangle's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Triangle<T, NDIMS>& tri)
{
  tri.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_TRIANGLE_HPP_
