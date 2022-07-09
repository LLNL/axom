// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_TETRAHEDRON_HPP_
#define AXOM_PRIMAL_TETRAHEDRON_HPP_

#include "axom/core.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Sphere.hpp"

#include "axom/slic/interface/slic.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
/*!
 * \class Tetrahedron
 *
 * \brief Represents a tetrahedral geometric shape defined by four points.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 */
template <typename T, int NDIMS>
class Tetrahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using SphereType = Sphere<T, NDIMS>;

  static constexpr int NUM_TET_VERTS = 4;

public:
  /// \brief Default constructor. Creates a degenerate tetrahedron.
  AXOM_HOST_DEVICE Tetrahedron()
    : m_points {PointType(), PointType(), PointType(), PointType()}
  { }

  /*!
   * \brief Tetrahedron constructor from 4 points A,B,C,D.
   * \param [in] A Vertex A of the tetrahedron
   * \param [in] B Vertex B of the tetrahedron
   * \param [in] C Vertex C of the tetrahedron
   * \param [in] D Vertex D of the tetrahedron
   */
  AXOM_HOST_DEVICE
  Tetrahedron(const PointType& A,
              const PointType& B,
              const PointType& C,
              const PointType& D)
    : m_points {A, B, C, D}
  { }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, or 3
   */
  AXOM_HOST_DEVICE PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_TET_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, or 3
   */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_TET_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Returns whether the tetrahedron is degenerate
   * \return true iff the tetrahedron is degenerate (has near zero volume)
   */
  AXOM_HOST_DEVICE
  bool degenerate(double eps = 1.0e-12) const
  {
    return axom::utilities::isNearlyEqual(ppedVolume(), 0.0, eps);
  }

  /*!
   * \brief Returns the barycentric coordinates of a point within a tetrahedron
   *
   * \param [in] p The point at which we want to compute barycentric coordinates
   * \param [in] skipNormalization Determines if the result should be normalized
   * by the volume of the tetrahedron. One might want to skip this if they were
   * only interested in the relative weights rather than their actual amounts
   * 
   * \post The barycentric coordinates sum to 1 when \a skipNormalization is false
   * Otherwise, the sum of coordinates will be proportional to volume of the tetrahedron
   * (Specifically, they should sum to the parallelpiped volume, which is 6x the volume).
   */
  Point<double, 4> physToBarycentric(const PointType& p,
                                     bool skipNormalization = false) const
  {
    Point<double, 4> bary;

    const PointType& A = m_points[0];
    const PointType& B = m_points[1];
    const PointType& C = m_points[2];
    const PointType& D = m_points[3];

    const auto pA = A - p;
    const auto pB = B - p;
    const auto pC = C - p;
    const auto pD = D - p;

    const double detA = VectorType::scalar_triple_product(pB, pC, pD);
    const double detB = -VectorType::scalar_triple_product(pC, pD, pA);
    const double detC = VectorType::scalar_triple_product(pD, pA, pB);
    const double detD = -VectorType::scalar_triple_product(pA, pB, pC);

    if(!skipNormalization)
    {
      constexpr double EPS = 1.0e-50;
      const double vol = VectorType::scalar_triple_product(B - A, C - A, D - A);
      // Compute one over denomenator; add a tiny amount to avoid division by zero
      const double ood = 1. / (vol - EPS);

      bary[0] = detA * ood;
      bary[1] = detB * ood;
      bary[2] = detC * ood;
      bary[3] = detD * ood;

      // Replace the smallest entry with the difference of 1 from the sum of the others
      const int amin = primal::abs(bary.array()).argMin();
      bary[amin] = 0.;
      bary[amin] = 1. - bary.array().sum();
    }
    else
    {
      bary[0] = detA;
      bary[1] = detB;
      bary[2] = detC;
      bary[3] = detD;
    }

    return bary;
  }

  /*!
   * \brief Returns the physical coordinates of a barycentric point
   * \param [in] bary Barycentric coordinates relative to this tetrahedron
   */
  PointType baryToPhysical(const Point<double, 4>& bary) const
  {
    SLIC_CHECK_MSG(
      axom::utilities::isNearlyEqual(1., bary[0] + bary[1] + bary[2] + bary[3]),
      "Barycentric coordinates must sum to (near) one.");

    PointType res;
    for(int i = 0; i < NUM_TET_VERTS; ++i)
    {
      res.array() += bary[i] * m_points[i].array();
    }

    return res;
  }

  /*!
   * \brief Simple formatted print of a tetrahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " "
       << m_points[3] << "}";

    return os;
  }

  /*!
   * \brief Returns the signed volume of the tetrahedron
   * \sa volume()
   */
  AXOM_HOST_DEVICE
  double signedVolume() const
  {
    constexpr double scale = 1. / 6.;
    return scale * ppedVolume();
  }

  /*!
   * \brief Returns the absolute (unsigned) volume of the tetrahedron
   * \sa signedVolume()
   */
  AXOM_HOST_DEVICE
  double volume() const { return axom::utilities::abs(signedVolume()); }

  /**
   * \brief Returns the circumsphere of the tetrahedron
   *
   * Implements formula from https://mathworld.wolfram.com/Circumsphere.html
   * \note This function is only available for 3D tetrahedra in 3D
   */
  template <int TDIM = NDIMS>
  typename std::enable_if<TDIM == 3, SphereType>::type circumsphere() const
  {
    using axom::numerics::determinant;
    using axom::numerics::dot_product;
    using axom::utilities::abs;
    using NumericArrayType = primal::NumericArray<T, NDIMS>;

    const PointType& p0 = m_points[0];
    const PointType& p1 = m_points[1];
    const PointType& p2 = m_points[2];
    const PointType& p3 = m_points[3];

    // clang-format off
    const Point<T, 4> sq { dot_product(p0.data(), p0.data(), NDIMS),
                           dot_product(p1.data(), p1.data(), NDIMS),
                           dot_product(p2.data(), p2.data(), NDIMS),
                           dot_product(p3.data(), p3.data(), NDIMS)};

    const double  a =  determinant(p0[0], p0[1], p0[2], 1.,
                                   p1[0], p1[1], p1[2], 1.,
                                   p2[0], p2[1], p2[2], 1.,
                                   p3[0], p3[1], p3[2], 1.);

    const double dx =  determinant(sq[0], p0[1], p0[2], 1.,
                                   sq[1], p1[1], p1[2], 1.,
                                   sq[2], p2[1], p2[2], 1.,
                                   sq[3], p3[1], p3[2], 1.);

    const double dy = -determinant(sq[0], p0[0], p0[2], 1.,
                                   sq[1], p1[0], p1[2], 1.,
                                   sq[2], p2[0], p2[2], 1.,
                                   sq[3], p3[0], p3[2], 1.);

    const double dz =  determinant(sq[0], p0[0], p0[1], 1.,
                                   sq[1], p1[0], p1[1], 1.,
                                   sq[2], p2[0], p2[1], 1.,
                                   sq[3], p3[0], p3[1], 1.);

    const double  c =  determinant(sq[0], p0[0], p0[1], p0[2],
                                   sq[1], p1[0], p1[1], p1[2],
                                   sq[2], p2[0], p2[1], p2[2],
                                   sq[3], p3[0], p3[1], p3[2]);
    // clang-format on

    const auto center = NumericArrayType {dx, dy, dz} / (2 * a);
    const T radius = sqrt(dx * dx + dy * dy + dz * dz - 4 * a * c) / (2 * abs(a));
    return SphereType(center.data(), radius);
  }

private:
  /*!
   * \brief Computes the signed volume of a parallelepiped defined by the
   * three edges of the tetrahedron incident to its first vertex
   *
   * \note The ppedVolume is a factor of 6 greater than that of the tetrahedron
   * \return The signed parallelepiped volume
   * \sa signedVolume(), volume()
   */
  AXOM_HOST_DEVICE
  double ppedVolume() const
  {
    return NDIMS != 3
      ? 0.
      : VectorType::scalar_triple_product(m_points[1] - m_points[0],
                                          m_points[2] - m_points[0],
                                          m_points[3] - m_points[0]);
  }

private:
  PointType m_points[4];
};  // namespace primal

//------------------------------------------------------------------------------
/// Free functions implementing Tetrahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Tetrahedron<T, NDIMS>& tet)
{
  tet.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_TETRAHEDRON_HPP_
