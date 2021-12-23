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

  static constexpr int NUM_TET_VERTS = 4;

public:
  /*!
   * \brief Default constructor. Creates a degenerate tetrahedron.
   */
  AXOM_HOST_DEVICE Tetrahedron() { }

  /*!
   * \brief Custom Constructor. Creates a tetrahedron from the 4 points A,B,C,D.
   * \param [in] A point instance corresponding to vertex A of the tetrahedron.
   * \param [in] B point instance corresponding to vertex B of the tetrahedron.
   * \param [in] C point instance corresponding to vertex C of the tetrahedron.
   * \param [in] D point instance corresponding to vertex D of the tetrahedron.
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
  bool degenerate(double eps = 1.0e-12) const
  {
    return axom::utilities::isNearlyEqual(ppedVolume(), 0.0, eps);
  }

  /*!
   * \brief Returns the barycentric coordinates of a point within a tetrahedron
   * \return The barycentric coordinates of the tetrahedron
   * \post The barycentric coordinates sum to 1.
   */
  Point<double, 4> physToBarycentric(const PointType& p) const
  {
    using axom::numerics::determinant;

    Point<double, 4> bary;

    const PointType& p0 = m_points[0];
    const PointType& p1 = m_points[1];
    const PointType& p2 = m_points[2];
    const PointType& p3 = m_points[3];

    double det0 = ppedVolume();

    SLIC_CHECK_MSG(
      !axom::utilities::isNearlyEqual(det0, 0.),
      "Attempting to find barycentric coordinates of degenerate tetrahedron");

    const double detScale = 1. / det0;

    // clang-format off
    const double det1 = determinant(1.0,  p[0],  p[1],  p[2],
                                    1.0, p1[0], p1[1], p1[2],
                                    1.0, p2[0], p2[1], p2[2],
                                    1.0, p3[0], p3[1], p3[2]);

    const double det2 = determinant(1.0, p0[0], p0[1], p0[2],
                                    1.0,  p[0],  p[1],  p[2],
                                    1.0, p2[0], p2[1], p2[2],
                                    1.0, p3[0], p3[1], p3[2]);

    const double det3 = determinant(1.0, p0[0], p0[1], p0[2],
                                    1.0, p1[0], p1[1], p1[2],
                                    1.0,  p[0],  p[1],  p[2],
                                    1.0, p3[0], p3[1], p3[2]);

    const double det4 = determinant(1.0, p0[0], p0[1], p0[2],
                                    1.0, p1[0], p1[1], p1[2],
                                    1.0, p2[0], p2[1], p2[2],
                                    1.0,  p[0],  p[1],  p[2]);
    // clang-format on

    bary[0] = det1 * detScale;
    bary[1] = det2 * detScale;
    bary[2] = det3 * detScale;
    bary[3] = det4 * detScale;

    SLIC_CHECK_MSG(
      axom::utilities::isNearlyEqual(bary[0] + bary[1] + bary[2] + bary[3], 1.),
      "Barycentric coordinates should sum to 1. rather than "
        << (bary[0] + bary[1] + bary[2] + bary[3]));

    return bary;
  }

  /*!
   * \brief Returns the physical coordinates of a barycentric point
   * \param [in] bary Barycentric coordinates relative to this tetrahedron
   * \return Physical point represented by bary
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
  double signedVolume() const
  {
    const double scale = 1. / 6.;
    return scale * ppedVolume();
  }

  /*!
   * \brief Returns the absolute (unsigned) volume of the tetrahedron
   * \sa signedVolume()
   */
  double volume() const { return axom::utilities::abs(signedVolume()); }

  /**
   * \brief Returns the circumsphere of the tetrahedron
   *
   * Implements formula from https://mathworld.wolfram.com/Circumsphere.html
   * \note This function is only available for 3D tetrahedra in 3D
   */
  template <int TDIM = NDIMS>
  typename std::enable_if<TDIM == 3, primal::Sphere<T, 3>>::type circumsphere() const
  {
    using axom::numerics::determinant;
    using axom::numerics::dot_product;
    using axom::utilities::abs;
    using NumericArrayType = primal::NumericArray<T, NDIMS>;
    using SphereType = primal::Sphere<T, NDIMS>;

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
    if(NDIMS != 3)
    {
      return 0.;
    }
    else
    {
      const PointType& p0 = m_points[0];
      const PointType& p1 = m_points[1];
      const PointType& p2 = m_points[2];
      const PointType& p3 = m_points[3];

      // clang-format off
      return axom::numerics::determinant<double>(
        1.0, p0[0], p0[1], p0[2],
        1.0, p1[0], p1[1], p1[2],
        1.0, p2[0], p2[1], p2[2],
        1.0, p3[0], p3[1], p3[2]);
      // clang-format on
    }
  }

private:
  PointType m_points[4];
};

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
