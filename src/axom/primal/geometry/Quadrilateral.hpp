// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_QUADRILATERAL_HPP_
#define AXOM_PRIMAL_QUADRILATERAL_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/fmt.hpp"

#include "axom/primal/constants.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Triangle.hpp"

#include <cmath>
#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Quadrilateral;

/**
 * \brief Overloaded output operator for quadrilaterals
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Quadrilateral<T, NDIMS>& quad);

/*!
 * \accelerated
 * \class Quadrilateral
 *
 * \brief Represents a quadrilateral geometric shape defined by four points.
 *
 * \accelerated
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * There are four vertices in the quadrilateral, labeled P through S as the
 * constructor's arguments.  They are accessible using the square-brackets
 * operator, with P being index 0, Q index 1, R index 2, and S index 3.
 *
 * Here's a diagram showing a square with the labeled vertices.
 *
 * \verbatim
 *
 * Q +---------+ R       +y
 *   |         |         
 *   |         |          ^
 *   |         |          |
 *   |         |          |
 * P +---------+ S        -----> +x
 *
 * \endverbatim
 */
template <typename T, int NDIMS>
class Quadrilateral
{
public:
  using PointType = Point<T, NDIMS>;
  using TriangleType = Triangle<T, NDIMS>;

  static constexpr int DIM = NDIMS;
  static constexpr int NUM_QUAD_VERTS = 4;

public:
  /// \brief Default constructor. Creates a degenerate quadrilateral.
  AXOM_HOST_DEVICE
  Quadrilateral() = default;

  /*!
   * \brief Custom Constructor. Creates a quadrilateral from the 4 points p, q, r, s.
   * \param [in] p point instance corresponding to vertex P of the quadrilateral.
   * \param [in] q point instance corresponding to vertex Q of the quadrilateral.
   * \param [in] r point instance corresponding to vertex R of the quadrilateral.
   * \param [in] s point instance corresponding to vertex S of the quadrilateral.
   */
  AXOM_HOST_DEVICE
  Quadrilateral(const PointType& p,
                const PointType& q,
                const PointType& r,
                const PointType& s)
    : m_points {p, q, r, s}
  { }

  /*!
   * \brief Quadrilateral constructor from an array of Points
   *
   * \param [in] pts An array containing at least 4 Points.
   *
   * \note It is the responsiblity of the caller to pass
   *       an array with at least 4 Points
   */
  AXOM_HOST_DEVICE
  explicit Quadrilateral(const PointType* pts)
  {
    for (int i = 0; i < NUM_QUAD_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Quadrilateral constructor from an Array of Points.
   *
   * \param [in] pts An ArrayView containing 4 Points.
   */
  AXOM_HOST_DEVICE
  explicit Quadrilateral(const axom::ArrayView<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_QUAD_VERTS);

    for(int i = 0; i < NUM_QUAD_VERTS; i++)
    {
      m_points[i] = pts[i];
    }
  }

  /*!
   * \brief Quadrilateral constructor from an initializer list of Points
   *
   * \param [in] pts an initializer list containing 4 Points
   */
  AXOM_HOST_DEVICE
  explicit Quadrilateral(std::initializer_list<PointType> pts)
  {
    SLIC_ASSERT(pts.size() == NUM_QUAD_VERTS);

    int i = 0;
    for(const auto& pt : pts)
    {
      m_points[i] = pt;
      i++;
    }
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1 or 2
   */
  AXOM_HOST_DEVICE
  PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_QUAD_VERTS);
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
    SLIC_ASSERT(idx >= 0 && idx < NUM_QUAD_VERTS);
    return m_points[idx];
  }

  /// \brief Returns the area of the quadrilateral (2D specialization)
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type area() const
  {
    return axom::utilities::abs(signedArea());
  }

  /**
   * \brief Returns the signed area of a 2D quadrilateral
   *
   * The area is positive when the vertices are oriented counter-clockwise.
   * \note This function is only available for quadrilaterals in 2D since
   * signed areas don't make sense for 3D quadrilaterals
   *
   * Compute the signed area by dividing the quadrilateral into two triangles
   * as shown below and summing their signed area
   *
   * \verbatim
   *
   * q +----+ r       +y
   *   |   /|         
   *   |  / |          ^
   *   | /  |          |
   *   |/   |          |
   * p +----+ s        -----> +x
   *
   * \endverbatim
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedArea() const
  {
    const PointType& p = m_points[0];
    const PointType& q = m_points[1];
    const PointType& r = m_points[2];
    const PointType& s = m_points[3];

    // clang-format off
    return TriangleType{p, q, r}.signedArea() +
           TriangleType{r, s, p}.signedArea();
    // clang-format on
  }

  /*!
   * \brief Returns the volume of the quadrilateral (synonym for area())
   * \sa area()
   */
  AXOM_HOST_DEVICE double volume() const { return area(); }

  /*!
   * \brief Returns the signed volume of a 2D quadrilateral (synonym for signedArea())
   * \sa signedArea()
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedVolume() const
  {
    return signedArea();
  }

  /*!
   * \brief Simple formatted print of a triangle instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " " << m_points[3] << "}";

    return os;
  }

private:
  PointType m_points[NUM_QUAD_VERTS];
};

//------------------------------------------------------------------------------
/// Free functions implementing Quadrilateral's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Quadrilateral<T, NDIMS>& quad)
{
  quad.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Quadrilateral using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Quadrilateral<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_QUADRILATERAL_HPP_
