// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Polygon.hpp
 *
 * \brief A Polygon primitive for primal
 */

#ifndef AXOM_PRIMAL_POLYGON_HPP_
#define AXOM_PRIMAL_POLYGON_HPP_

#include "axom/core/Array.hpp"
#include "axom/primal/geometry/Point.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Polygon;

/// \brief Overloaded output operator for polygons
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polygon<T, NDIMS>& poly);

/*!
 * \class Polygon
 *
 * \brief Represents a polygon defined by an array of points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 * \note The polygon vertices should be ordered in a counter clockwise
 *       orientation with respect to the polygon's desired normal vector
 */
template <typename T, int NDIMS>
class Polygon
{
public:
  using PointType = Point<T, NDIMS>;

public:
  /// Default constructor for an empty polygon
  Polygon() { }

  /*!
   * \brief Constructor for an empty polygon that reserves space for
   *  the given number of vertices
   *
   * \param [in] numExpectedVerts number of vertices for which to reserve space
   * \pre numExpectedVerts is not negative
   *
   */
  explicit Polygon(int numExpectedVerts)
  {
    SLIC_ASSERT(numExpectedVerts >= 0);
    m_vertices.reserve(numExpectedVerts);
  }

  /// \brief Constructor for a polygon with the given vertices
  Polygon(const axom::Array<PointType>& vertices) { m_vertices = vertices; }

  /// Return the number of vertices in the polygon
  int numVertices() const { return static_cast<int>(m_vertices.size()); }

  /// Appends a vertex to the list of vertices
  void addVertex(const PointType& pt) { m_vertices.push_back(pt); }

  /// Clears the list of vertices
  void clear() { m_vertices.clear(); }

  /// Retrieves the vertex at index idx
  PointType& operator[](int idx) { return m_vertices[idx]; }
  /// Retrieves the vertex at index idx
  const PointType& operator[](int idx) const { return m_vertices[idx]; }

  /*!
   * \brief Computes the average of the polygon's vertex positions
   *
   * \return A point at the mean of the polygon's vertices
   * \pre  polygon.isValid() is true
   */
  PointType vertexMean() const
  {
    SLIC_ASSERT(isValid());

    PointType sum;

    const int sz = numVertices();
    for(int i = 0; i < sz; ++i)
    {
      sum.array() += m_vertices[i].array();
    }
    sum.array() /= sz;

    return sum;
  }

  /*!
   * \brief Simple formatted print of a polygon instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numVertices();

    os << "{" << sz << "-gon:";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << m_vertices[i] << ",";
    }
    if(sz >= 2)
    {
      os << m_vertices[sz - 1];
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Simple check for validity of a polygon
   *
   * Initial check is that the polygon has three or more vertices
   * \return True, if the polygon is valid, False otherwise
   */
  bool isValid() const { return m_vertices.size() >= 3; }

private:
  axom::Array<PointType> m_vertices;
};

//------------------------------------------------------------------------------
/// Free functions implementing Polygon's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polygon<T, NDIMS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_POLYGON_HPP_
