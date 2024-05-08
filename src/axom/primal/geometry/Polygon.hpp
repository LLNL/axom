// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
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
#include "axom/primal/geometry/Vector.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
namespace
{
static constexpr int DEFAULT_MAX_NUM_VERTICES = 20;
} /* end anonymous namespace */

// Forward declare the templated classes and operator functions
template <typename T, int NDIMS, int MAX_VERTS>
class Polygon;

/// \brief Overloaded output operator for polygons
template <typename T, int NDIMS, int MAX_VERTS>
std::ostream& operator<<(std::ostream& os,
                         const Polygon<T, NDIMS, MAX_VERTS>& poly);

/*!
 * \class Polygon
 *
 * \brief Represents a polygon defined by an array of points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 * \tparam MAX_VERTS the max number of vertices to preallocate space for
 *         (default max is 20 vertices)
 * \note The polygon vertices should be ordered in a counter clockwise
 *       orientation with respect to the polygon's desired normal vector
 */
template <typename T, int NDIMS, int MAX_VERTS = DEFAULT_MAX_NUM_VERTICES>
class Polygon
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

public:
  /// Default constructor for an empty polygon
  AXOM_HOST_DEVICE
  Polygon() { }

  /// \brief Constructor for a polygon with the given vertices
  explicit Polygon(const axom::Array<PointType>& vertices)
  {
    SLIC_ASSERT(static_cast<int>(vertices.size()) <= MAX_VERTS);

    for(int i = 0; i < vertices.size(); i++)
    {
      m_vertices[i] = vertices[i];
    }
    m_num_vertices = vertices.size();
  }

  // \brief Constructor for a polygon with an initializer list of Points
  AXOM_HOST_DEVICE
  explicit Polygon(std::initializer_list<PointType> vertices)
  {
    SLIC_ASSERT(static_cast<int>(vertices.size()) <= MAX_VERTS);

    int i = 0;
    for(const auto& vertex : vertices)
    {
      m_vertices[i] = vertex;
      i++;
    }
    m_num_vertices = vertices.size();
  }

  /// Return the number of vertices in the polygon
  AXOM_HOST_DEVICE
  int numVertices() const { return m_num_vertices; }

  /// Appends a vertex to the list of vertices
  AXOM_HOST_DEVICE
  void addVertex(const PointType& pt)
  {
    SLIC_ASSERT(m_num_vertices + 1 < MAX_VERTS);
    m_vertices[m_num_vertices] = pt;
    m_num_vertices++;
  }

  /// Clears the list of vertices
  AXOM_HOST_DEVICE
  void clear()
  {
    for(int i = 0; i < MAX_VERTS; i++)
    {
      m_vertices[i] = PointType();
    }
    m_num_vertices = 0;
  }

  /// Retrieves the vertex at index idx
  AXOM_HOST_DEVICE
  PointType& operator[](int idx) { return m_vertices[idx]; }

  /// Retrieves the vertex at index idx
  AXOM_HOST_DEVICE
  const PointType& operator[](int idx) const { return m_vertices[idx]; }

  /*! 
   * \brief Robustly returns the normal of the polygon (not normalized)
   * \pre This function is only valid when NDIMS = 3 and the polygon is valid
   * \return normal polygon normal vector
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, VectorType>::type normal() const
  {
    SLIC_ASSERT(isValid());
    const int nverts = numVertices();

    VectorType v0(m_vertices[nverts - 1]), v1(m_vertices[0]);

    VectorType normal = VectorType::cross_product(v0, v1);

    // Iterate over each pair of vertices
    for(int i = 1; i < nverts; ++i)
    {
      v0 = v1;
      v1 = VectorType(m_vertices[i]);
      normal += VectorType::cross_product(v0, v1);
    }

    return normal;
  }

  /*!
   * \brief Computes the average of the polygon's vertex positions
   *
   * \return A point at the mean of the polygon's vertices
   * \pre  polygon.isValid() is true
   */
  AXOM_HOST_DEVICE
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

  /**
   * \brief Returns the area of the polygon (3D specialization)
   *
   * The algorithm sums up contributions from triangles defined by each edge 
   * of the polygon and a local origin (at the Polygon's \a vertexMean()).
   * The algorithm assumes a planar polygon embedded in 3D, however it will return 
   * a value for non-planar polygons.
   *
   * The area is always non-negative since 3D polygons do not have a unique orientation.
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 3, double>::type area() const
  {
    const int nVerts = numVertices();

    // check for early return
    if(nVerts < 3)
    {
      return 0.0;
    }

    // Add up areas of triangles connecting polygon edges the vertex average
    VectorType sum;
    const auto O = vertexMean();  // 'O' for (local) origin
    for(int curr = 0, prev = nVerts - 1; curr < nVerts; prev = curr++)
    {
      sum +=
        VectorType::cross_product(m_vertices[prev] - O, m_vertices[curr] - O);
    }

    return 0.5 * axom::utilities::abs(sum.norm());
  }

  /**
   * \brief Returns the area of a planar polygon (2D specialization)
   *
   * The area is always non-negative.
   * \sa signedArea()
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type area() const
  {
    return axom::utilities::abs(signedArea());
  }

  /**
   * \brief Returns the signed area of a 2D polygon
   *
   * The signed area accounts for the orientation of the polygon.
   * It is positive when the vertices are oriented counter-clockwise
   * and negative when the vertices are oriented clockwise.
   * \note Signed area is only defined in 2D.
   * \sa area()
   */
  template <int TDIM = NDIMS>
  AXOM_HOST_DEVICE typename std::enable_if<TDIM == 2, double>::type signedArea() const
  {
    const int nVerts = numVertices();
    double sum = 0.;

    // check for early return
    if(nVerts < 3)
    {
      return sum;
    }

    // use shoelace algorithm
    for(int curr = 0, prev = nVerts - 1; curr < nVerts; prev = curr++)
    {
      sum += m_vertices[prev][0] * m_vertices[curr][1]  //
        - m_vertices[curr][0] * m_vertices[prev][1];
    }

    return 0.5 * sum;
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
  AXOM_HOST_DEVICE
  bool isValid() const { return m_num_vertices >= 3; }

private:
  PointType m_vertices[MAX_VERTS];
  int m_num_vertices {0};
};

//------------------------------------------------------------------------------
/// Free functions implementing Polygon's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS, int MAX_VERTS>
std::ostream& operator<<(std::ostream& os,
                         const Polygon<T, NDIMS, MAX_VERTS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::Polygon using fmt
template <typename T, int NDIMS, int MAX_VERTS>
struct axom::fmt::formatter<axom::primal::Polygon<T, NDIMS, MAX_VERTS>>
  : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_POLYGON_HPP_
