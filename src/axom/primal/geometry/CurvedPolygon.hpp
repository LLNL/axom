// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file CurvedPolygon.hpp
 *
 * \brief A CurvedPolygon primitive whose edges are Bezier Curves
 */

#ifndef AXOM_PRIMAL_CURVEDPOLYGON_HPP_
#define AXOM_PRIMAL_CURVEDPOLYGON_HPP_

#include "axom/slic.hpp"

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"

#include <vector>
#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class CurvedPolygon;

/*! \brief Overloaded output operator for polygons */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const CurvedPolygon<T, NDIMS>& poly);

/*!
 * \class CurvedPolygon
 *
 * \brief Represents a polygon with curved edges defined by BezierCurves
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 * \note The component curves should be ordered in a counter clockwise
 *       orientation with respect to the polygon's normal vector
 */
template <typename T, int NDIMS>
class CurvedPolygon
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = NumericArray<T, NDIMS>;
  using BezierCurveType = BezierCurve<T, NDIMS>;

public:
  /*! Default constructor for an empty polygon   */
  CurvedPolygon() = default;

  /*!
   * \brief Constructor for an empty CurvedPolygon that reserves space for
   *  the given number of Edges
   *
   * \param [in] numExpectedEdges number of edges for which to reserve space
   * \pre numExpectedEdges is at least 1
   *
   */
  explicit CurvedPolygon(int nEdges)
  {
    SLIC_ASSERT(nEdges >= 1);
    m_edges.reserve(nEdges);
    m_edges.resize(nEdges);
  }

  CurvedPolygon(BezierCurveType* curves, int nEdges)
  {
    SLIC_ASSERT(curves != nullptr);
    SLIC_ASSERT(nEdges >= 1);

    m_edges.reserve(nEdges);

    for(int e = 0; e < nEdges; ++e)
    {
      this->addEdge(curves[e]);
    }
  }

  /*! Return the number of edges in the polygon */
  int numEdges() const { return m_edges.size(); }

  void setNumEdges(int ngon)
  {
    SLIC_ASSERT(ngon >= 0);
    m_edges.resize(ngon);
  }

  /* Checks equality of two Bezier Curve */
  friend inline bool operator==(const CurvedPolygon<T, NDIMS>& lhs,
                                const CurvedPolygon<T, NDIMS>& rhs)
  {
    return lhs.m_edges == rhs.m_edges;
  }

  friend inline bool operator!=(const CurvedPolygon<T, NDIMS>& lhs,
                                const CurvedPolygon<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /*! Appends a BezierCurve to the list of edges */
  void addEdge(const BezierCurveType& c1) { m_edges.push_back(c1); }

  /*! Splits an edge "in place" */
  void splitEdge(int idx, T t)
  {
    SLIC_ASSERT(idx < static_cast<int>(m_edges.size()));
    m_edges.insert(m_edges.begin() + idx + 1, 1, m_edges[idx]);
    auto& csplit = m_edges[idx];
    csplit.split(t, m_edges[idx], m_edges[idx + 1]);
  }

  /*! Clears the list of edges */
  void clear() { m_edges.clear(); }

  std::vector<BezierCurve<T, NDIMS>> getEdges() const { return m_edges; }

  /*! Retrieves the Bezier Curve at index idx */
  BezierCurveType& operator[](int idx) { return m_edges[idx]; }
  /*! Retrieves the vertex at index idx */
  const BezierCurveType& operator[](int idx) const { return m_edges[idx]; }

  /*!
   * \brief Simple formatted print of a CurvedPolygon instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numEdges();

    os << "{" << sz << "-sided Bezier polygon:";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << m_edges[i] << ",";
    }
    if(sz >= 2)
    {
      os << m_edges[sz - 1];
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Check closedness of a CurvedPolygon
   *
   * Check is that the endpoint of each edge coincides with startpoint of next edge
   * \return True, if the polygon is closed, False otherwise
   */
  bool isClosed(double tol = 1e-5) const
  {
    const int ngon = numEdges();
    if(ngon <= 2)
    {
      return false;
    }
    else
    {
      for(int p = 0; p < NDIMS; ++p)
      {
        for(int i = 0; i < (ngon - 1); ++i)
        {
          if(!axom::utilities::isNearlyEqual(m_edges[i][m_edges[i].getOrder()][p],
                                             m_edges[i + 1][0][p],
                                             tol))
          {
            return false;
          }
        }
        if(!axom::utilities::isNearlyEqual(
             m_edges[ngon - 1][m_edges[ngon - 1].getOrder()][p],
             m_edges[0][0][p],
             tol))
        {
          return false;
        }
      }
    }
    return true;
  }

  /*!
   * \brief Check closedness of a CurvedPolygon
   *
   * Check is that the endpoint of each edge coincides with startpoint of next edge
   * \return True, if the polygon is closed, False otherwise
   */
  T area(double tol = 1e-8) const
  {
    const int ngon = numEdges();
    T A = 0.0;
    if(!isClosed(1e3 * tol))
    {
      SLIC_DEBUG("Warning! The area is 0 because the element is not closed.");
      return A;
    }
    else
    {
      for(int ed = 0; ed < ngon; ++ed)
      {
        A += m_edges[ed].sectorArea();
      }
      return A;
    }
  }

  PointType centroid(double tol = 1e-8) const
  {
    const int ngon = numEdges();
    PointType M = PointType::make_point(0.0, 0.0);
    if(!isClosed(1e3 * tol))
    {
      SLIC_DEBUG(
        "Warning! The moments are 0 because the element is not closed.");
      return M;
    }
    else
    {
      const T A = area();
      if(A != 0.)
      {
        for(int ed = 0; ed < ngon; ++ed)
        {
          PointType Mc = m_edges[ed].sectorCentroid();
          M[0] += (Mc[0]);
          M[1] += (Mc[1]);
        }
        M.array() /= A;
      }
      return M;
    }
  }

  /*!
   * \brief Reverses orientation of a CurvedPolygon
   */
  void reverseOrientation()
  {
    const int ngon = numEdges();
    std::vector<BezierCurve<T, NDIMS>> old_edges = m_edges;
    for(int i = 0; i < ngon; ++i)
    {
      old_edges[ngon - 1 - i].reverseOrientation();
      m_edges[i] = old_edges[ngon - 1 - i];
    }
  }

private:
  std::vector<BezierCurve<T, NDIMS>> m_edges;
};

//------------------------------------------------------------------------------
/// Free functions implementing Polygon's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const CurvedPolygon<T, NDIMS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_CURVEDPOLYGON_HPP_
