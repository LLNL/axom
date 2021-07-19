// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file Polyhedron.hpp
 *
 * \brief A Polyhedron primitive for primal
 */

#ifndef AXOM_PRIMAL_POLYHEDRON_HPP_
#define AXOM_PRIMAL_POLYHEDRON_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/NumericArray.hpp"

#include <vector>
#include <ostream>  // for std::ostream

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class Polyhedron;

/*! \brief Overloaded output operator for polyhedrons */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polyhedron<T, NDIMS>& poly);

/*!
 * \class Polyhedron
 *
 * \brief Represents a polyhedron defined by an array of points and optionally
 *        their neighbors (a point is a neighbor of another point if there is
 *        an edge between them)
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \note The Polyhedron functions do not check that points defining a face are
 *       coplanar. It is the responsibility of the caller to pass a
 *       valid set of points and neighbors.
 *
 * \note The polyhedron neighbors should be ordered counter clockwise
 *       as when viewing the polyhedron externally.
 *
 *          4--------7
 *         /|       /|
 *        / |      / |
 *       5--------6  |
 *       |  0-----|--3
 *       | /      | /
 *       |/       |/
 *       1--------2
 *
 *       For example, vertex 5 above should have neighbors (4,1,6).
 *
 * \note The Polyhedron functions do not check that neighbors are ordered
 *       counter clockwise. It is the responsibility of the caller to pass a
 *       valid neighbors ordering.
 */
template <typename T, int NDIMS = 3>
class Polyhedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = NumericArray<T, NDIMS>;

private:
  using Coords = std::vector<PointType>;

public:
  /*! Default constructor for an empty polyhedron   */
  Polyhedron() { }

  /*!
   * \brief Constructor for an empty polyhedron that reserves space for
   *  the given number of vertices and optionally neighbors.
   *
   * \param [in] numExpectedVerts number of vertices for which to reserve space
   * \param [in] defineNeighbors optionally reserves space for neighbors
   *             (false by default)
   * \pre numExpectedVerts is not negative
   *
   */
  Polyhedron(int numExpectedVerts, bool defineNeighbors = false)
  {
    SLIC_ASSERT(numExpectedVerts >= 0);
    m_vertices.reserve(numExpectedVerts);

    if(defineNeighbors)
    {
      m_neighbors.reserve(numExpectedVerts);
    }
  }

  /*! Return the number of vertices in the polyhedron */
  int numVertices() const { return static_cast<int>(m_vertices.size()); }

  /*! Appends a vertex to the list of vertices */
  void addVertex(const PointType& pt) { m_vertices.push_back(pt); }

  /*!
   * \brief Stores nbrs as the neighbors of vertex pt.
   *        If vertex pt is not in the polyhedron, neighbors are not added.
   *        Resizes neighbors list to number of vertices if neighbors list
   *        size is less than pt's index.
   *
   * \note pt should be a vertex of this polyhedron. Otherwise, pt may not
   *       be found due to floating point precision differences.
   *
   * \param [in] pt The vertex to add neighbors
   * \param [in] nbrs The neighbors to add to the list of neighbors
   */
  void addNeighbors(const PointType& pt, const std::vector<int>& nbrs)
  {
    for(unsigned int i = 0; i < m_vertices.size(); i++)
    {
      if(m_vertices[i] == pt)
      {
        if(m_neighbors.size() < i + 1)
        {
          m_neighbors.resize(numVertices());
        }
        m_neighbors[i] = nbrs;
      }
    }
  }

  /*! Clears the list of vertices and neighbors */
  void clear()
  {
    m_vertices.clear();
    m_neighbors.clear();
  }

  /*! Retrieves the vertices */
  Coords& getVertices() { return m_vertices; }

  /*! Retrieves the vertex at index idx */
  PointType& operator[](int idx) { return m_vertices[idx]; }
  /*! Retrieves the vertex at index idx */
  const PointType& operator[](int idx) const { return m_vertices[idx]; }

  /*! Retrieves the neighbors */
  std::vector<std::vector<int>>& getNeighbors() { return m_neighbors; }

  /*! Retrieves the neighbors for vertex i */
  std::vector<int>& getNeighbors(int i) { return m_neighbors[i]; }
  /*! Retrieves the neighbors for vertex i */
  const std::vector<int>& getNeighbors(int i) const { return m_neighbors[i]; }

  /*!
   * \brief Computes the centroid as the average of the polyhedron's vertex
   *  positions
   *
   * \return The centroid of the polyhedron's vertices
   *
   * \pre  polyhedron.isValid() is true
   */
  PointType centroid() const
  {
    SLIC_ASSERT(isValid());

    NumArrayType sum;

    for(int i = 0; i < numVertices(); ++i)
    {
      sum += m_vertices[i].array();
    }
    sum /= numVertices();

    return sum;
  }

  /*!
   * \brief Finds the faces of the Polyhedron, assuming the vertex neighbors
   *        are in counter-clockwise ordering.
   *
   * \return The polyhedron faces as vertex indices for each face
   *
   * \note Function is based off extractFaces() in Mike Owen's PolyClipper.
   *
   * \pre polyhedron vertex neighbors are defined
   */
  std::vector<std::vector<int>> getFaces() const
  {
    SLIC_CHECK_MSG(
      hasNeighbors(),
      "Polyhedron::getFaces() is only valid with vertex neighbors.");

    std::vector<std::vector<int>> faces;
    std::vector<std::vector<int>> checkedEdges;

    // Check each vertex
    for(int i = 0; i < numVertices(); i++)
    {
      // Check each neighbor index (edge) of the vertex
      for(int ni : m_neighbors[i])
      {
        // Check if edge has not been visited
        std::vector<int> edgeToCheck {ni, i};
        if(std::find(checkedEdges.begin(), checkedEdges.end(), edgeToCheck) ==
           checkedEdges.end())
        {
          std::vector<int> face {ni};
          int vstart = ni;
          int vnext = i;
          int vprev = ni;

          // Add neighboring vertices until we reach the starting vertex.
          while(vnext != vstart)
          {
            face.push_back(vnext);
            checkedEdges.push_back({vprev, vnext});
            auto itr = std::find(m_neighbors[vnext].begin(),
                                 m_neighbors[vnext].end(),
                                 vprev);
            vprev = vnext;
            if(itr == m_neighbors[vnext].begin())
            {
              vnext = m_neighbors[vnext].back();
            }
            else
            {
              vnext = *(itr - 1);
            }
          }  // end of while loop

          // Add last edge connecting the face
          checkedEdges.push_back({vprev, vnext});
          faces.push_back(face);
        }
      }
    }

    return faces;
  }

  /*!
   * \brief Finds the volume of the polyhedron.
   *
   * \return The volume of the polyhedron
   *
   * \note Function is based off moments() in Mike Owen's PolyClipper.
   *
   * \pre polyhedron vertex neighbors are defined, and polyhedron is 3D
   */
  double volume() const
  {
    SLIC_CHECK_MSG(hasNeighbors(),
                   "Polyhedron::volume() is only valid with vertex neighbors.");

    double retVol = 0.0;

    // 0 if less than 4 vertices
    if(numVertices() < 4)
    {
      return retVol;
    }

    // Finds the volume of tetrahedrons formed from vertices of the Polyhedron
    // faces and an arbitrary origin (the first vertex)
    else
    {
      std::vector<std::vector<int>> faces = getFaces();
      VectorType origin(m_vertices[0].data());
      for(std::vector<int>& face : faces)
      {
        int n = face.size();
        VectorType v0(m_vertices[face[0]].data());
        v0 -= origin;

        for(int i = 1; i < n - 1; ++i)
        {
          VectorType v1(m_vertices[face[i]].data());
          v1 -= origin;
          VectorType v2(m_vertices[face[(i + 1) % n]].data());
          v2 -= origin;
          double partialVol = v0.dot(VectorType::cross_product(v1, v2));
          retVol += partialVol;
        }
      }
    }

    return retVol / 6.0;
  }

  /*!
   * \brief Simple formatted print of a polyhedron instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int sz = numVertices();

    os << "{" << sz << "-vertex polyhedron:\n";
    for(int i = 0; i < sz - 1; ++i)
    {
      os << "Vertex " << i << ": " << m_vertices[i] << "\n";
      if(hasNeighbors() && !m_neighbors[i].empty())
      {
        os << "Neighbors of vertex " << i << ": ";
        int nz = m_neighbors[i].size();
        for(int j = 0; j < nz; j++)
        {
          os << m_neighbors[i][j] << " ";
        }
        os << "\n";
      }
    }
    if(sz >= 2)
    {
      os << "Vertex " << sz - 1 << ": " << m_vertices[sz - 1] << "\n";
      if(hasNeighbors() && !m_neighbors[sz - 1].empty())
      {
        os << "Neighbors of vertex " << sz - 1 << ": ";
        int nz = m_neighbors[sz - 1].size();
        for(int j = 0; j < nz; j++)
        {
          os << m_neighbors[sz - 1][j] << " ";
        }
        os << "\n";
      }
    }
    os << "}";

    return os;
  }

  /*!
   * \brief Simple check for validity of a polyhedron
   *
   * Initial check is that the polyhedron has four or more vertices
   * \return True, if the polyhedron is valid, False otherwise
   */
  bool isValid() const { return m_vertices.size() >= 4; }

  /*!
   * \brief Check if vertex neighbors information is available
   *
   * Checks if neighbors is not empty
   * \return True, if the polyhedron has neighbors info for each vertex,
   *         False otherwise
   */
  bool hasNeighbors() const { return !m_neighbors.empty(); }

private:
  Coords m_vertices;
  std::vector<std::vector<int>> m_neighbors;
};

//------------------------------------------------------------------------------
/// Free functions implementing Polyhedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Polyhedron<T, NDIMS>& poly)
{
  poly.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_POLYHEDRON_HPP_
