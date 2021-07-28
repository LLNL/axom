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

#include "axom/core/StackArray.hpp"

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
 * \class NeighborCollection
 *
 * \brief Represents a collection of neighbor relations between vertices.
 */
class NeighborCollection
{
public:
  static constexpr int MAX_VERTS = 32;
  static constexpr int MAX_NBRS_PER_VERT = 8;

  using VertexNbrs = axom::StackArray<axom::int8, MAX_NBRS_PER_VERT>;

public:
  AXOM_HOST_DEVICE NeighborCollection() : num_nbrs {0} { }

  AXOM_HOST_DEVICE void clear()
  {
    for(int i = 0; i < MAX_VERTS; i++)
    {
      num_nbrs[i] = 0;
    }
  }

  AXOM_HOST_DEVICE int getNumNeighbors(int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return num_nbrs[vtx];
  }

  AXOM_HOST_DEVICE const VertexNbrs& operator[](int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  AXOM_HOST_DEVICE VertexNbrs& operator[](int vtx)
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  AXOM_HOST_DEVICE const VertexNbrs& getNeighbors(int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  AXOM_HOST_DEVICE VertexNbrs& getNeighbors(int vtx)
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  AXOM_HOST_DEVICE void addNeighbors(axom::int8 idx1,
                                     std::initializer_list<axom::int8> idx2)
  {
    SLIC_ASSERT(num_nbrs[idx1] + idx2.size() <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(idx1 >= 0 && idx1 < MAX_VERTS);
    for(axom::int8 nbr : idx2)
    {
      axom::int8 idx_insert = num_nbrs[idx1];
      nbrs[idx1][idx_insert] = nbr;
      num_nbrs[idx1]++;
    }
  };

  AXOM_HOST_DEVICE void insertNeighborAtPos(axom::int8 idx1,
                                            axom::int8 idx2,
                                            axom::int8 pos)
  {
    SLIC_ASSERT(num_nbrs[idx1] + 1 <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(idx1 >= 0 && idx1 < MAX_VERTS);
    axom::uint8 old_nbrs[MAX_NBRS_PER_VERT];
    // copy elements from [pos, nnbrs)
    for(int ip = pos; ip < num_nbrs[idx1]; ip++)
    {
      old_nbrs[ip - pos] = nbrs[idx1][ip];
    }
    nbrs[idx1][pos] = idx2;
    for(int ip = pos; ip < num_nbrs[idx1]; ip++)
    {
      nbrs[idx1][ip + 1] = old_nbrs[ip - pos];
    }
    num_nbrs[idx1]++;
  }

  AXOM_HOST_DEVICE void pruneNeighbors()
  {
    for(int iv = 0; iv < MAX_VERTS; iv++)
    {
      int curr_idx = 0;
      for(int inbr = 0; inbr < num_nbrs[iv]; inbr++)
      {
        if(nbrs[iv][inbr] == -1) continue;
        nbrs[iv][curr_idx] = nbrs[iv][inbr];
        curr_idx++;
      }
      num_nbrs[iv] = curr_idx;
    }
  }

private:
  axom::int8 num_nbrs[MAX_VERTS];
  VertexNbrs nbrs[MAX_VERTS];
};

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

  constexpr static int MAX_VERTS = 32;

private:
  using Coords = StackArray<PointType, MAX_VERTS>;
  using Neighbors = NeighborCollection;

public:
  /*! Default constructor for an empty polyhedron   */
  AXOM_HOST_DEVICE Polyhedron() : m_num_vertices(0) { }

  /*! Return the number of vertices in the polyhedron */
  AXOM_HOST_DEVICE int numVertices() const { return m_num_vertices; }

  /*!
   * \brief Appends a vertex to the list of vertices.
   *
   * \return The index where the vertex was inserted into.
   */
  AXOM_HOST_DEVICE int addVertex(const PointType& pt)
  {
    SLIC_ASSERT(m_num_vertices < MAX_VERTS);
    m_vertices[m_num_vertices] = pt;
    m_num_vertices++;
    return m_num_vertices - 1;
  }

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
  AXOM_HOST_DEVICE
  void addNeighbors(const PointType& pt, std::initializer_list<axom::int8> nbrs)
  {
    for(unsigned int i = 0; i < m_num_vertices; i++)
    {
      if(m_vertices[i] == pt)
      {
        m_neighbors.addNeighbors(i, nbrs);
      }
    }
  }

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
  AXOM_HOST_DEVICE
  void addNeighbors(int idx, std::initializer_list<axom::int8> nbrs)
  {
    m_neighbors.addNeighbors(idx, nbrs);
  }

  /*! Clears the list of vertices and neighbors */
  AXOM_HOST_DEVICE void clear()
  {
    m_num_vertices = 0;
    m_neighbors.clear();
  }

  /*! Retrieves the vertices */
  Coords& getVertices() { return m_vertices; }

  /*! Retrieves the vertex at index idx */
  AXOM_HOST_DEVICE PointType& operator[](int idx) { return m_vertices[idx]; }
  /*! Retrieves the vertex at index idx */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    return m_vertices[idx];
  }

  /*! Retrieves the neighbors */
  AXOM_HOST_DEVICE Neighbors& getNeighbors() { return m_neighbors; }

  AXOM_HOST_DEVICE int getNumNeighbors(int i) const
  {
    return m_neighbors.getNumNeighbors(i);
  }

  /*! Retrieves the neighbors for vertex i */
  AXOM_HOST_DEVICE
  typename Neighbors::VertexNbrs& getNeighbors(int i) { return m_neighbors[i]; }
  /*! Retrieves the neighbors for vertex i */
  AXOM_HOST_DEVICE
  const typename Neighbors::VertexNbrs& getNeighbors(int i) const
  {
    return m_neighbors[i];
  }

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
      for(int nidx = 0; nidx < m_neighbors.getNumNeighbors(i); nidx++)
      {
        int ni = m_neighbors[i][nidx];
        // Check if edge has not been visited
        std::vector<int> edgeToCheck {ni, i};
        if(std::find(checkedEdges.begin(), checkedEdges.end(), edgeToCheck) ==
           checkedEdges.end())
        {
          std::vector<int> face {ni};
          axom::int8 vstart = ni;
          axom::int8 vnext = i;
          axom::int8 vprev = ni;

          // Add neighboring vertices until we reach the starting vertex.
          while(vnext != vstart)
          {
            face.push_back(vnext);
            checkedEdges.push_back({vprev, vnext});
            int numNeighbors = m_neighbors.getNumNeighbors(vnext);
            auto itr = std::find(m_neighbors[vnext] + 0,
                                 m_neighbors[vnext] + numNeighbors,
                                 vprev);
            vprev = vnext;
            if(itr == m_neighbors[vnext] + 0)
            {
              vnext = m_neighbors[vnext][numNeighbors - 1];
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
      if(hasNeighbors() && m_neighbors.getNumNeighbors(i) > 0)
      {
        os << "Neighbors of vertex " << i << ": ";
        int nz = m_neighbors.getNumNeighbors(i);
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
      if(hasNeighbors() && m_neighbors.getNumNeighbors(sz - 1) > 0)
      {
        os << "Neighbors of vertex " << sz - 1 << ": ";
        int nz = m_neighbors.getNumNeighbors(sz - 1);
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
  bool isValid() const { return m_num_vertices >= 4; }

  /*!
   * \brief Check if vertex neighbors information is available
   *
   * Checks if neighbors is not empty
   * \return True, if the polyhedron has neighbors info for each vertex,
   *         False otherwise
   */
  bool hasNeighbors() const
  {
    bool has_nbrs = (m_num_vertices > 0);
    for(int i = 0; i < m_num_vertices; i++)
    {
      has_nbrs = has_nbrs && (m_neighbors.getNumNeighbors(i) > 0);
    }
    return has_nbrs;
  }

private:
  int m_num_vertices;
  Coords m_vertices;
  Neighbors m_neighbors;
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
