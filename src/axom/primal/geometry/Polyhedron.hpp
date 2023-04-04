// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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

#include <ostream>

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

  using VertexNbrs = axom::StackArray<std::int8_t, MAX_NBRS_PER_VERT>;

public:
  /*!
   * \brief Constructs an empty NeighborCollection.
   */
  AXOM_HOST_DEVICE NeighborCollection() : num_nbrs {0} { }

  /*!
   * \brief Clears the set of neighbors.
   */
  AXOM_HOST_DEVICE void clear()
  {
    for(int i = 0; i < MAX_VERTS; i++)
    {
      num_nbrs[i] = 0;
    }
  }

  /*!
   * \brief Returns the number of neighbors of a given vertex.
   *
   * \param [in] vtx index of the vertex to check.
   * \return num_nbrs the number of neighbors the vertex has.
   */
  AXOM_HOST_DEVICE int getNumNeighbors(int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return num_nbrs[vtx];
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE const VertexNbrs& operator[](int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return getNeighbors(vtx);
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE VertexNbrs& operator[](int vtx)
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return getNeighbors(vtx);
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE const VertexNbrs& getNeighbors(int vtx) const
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  /*!
   * \brief Gets the array of neighbors for a given vertex index vtx.
   */
  AXOM_HOST_DEVICE VertexNbrs& getNeighbors(int vtx)
  {
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    return nbrs[vtx];
  }

  /*!
   * \brief Adds a set of neighbors to a given vertex.
   *
   * \param [in] vtx the index of the vertex to add new neighbors to
   * \param [in] nbr a set of indices of neighboring vertices
   *
   * \pre vtx < MAX_VERTS
   */
  AXOM_HOST_DEVICE void addNeighbors(std::int8_t vtx,
                                     std::initializer_list<std::int8_t> nbrIds)
  {
    SLIC_ASSERT(num_nbrs[vtx] + nbrIds.size() <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    for(std::int8_t nbr : nbrIds)
    {
      std::int8_t idx_insert = num_nbrs[vtx];
      nbrs[vtx][idx_insert] = nbr;
      num_nbrs[vtx]++;
    }
  };

  /*!
   * \brief Adds a single neighbor to a given vertex.
   *
   * \param [in] vtx the index of the vertex to add the new neighbor to
   * \param [in] nbr an index of a neighboring vertex
   *
   * \pre vtx < MAX_VERTS
   */
  AXOM_HOST_DEVICE void addNeighbors(std::int8_t vtx, std::int8_t nbrId)
  {
    SLIC_ASSERT(num_nbrs[vtx] + 1 <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    std::int8_t idx_insert = num_nbrs[vtx];
    nbrs[vtx][idx_insert] = nbrId;
    num_nbrs[vtx]++;
  };

  /*!
   * \brief Inserts a new neighbor for a vertex at a given index in the
   *  neighbor list.
   *
   * \param [in] vtx the index of the vertex to add a new neighbor to
   * \param [in] nbr the index of the neighboring vertex
   * \param [in] pos the position in the neighbor list to insert into
   *
   * \pre vtx < MAX_VERTS
   * \pre pos <= num_nbrs[vtx]
   */
  AXOM_HOST_DEVICE void insertNeighborAtPos(std::int8_t vtx,
                                            std::int8_t nbr,
                                            std::int8_t pos)
  {
    SLIC_ASSERT(num_nbrs[vtx] + 1 <= MAX_NBRS_PER_VERT);
    SLIC_ASSERT(vtx >= 0 && vtx < MAX_VERTS);
    SLIC_ASSERT(pos <= num_nbrs[vtx]);
    std::uint8_t old_nbrs[MAX_NBRS_PER_VERT];
    // copy elements from [pos, nnbrs)
    for(int ip = pos; ip < num_nbrs[vtx]; ip++)
    {
      old_nbrs[ip - pos] = nbrs[vtx][ip];
    }
    nbrs[vtx][pos] = nbr;
    for(int ip = pos; ip < num_nbrs[vtx]; ip++)
    {
      nbrs[vtx][ip + 1] = old_nbrs[ip - pos];
    }
    num_nbrs[vtx]++;
  }

  /*!
   * \brief Compacts all neighbor sets by removing neighbors with index -1.
   */
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
  std::int8_t num_nbrs[MAX_VERTS];
  VertexNbrs nbrs[MAX_VERTS];
};

/*!
 * \class Polyhedron
 *
 * \brief Represents a convex polyhedron defined by an array of points and
 *        optionally their neighbors (a point is a neighbor of another point if
 *        there is an edge between them)
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \note The Polyhedron functions do not check that points defining a face are
 *       coplanar. It is the responsibility of the caller to pass a
 *       valid set of points and neighbors representing a convex polyhedron.
 *
 * \note The polyhedron neighbors should be ordered counter clockwise
 *       as when viewing the polyhedron externally.
 *
 *       <pre>
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
 *       </pre>
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

  constexpr static int MAX_VERTS = NeighborCollection::MAX_VERTS;

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
    SLIC_ASSERT(m_num_vertices + 1 < MAX_VERTS);
    m_vertices[m_num_vertices] = pt;
    m_num_vertices++;
    return m_num_vertices - 1;
  }

  /*!
   * \brief Stores nbrs as the neighbors of vertex pt.
   *        If vertex pt is not in the polyhedron, neighbors are not added.
   *
   * \note pt should be a vertex of this polyhedron. Otherwise, pt may not
   *       be found due to floating point precision differences.
   *
   * \param [in] pt The vertex to add neighbors
   * \param [in] nbrs The neighbors to add to the list of neighbors
   */
  AXOM_HOST_DEVICE
  void addNeighbors(const PointType& pt, std::initializer_list<std::int8_t> nbrs)
  {
    for(int i = 0; i < m_num_vertices; i++)
    {
      if(m_vertices[i] == pt)
      {
        m_neighbors.addNeighbors(i, nbrs);
      }
    }
  }

  /*!
   * \brief Stores nbrs as the neighbors of vertex at a given index.
   *
   * \note Caller is responsible for ensuring a given vertex has enough space
   *       remaining to fit the new neighbors.
   *
   * \param [in] vtxId The vertex id to add neighbors
   * \param [in] nbrs The neighbors to add to the list of neighbors
   *
   * \pre vtxId < getVertices()
   */
  AXOM_HOST_DEVICE
  void addNeighbors(int vtxId, std::initializer_list<std::int8_t> nbrs)
  {
    m_neighbors.addNeighbors(vtxId, nbrs);
  }

  /*!
   * \brief Adds a single nbr as the neighbor of vertex at a given index.
   *
   * \note Caller is responsible for ensuring a given vertex has enough space
   *       remaining to fit the new neighbor.
   *
   * \param [in] vtxId The vertex id to add the neighbor
   * \param [in] nbr The neighbor to add to the list of neighbors
   *
   * \pre vtxId < getVertices()
   */
  AXOM_HOST_DEVICE
  void addNeighbors(int vtxId, int nbr)
  {
    m_neighbors.addNeighbors(vtxId, nbr);
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
  AXOM_HOST_DEVICE
  PointType centroid() const
  {
    SLIC_ASSERT(isValid());

    NumArrayType sum;

    const int nVerts = numVertices();
    if(nVerts > 0)
    {
      for(int i = 0; i < nVerts; ++i)
      {
        sum += m_vertices[i].array();
      }
      sum /= nVerts;
    }
    return PointType(sum);
  }

public:
  /*!
   * \brief Helper function to find the faces of the Polyhedron, assuming the
   *        vertex neighbors are in counter-clockwise ordering.
   *
   * \param [out] faces is the vertex indices for faces
   * \param [out] face_offset is the offset for each face
   * \param [out] face_size is the number of vertices for each face
   * \param [out] face_count is the number of faces
   *
   * \warning Function is experimental, input parameters and/or output may
   *          change in the future.
   *
   * \note This function does not check that the provided buffers are
   *       large enough to hold all values. It is the responsibility of the
   *       caller to know ahead of time and pass buffers with appropriate size.
   *
   * \note Function is based off extractFaces() in Mike Owen's PolyClipper.
   *
   * \pre polyhedron vertex neighbors are defined
   */
  AXOM_HOST_DEVICE
  void getFaces(int* faces, int* face_size, int* face_offset, int& face_count) const
  {
    std::int8_t curFaceIndex = 0;
    std::int8_t checkedSize = 0;
    std::int8_t facesAdded = 0;
    // # edges * (# vertices per edge) * (# orientation per edge)
    std::int8_t checkedEdges[MAX_VERTS * 2 * 2] = {0};

    // Check each vertex
    for(int i = 0; i < numVertices(); ++i)
    {
      // Check each neighbor index (edge) of the vertex
      for(int j = 0; j < m_neighbors.getNumNeighbors(i); j++)
      {
        // Check if edge has not been visited
        int ni = m_neighbors[i][j];
        int edgeToCheck[2] = {i, ni};

        bool alreadyChecked = false;

        for(int edge = 0; edge < checkedSize; edge++)
        {
          if(edgeToCheck[0] == checkedEdges[edge * 2] &&
             edgeToCheck[1] == checkedEdges[(edge * 2) + 1])
          {
            alreadyChecked = true;
            break;
          }
        }

        if(!alreadyChecked)
        {
          face_offset[facesAdded] = curFaceIndex;
          faces[curFaceIndex++] = i;
          std::int8_t curFaceSize = 1;
          std::int8_t vstart = i;
          std::int8_t vnext = ni;
          std::int8_t vprev = i;

          // Add neighboring vertices until we reach the starting vertex.
          while(vnext != vstart)
          {
            faces[curFaceIndex++] = vnext;
            curFaceSize++;
            checkedEdges[checkedSize * 2] = vprev;
            checkedEdges[(checkedSize * 2) + 1] = vnext;
            checkedSize++;

            int numNeighbors = m_neighbors.getNumNeighbors(vnext);
            int itr = -1;
            for(int k = 0; k < m_neighbors.getNumNeighbors(vnext); k++)
            {
              if(m_neighbors[vnext][k] == vprev)
              {
                itr = k;
                break;
              }
            }
            vprev = vnext;
            if(itr == 0)
            {
              vnext = m_neighbors[vnext][numNeighbors - 1];
            }
            else
            {
              vnext = m_neighbors[vnext][itr - 1];
            }
          }  // end of while loop

          // Add last edge connecting the face
          checkedEdges[checkedSize * 2] = vprev;
          checkedEdges[(checkedSize * 2) + 1] = vnext;
          checkedSize++;

          face_size[facesAdded++] = curFaceSize++;
        }
      }
    }
    face_count = facesAdded;
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
  AXOM_HOST_DEVICE
  double volume() const
  {
    double retVol = 0.0;

    if(!isValid())
    {
      return retVol;
    }

    // Finds the volume of tetrahedrons formed from vertices of the Polyhedron
    // faces and an arbitrary origin (the first vertex)
    else
    {
      SLIC_CHECK_MSG(
        hasNeighbors(),
        "Polyhedron::volume() is only valid with vertex neighbors.");

      // faces is an overestimation
      int faces[MAX_VERTS * MAX_VERTS];
      int face_size[MAX_VERTS * 2];
      int face_offset[MAX_VERTS * 2];
      int face_count;
      getFaces(faces, face_size, face_offset, face_count);

      const PointType& origin = m_vertices[0];

      for(int i = 0; i < face_count; ++i)
      {
        const int N = face_size[i];
        const int i_offset = face_offset[i];
        const VectorType v0 = m_vertices[faces[i_offset]] - origin;

        for(int j = 1, k = 2; j < N - 1; ++j, ++k)
        {
          retVol += VectorType::scalar_triple_product(
            v0,
            m_vertices[faces[i_offset + j]] - origin,
            m_vertices[faces[i_offset + k]] - origin);
        }
      }
    }

    return retVol / 6.;
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
          int curN = getNeighbors(i)[j];
          os << curN << " ";
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
          int curN = getNeighbors(sz - 1)[j];
          os << curN << " ";
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
  AXOM_HOST_DEVICE
  bool isValid() const { return m_num_vertices >= 4; }

  /*!
   * \brief Check if vertex neighbors information is available
   *
   * Checks if neighbors is not empty
   * \return True, if the polyhedron has neighbors info for each vertex,
   *         False otherwise
   */
  AXOM_HOST_DEVICE
  bool hasNeighbors() const
  {
    if(m_num_vertices <= 1)
    {
      return false;
    }
    bool has_nbrs = true;
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

/// Overload to format a primal::Polyhedron using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::Polyhedron<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_POLYHEDRON_HPP_
