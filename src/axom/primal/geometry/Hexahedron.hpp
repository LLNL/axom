// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef HEXAHEDRON_HPP_
#define HEXAHEDRON_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/slic/interface/slic.hpp"

#include <ostream>  // for std::ostream

namespace axom
{
namespace primal
{
/*!
 * \class Hexahedron
 *
 * \brief Represents an hexahedral geometric shape defined by eight points.
 *
 * \accelerated
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 *
 * There are eight vertices in the hexahedron, labelled P through W as the
 * constructor's arguments.  They are accessible using the square-brackets
 * operator, with P being index 0, Q index 1, through W as index 7.
 *
 * Here's a diagram showing
 * a cube with the labeled vertices.
 *
 * \verbatim
 *
 * W +---------+ V
 *   |\        |\
 *   |  \      |  \
 *   | T + --------+ U
 * S +---|-----+ R |
 *   \   |     \   |
 *    \  |      \  |
 *     \ |       \ |
 *   P  +----------+ Q
 *
 * \endverbatim
 *
 */
template <typename T, int NDIMS = 3>
class Hexahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  enum
  {
    NUM_HEX_VERTS = 8
  };

public:
  /*!
   * \brief Default constructor. Creates a degenerate hexahedron.
   */
  AXOM_HOST_DEVICE Hexahedron() { }

  /*!
   * \brief Creates an hexahedron from the 8 points p,q,r,s,t,u,v,w.
   * \param [in] Point corresponding to vertex p of the hexahedron.
   * \param [in] Point corresponding to vertex q of the hexahedron.
   * \param [in] Point corresponding to vertex r of the hexahedron.
   * \param [in] Point corresponding to vertex s of the hexahedron.
   * \param [in] Point corresponding to vertex t of the hexahedron.
   * \param [in] Point corresponding to vertex u of the hexahedron.
   * \param [in] Point corresponding to vertex v of the hexahedron.
   * \param [in] Point corresponding to vertex w of the hexahedron.
   *
   */
  AXOM_HOST_DEVICE
  Hexahedron(const PointType& p,
             const PointType& q,
             const PointType& r,
             const PointType& s,
             const PointType& t,
             const PointType& u,
             const PointType& v,
             const PointType& w)
  {
    m_points[0] = p;
    m_points[1] = q;
    m_points[2] = r;
    m_points[3] = s;
    m_points[4] = t;
    m_points[5] = u;
    m_points[6] = v;
    m_points[7] = w;
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, 5, 6, or 7
   */
  AXOM_HOST_DEVICE PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_HEX_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, 5, 6, or 7
   */
  AXOM_HOST_DEVICE const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_HEX_VERTS);
    return m_points[idx];
  }

  /*!
    * \brief Test if this Hexahedron is equal to another, within a tolerance
    */
  AXOM_HOST_DEVICE
  bool equals(const Hexahedron& other, double eps = 1.e-24) const
  {
    // Two hexes are equal if each vertex is closer than eps to a vertex of the other.
    int matched[NUM_HEX_VERTS];
    for(int i = 0; i < NUM_HEX_VERTS; ++i)
    {
      matched[i] = 0;
    }
    for(int ourvert = 0; ourvert < NUM_HEX_VERTS; ++ourvert)
    {
      for(int theirvert = 0; theirvert < NUM_HEX_VERTS; ++theirvert)
      {
        if(!matched[theirvert] &&
           squared_distance(m_points[ourvert], other[theirvert]) < eps)
        {
          matched[theirvert] = 1;
        }
      }
    }
    int matchedcount = 0;
    for(int i = 0; i < NUM_HEX_VERTS; ++i)
    {
      matchedcount += matched[i];
    }

    return (matchedcount == NUM_HEX_VERTS);
  }

  /*!
   * \brief Simple formatted print of an hexahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " "
       << m_points[3] << " " << m_points[4] << " " << m_points[5] << " "
       << m_points[6] << " " << m_points[7] << "}";

    return os;
  }

private:
  PointType m_points[NUM_HEX_VERTS];
};

//------------------------------------------------------------------------------
/// Free functions implementing Hexahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Hexahedron<T, NDIMS>& hex)
{
  hex.print(os);
  return os;
}

} /* namespace primal */
} /* namespace axom */

#endif /* HEXAHEDRON_HPP_ */
