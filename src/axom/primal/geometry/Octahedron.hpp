// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef OCTAHEDRON_HPP_
#define OCTAHEDRON_HPP_

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
 * \class Octahedron
 *
 * \brief Represents an octahedral geometric shape defined by six points.
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of spatial dimensions
 */
template <typename T, int NDIMS=3>
class Octahedron
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;

  enum
  {
    NUM_OCT_VERTS = 6
  };

public:
  /*!
   * \brief Default constructor. Creates a degenerate octahedron.
   */
  Octahedron() { }

  /*!
   * \brief Creates an octahedron from the 6 points p,q,r,s,t,u.
   * \param [in] Point corresponding to vertex p of the octahedron.
   * \param [in] Point corresponding to vertex q of the octahedron.
   * \param [in] Point corresponding to vertex r of the octahedron.
   * \param [in] Point corresponding to vertex s of the octahedron.
   * \param [in] Point corresponding to vertex t of the octahedron.
   * \param [in] Point corresponding to vertex u of the octahedron.
   *
   * p is opposite s, q is opposite t, r is opposite u.
   */
  Octahedron(const PointType& p,
             const PointType& q,
             const PointType& r,
             const PointType& s,
             const PointType& t,
             const PointType& u)
  {
    m_points[0] = p;
    m_points[1] = q;
    m_points[2] = r;
    m_points[3] = s;
    m_points[4] = t;
    m_points[5] = u;
  }

  /*!
   * \brief Destructor
   */
  ~Octahedron() { }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, or 5
   */
  PointType& operator[](int idx)
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_OCT_VERTS);
    return m_points[idx];
  }

  /*!
   * \brief Index operator to get the i^th vertex
   * \param idx The index of the desired vertex
   * \pre idx is 0, 1, 2, 3, 4, or 5
   */
  const PointType& operator[](int idx) const
  {
    SLIC_ASSERT(idx >= 0 && idx < NUM_OCT_VERTS);
    return m_points[idx];
  }

   /*!
    * \brief Test if this Octahedron is equal to another, within a tolerance
    */
   bool equals(const Octahedron & other, double eps=1.e-24) const
   {
      // Two octs are equal if each vertex is closer than eps to a vertex of the other.
      int matched[NUM_OCT_VERTS];
      for (int i = 0; i < NUM_OCT_VERTS; ++i) { matched[i] = 0; }
      for (int ourvert = 0; ourvert < NUM_OCT_VERTS; ++ourvert)
      {
         for (int theirvert = 0; theirvert < NUM_OCT_VERTS; ++theirvert)
         {
            if (!matched[theirvert] && squared_distance(m_points[ourvert], other[theirvert]) < eps)
            {
               matched[theirvert] = 1;
            }
         }
      }
      int matchedcount = 0;
      for (int i = 0; i < NUM_OCT_VERTS; ++i)
      {
         matchedcount += matched[i];
      }
      
      return (matchedcount == NUM_OCT_VERTS);
   }

  /*!
   * \brief Simple formatted print of an octahedron instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{" << m_points[0] << " " << m_points[1] << " " << m_points[2] << " "
       << m_points[3] << " " << m_points[4] << " " << m_points[5] << "}";

    return os;
  }

private:
  PointType m_points[NUM_OCT_VERTS];
};

//------------------------------------------------------------------------------
/// Free functions implementing Octahedron's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Octahedron<T, NDIMS>& oct)
{
  oct.print(os);
  return os;
}

} /* namespace primal */
} /* namespace axom */

#endif /* OCTAHEDRON_HPP_ */
