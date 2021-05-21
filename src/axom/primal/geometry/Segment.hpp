// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SEGMENT_HPP_
#define SEGMENT_HPP_

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int DIM>
class Segment;

/*!
 * \brief Overloaded output operator for Segment
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Segment<T, NDIMS>& seg);

/*!
 * \class
 *
 * \brief Represents a directed straight line connecting two points
 *  \f$ A,B \in \mathcal{R}^n \f$ where point A is referred to as the source
 *  point and point B is the target.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Segment
{
public:
  typedef Point<T, NDIMS> PointType;

public:
  /*!
   * \brief Creates a segment instance from point A to point B.
   * \param A user-supplied source point
   * \param B user-supplied target point
   */
  Segment(const PointType& A, const PointType& B);

  /*!
   * \brief Destructor.
   */
  ~Segment();

  /*!
   * \brief Returns the source point of the segment.
   * \return s the source point of the segment.
   */
  const PointType& source() const { return m_source; };

  /*!
   * \brief Returns the target point of the segment.
   * \return t the target point of the segment.
   */
  const PointType& target() const { return m_target; };

  /*!
   * \brief Returns a point \f$ (1 - t)A + tB \f$
   * \param [in] t user-supplied parameter for blending A and B
   * \post Return point P will always be collinear with the segment's
   *       start point A and end point B.
   * \post If \f$ t = 0, \f$ the return point \f$ P = A. \f$
   * \post If \f$ t = 1, \f$ the return point \f$ P = B. \f$
   */
  PointType at(const T& t) const
  {
    return PointType::lerp(m_source, m_target, t);
  }

  /*!
   * \brief Returns the length of the segment
   */
  double length() const
  {
    typedef Vector<T, NDIMS> VectorType;
    return VectorType(m_source, m_target).norm();
  }

  /*!
   * \brief Simple formatted print of a segment instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    os << "{source:" << m_source << "; target:" << m_target << "}";

    return os;
  }

private:
  /*!
   * \brief Default Constructor. Does nothing.
   * \note Made private to prevent its use in application code.
   */
  Segment() {};

  PointType m_source;
  PointType m_target;
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
//  Segment Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
template <typename T, int NDIMS>
Segment<T, NDIMS>::Segment(const PointType& A, const PointType& B)
  : m_source(A)
  , m_target(B)
{ }

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Segment<T, NDIMS>::~Segment()
{ }

//------------------------------------------------------------------------------
/// Free functions implementing Segments's operators
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Segment<T, NDIMS>& seg)
{
  seg.print(os);
  return os;
}

} /* namespace primal */
} /* namespace axom */

#endif /* SEGMENT_HPP_ */
