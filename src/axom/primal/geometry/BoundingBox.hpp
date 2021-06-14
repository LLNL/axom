// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_BOUNDINGBOX_HPP_
#define AXOM_PRIMAL_BOUNDINGBOX_HPP_

#include <limits>

#include "axom/config.hpp"

#include "axom/core/Macros.hpp"  // for AXOM_HOST__DEVICE

#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

#include "axom/primal/operators/detail/intersect_bounding_box_impl.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class BoundingBox;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Equality comparison operator for bounding boxes.
 * Two bounding boxes are equal when they have the same bounds
 */
template <typename T, int NDIMS>
bool operator==(const BoundingBox<T, NDIMS>& lhs,
                const BoundingBox<T, NDIMS>& rhs);

/*!
 * \brief Inequality comparison operator for bounding boxes.
 * Two bounding boxes are unequal when they have different bounds
 */
template <typename T, int NDIMS>
bool operator!=(const BoundingBox<T, NDIMS>& lhs,
                const BoundingBox<T, NDIMS>& rhs);

/*!
 * \brief Overloaded output operator for bounding boxes
 */
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BoundingBox<T, NDIMS>& bb);

///@}

/*!
 * \class
 *
 * \brief BoundingBox represents and axis-aligned bounding box defined by
 * its min and max coordinates.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class BoundingBox
{
public:
  typedef T CoordType;
  typedef Point<T, NDIMS> PointType;
  typedef Vector<T, NDIMS> VectorType;
  typedef BoundingBox<T, NDIMS> BoxType;

  static constexpr T InvalidMin = std::numeric_limits<T>::max();
  static constexpr T InvalidMax = std::numeric_limits<T>::lowest();

public:
  /*!
   * \brief Constructor. Creates a bounding box with an invalid bound
   * The lower bound is set to the greatest possible point and the upper bound
   * is set to the smallest possible point.  This way adding any point resets
   * the bounds to a valid range.
   */
  AXOM_HOST_DEVICE
  BoundingBox() : m_min(PointType(InvalidMin)), m_max(PointType(InvalidMax)) { }

  /*!
   * \brief Constructor. Creates a bounding box containing a single point
   */
  AXOM_HOST_DEVICE
  BoundingBox(const PointType& pt) : m_min(pt), m_max(pt) { }

  /*!
   * \brief Constructor. Creates a bounding box containing the collection of
   * points.
   * \pre pt must point to at least n valid point
   * \note If n <= 0, defaults to default constructor values
   */
  BoundingBox(const PointType* pts, int n);

  /*!
   * \brief Constructor. Creates a bounding box with a given min and max point
   *  The code ensures that the bounds are valid.
   */
  AXOM_HOST_DEVICE
  BoundingBox(const PointType& lowerPt, const PointType& upperPt)
    : m_min(lowerPt)
    , m_max(upperPt)
  {
    this->checkAndFixBounds();
  }

  /*!
   * \brief Resets the bounds to those of the default constructor
   * \note This invalidates the bounding box (i.e. isValid() will be false)
   */
  void clear();

  /*!
   * \brief Returns const reference to the min corner of the bounding box.
   * \return const reference to the min corner of the bounding box.
   */
  AXOM_HOST_DEVICE
  const PointType& getMin() const { return m_min; };

  /*!
   * \brief Returns const reference to the max corner of the bounding box.
   * \return const reference to the max corner of the bounding box.
   */
  AXOM_HOST_DEVICE
  const PointType& getMax() const { return m_max; };

  /*!
   * \brief Returns the centroid (midpoint) of the bounding box.
   * \return Point at the bounding box centroid.
   */
  AXOM_HOST_DEVICE
  PointType getCentroid() const { return PointType::midpoint(m_min, m_max); }

  /*!
   * \brief Returns a vector from the min to the max points of the bounding box
   * \return Vector from min point to max point of bounding box.
   */
  AXOM_HOST_DEVICE
  VectorType range() const { return VectorType(m_min, m_max); };

  /*!
   * \brief Updates bounds to include the provided point.
   * \param [in] pt to include.
   */
  template <typename OtherType>
  AXOM_HOST_DEVICE void addPoint(const Point<OtherType, NDIMS>& pt);

  /*!
   * \brief Updates bounds to include the provided bounding box.
   * Convenience function -- equivalent to adding the min and max point of bbox
   * \param [in] bbox to include.
   */
  template <typename OtherType>
  AXOM_HOST_DEVICE void addBox(const BoundingBox<OtherType, NDIMS>& bbox);

  /*!
   * \brief Returns the dimension of the ambient space for this bounding box.
   * \return d the dimension of this bounding box instance.
   * \post d >= 1.
   */
  int dimension() const { return NDIMS; };

  /*!
   * \brief Finds the longest dimension of the bounding box
   * \return idx the index of the longest dimension.
   * \post idx >= 0 < NDIMS
   * \note In the case of ties, where the bounding box has more than one side
   *  with the same length, the code picks the first dimension as the longest
   *  dimension.
   */
  int getLongestDimension() const;

  /*!
   * \brief Intersects the current bounding box with another bounding box
   * \param [in] otherBox The other box to intersect
   * \note If the intersection is empty, the bounding box will be cleared
   * \return A reference to the bounding box after it has been intersected
   */
  BoundingBox& intersect(const BoundingBox& otherBox);

  /*!
   * \brief Expands the lower and upper bounds by the given amount.
   * \param [in] expansionAmount an absolute amount to expand
   * \note Moves min & max point expansionAmount away from the center.
   * \note This function checks to ensure the bounding box is valid afterwards.
   * \note If expansionAmount is negative, the bounding box will contract
   * \return A reference to the bounding box after it has been expanded
   */
  BoundingBox& expand(T expansionAmount);

  /*!
   * \brief Scales the bounding box about its center by a given amount.
   * \param [in] scaleFactor the multiplicative factor by which to scale
   * \note Checks to ensure that the bounding box is valid after inflation.
   * \note If scaleFactor is less than 1, the bounding box will shrink.
   * \note If scaleFactor is 0, the bounding box will shrink to its midpoint
   * \note The sign of the shrinkFactor has no effect since we are shrinking
   *  towards the center, and we fix the bounds after shrinking
   *  \return A reference to the bounding box after it has been scaled
   */
  AXOM_HOST_DEVICE BoundingBox& scale(double scaleFactor);

  /*!
   * \brief Shifts the bounding box by a fixed displacement.
   * \param [in] displacement the amount with which to move the bounding box
   * \return A reference to the bounding box after it has been shifted
   */
  BoundingBox& shift(const VectorType& displacement);

  /*!
   * \brief Checks whether the box contains the point
   * \param [in] otherPt the point that we are checking
   * \return status true if point inside the box, else false.
   *
   * \note This function assumes all intervals are closed
   * (i.e. contain their boundaries).  We may need to deal with open
   * and half open boundaries in the future.
   */
  template <typename OtherType>
  AXOM_HOST_DEVICE bool contains(const Point<OtherType, NDIMS>& otherPt) const;

  /*!
   * \brief Checks whether the box fully contains another bounding box
   * \param [in] otherBB the bounding box that we are checking
   * \return status true if bb is inside the box, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two Ts are comparable with
   *  operator<().
   */
  template <typename OtherType>
  bool contains(const BoundingBox<OtherType, NDIMS>& otherBB) const;

  /*!
   * \param [in] otherBB the bounding box that we are checking.
   * \return status true if bb intersects otherBB, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two Ts are comparable with
   *  operator<().
   */
  template <typename OtherType>
  AXOM_HOST_DEVICE bool intersectsWith(
    const BoundingBox<OtherType, NDIMS>& otherBB) const;

  /*!
   * \brief Checks that we have a valid bounding box.
   * \note A bounding box is valid when the length of each dimension is greater
   *  than or equal to zero.
   * \return status true if point inside the box, else false.
   */
  bool isValid() const;

  /*!
   * \brief Subdivides this bounding box instance into two sub-boxes by
   *  splitting along the given dimension. If a dimension is not provided, this
   *  method will split the bounding box along the longest dimension.
   * \param [in,out] right the right sub-box.
   * \param [in,out] left  the left sub-box.
   * \param [in] dimension the dimension to split along (optional)
   * \pre dimension >= -1 && dimension < NDIMS
   * \note if dimension==-1, the bounding box is split along its longest edge.
   */
  void bisect(BoxType& right, BoxType& left, int dimension = -1) const;

  /*!
   * \brief Simple formatted print of a bounding box instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /// \name Static methods
  /// @{

  /*!
   * \brief Returns the list of points of a 2-D BoundingBox instance.
   * \param [in]  bb user-supplied instance of the bounding box.
   * \param [out] pnts the list of points
   * \post pnts.size() == 4
   * \note The ordering of the points has as follows:
   * \verbatim
   *
   *     3      2
   *     +------+
   *     |      |
   *     |      |
   *     +------+
   *     0      1
   *
   * \endverbatim
   */
  static void getPoints(const BoundingBox<T, 2>& bb,
                        std::vector<Point<T, 2>>& pnts);

  /*!
   * \brief Returns the list of points of a 3-D BoundingBox instance.
   * \param [in]  bb user-supplied instance of the bounding box.
   * \param [out] pnts the list of points
   * \post pnts.size() == 8
   * \note The ordering of the points has as follows:
   * \verbatim
   *
   *         4          7
   *         +----------+
   *        /|         /|
   *     5 / |      6 / |
   *      +--|-------+  |
   *      |  +-------|--+
   *      | / 0      | / 3
   *      |/         |/
   *      +----------+
   *      1          2
   *
   * \endverbatim
   */
  static void getPoints(const BoundingBox<T, 3>& bb,
                        std::vector<Point<T, 3>>& pnts);

  /// @}

private:
  /*!
   * \brief Sets the min point for this bounding box instance.
   * \param [in] newMin the new min point.
   */
  inline void setMin(const PointType& newMin) { m_min = newMin; };

  /*!
   * \brief Sets the max point for this bounding box instance.
   * \param [in] newMax the new max point.
   */
  inline void setMax(const PointType& newMax) { m_max = newMax; };

  /*!
   * \brief Ensures that the bounds are valid.
   * A bounding box is valid when its extent in each dimension
   * (max coordinate minus min coordinate) is greater than or equal to zero
   */
  AXOM_HOST_DEVICE void checkAndFixBounds();

private:
  PointType m_min;
  PointType m_max;
};

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  BoundingBox implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------

template <typename T, int NDIMS>
constexpr T BoundingBox<T, NDIMS>::InvalidMin;

template <typename T, int NDIMS>
constexpr T BoundingBox<T, NDIMS>::InvalidMax;

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherT>
bool BoundingBox<T, NDIMS>::contains(const Point<OtherT, NDIMS>& otherPt) const
{
  for(int dim = 0; dim < NDIMS; ++dim)
  {
    if(otherPt[dim] < m_min[dim] || otherPt[dim] > m_max[dim])
    {
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BoundingBox<T, NDIMS>::BoundingBox(const PointType* pts, int n)
{
  if(n <= 0)
  {
    clear();
  }
  // sanity check:
  SLIC_ASSERT(pts != nullptr);

  this->m_min = this->m_max = pts[0];

  for(int i = 1; i < n; i++)
  {
    this->addPoint(pts[i]);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherT>
bool BoundingBox<T, NDIMS>::contains(const BoundingBox<OtherT, NDIMS>& otherBB) const
{
  return this->contains(otherBB.getMin()) && this->contains(otherBB.getMax());
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherType>
bool BoundingBox<T, NDIMS>::intersectsWith(
  const BoundingBox<OtherType, NDIMS>& otherBB) const
{
  bool status = true;

  // AABBs cannot intersect if they are separated along any dimension
  for(int i = 0; i < NDIMS; ++i)
  {
    status &= detail::intersect_bbox_bbox(m_min[i],
                                          m_max[i],
                                          otherBB.m_min[i],
                                          otherBB.m_max[i]);
  }  // END for all dimensions

  return status;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool BoundingBox<T, NDIMS>::isValid() const
{
  for(int dim = 0; dim < NDIMS; ++dim)
  {
    if(m_min[dim] > m_max[dim])
    {
      return false;
    }
  }
  return true;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherT>
void BoundingBox<T, NDIMS>::addPoint(const Point<OtherT, NDIMS>& pt)
{
  for(int dim = 0; dim < NDIMS; ++dim)
  {
    T coord = static_cast<T>(pt[dim]);

    if(coord < m_min[dim])
    {
      m_min[dim] = coord;
    }

    if(coord > m_max[dim])
    {
      m_max[dim] = coord;
    }
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
template <typename OtherT>
void BoundingBox<T, NDIMS>::addBox(const BoundingBox<OtherT, NDIMS>& bbox)
{
  this->addPoint(bbox.getMin());
  this->addPoint(bbox.getMax());
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
int BoundingBox<T, NDIMS>::getLongestDimension() const
{
  SLIC_ASSERT(this->isValid());

  int maxDim = 0;
  T max = std::numeric_limits<T>::min();
  for(int i = 0; i < NDIMS; ++i)
  {
    T dx = m_max[i] - m_min[i];
    if(dx > max)
    {
      max = dx;
      maxDim = i;
    }

  }  // END for all dimensions

  return maxDim;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BoundingBox<T, NDIMS>& BoundingBox<T, NDIMS>::expand(T expansionAmount)
{
  for(int dim = 0; dim < NDIMS; ++dim)
  {
    m_min[dim] -= expansionAmount;
    m_max[dim] += expansionAmount;
  }

  this->checkAndFixBounds();
  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BoundingBox<T, NDIMS>& BoundingBox<T, NDIMS>::scale(double scaleFactor)
{
  const PointType midpoint = getCentroid();
  const VectorType r = static_cast<T>(scaleFactor * 0.5) * range();

  m_min = PointType(midpoint.array() - r.array());
  m_max = PointType(midpoint.array() + r.array());

  this->checkAndFixBounds();

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BoundingBox<T, NDIMS>& BoundingBox<T, NDIMS>::shift(const VectorType& disp)
{
  m_min.array() += disp.array();
  m_max.array() += disp.array();

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BoundingBox<T, NDIMS>::checkAndFixBounds()
{
  for(int dim = 0; dim < NDIMS; ++dim)
  {
    if(m_min[dim] > m_max[dim])
    {
      utilities::swap(m_min[dim], m_max[dim]);
    }
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BoundingBox<T, NDIMS>::clear()
{
  m_min = PointType(std::numeric_limits<T>::max());
  m_max = PointType(std::numeric_limits<T>::lowest());
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& BoundingBox<T, NDIMS>::print(std::ostream& os) const
{
  os << "{ min:" << m_min << "; max:" << m_max << "; range:" << range() << " }";
  return os;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BoundingBox<T, NDIMS>& BoundingBox<T, NDIMS>::intersect(const BoundingBox& otherBox)
{
  for(int i = 0; i < NDIMS; ++i)
  {
    m_min[i] = std::max(m_min[i], otherBox.m_min[i]);
    m_max[i] = std::min(m_max[i], otherBox.m_max[i]);
  }

  if(!isValid())
  {
    clear();
  }

  return *this;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BoundingBox<T, NDIMS>::bisect(BoxType& right, BoxType& left, int dim) const
{
  SLIC_ASSERT(this->isValid());

  if(dim < 0)
  {
    dim = this->getLongestDimension();
  }
  SLIC_ASSERT(dim >= 0 && dim < NDIMS);

  // calculate mid along the given dimension
  T mid = 0.5 * (m_max[dim] + m_min[dim]);

  // update right
  right.setMin(this->getMin());
  PointType new_right_max = this->getMax();
  new_right_max[dim] = mid;
  right.setMax(new_right_max);
  SLIC_ASSERT(right.isValid());

  // update left
  left.setMax(this->getMax());
  PointType new_left_min = this->getMin();
  new_left_min[dim] = mid;
  left.setMin(new_left_min);
  SLIC_ASSERT(left.isValid());
}

//------------------------------------------------------------------------------
//    Implementation of static methods
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline void BoundingBox<T, NDIMS>::getPoints(const BoundingBox<T, 2>& bb,
                                             std::vector<Point<T, 2>>& pnts)
{
  pnts.resize(4);
  const Point<T, 2>& min = bb.getMin();
  const Point<T, 2>& max = bb.getMax();

  pnts[0] = Point<T, 2>::make_point(min[0], min[1]);
  pnts[1] = Point<T, 2>::make_point(max[0], min[1]);
  pnts[2] = Point<T, 2>::make_point(max[0], max[1]);
  pnts[3] = Point<T, 2>::make_point(min[0], max[1]);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline void BoundingBox<T, NDIMS>::getPoints(const BoundingBox<T, 3>& bb,
                                             std::vector<Point<T, 3>>& pnts)
{
  pnts.resize(8);
  const Point<T, 3>& min = bb.getMin();
  const Point<T, 3>& max = bb.getMax();

  pnts[0] = Point<T, 3>::make_point(min[0], min[1], min[2]);
  pnts[1] = Point<T, 3>::make_point(max[0], min[1], min[2]);
  pnts[2] = Point<T, 3>::make_point(max[0], max[1], min[2]);
  pnts[3] = Point<T, 3>::make_point(min[0], max[1], min[2]);

  pnts[4] = Point<T, 3>::make_point(min[0], min[1], max[2]);
  pnts[5] = Point<T, 3>::make_point(max[0], min[1], max[2]);
  pnts[6] = Point<T, 3>::make_point(max[0], max[1], max[2]);
  pnts[7] = Point<T, 3>::make_point(min[0], max[1], max[2]);
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template <typename T, int NDIMS>
bool operator==(const BoundingBox<T, NDIMS>& lhs, const BoundingBox<T, NDIMS>& rhs)
{
  return lhs.getMin() == rhs.getMin() && lhs.getMax() == rhs.getMax();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool operator!=(const BoundingBox<T, NDIMS>& lhs, const BoundingBox<T, NDIMS>& rhs)
{
  return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BoundingBox<T, NDIMS>& bb)
{
  return bb.print(os);
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BOUNDINGBOX_HPP_
