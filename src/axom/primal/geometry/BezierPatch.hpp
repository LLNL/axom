// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BezierPatch.hpp
 *
 * \brief A BezierPatch primitive
 */

#ifndef AXOM_PRIMAL_BEZIERPATCH_HPP_
#define AXOM_PRIMAL_BEZIERPATCH_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T>
class BezierPatch;

/*! \brief Overloaded output operator for Bezier Patches*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T>& bPatch);

/*!
 * \class BezierPatch
 *
 * \brief Represents a 3D Bezier patch defined by a 2D array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * The order of a Bezier patch with (N+1)(M+1) control points is (N, M).
 * The patch is approximated by the control points,
 * parametrized from u=0 to u=1 and v=0 to v=1.
 * 
 * Contains a 2D array of positive weights to represent a rational Bezier patch.
 * Nonrational Bezier patches are identified by an empty weights array.
 * Algorithms for Rational Bezier curves derived from 
 * Gerald Farin, "Algorithms for rational Bezier curves"
 * Computer-Aided Design, Volume 15, Number 2, 1983,
 */
template <typename T>
class BezierPatch
{
public:
  using PointType = Point<T, 3>;
  using VectorType = Vector<T, 3>;
  using PlaneType = Plane<T, 3>;

  using CoordsVec = axom::Array<PointType, 1>;
  using CoordsMat = axom::Array<PointType, 2>;
  using WeightsVec = axom::Array<T, 1>;
  using WeightsMat = axom::Array<T, 2>;

  using BoundingBoxType = BoundingBox<T, 3>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, 3>;
  using BezierCurveType = primal::BezierCurve<T, 3>;

  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A Bezier Patch must be defined using an arithmetic type");

public:
  /*!
   * \brief Constructor for a nonrational Bezier Patch that reserves 
   *  space for the given order of the surface
   *
   * Constructs an empty patch by default (no nodes/weights on either axis)
   * 
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre ord_u, ord_v greater than or equal to -1.
   */
  BezierPatch(int ord_u = -1, int ord_v = -1)
  {
    SLIC_ASSERT(ord_u >= -1 && ord_v >= -1);

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Bezier Patch from an array of coordinates
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order in both directions is greater than or equal to zero
   *
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[0]               -> nodes[0, 0],     ..., pts[ord_v]           -> nodes[0, ord_v]
   * pts[ord_v+1]         -> nodes[1, 0],     ..., pts[2*ord_v]         -> nodes[1, ord_v]
   *                                          ...
   * pts[ord_u*(ord_v-1)] -> nodes[ord_u, 0], ..., pts[ord_u*ord_v] -> nodes[ord_u, ord_v]
   * 
   */
  BezierPatch(PointType* pts, int ord_u, int ord_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord_u >= 0 && ord_v >= 0);

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);

    for(int t = 0; t < sz_u * sz_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Rational Bezier Patch from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D C-style array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order is greater than or equal to zero in each direction
   * 
   * Elements of pts and weights are mapped to control nodes (p, q) lexicographically
   * 
   * If \p weights is the null pointer, creates a nonrational curve
   */
  BezierPatch(PointType* pts, T* weights, int ord_u, int ord_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord_u >= 0 && ord_v >= 0);

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);

    for(int t = 0; t < sz_u * sz_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    if(weights == nullptr)
    {
      makeNonrational();
    }
    else
    {
      m_weights.resize(sz_u, sz_v);

      for(int t = 0; t < sz_u * sz_v; ++t)
      {
        m_weights.flatIndex(t) = weights[t];
      }
    }

    SLIC_ASSERT(isValidRational());
  }

  /*!
   * \brief Constructor for a Bezier Patch from a 1D Axom array of coordinates
   *
   * \param [in] pts A 1D Axom array of (ord_u+1)*(ord_v+1) control points
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order in both directions is greater than or equal to zero
   * 
   * Elements of pts are mapped to control nodes (p, q) lexicographically
   */
  BezierPatch(const CoordsVec& pts, int ord_u, int ord_v)
  {
    SLIC_ASSERT(ord_u >= 0 && ord_v >= 0);

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);

    for(int t = 0; t < sz_u * sz_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts.flatIndex(t);
    }

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Rational Bezier Patch from 1D Axom arrays of coordinates and weights
   *
   * \param [in] pts A 1D Axom array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D Axom array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order is greater than or equal to zero in each direction
   * 
   * Elements of pts and weights are mapped to control nodes (p, q) lexicographically.
   */
  BezierPatch(const CoordsVec& pts, const WeightsVec& weights, int ord_u, int ord_v)
  {
    SLIC_ASSERT(ord_u >= 0 && ord_v >= 0);
    SLIC_ASSERT(weights.size() == pts.size());

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);
    for(int t = 0; t < sz_u * sz_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts.flatIndex(t);
    }

    m_weights.resize(sz_u, sz_v);
    for(int t = 0; t < sz_u * sz_v; ++t)
    {
      m_weights.flatIndex(t) = weights.flatIndex(t);
    }

    SLIC_ASSERT(isValidRational());
  }

  /*!
   * \brief Constructor for a Bezier Patch from an Axom array of coordinates
   *
   * \param [in] pts A 2D Axom array with (ord_u+1, ord_v+1) control points
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order is greater than or equal to zero in each direction
   *
   */
  BezierPatch(const CoordsMat& pts, int ord_u, int ord_v)
  {
    SLIC_ASSERT((ord_u >= 0) && (ord_v >= 0));

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);
    m_controlPoints = pts;

    makeNonrational();
  }

  /*!
   * \brief Constructor for a rational Bezier Patch from a 2D Axom array 
   *   of weights and coordinates
   *
   * \param [in] pts A 2D Axom array with (ord_u+1, ord_v+1) control points
   * \param [in] weights A 2D Axom array with (ord_u+1, ord_v+1) weights
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   * \pre order is greater than or equal to zero in each direction
   */
  BezierPatch(const CoordsMat& pts, const WeightsMat& weights, int ord_u, int ord_v)
  {
    SLIC_ASSERT(ord_u >= 0 && ord_v >= 0);
    SLIC_ASSERT(pts.shape()[0] == weights.shape()[0]);
    SLIC_ASSERT(pts.shape()[1] == weights.shape()[1]);

    const int sz_u = utilities::max(0, ord_u + 1);
    const int sz_v = utilities::max(0, ord_v + 1);

    m_controlPoints.resize(sz_u, sz_v);
    m_controlPoints = pts;

    m_weights.resize(sz_u, sz_v);
    m_weights = weights;

    SLIC_ASSERT(isValidRational());
  }

  /*!
   * \brief Sets the order of Bezier patch
   *
   * \param [in] ord_u The patch's polynomial order on the first axis
   * \param [in] ord_v The patch's polynomial order on the second axis
   *
   * \note Will only resize the arrays, and likely make the patch invalid
   */
  void setOrder(int ord_u, int ord_v)
  {
    m_controlPoints.resize(ord_u + 1, ord_v + 1);
    if(isRational())
    {
      m_weights.resize(ord_u + 1, ord_v + 1);
    }
  }

  /// Returns the order of the Bezier Patch on the first axis
  int getOrder_u() const
  {
    return static_cast<int>(m_controlPoints.shape()[0]) - 1;
  }

  /// Returns the order of the Bezier Patch on the second axis
  int getOrder_v() const
  {
    return static_cast<int>(m_controlPoints.shape()[1]) - 1;
  }

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      const int ord_u = getOrder_u();
      const int ord_v = getOrder_v();

      m_weights.resize(ord_u + 1, ord_v + 1);
      m_weights.fill(1.0);
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.clear(); }

  /// Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    m_controlPoints.clear();
    makeNonrational();
  }

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int ui, int vi) { return m_controlPoints(ui, vi); }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int ui, int vi) const
  {
    return m_controlPoints(ui, vi);
  }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] ui The index of the weight on the first axis
   * \param [in] vi The index of the weight on the second axis
   * \pre Requires that the surface be rational
   */
  const T& getWeight(int ui, int vi) const
  {
    SLIC_ASSERT(isRational());
    return m_weights(ui, vi);
  }

  /*!
   * \brief Set the weight at a specific index
   *
   * \param [in] ui The index of the weight in on the first axis
   * \param [in] vi The index of the weight in on the second axis
   * \param [in] weight The updated value of the weight
   * \pre Requires that the surface be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int ui, int vi, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights(ui, vi) = weight;
  };

  /// Checks equality of two Bezier Patches
  friend inline bool operator==(const BezierPatch<T>& lhs,
                                const BezierPatch<T>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_weights == rhs.m_weights);
  }

  friend inline bool operator!=(const BezierPatch<T>& lhs,
                                const BezierPatch<T>& rhs)
  {
    return !(lhs == rhs);
  }

  /// Returns a copy of the Bezier patch's control points
  CoordsMat getControlPoints() const { return m_controlPoints; }

  /// Returns a copy of the Bezier patch's weights
  WeightsMat getWeights() const { return m_weights; }

  /*!
   * \brief Reverses the order of one direction of the Bezier patch's control points and weights
   *
   * \param [in] axis orientation of patch. 0 to reverse in u, 1 for reverse in v
   */
  void reverseOrientation(int axis)
  {
    if(axis == 0)
    {
      reverseOrientation_u();
    }
    else
    {
      reverseOrientation_v();
    }
  }

  /// Reverses the order of the Bezier patch's control points and weights on the first axis
  void reverseOrientation_u()
  {
    const int ord_u = getOrder_u();
    const int mid_u = (ord_u + 1) / 2;

    const int ord_v = getOrder_v();

    for(int q = 0; q <= ord_v; ++q)
    {
      for(int i = 0; i < mid_u; ++i)
      {
        axom::utilities::swap(m_controlPoints(i, q),
                              m_controlPoints(ord_u - i, q));
      }

      if(isRational())
      {
        for(int i = 0; i < mid_u; ++i)
        {
          axom::utilities::swap(m_weights(i, q), m_weights(ord_u - i, q));
        }
      }
    }
  }

  /// Reverses the order of the Bezier patch's control points and weights on the second axis
  void reverseOrientation_v()
  {
    const int ord_u = getOrder_u();

    const int ord_v = getOrder_v();
    const int mid_v = (ord_v + 1) / 2;

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int i = 0; i < mid_v; ++i)
      {
        axom::utilities::swap(m_controlPoints(p, i),
                              m_controlPoints(p, ord_v - i));
      }

      if(isRational())
      {
        for(int i = 0; i < mid_v; ++i)
        {
          axom::utilities::swap(m_weights(p, i), m_weights(p, ord_v - i));
        }
      }
    }
  }

  /// Swap the axes such that s(u, v) becomes s(v, u)
  void swapAxes()
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    CoordsMat new_controlPoints(ord_v + 1, ord_u + 1);

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        new_controlPoints(q, p) = m_controlPoints(p, q);
      }
    }

    m_controlPoints = new_controlPoints;

    if(isRational())
    {
      WeightsMat new_weights(ord_v + 1, ord_u + 1);
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          new_weights(q, p) = m_weights(p, q);
        }
      }

      m_weights = new_weights;
    }
  }

  /// Returns an axis-aligned bounding box containing the Bezier patch
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /// Returns an oriented bounding box containing the Bezier patch
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(),
                                   static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Evaluates a slice Bezier patch for a fixed parameter value of \a u or \a v
   *
   * \param [in] u parameter value at which to evaluate the first axis
   * \param [in] v parameter value at which to evaluate the second axis
   * \param [in] axis orientation of curve. 0 for fixed u, 1 for fixed v
   * \return p the value of the Bezier patch at (u, v)
   *
   * \note We typically evaluate the patch at \a u or \a v between 0 and 1
   */
  BezierCurveType isocurve(T uv, int axis) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    if(axis == 0)
    {
      return isocurve_u(uv);
    }
    else
    {
      return isocurve_v(uv);
    }
  }

  /// Return an isocurve for a fixed value of u
  BezierCurveType isocurve_u(T u) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    BezierCurveType c(ord_v);
    axom::Array<T> dCarray(ord_u + 1);

    if(isRational())
    {
      c.makeRational();
      // If looking for an edge, don't need to use De Casteljau
      if(u == 0.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(0, q);
          c.setWeight(q, m_weights(0, q));
        }
        return c;
      }
      if(u == 1.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(ord_u, q);
          c.setWeight(q, m_weights(ord_u, q));
        }
        return c;
      }

      // Otherwise, do De Casteljau along the v-axis
      axom::Array<T> dWarray(ord_u + 1);
      for(int q = 0; q <= ord_v; ++q)
      {
        for(int i = 0; i < 3; ++i)
        {
          for(int p = 0; p <= ord_u; ++p)
          {
            dCarray[p] = m_controlPoints(p, q)[i] * m_weights(p, q);
            dWarray[p] = m_weights(p, q);
          }

          for(int p = 1; p <= ord_u; ++p)
          {
            const int end = ord_u - p;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], u);
              dWarray[k] = lerp(dWarray[k], dWarray[k + 1], u);
            }
          }

          c[q][i] = dCarray[0] / dWarray[0];
          c.setWeight(q, dWarray[0]);
        }
      }
    }
    else
    {
      // If looking for an edge, don't need to use De Casteljau
      if(u == 0.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(0, q);
        }
        return c;
      }
      if(u == 1.0)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          c[q] = m_controlPoints(ord_u, q);
        }
        return c;
      }

      // Otherwise, do De Casteljau along the v-axis
      for(int q = 0; q <= ord_v; ++q)
      {
        for(int i = 0; i < 3; ++i)
        {
          for(int p = 0; p <= ord_u; ++p)
          {
            dCarray[p] = m_controlPoints(p, q)[i];
          }

          for(int p = 1; p <= ord_u; ++p)
          {
            const int end = ord_u - p;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], u);
            }
          }
          c[q][i] = dCarray[0];
        }
      }
    }

    return c;
  }

  /// Return an isocurve for a fixed value of v
  BezierCurveType isocurve_v(T v) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    BezierCurveType c(ord_u);
    axom::Array<T> dCarray(ord_v + 1);

    if(isRational())
    {
      c.makeRational();

      // If looking for an edge, don't need to use De Casteljau
      if(v == 0.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, 0);
          c.setWeight(p, m_weights(p, 0));
        }
        return c;
      }
      if(v == 1.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, ord_v);
          c.setWeight(p, m_weights(p, ord_v));
        }
        return c;
      }

      // Otherwise, do De Casteljau along the u-axis
      axom::Array<T> dWarray(ord_v + 1);
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int i = 0; i < 3; ++i)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCarray[q] = m_controlPoints(p, q)[i] * m_weights(p, q);
            dWarray[q] = m_weights(p, q);
          }

          for(int q = 1; q <= ord_v; ++q)
          {
            const int end = ord_v - q;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], v);
              dWarray[k] = lerp(dWarray[k], dWarray[k + 1], v);
            }
          }
          c[p][i] = dCarray[0] / dWarray[0];
          c.setWeight(p, dWarray[0]);
        }
      }
    }
    else
    {
      // If looking for an edge, don't need to use De Casteljau
      if(v == 0.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, 0);
        }
        return c;
      }
      if(v == 1.0)
      {
        for(int p = 0; p <= ord_u; ++p)
        {
          c[p] = m_controlPoints(p, ord_v);
        }
        return c;
      }

      // Otherwise, do De Casteljau along the u-axis
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int i = 0; i < 3; ++i)
        {
          for(int q = 0; q <= ord_v; ++q)
          {
            dCarray[q] = m_controlPoints(p, q)[i];
          }

          for(int q = 1; q <= ord_v; ++q)
          {
            const int end = ord_v - q;
            for(int k = 0; k <= end; ++k)
            {
              dCarray[k] = lerp(dCarray[k], dCarray[k + 1], v);
            }
          }
          c[p][i] = dCarray[0];
        }
      }
    }

    return c;
  }

  /*!
   * \brief Evaluates a Bezier patch at a particular parameter value (\a u, \a v)
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return p the value of the Bezier patch at (u, v)
   *
   * \note We typically evaluate the patch at \a u and \a v between 0 and 1
   */
  PointType evaluate(T u, T v) const
  {
    if(getOrder_u() >= getOrder_v())
    {
      return isocurve_u(u).evaluate(v);
    }
    else
    {
      return isocurve_v(v).evaluate(u);
    }
  }

  /*!
   * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the u axis
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec a tangent vector of the Bezier patch at (u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType du(T u, T v) const { return isocurve_v(v).dt(u); }

  /*!
   * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the v axis
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec a tangent vector of the Bezier patch at (u, v)
   *
   * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
   */
  VectorType dv(T u, T v) const { return isocurve_u(u).dt(v); }

  /*!
   * \brief Computes the normal vector of a Bezier patch at a particular parameter value (\a u, \a v)
   *
   * \param [in] u parameter value at which to evaluate on the first axis
   * \param [in] v parameter value at which to evaluate on the second axis
   * \return vec the normal vector of the Bezier patch at (u, v)
   *
   * \note We typically find the normal of the patch at \a u and \a v between 0 and 1
   */
  VectorType normal(T u, T v) const
  {
    return VectorType::cross_product(du(u, v), dv(u, v));
  }

  /*!
   * \brief Splits a Bezier patch into two Bezier patches
   *
   * \param [in] uv parameter value between 0 and 1 at which to bisect the patch
   * \param [in] axis orientation of split. 0 for fixed u, 1 for fixed v
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   *
   * \pre Parameter \a uv must be between 0 and 1
   */
  void split(T uv, int axis, BezierPatch& p1, BezierPatch& p2) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    if(axis == 0)
    {
      split_u(uv, p1, p2);
    }
    else
    {
      split_v(uv, p1, p2);
    }
  }

  /// Split the patch along a fixed value of u
  void split_u(T u, BezierPatch& p1, BezierPatch& p2) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    p1.setOrder(ord_u, ord_v);

    // Note: The second patch's control points are computed inline
    //       as we find the first patch's control points
    p2 = *this;

    if(isRational())
    {
      p1.makeRational();  // p2 already rational

      for(int q = 0; q <= ord_v; ++q)
      {
        // Do the rational de Casteljau algorithm
        p1(0, q) = p2(0, q);
        p1.setWeight(0, q, p2.getWeight(0, q));

        for(int p = 1; p <= ord_u; ++p)
        {
          const int end = ord_u - p;

          for(int k = 0; k <= end; ++k)
          {
            double temp_weight =
              lerp(p2.getWeight(k, q), p2.getWeight(k + 1, q), u);

            for(int i = 0; i < 3; ++i)
            {
              p2(k, q)[i] = lerp(p2.getWeight(k, q) * p2(k, q)[i],
                                 p2.getWeight(k + 1, q) * p2(k + 1, q)[i],
                                 u) /
                temp_weight;
            }

            p2.setWeight(k, q, temp_weight);
          }

          p1(p, q) = p2(0, q);
          p1.setWeight(p, q, p2.getWeight(0, q));
        }
      }
    }
    else
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        p1(0, q) = m_controlPoints(0, q);

        for(int i = 0; i < 3; ++i)
        {
          for(int p = 1; p <= ord_u; ++p)
          {
            const int end = ord_u - p;
            for(int k = 0; k <= end; ++k)
            {
              p2(k, q)[i] = lerp(p2(k, q)[i], p2(k + 1, q)[i], u);
            }

            p1(p, q)[i] = p2(0, q)[i];
          }
        }
      }
    }
  }

  void split_v(T v, BezierPatch& p1, BezierPatch& p2) const
  {
    using axom::utilities::lerp;

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    p1.setOrder(ord_u, ord_v);

    // Note: The second patch's control points are computed inline
    //       as we find the first patch's control points
    p2 = *this;

    if(isRational())
    {
      p1.makeRational();  // p2 already rational

      for(int p = 0; p <= ord_u; ++p)
      {
        // Do the rational de Casteljau algorithm
        p1(p, 0) = p2(p, 0);
        p1.setWeight(p, 0, p2.getWeight(p, 0));

        for(int q = 1; q <= ord_v; ++q)
        {
          const int end = ord_v - q;

          for(int k = 0; k <= end; ++k)
          {
            double temp_weight =
              lerp(p2.getWeight(p, k), p2.getWeight(p, k + 1), v);

            for(int i = 0; i < 3; ++i)
            {
              p2(p, k)[i] = lerp(p2.getWeight(p, k) * p2(p, k)[i],
                                 p2.getWeight(p, k + 1) * p2(p, k + 1)[i],
                                 v) /
                temp_weight;
            }

            p2.setWeight(p, k, temp_weight);
          }

          p1(p, q) = p2(p, 0);
          p1.setWeight(p, q, p2.getWeight(p, 0));
        }
      }
    }
    else
    {
      for(int p = 0; p <= ord_u; ++p)
      {
        p1(p, 0) = m_controlPoints(p, 0);

        for(int i = 0; i < 3; ++i)
        {
          for(int q = 1; q <= ord_v; ++q)
          {
            const int end = ord_v - q;
            for(int k = 0; k <= end; ++k)
            {
              p2(p, k)[i] = lerp(p2(p, k)[i], p2(p, k + 1)[i], v);
            }

            p1(p, q)[i] = p2(p, 0)[i];
          }
        }
      }
    }
  }

  /*!
   * \brief Splits a Bezier patch into four Bezier patches
   *
   * \param [in] u parameter value between 0 and 1 at which to bisect on the first axis
   * \param [in] v parameter value between 0 and 1 at which to bisect on the second axis
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   * \param [out] p3 Third output Bezier patch
   * \param [out] p4 Fourth output Bezier patch
   *
   *   v = 1
   *   ----------------------
   *   |         |          |
   *   |   p3    |    p4    |
   *   |         |          |
   *   --------(u,v)---------
   *   |         |          |
   *   |   p1    |    p2    |
   *   |         |          |
   *   ---------------------- u = 1
   *
   * \pre Parameter \a u and \a v must be between 0 and 1
   */
  void split(T u,
             T v,
             BezierPatch& p1,
             BezierPatch& p2,
             BezierPatch& p3,
             BezierPatch& p4) const
  {
    // Bisect the patch along the u direction
    split_u(u, p1, p2);

    // Temporarily store the result in each half and split again
    BezierPatch p0(p1);
    p0.split_v(v, p1, p3);

    p0 = p2;
    p0.split_v(v, p2, p4);
  }

  /*!
   * \brief Predicate to check if the Bezier patch is approximately planar
   *
   * This function checks if all control points of the BezierPatch
   * are approximately on the plane defined by its four corners
   *
   * \param [in] tol Threshold for sum of squared distances
   * \return True if c1 is near-planar
   */
  bool isPlanar(double tol = 1E-8) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(ord_u <= 0 && ord_v <= 0)
    {
      return true;
    }
    if(ord_u == 1 && ord_v == 0)
    {
      return true;
    }
    if(ord_u == 0 && ord_v == 1)
    {
      return true;
    }

    // Check that the four corners aren't coplanar
    VectorType v1(m_controlPoints(0, 0), m_controlPoints(0, ord_v));
    VectorType v2(m_controlPoints(0, 0), m_controlPoints(ord_u, 0));
    VectorType v3(m_controlPoints(0, 0), m_controlPoints(ord_u, ord_v));
    if(!axom::utilities::isNearlyEqual(
         VectorType::scalar_triple_product(v1, v2, v3),
         0.0,
         tol))
      return false;

    // Find three points that produce a nonzero normal
    Vector3D plane_normal = VectorType::cross_product(v1, v2);
    if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, tol))
      plane_normal = VectorType::cross_product(v1, v3);
    if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, tol))
      plane_normal = VectorType::cross_product(v2, v3);
    plane_normal = plane_normal.unitVector();

    double sqDist = 0.0;

    // Check all control points for simplicity
    for(int p = 0; p <= ord_u && sqDist <= tol; ++p)
    {
      for(int q = ((p == 0) ? 1 : 0); q <= ord_v && sqDist <= tol; ++q)
      {
        double signedDist = VectorType::dot_product(
          plane_normal,
          Vector3D(m_controlPoints(0, 0), m_controlPoints(p, q)));
        sqDist += signedDist * signedDist;
      }
    }

    return (sqDist <= tol);
  }

  /*!
   * \brief Predicate to check if the Bezier patch is approximately planar-polygonal
   *
   * This function checks if the edges of a BezierPatch are approximately linear
   *
   * \param [in] tol Threshold for sum of squared distances
   * \return True if c1 is near-planar-polygonal
   */
  bool isPolygonal(double tol = 1E-8) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(ord_u <= 0 && ord_v <= 0) return true;
    if(ord_u == 1 && ord_v == 0) return true;
    if(ord_u == 0 && ord_v == 1) return true;

    // Check if the patch is planar
    if(!isPlanar(tol)) return false;

    // Check if each bounding curve is linear
    if(!isocurve_u(0).isLinear(tol)) return false;
    if(!isocurve_v(0).isLinear(tol)) return false;
    if(!isocurve_u(1).isLinear(tol)) return false;
    if(!isocurve_v(1).isLinear(tol)) return false;

    return true;
  }

  /*!
   * \brief Simple formatted print of a Bezier Patch instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    os << "{ order (" << ord_u << ',' << ord_v << ") Bezier Patch ";

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        os << m_controlPoints(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
      }
    }

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p <= ord_u; ++p)
      {
        for(int q = 0; q <= ord_v; ++q)
        {
          os << m_weights(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
        }
      }
    }
    os << "}";

    return os;
  }

private:
  /// Check that the weights used are positive, and
  ///  that there is one for each control node
  bool isValidRational() const
  {
    if(!isRational())
    {
      return true;
    }

    const int ord_u = getOrder_u();
    const int ord_v = getOrder_v();

    if(m_weights.shape()[0] != (ord_u + 1) || m_weights.shape()[1] != (ord_v + 1))
    {
      return false;
    }

    for(int p = 0; p <= ord_u; ++p)
    {
      for(int q = 0; q <= ord_v; ++q)
      {
        if(m_weights(p, q) <= 0)
        {
          return false;
        }
      }
    }

    return true;
  }

  CoordsMat m_controlPoints;
  WeightsMat m_weights;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierPatch
//------------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T>& bPatch)
{
  bPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BEZIERPATCH_HPP_
