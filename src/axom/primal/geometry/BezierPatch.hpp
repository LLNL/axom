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

#include <vector>
#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class BezierPatch;

/*! \brief Overloaded output operator for Bezier Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T, NDIMS>& bCurve);

/*!
 * \class BezierPatch
 *
 * \brief Represents a 3D Bezier patch defined by a 2D array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * The order of a Bezier patch with (N+1)(M+1) control points is (N, M).
 * The patch is approximated by the control points,
 * parametrized from t=0 to t=1 and s=0 to s=1.
 * 
 * Contains a 2D array of positive weights to represent a rational Bezier patch.
 * Nonrational Bezier patches are identified by an empty weights array.
 * Algorithms for Rational Bezier curves derived from 
 * < Not yet known. Probably the same Farin book? >
 */
template <typename T, int NDIMS = 3>
class BezierPatch
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = NumericArray<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;
  using CoordsVec = axom::Array<PointType>;
  using CoordsMat = axom::Array<PointType, 2>;
  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using BezierCurveType = primal::BezierCurve<T, NDIMS>;

  // TODO: Unsure if this behavior is captured by the templating statement above
  AXOM_STATIC_ASSERT_MSG(NDIMS == 3,
                         "A Bezier Patch object must be defined in 3-D");
  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A Bezier Curve must be defined using an arithmetic type");

public:
  /*!
   * \brief Constructor for a Bezier Patch that reserves space for
   *  the given order of the surface
   *
   * \param [in] order the order of the resulting Bezier curve
   * \pre order is greater than or equal to -1.
   */
  BezierPatch(int ord_t = -1, int ord_s = -1)
  {
    SLIC_ASSERT(ord_t >= -1, ord_s >= -1);

    const int sz_t = utilities::max(0, ord_t + 1);
    const int sz_s = utilities::max(0, ord_s + 1);

    m_controlPoints.resize(sz_t, sz_s);

    makeNonrational();
  }

  /*!
   * TODO: Make this constructor work 
   * (or more likely, don't. That's so many points to initialize in an array...)
   * \brief Constructor for an order \a ord= n Bezier curve
   * from a list of coordinates:
   * \verbatim {x_0, x_1, x_2,...,x_n,
   *            y_0, y_1, y_2,...,y_n,
   *            z_0, z_1, z_2,...,z_n}
   *
   * \param [in] pts an array with (p+1)*NDIMS entries, ordered by coordinate
   * then by polynomial order
   * \param [in] ord Polynomial order of the curve
   * \pre order is greater than or equal to zero
   */
  //BezierCurve(T* pts, int ord)
  //{
  //  SLIC_ASSERT(pts != nullptr);
  //  SLIC_ASSERT(ord >= 0);

  //  const int sz = utilities::max(0, ord + 1);
  //  m_controlPoints.resize(sz);

  //  for(int p = 0; p <= ord; p++)
  //  {
  //    auto& pt = m_controlPoints[p];
  //    for(int j = 0; j < NDIMS; j++)
  //    {
  //      pt[j] = pts[j * (ord + 1) + p];
  //    }
  //  }

  //  makeNonrational();
  //}

  /*!
   * \brief Constructor for a Bezier Patch from a 2D array of coordinates
   *
   * \param [in] pts an array with ord_t+1 arrays with ord_m+1 control points
   * \param [in] ord_t The patches's polynomial order in p
   * \param [in] ord_s The patches's polynomial order in q
   * \pre order in both directions is greater than or equal to zero
   *
   */
  BezierPatch(PointType* pts, int ord_t, int ord_s)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord_t >= 0, ord_s >= 0);

    const int sz_t = utilities::max(0, ord_t + 1);
    const int sz_s = utilities::max(0, ord_s + 1);

    m_controlPoints.resize(sz_t, sz_s);

    for(int t = 0; t < sz_t * sz_s; ++t)
      m_controlPoints(t / sz_s, t % sz_s) = pts[t];

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Rational Bezier Path from an array of coordinates and weights
   *
   * \param [in] pts a matrix with (ord_t+1) * (ord_s+1) control points
   * \param [in] weights a matrix with (ord_t+1) * (ord_s+1) positive weights
   * \param [in] ord_t The patches's polynomial order in p
   * \param [in] ord_s The patches's polynomial order in q
   * \pre order is greater than or equal to zero in each direction
   *
   */
  BezierPatch(PointType* pts, T* weights, int ord_t, int ord_s)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord_t >= 0, ord_s >= 0);

    const int sz_t = utilities::max(0, ord_t + 1);
    const int sz_s = utilities::max(0, ord_s + 1);

    m_controlPoints.resize(sz_t, sz_s);

    for(int t = 0; t < sz_t * sz_s; ++t)
      m_controlPoints(t / sz_s, t % sz_s) = pts[t];

    if(weights == nullptr)
    {
      m_weights.resize(0, 0);
    }
    else
    {
      m_weights.resize(sz_t, sz_s);

      for(int t = 0; t < sz_t * sz_s; ++t)
        m_weights(t / sz_s, t % sz_s) = weights[t];

      SLIC_ASSERT(isValidRational());
    }
  }

  /*!
   * \brief Constructor for a Bezier Patch from an matrix of coordinates
   *
   * \param [in] pts a vector with (ord_t+1) * (ord_s+1) control points
   * \param [in] ord_t The patch's polynomial order in p
   * \param [in] ord_s The patch's polynomial order in q
   * \pre order is greater than or equal to zero in each direction
   *
   */
  BezierPatch(const CoordsMat& pts, int ord_t, int ord_s)
  {
    SLIC_ASSERT((ord_t >= 0) && (ord_s >= 0));

    const int sz_t = utilities::max(0, ord_t + 1);
    const int sz_s = utilities::max(0, ord_s + 1);

    m_controlPoints.resize(sz_t, sz_s);
    m_controlPoints = pts;

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Rational Bezier Patch from a matrix
   * of coordinates and weights
   *
   * \param [in] pts a vector with (ord_t+1) * (ord_s+1) control points
   * \param [in] ord_t The patch's polynomial order in p
   * \param [in] ord_s The patch's polynomial order in q
   * \pre order is greater than or equal to zero in each direction
   *
   * TODO: Can we just do this entire thing with initalizer lists? 
   */
  BezierPatch(const CoordsMat& pts,
              const axom::Array<axom::Array<T>>& weights,
              int ord_t,
              int ord_s)
  {
    SLIC_ASSERT(ord >= 0);
    SLIC_ASSERT(pts.shape()[0] == weights.shape()[0]);
    SLIC_ASSERT(pts.shape()[1] == weights.shape()[1]);

    const int sz_t = utilities::max(0, ord_t + 1);
    const int sz_s = utilities::max(0, ord_s + 1);

    m_controlPoints.resize(sz_t, sz_s);
    m_controlPoints = pts;

    m_weights.resize(sz_t, sz_s);
    m_weights = weights;

    SLIC_ASSERT(isValidRational());
  }

  /// Sets the order of the Bezier Curve
  void setOrder(int ord_t, int ord_s)
  {
    m_controlPoints.resize(ord_t + 1, ord_s + 1);
  }

  /// Returns the order of the Bezier Curve on the first axis
  int getOrder_t() const
  {
    return static_cast<int>(m_controlPoints.shape()[0]) - 1;
  }

  /// Returns the order of the Bezier Curve
  int getOrder_s() const
  {
    return static_cast<int>(m_controlPoints.shape()[1]) - 1;
  }

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      const int ord_t = getOrder_t();
      const int ord_s = getOrder_s();

      m_weights.resize(ord_t + 1, ord_s + 1);

      for(int i = 0; i <= ord_t; ++i)
        for(int j = 0; j <= ord_s; ++j) m_weights[i][j] = 1.0;
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.resize(0, 0); }

  /// Use array size as flag for rationality
  bool isRational() const { return (m_weights.size() != 0); }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    for(int p = 0; p <= ord; ++p)
      for(int q = 0; q <= getOrder_s; ++q) m_controlPoints(p, q) = PointType();

    makeNonrational();
  }

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int idx_p, int idx_q)
  {
    return m_controlPoints(idx_p, idx_q);
  }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int idx_p, int idx_q) const
  {
    return m_controlPoints(idx_p, idx_q);
  }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] idx_p The index of the weight in p
   * \param [in] idx_q The index of the weight in q
   * \pre Requires that the surface be rational
   */
  const T& getWeight(int idx_p, int idx_q) const
  {
    SLIC_ASSERT(isRational());
    return m_weights(idx_p, idx_q);
  }

  /*!
   * \brief Set the weight at a specific index
   *
   * \param [in] idx_p The index of the weight in p
   * \param [in] idx_q The index of the weight in q
   * \param [in] weight The updated value of the weight
   * \pre Requires that the surface be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int idx_p, int idx_q, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights(idx_p, idx_q) = weight;
  };

  /// Checks equality of two Bezier Patches
  friend inline bool operator==(const BezierPatch<T, NDIMS>& lhs,
                                const BezierPatch<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_weights == rhs.m_weights);
  }

  friend inline bool operator!=(const BezierPatch<T, NDIMS>& lhs,
                                const BezierPatch<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /// Returns a copy of the Bezier curve's control points
  CoordsMat getControlPoints() const { return m_controlPoints; }

  /// Returns a copy of the Bezier curve's weights
  axom::Array<T, 2> getWeights() const { return m_weights; }

  /// Reverses the order of the Bezier curve's control points and weights
  // void reverseOrientation() {}

  /// Returns an axis-aligned bounding box containing the Bezier curve
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /// Returns an oriented bounding box containing the Bezier curve
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(),
                                   static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Evaluates a slice Bezier patch for a fixed parameter value of \a t or \a s
   *
   * \param [in] t parameter value at which to evaluate
   * \param [in] s parameter value at which to evaluate
   * \param [in] axis orientation of curve. 0 for fixed t, 1 for fixed s
   * \return p the value of the Bezier patch at (t, s)
   *
   * \note We typically evaluate the curve at \a t or \a s between 0 and 1
   */
  BezierCurveType slice(T ts, int axis = 0) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    using axom::utilities::lerp;

    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    BezierCurveType c(-1);

    if(isRational())
    {
      if(axis == 0)  // Keeping a fixed value of t
      {
        c.setOrder(ord_s);
        c.makeRational();

        axom::Array<T> dCarray(ord_t + 1);
        axom::Array<T> dWarray(ord_t + 1);

        // Run de Casteljau algorithm on each row of control nodes and each dimension
        for(int q = 0; q <= ord_s; ++q)
          for(int i = 0; i < NDIMS; ++i)
          {
            for(int p = 0; p <= ord_t; ++p)
            {
              dCarray[p] = m_controlPoints(p, q)[i] * m_weights(p, q);
              dWarray[p] = m_weights(p, q);
            }

            for(int p = 1; p <= ord_t; ++p)
            {
              const int end = ord_t - p;
              for(int k = 0; k <= end; ++k)
              {
                dCarray[k] = lerp(dCarray[k], dCarray[k + 1], ts);
                dWarray[k] = lerp(dWarray[k], dWarray[k + 1], ts);
              }
            }

            c[q][i] = dCarray[0] / dWarray[0];
            c.setWeight(q, dWarray[0]);
          }
      }
      // Run de Casteljau algorithm on each column of control nodes and each dimension
      else  // Keeping a fixed value of s
      {
        c.setOrder(ord_t);
        c.makeRational();

        axom::Array<T> dCarray(ord_s + 1);
        axom::Array<T> dWarray(ord_s + 1);

        for(int p = 0; p <= ord_t; ++p)
          for(int i = 0; i < NDIMS; ++i)
          {
            for(int q = 0; q <= ord_s; ++q)
            {
              dCarray[q] = m_controlPoints(p, q)[i] * m_weights(p, q);
              dWarray[q] = m_weights(p, q);
            }

            for(int q = 1; q <= ord_s; ++q)
            {
              const int end = ord_s - q;
              for(int k = 0; k <= end; ++k)
              {
                dCarray[k] = lerp(dCarray[k], dCarray[k + 1], ts);
                dWarray[k] = lerp(dWarray[k], dWarray[k + 1], ts);
              }
            }
            c[p][i] = dCarray[0] / dWarray[0];
            c.setWeight(p, dWarray[0]);
          }
      }
    }
    else
    {
      if(axis == 0)  // Keeping a fixed value of t
      {
        c.setOrder(ord_s);
        axom::Array<T> dCarray(ord_t + 1);  // Temp array

        // Run de Casteljau algorithm on each row of control nodes and each dimension
        for(int q = 0; q <= ord_s; ++q)
          for(int i = 0; i < NDIMS; ++i)
          {
            for(int p = 0; p <= ord_t; ++p)
              dCarray[p] = m_controlPoints(p, q)[i];

            for(int p = 1; p <= ord_t; ++p)
            {
              const int end = ord_t - p;
              for(int k = 0; k <= end; ++k)
                dCarray[k] = lerp(dCarray[k], dCarray[k + 1], ts);
            }
            c[q][i] = dCarray[0];
          }
      }
      // Run de Casteljau algorithm on each column of control nodes and each dimension
      else  // Keeping a fixed value of s
      {
        c.setOrder(ord_t);
        axom::Array<T> dCarray(ord_s + 1);  // Temp array

        for(int p = 0; p <= ord_t; ++p)
          for(int i = 0; i < NDIMS; ++i)
          {
            for(int q = 0; q <= ord_s; ++q)
              dCarray[q] = m_controlPoints(p, q)[i];

            for(int q = 1; q <= ord_s; ++q)
            {
              const int end = ord_s - q;
              for(int k = 0; k <= end; ++k)
              {
                dCarray[k] = lerp(dCarray[k], dCarray[k + 1], ts);
              }
            }
            c[p][i] = dCarray[0];
          }
      }
    }

    return c;
  }

  /*!
   * \brief Evaluates a Bezier patch at a particular parameter value (\a t, \a s)
   *
   * \param [in] t parameter value at which to evaluate
   * \param [in] s parameter value at which to evaluate
   * \return p the value of the Bezier patch at (t, s)
   *
   * \note We typically evaluate the curve at \a t and \a s between 0 and 1
   */
  PointType evaluate(T t, T s) const { return slice(t).evaluate(s); }

  /*!
   * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a t, \a s) oriented in \a axis
   *
   * \param [in] t parameter value at which to evaluate
   * \param [in] s parameter value at which to evaluate
   * \param [in] axis orientation of vector. 0 for fixed t, 1 for fixed s
   * \return v a tangent vector of the Bezier curve at (t, s)
   *
   * \note We typically find the tangent of the curve at \a t and \a s between 0 and 1
   */
  VectorType dt(T t, T s, int axis) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));
    if(axis == 0)  // Get slice at fixed s
      return slice(s, 1).dt(t);
    else  // Get slice at fixed t
      return slice(t, 0).dt(s);
  }

  /*!
   * \brief Computes the normal vector of a Bezier patch at a particular parameter value (\a t, \a s)
   *
   * \param [in] t parameter value at which to evaluate
   * \param [in] s parameter value at which to evaluate
   * \return v the normal vector of the Bezier curve at (t, s)
   *
   * \note We typically find the normal of the curve at \a t and \a s between 0 and 1
   */
  VectorType normal(T t, T s) const
  {
    VectorType tangent_t = dt(t, s, 0);
    VectorType tangent_s = dt(t, s, 1);
    return VectorType::cross_product(tangent_t, tangent_s);
  }

  /*!
   * \brief Splits a Bezier patch into two Bezier patches
   *
   * \param [in] ts parameter value between 0 and 1 at which to bisect the patch
   * \param [in] axis orientation of split. 0 for fixed t, 1 for fixed s
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   *
   * \pre Parameter \a ts must be between 0 and 1
   */
  void split(T ts, int axis, BezierPatch& p1, BezierPatch& p2) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    using axom::utilities::lerp;

    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    SLIC_ASSERT((ord_t >= 0) && (ord_s >= 0));

    // Note: The second patch's control points are computed inline
    //       as we find the first patch's control points
    p2 = *this;

    p1.setOrder(ord_t, ord_s);
    if(isRational())
    {
      p1.makeRational();  // p2 already rational
      if(axis == 0)       // Split across a fixed value of t
      {
        // Run algorithm across each row of control nodes
        for(int q = 0; q <= ord_s; ++q)
        {
          // Do the rational de Casteljau algorithm
          p1(0, q) = p2(0, q);
          p1.setWeight(0, q, p2.getWeight(0, q));

          for(int p = 1; p <= ord_t; ++p)
          {
            const int end = ord_t - p;

            for(int k = 0; k <= end; ++k)
            {
              double temp_weight =
                lerp(p2.getWeight(k, q), p2.getWeight(k + 1, q), ts);

              for(int i = 0; i < NDIMS; ++i)
              {
                p2(k, q)[i] = lerp(p2.getWeight(k, q) * p2(k, q)[i],
                                   p2.getWeight(k + 1, q) * p2(k + 1, q)[i],
                                   ts) /
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
        // Run algorithm across each column of control nodes
        for(int p = 0; p <= ord_t; ++p)
        {
          // Do the rational de Casteljau algorithm
          p1(p, 0) = p2(p, 0);
          p1.setWeight(p, 0, p2.getWeight(p, 0));

          for(int q = 1; q <= ord_s; ++q)
          {
            const int end = ord_s - q;

            for(int k = 0; k <= end; ++k)
            {
              double temp_weight =
                lerp(p2.getWeight(p, k), p2.getWeight(p, k + 1), ts);

              for(int i = 0; i < NDIMS; ++i)
              {
                p2(p, k)[i] = lerp(p2.getWeight(p, k) * p2(p, k)[i],
                                   p2.getWeight(p, k + 1) * p2(p, k + 1)[i],
                                   ts) /
                  temp_weight;
              }

              p2.setWeight(p, k, temp_weight);
            }

            p1(p, q) = p2(p, 0);
            p1.setWeight(p, q, p2.getWeight(p, 0));
          }
        }
      }
    }
    else
    {
      if(axis == 0)  // Split across a fixed value of t
      {
        // Do the split for each row of control nodes
        for(int q = 0; q <= ord_s; ++q)
        {
          p1(0, q) = m_controlPoints(0, q);

          for(int i = 0; i < NDIMS; ++i)
            for(int p = 1; p <= ord_t; ++p)
            {
              const int end = ord_t - p;
              for(int k = 0; k <= end; ++k)
                p2(k, q)[i] = lerp(p2(k, q)[i], p2(k + 1, q)[i], ts);

              p1(p, q)[i] = p2(0, q)[i];
            }
        }
      }
      else
      {
        // Do the split for each column of control nodes
        for(int p = 0; p <= ord_t; ++p)
        {
          p1(p, 0) = m_controlPoints(p, 0);

          for(int i = 0; i < NDIMS; ++i)
            for(int q = 1; q <= ord_s; ++q)
            {
              const int end = ord_s - q;
              for(int k = 0; k <= end; ++k)
                p2(p, k)[i] = lerp(p2(p, k)[i], p2(p, k + 1)[i], ts);

              p1(p, q)[i] = p2(p, 0)[i];
            }
        }
      }
    }
  }

  /*!
   * \brief Splits a Bezier patch into four Bezier patches
   *
   * \param [in] t parameter value between 0 and 1 at which to bisect the patch in t
   * \param [in] s parameter value between 0 and 1 at which to bisect the patch in s
   * \param [out] p1 First output Bezier patch
   * \param [out] p2 Second output Bezier patch
   * \param [out] p3 Third output Bezier patch
   * \param [out] p4 Fourth output Bezier patch
   *
   *   s = 1
   *   ----------------------
   *   |         |          |
   *   |   p3    |    p4    |
   *   |         |          |
   *   --------(t,s)---------
   *   |         |          |
   *   |   p1    |    p2    |
   *   |         |          |
   *   ---------------------- t = 1
   *
   * \pre Parameter \a t and \a s must be between 0 and 1
   */
  void split(T t,
             T s,
             BezierPatch& p1,
             BezierPatch& p2,
             BezierPatch& p3,
             BezierPatch& p4) const
  {
    // Bisect the patch along the t direction
    split(t, 0, p1, p2);

    // Temporarily store the result in each half and split again
    BezierPatch p0(p1);
    p0.split(s, 1, p1, p3);

    p0 = p2;
    p0.split(s, 1, p2, p4);
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
    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    if(ord_t <= 1 && ord_s <= 1) return true;

    SegmentType the_plane =
      make_plane(m_controlPoints[0][0],
                 m_controlPoints[0][ord_s] m_controlPoints[ord_t][0]);

    double sqDist = 0.0;

    // Check all control points for simplicity
    for(int p = 1; p < ord_t && sqDist <= tol; ++p)
      for(int q = 1; p < ord_s && sqDist <= tol; ++q)
        sqDist += squared_distance(m_controlPoints(p, q), the_plane);

    return (sqDist <= tol);
  }

  /*!
   * \brief Simple formatted print of a Bezier Curve instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    os << "{ order (" << ord_t << ',' << ord_s << ") Bezier Patch ";

    for(int p = 0; p <= ord_t; ++p)
      for(int q = 0; q <= ord_s; ++q)
        os << m_controlPoints(p, q) << ((p < ord_t || q < ord_s) ? "," : "");

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p <= ord_t; ++p)
        for(int q = 0; q <= ord_s; ++q)
          os << m_weights(p, q) << ((p < ord_t || q < ord_s) ? "," : "");
    }
    os << "}";

    return os;
  }

private:
  /// Check that the weights used are positive, and
  ///  that there is one for each control node
  bool isValidRational() const
  {
    if(!isRational()) return true;

    const int ord_t = getOrder_t();
    const int ord_s = getOrder_s();

    if(m_weights.shape()[0] != (ord_t + 1) || m_weights.shape()[1] != (ord_s + 1))
      return false;

    for(int p = 0; p <= ord_t; ++p)
      for(int q = 0; q <= ord_s; ++q)
        if(m_weights(p, q) <= 0) return false;

    return true;
  }

  CoordsMat m_controlPoints;
  axom::Array<T, 2> m_weights;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierCurve
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierPatch<T, NDIMS>& bPatch)
{
  bPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BEZIERPATCH_HPP_
