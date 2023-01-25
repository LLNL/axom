// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file BezierCurve.hpp
 *
 * \brief A BezierCurve primitive
 */

#ifndef AXOM_PRIMAL_BEZIERCURVE_HPP_
#define AXOM_PRIMAL_BEZIERCURVE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
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
class BezierCurve;

/*! \brief Overloaded output operator for Bezier Curves*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierCurve<T, NDIMS>& bCurve);

/*!
 * \class BezierCurve
 *
 * \brief Represents a Bezier curve defined by an array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * The order of a Bezier curve with N+1 control points is N.
 * The curve is approximated by the control points,
 * parametrized from t=0 to t=1.
 * 
 * Contains an array of positive weights to represent a rational Bezier curve.
 * Nonrational Bezier curves are identified by an empty weights array.
 * Algorithms for Rational Bezier curves derived from 
 * Gerald Farin, "Algorithms for rational Bezier curves"
 * Computer-Aided Design, Volume 15, Number 2, 1983,
 */
template <typename T, int NDIMS>
class BezierCurve
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = NumericArray<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;
  using CoordsVec = axom::Array<PointType>;
  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A Bezier Curve object may be defined in 2-D or 3-D");
  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A Bezier Curve must be defined using an arithmetic type");

public:
  /*!
   * \brief Constructor for a Bezier Curve that reserves space for
   *  the given order of the curve
   *
   * \param [in] order the order of the resulting Bezier curve
   * \pre order is greater than or equal to -1.
   */
  explicit BezierCurve(int ord = -1)
  {
    SLIC_ASSERT(ord >= -1);
    const int sz = utilities::max(-1, ord + 1);
    m_controlPoints.resize(sz);

    makeNonrational();
  }

  /*!
   * \brief Constructor for an order \a ord= n Bezier curve
   * from a list of coordinates:
   * \verbatim {x_0, x_1, x_2,...,x_n,
   *            y_0, y_1, y_2,...,y_n,
   *            z_0, z_1, z_2,...,z_n}
   *
   * \param [in] pts an array with (n+1)*NDIMS entries, ordered by coordinate
   * then by polynomial order
   * \param [in] ord Polynomial order of the curve
   * \pre order is greater than or equal to zero
   */
  BezierCurve(T* pts, int ord)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord >= 0);

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);

    for(int p = 0; p <= ord; p++)
    {
      auto& pt = m_controlPoints[p];
      for(int j = 0; j < NDIMS; j++)
      {
        pt[j] = pts[j * (ord + 1) + p];
      }
    }

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Bezier Curve from an array of coordinates
   *
   * \param [in] pts a vector with ord+1 control points
   * \param [in] ord The Curve's polynomial order
   * \pre order is greater than or equal to zero
   *
   */
  BezierCurve(PointType* pts, int ord)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord >= 0);

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);

    for(int p = 0; p <= ord; ++p)
    {
      m_controlPoints[p] = pts[p];
    }

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Bezier Curve from an array of coordinates
   *
   * \param [in] pts a vector with ord+1 control points
   * \param [in] weights a vector with ord+1 positive weights
   * \param [in] ord The Curve's polynomial order
   * \pre order is greater than or equal to zero
   *
   */
  BezierCurve(PointType* pts, T* weights, int ord)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord >= 0);

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);
    for(int p = 0; p <= ord; ++p) m_controlPoints[p] = pts[p];

    if(weights == nullptr)
      m_weights.resize(0);
    else
    {
      m_weights.resize(sz);
      for(int p = 0; p <= ord; ++p) m_weights[p] = weights[p];
      SLIC_ASSERT(isValidRational());
    }
  }

  /*!
   * \brief Constructor for a Bezier Curve from an vector of coordinates
   *
   * \param [in] pts a vector with ord+1 control points
   * \param [in] ord The Curve's polynomial order
   * \pre order is greater than or equal to zero
   *
   */
  BezierCurve(const axom::Array<PointType>& pts, int ord)
  {
    SLIC_ASSERT(ord >= 0);

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);
    m_controlPoints = pts;

    makeNonrational();
  }

  /*!
   * \brief Constructor for a Rational Bezier Curve from an vector 
   * of coordinates and weights
   *
   * \param [in] pts a vector with ord+1 control points
   * \param [in] weights a vector with ord+1 positive weights
   * \param [in] ord The Curve's polynomial order
   * \pre order is greater than or equal to zero
   *
   */
  BezierCurve(const axom::Array<PointType>& pts,
              const axom::Array<T>& weights,
              int ord)
  {
    SLIC_ASSERT(ord >= 0);
    SLIC_ASSERT(pts.size() == weights.size());

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);
    m_weights.resize(sz);

    m_controlPoints = pts;
    m_weights = weights;

    SLIC_ASSERT(isValidRational());
  }

  /// Sets the order of the Bezier Curve
  void setOrder(int ord) { m_controlPoints.resize(ord + 1); }

  /// Returns the order of the Bezier Curve
  int getOrder() const { return static_cast<int>(m_controlPoints.size()) - 1; }

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      const int ord = getOrder();
      m_weights.resize(ord + 1);
      for(int i = 0; i <= ord; i++) m_weights[i] = 1.0;
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.resize(0); }

  /// Use array size as flag for rationality
  bool isRational() const { return (m_weights.size() != 0); }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    const int ord = getOrder();
    for(int p = 0; p <= ord; ++p)
    {
      m_controlPoints[p] = PointType();
    }

    makeNonrational();
  }

  /// Retrieves the control point at index \a idx
  PointType& operator[](int idx) { return m_controlPoints[idx]; }

  /// Retrieves the control point at index \a idx
  const PointType& operator[](int idx) const { return m_controlPoints[idx]; }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] idx The index of the weight
   * \pre Requires that the curve be rational
   */
  const T& getWeight(int idx) const
  {
    SLIC_ASSERT(isRational());
    return m_weights[idx];
  }

  /*!
   * \brief Set the weight at a specific index
   *
   * \param [in] idx The index of the weight
   * \param [in] weight The updated value of the weight
   * \pre Requires that the curve be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int idx, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights[idx] = weight;
  };

  /// Checks equality of two Bezier Curve
  friend inline bool operator==(const BezierCurve<T, NDIMS>& lhs,
                                const BezierCurve<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_weights == rhs.m_weights);
  }

  friend inline bool operator!=(const BezierCurve<T, NDIMS>& lhs,
                                const BezierCurve<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /// Returns a copy of the Bezier curve's control points
  CoordsVec getControlPoints() const { return m_controlPoints; }

  /// Returns a copy of the Bezier curve's weights
  axom::Array<T> getWeights() const { return m_weights; }

  /// Reverses the order of the Bezier curve's control points and weights
  void reverseOrientation()
  {
    const int ord = getOrder();
    const int mid = (ord + 1) / 2;
    for(int i = 0; i < mid; ++i)
    {
      axom::utilities::swap(m_controlPoints[i], m_controlPoints[ord - i]);
    }

    if(isRational())
      for(int i = 0; i < mid; ++i)
      {
        axom::utilities::swap(m_weights[i], m_weights[ord - i]);
      }
  }

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
   * \brief Evaluates a Bezier curve at a particular parameter value \a t
   *
   * \param [in] t parameter value at which to evaluate
   * \return p the value of the Bezier curve at t
   *
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  PointType evaluate(T t) const
  {
    using axom::utilities::lerp;

    PointType ptval;

    const int ord = getOrder();
    axom::Array<T> dCarray(ord + 1);

    if(isRational())
    {
      axom::Array<T> dWarray(ord + 1);

      // Run algorithm from Farin '83 on each dimension
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord; ++p)
        {
          dCarray[p] = m_controlPoints[p][i] * m_weights[p];
          dWarray[p] = m_weights[p];
        }

        for(int p = 1; p <= ord; ++p)
        {
          const int end = ord - p;
          for(int k = 0; k <= end; ++k)
          {
            dCarray[k] = lerp(dCarray[k], dCarray[k + 1], t);
            dWarray[k] = lerp(dWarray[k], dWarray[k + 1], t);
          }
        }
        ptval[i] = dCarray[0] / dWarray[0];
      }

      return ptval;
    }
    else  // If ordinary Bezier curve
    {
      // Run de Casteljau algorithm on each dimension
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord; ++p)
        {
          dCarray[p] = m_controlPoints[p][i];
        }

        for(int p = 1; p <= ord; ++p)
        {
          const int end = ord - p;
          for(int k = 0; k <= end; ++k)
          {
            dCarray[k] = lerp(dCarray[k], dCarray[k + 1], t);
          }
        }
        ptval[i] = dCarray[0];
      }

      return ptval;
    }
  }

  /*!
   * \brief Computes the tangent of  a Bezier curve at a particular parameter value \a t
   *
   * \param [in] t parameter value at which to compute tangent 
   * \return p the tangent vector of the Bezier curve at t
   *
   * \note We typically find the tangent of the curve at \a t between 0 and 1
   */
  VectorType dt(T t) const
  {
    using axom::utilities::lerp;
    VectorType val;

    const int ord = getOrder();
    axom::Array<T> dCarray(ord + 1);

    if(isRational())
    {
      axom::Array<T> dWarray(ord + 1);

      // Run algorithm from Farin '83 on each dimension
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord; ++p)
        {
          dCarray[p] = m_controlPoints[p][i] * m_weights[p];
          dWarray[p] = m_weights[p];
        }

        // stop one step early and do calculation on last two values
        for(int p = 1; p <= ord - 1; ++p)
        {
          const int end = ord - p;
          for(int k = 0; k <= end; ++k)
          {
            dCarray[k] = lerp(dCarray[k], dCarray[k + 1], t);
            dWarray[k] = lerp(dWarray[k], dWarray[k + 1], t);
          }
        }

        val[i] = ord * (dWarray[0] * dCarray[1] - dWarray[1] * dCarray[0]);
        dWarray[0] = lerp(dWarray[0], dWarray[1], t);
        val[i] /= dWarray[0] * dWarray[0];
      }

      return val;
    }
    else
    {
      // Run de Casteljau algorithm on each dimension
      for(int i = 0; i < NDIMS; ++i)
      {
        for(int p = 0; p <= ord; ++p)
        {
          dCarray[p] = m_controlPoints[p][i];
        }

        // stop one step early and take difference of last two values
        for(int p = 1; p <= ord - 1; ++p)
        {
          const int end = ord - p;
          for(int k = 0; k <= end; ++k)
          {
            dCarray[k] = lerp(dCarray[k], dCarray[k + 1], t);
          }
        }
        val[i] = ord * (dCarray[1] - dCarray[0]);
      }

      return val;
    }
  }

  /*!
   * \brief Splits a Bezier curve into two Bezier curves at a given parameter value
   *
   * \param [in] t parameter value between 0 and 1 at which to evaluate
   * \param [out] c1 First output Bezier curve
   * \param [out] c2 Second output Bezier curve
   *
   * \pre Parameter \a t must be between 0 and 1
   */
  void split(T t, BezierCurve& c1, BezierCurve& c2) const
  {
    int ord = getOrder();
    SLIC_ASSERT(ord >= 0);

    // Note: the second curve's control points are computed inline
    //       as we find the first curve's control points
    c2 = *this;

    c1.setOrder(ord);
    c1[0] = m_controlPoints[0];

    if(isRational())
    {
      c1.makeRational();
      c1.setWeight(0, c2.getWeight(0));

      // After each iteration, save the first control point and weight into c1
      for(int p = 1; p <= ord; ++p)
      {
        const int end = ord - p;
        for(int k = 0; k <= end; ++k)
        {
          // Do linear interpolation on weights
          double temp_weight =
            axom::utilities::lerp(c2.getWeight(k), c2.getWeight(k + 1), t);

          for(int i = 0; i < NDIMS; ++i)
          {
            // Do weighted interpolation on nodes
            c2[k][i] = axom::utilities::lerp(c2.getWeight(k) * c2[k][i],
                                             c2.getWeight(k + 1) * c2[k + 1][i],
                                             t) /
              temp_weight;
          }

          c2.setWeight(k, temp_weight);

          c1[p] = c2[0];
          c1.setWeight(p, c2.getWeight(0));
        }
      }
    }
    else  // Code can be simpler if not rational Bezier curves
    {
      // Run de Casteljau algorithm
      // After each iteration, save the first control point into c1
      for(int p = 1; p <= ord; ++p)
      {
        const int end = ord - p;
        for(int k = 0; k <= end; ++k)
        {
          c2[k] = PointType::lerp(c2[k], c2[k + 1], t);
        }
        c1[p] = c2[0];
      }
    }

    return;
  }

  /*!
   * \brief Predicate to check if the Bezier curve is approximately linear
   *
   * This function checks if the internal control points of the BezierCurve
   * are approximately on the line defined by its two endpoints
   *
   * \param [in] tol Threshold for sum of squared distances
   * \return True if c1 is near-linear
   */
  bool isLinear(double tol = 1E-8) const
  {
    const int ord = getOrder();
    if(ord <= 1)
    {
      return true;
    }

    SegmentType seg(m_controlPoints[0], m_controlPoints[ord]);
    double sqDist = 0.0;
    for(int p = 1; p < ord && sqDist <= tol; ++p)  // check interior control points
    {
      sqDist += squared_distance(m_controlPoints[p], seg);
    }
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
    const int ord = getOrder();

    os << "{ order " << ord << " Bezier Curve ";
    for(int p = 0; p <= ord; ++p)
    {
      os << m_controlPoints[p] << (p < ord ? "," : "");
    }

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p <= ord; ++p)
        os << m_weights[p] << (p < ord ? ", " : "]");
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

    const int ord = getOrder();

    if(m_weights.size() != (ord + 1)) return false;

    for(int i = 0; i <= ord; ++i)
      if(m_weights[i] <= 0) return false;

    return true;
  }

  CoordsVec m_controlPoints;
  axom::Array<T> m_weights;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierCurve
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const BezierCurve<T, NDIMS>& bCurve)
{
  bCurve.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_BEZIERCURVE_HPP_
