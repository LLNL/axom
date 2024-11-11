// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSCurve.hpp
 *
 * \brief A NURBS curve primitive
 */

#ifndef AXOM_PRIMAL_NURBSCURVE_HPP_
#define AXOM_PRIMAL_NURBSCURVE_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/KnotVector.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <vector>
#include <ostream>
#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class NURBSCurve;

/*! \brief Overloaded output operator for NURBS Curves*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSCurve<T, NDIMS>& bCurve);

/*!
 * \class NURBSCurve
 *
 * \brief Represents a NURBS curve defined by an array of control points, weights and knots
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * A NURBS curve has degree `p`, `n+1` control points, optionally `n+1` weights, 
 * and a knot vector of length `k+1`. A valid curve has k+1 = n+p+2
 * The curve must be open (clamped on each end) and continuous (unless p = 0)
 * 
 * Nonrational NURBS curves are identified by an empty weights array.
 */
template <typename T, int NDIMS>
class NURBSCurve
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;
  using WeightsVec = axom::Array<T>;
  using CoordsVec = axom::Array<PointType>;
  using KnotVectorType = KnotVector<T>;
  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG(
    (NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
    "A NURBS Curve object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A NURBS Curve must be defined using an arithmetic type");

public:
  /*!
   * \brief Default constructor for a NURBS Curve
   *
   * \note An empty NURBS curve is not valid
   */
  explicit NURBSCurve()
  {
    m_controlPoints.resize(0);
    m_knotvec = KnotVectorType();
    makeNonrational();
  }

  /*!
   * \brief Constructor for a simple NURBS Curve that reserves space for
   *  the minimum (sensible) number of points for the given degree
   *
   * \param [in] deg the degree of the resulting curve
   * 
   * A uniform knot vector is constructed such that the curve is continuous
   * 
   * \pre degree is greater than or equal to 0.
   */
  explicit NURBSCurve(int degree)
  {
    SLIC_ASSERT(degree >= 0);
    m_controlPoints.resize(degree + 1);
    m_knotvec = KnotVectorType(degree + 1, degree);

    makeNonrational();
  }

  /*!
   * \brief Constructor for a NURBS curve with given number of control points and degree
    *
    * \param [in] npts the number of control points
    * \param [in] degree the degree of the curve
    * 
    * A uniform knot vector is constructed such that the curve is continuous
    * 
    * \pre npts > degree, degree >= 0
    */
  NURBSCurve(int npts, int degree)
  {
    SLIC_ASSERT(npts > degree);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    m_knotvec = KnotVectorType(npts, degree);

    makeNonrational();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS curve from a Bezier curve
   *
   * \param [in] bezierCurve the Bezier curve to convert to a NURBS curve 
   */
  explicit NURBSCurve(const BezierCurve<T, NDIMS>& bezierCurve)
  {
    m_controlPoints = bezierCurve.getControlPoints();
    m_weights = bezierCurve.getWeights();

    int degree = bezierCurve.getOrder();
    m_knotvec = KnotVectorType(degree + 1, degree);
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points and degree
   *
   * \param [in] pts the control points of the curve
   * \param [in] npts the number of control points
   * \param [in] degree the degree of the curve
   * 
   * A uniform knot vector is constructed such that the curve is continuous
   * 
   * \pre Requires valid pointers, npts > degree, degree >= 0
   */
  NURBSCurve(const PointType* pts, int npts, int degree)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(npts >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    makeNonrational();

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
    }

    m_knotvec = KnotVectorType(npts, degree);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points, weights, and degree
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] npts the number of control points and weights
   * \param [in] degree the degree of the curve
   * 
   * A uniform knot vector is constructed such that the curve is continuous
   * 
   * \pre Requires valid pointers, npts > degree, npts == nwts, and degree >= 0
   */
  NURBSCurve(const PointType* pts, const T* weights, int npts, int degree)
  {
    SLIC_ASSERT(pts != nullptr && weights != nullptr);
    SLIC_ASSERT(npts >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    m_weights.resize(npts);

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
      m_weights[i] = weights[i];
    }

    // Knots for the clamped curve
    m_knotvec = KnotVectorType(npts, degree);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points and knots
   *
   * \param [in] pts the control points of the curve
   * \param [in] npts the number of control points and weights
   * \param [in] knots the weights of the control points
   * \param [in] nkts the length of the knot vector
   * 
   * For clamped and continuous curves, npts and nkts 
   *   uniquely determine the degree
   * 
   * \pre Requires valid pointers, a valid knot vector and npts > degree
   */
  NURBSCurve(const PointType* pts, int npts, const T* knots, int nkts)
  {
    SLIC_ASSERT(nkts >= 0);
    SLIC_ASSERT(npts >= 0);

    m_controlPoints.resize(npts);

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
    }
    makeNonrational();

    m_knotvec = KnotVectorType(knots, nkts, nkts - npts - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points, weights, and knots
   *
   * \param [in] pts the control points of the curve
   * \param [in] npts the number of control points and weights
   * \param [in] weights the weights of the control points
   * \param [in] knots the weights of the control points
   * \param [in] nkts the length of the knot vector
   * 
   * For clamped and continuous curves, npts and nkts 
   *   uniquely determine the degree
   * 
   * \pre Requires valid pointers, a valid knot vector, npts > degree, npts == nwts, wts > 0
   */
  NURBSCurve(const PointType* pts, const T* weights, int npts, const T* knots, int nkts)
  {
    SLIC_ASSERT(pts != nullptr && weights != nullptr);
    SLIC_ASSERT(nkts >= 0);
    SLIC_ASSERT(npts >= 0);

    m_controlPoints.resize(npts);
    m_weights.resize(npts);

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
      m_weights[i] = weights[i];
    }

    m_knotvec = KnotVectorType(knots, nkts, nkts - npts - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of data
   *
   * \param [in] pts the control points of the curve
   * \param [in] degree the degree of the curve
   *
   * A uniform knot vector is constructed such that the curve is continuous
   * 
   * \pre Requires npts > degree and degree >= 0
   */
  NURBSCurve(const axom::Array<PointType>& pts, int degree)
    : m_controlPoints(pts)
  {
    makeNonrational();
    m_knotvec = KnotVectorType(pts.size(), degree);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of data
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] degree the degree of the curve
   *
   * A uniform knot vector is constructed such that the curve is continuous
   * 
   * \pre Requires npts > degree, degree >= 0, and npts == nwts
   */
  NURBSCurve(const axom::Array<PointType>& pts,
             const axom::Array<T>& weights,
             int degree)
    : m_controlPoints(pts)
    , m_weights(weights)
  {
    m_knotvec = KnotVectorType(pts.size(), degree);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of nodes and knots
   *
   * \param [in] pts the control points of the curve
   * \param [in] knots the knot vector of the curve
   *
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector and npts > degree
   */
  NURBSCurve(const axom::Array<PointType>& pts, const axom::Array<T>& knots)
    : m_controlPoints(pts)
  {
    makeNonrational();

    m_knotvec = KnotVectorType(knots, knots.size() - pts.size() - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of nodes, weights, and knots
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] knots the knot vector of the curve
   *
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector, npts > degree, npts == nwts, wts > 0
   */
  NURBSCurve(const axom::Array<PointType>& pts,
             const axom::Array<T>& weights,
             const axom::Array<T>& knots)
    : m_controlPoints(pts)
    , m_weights(weights)
  {
    m_knotvec = KnotVectorType(knots, knots.size() - pts.size() - 1);

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor from axom array of nodes and KnotVector object
   *
   * \param [in] pts the control points of the curve
   * \param [in] knotVector A KnotVector object
   * 
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   *
   * \pre Requires a valid knot vector and npts > degree
   */
  NURBSCurve(const axom::Array<PointType>& pts, const KnotVectorType& knotVector)
    : m_controlPoints(pts)
    , m_knotvec(knotVector)
  {
    makeNonrational();

    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Constructor from axom array of nodes, weights, and KnotVector object
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] knotVector A KnotVector object
   *
   * For clamped and continuous curves, npts and the knot vector 
   *   uniquely determine the degree
   * 
   * \pre Requires a valid knot vector, npts > degree, npts == nwts, wts > 0
   */
  NURBSCurve(const axom::Array<PointType>& pts,
             const axom::Array<T>& weights,
             const KnotVectorType& knotVector)
    : m_controlPoints(pts)
    , m_weights(weights)
    , m_knotvec(knotVector)
  {
    SLIC_ASSERT(isValidNURBS());
  }

  /*!
   * \brief Evaluate a NURBS Curve at a particular parameter value \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \return p The point of the NURBS curve at t
   * 
   * Adapted from Algorithm A4.1 on page 124 of "The NURBS Book"
   * 
   * \pre Requires \a t in the span of the knots
   */
  PointType evaluate(T t) const
  {
    SLIC_ASSERT(t >= m_knotvec[0] && t <= m_knotvec[m_knotvec.getNumKnots() - 1]);

    const auto span = m_knotvec.findSpan(t);
    const auto N_evals = m_knotvec.calculateBasisFunctionsBySpan(span, t);
    const int degree = m_knotvec.getDegree();

    PointType P;

    if(isRational())
    {
      Point<T, NDIMS + 1> H;
      for(int j = 0; j <= degree; ++j)
      {
        const auto offset = span - degree + j;
        const T weight = m_weights[offset];
        const auto& controlPoint = m_controlPoints[offset];

        for(int i = 0; i < NDIMS; ++i)
        {
          H[i] += N_evals[j] * weight * controlPoint[i];
        }
        H[NDIMS] += N_evals[j] * weight;
      }

      for(int i = 0; i < NDIMS; ++i)
      {
        P[i] = H[i] / H[NDIMS];
      }

      return P;
    }
    else
    {
      for(int j = 0; j <= degree; ++j)
      {
        const auto offset = span - degree + j;
        const auto& controlPoint = m_controlPoints[offset];

        for(int i = 0; i < NDIMS; ++i)
        {
          P[i] += N_evals[j] * controlPoint[i];
        }
      }

      return P;
    }
  }

  /*!
   * \brief Evaluate the first derivative of NURBS Curve at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \return p The vector of NURBS curve's derivative at t
   * 
   * \pre Requires \a t in the span of the knots
   */
  VectorType dt(T t) const
  {
    PointType eval;
    axom::Array<VectorType> ders;

    evaluateDerivatives(t, 1, eval, ders);
    return ders[0];
  }

  /*!
   * \brief Evaluate the second derivative of NURBS Curve at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \return p The vector of NURBS curve's 2nd derivative at t
   * 
   * \pre Requires \a t in the span of the knots
   */
  VectorType dtdt(T t) const
  {
    PointType eval;
    axom::Array<VectorType> ders;

    evaluateDerivatives(t, 2, eval, ders);
    return ders[1];
  }

  /*!
   * \brief Evaluate the curve and the first derivative at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \param [out] eval The point on the curve at t
   * \param [out] Dt The first derivative of the curve at t
   * 
   * \pre Requires \a t in the span of the knots
   */
  void evaluate_first_derivative(T t, PointType& eval, VectorType& Dt) const
  {
    axom::Array<VectorType> ders;

    evaluateDerivatives(t, 1, eval, ders);
    Dt = ders[0];
  }

  /*!
   * \brief Evaluate the curve and the first two derivatives at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \param [out] eval The point on the curve at t
   * \param [out] Dt The first derivative of the curve at t
   * \param [out] DtDt The second derivative of the curve at t
   * 
   * \pre Requires \a t in the span of the knots
   */
  void evaluate_second_derivative(T t,
                                  PointType& eval,
                                  VectorType& Dt,
                                  VectorType& DtDt) const
  {
    axom::Array<VectorType> ders;

    evaluateDerivatives(t, 2, eval, ders);
    Dt = ders[0];
    DtDt = ders[1];
  }

  /*!
   * \brief Evaluate the curve and the first \a d derivatives at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \param [in] d The number of derivatives to evaluate
   * \param [out] eval The point on the curve at t
   * \param [out] ders An array of the first d derivatives at t
   * 
   * Implementation adapted from Algorithm A3.2 on p. 93 of "The NURBS Book".
   * Rational derivatives from Algorithm A4.2 on p. 127 of "The NURBS Book".
   * 
   * \pre Requires \a t in the span of the knots
   */
  void evaluateDerivatives(T t,
                           int d,
                           PointType& eval,
                           axom::Array<VectorType>& ders) const
  {
    SLIC_ASSERT(t >= m_knotvec[0] && t <= m_knotvec[m_knotvec.getNumKnots() - 1]);

    const int p = m_knotvec.getDegree();
    ders.resize(d);

    const bool isCurveRational = this->isRational();

    int du = std::min(d, p);
    const auto span = m_knotvec.findSpan(t);
    const auto N_evals = m_knotvec.derivativeBasisFunctionsBySpan(span, t, du);

    // Store w(u) in Awders[NDIMS][0], w'(u) in Awders[NDIMS][1], ...
    axom::Array<Point<T, NDIMS + 1>> Awders(d + 1);

    // Compute the homogenous point and its d derivatives
    for(int k = 0; k <= du; k++)
    {
      Point<T, NDIMS + 1> Pw(0.0);
      for(int j = 0; j <= p; j++)
      {
        auto offset = span - p + j;
        const double weight = isCurveRational ? m_weights[offset] : 1.0;

        // Compute the weighted point.
        for(int i = 0; i < NDIMS; ++i)
        {
          Pw[i] += N_evals[k][j] * weight * m_controlPoints[offset][i];
        }
        Pw[NDIMS] += N_evals[k][j] * weight;
      }

      Awders[k] = Pw;
    }

    // Do the rational part from A4.2. Note that Awders[0][0-NDIMS] is the point A(u)
    // and Awders[0][NDIMS] is w(u).

    // Zero out the points
    for(int i = 0; i < NDIMS; ++i)
    {
      eval[i] = 0.0;
      for(int k = 0; k < d; k++)
      {
        ders[k][i] = 0.0;
      }
    }

    // Separate k = 0 case
    Point<T, NDIMS + 1> v = Awders[0];
    for(int i = 0; i < NDIMS; ++i)
    {
      eval[i] = v[i] / Awders[0][NDIMS];
    }

    // Separate k = 1 case
    v = Awders[1];
    for(int j = 0; j < NDIMS; ++j)
    {
      ders[0][j] = (v[j] - Awders[1][NDIMS] * eval[j]) / Awders[0][NDIMS];
    }

    // Recursive formula for k >= 2
    for(int k = 2; k <= d; k++)
    {
      v = Awders[k];
      for(int j = 0; j < NDIMS; ++j)
      {
        v[j] = v[j] - Awders[k][NDIMS] * eval[j];
      }

      for(int i = 1; i < k; i++)
      {
        auto bin = axom::utilities::binomialCoefficient(k, i);
        for(int j = 0; j < NDIMS; ++j)
        {
          v[j] = v[j] - bin * Awders[i][NDIMS] * ders[k - i - 1][j];
        }
      }

      for(int j = 0; j < NDIMS; ++j)
      {
        ders[k - 1][j] = v[j] / Awders[0][NDIMS];
      }
    }
  }

  /*! 
   * \brief Insert a knot with given multiplicity
   *
   * \param [in] t The parameter value of the knot to insert
   * \param [in] target_multiplicity The multiplicity of the knot to insert
   * \return The index of the new knot
   * 
   * Algorithm A5.1 on p. 151 of "The NURBS Book"
   * 
   * \note If the knot is already present, it will be inserted
   *  up to the given multiplicity, or the maximum permitted by the degree
   * 
   * \pre Requires \a t in the span of the knots
   * 
   * \return The (maximum) index of the new knot
   */
  axom::IndexType insertKnot(T t, int target_multiplicity = 1)
  {
    SLIC_ASSERT(t >= m_knotvec[0] && t <= m_knotvec[m_knotvec.getNumKnots() - 1]);
    SLIC_ASSERT(target_multiplicity > 0);

    const bool isRational = this->isRational();

    const int n = getNumControlPoints() - 1;
    const int p = m_knotvec.getDegree();

    // Find the span and multiplicity of the knot
    int s = 0;
    const auto span = m_knotvec.findSpan(t, s);

    // Fix the maximum multiplicity of the knot
    int r = std::min(target_multiplicity - s, p - s);
    if(r <= 0)
    {
      return span;  // Early exit if no knots to add
    }

    // Save unaltered control points
    CoordsVec newControlPoints(m_controlPoints.size() + r);
    WeightsVec newWeights(isRational ? m_weights.size() + r : 0);
    for(int i = 0; i <= span - p; ++i)
    {
      newControlPoints[i] = m_controlPoints[i];
      if(isRational)
      {
        newWeights[i] = m_weights[i];
      }
    }

    for(auto i = span - s; i <= n; ++i)
    {
      newControlPoints[i + r] = m_controlPoints[i];
      if(isRational)
      {
        newWeights[i + r] = m_weights[i];
      }
    }

    // Insert the new control points
    CoordsVec tempControlPoints(p + 1);
    WeightsVec tempWeights(isRational ? p + 1 : 0);
    for(int i = 0; i <= p - s; ++i)
    {
      for(int N = 0; N < NDIMS; ++N)
      {
        tempControlPoints[i][N] = m_controlPoints[span - p + i][N] *
          (isRational ? m_weights[span - p + i] : 1.0);
      }

      if(isRational)
      {
        tempWeights[i] = m_weights[span - p + i];
      }
    }

    // Insert the knot r times
    axom::IndexType L;
    for(int j = 1; j <= r; ++j)
    {
      L = span - p + j;
      for(int i = 0; i <= p - j - s; ++i)
      {
        T alpha =
          (t - m_knotvec[L + i]) / (m_knotvec[i + span + 1] - m_knotvec[L + i]);

        tempControlPoints[i].array() =
          (1.0 - alpha) * tempControlPoints[i].array() +
          alpha * tempControlPoints[i + 1].array();

        if(isRational)
        {
          tempWeights[i] =
            (1.0 - alpha) * tempWeights[i] + alpha * tempWeights[i + 1];
        }
      }

      for(int N = 0; N < NDIMS; ++N)
      {
        newControlPoints[L][N] =
          tempControlPoints[0][N] / (isRational ? tempWeights[0] : 1.0);
        newControlPoints[span + r - j - s][N] = tempControlPoints[p - j - s][N] /
          (isRational ? tempWeights[p - j - s] : 1.0);
      }

      if(isRational)
      {
        newWeights[L] = tempWeights[0];
        newWeights[span + r - j - s] = tempWeights[p - j - s];
      }
    }

    for(auto i = L + 1; i < span - s; ++i)
    {
      for(int N = 0; N < NDIMS; ++N)
      {
        newControlPoints[i][N] =
          tempControlPoints[i - L][N] / (isRational ? tempWeights[i - L] : 1.0);
      }

      if(isRational)
      {
        newWeights[i] = tempWeights[i - L];
      }
    }

    // Update the knot vector and control points
    m_knotvec.insertKnotBySpan(span, t, r);
    m_controlPoints = newControlPoints;
    m_weights = newWeights;

    return span + r;
  }

  /*!
   * \brief Splits a NURBS curve into two curves at a given parameter value
   *
   * \param [in] t parameter value between 0 and 1 at which to evaluate
   * \param [out] n1 First output NURBS curve
   * \param [out] n2 Second output NURBS curve
   * \param [in] normalize Whether to normalize the output curves
   * 
   * If t is at the knot u_0 or u_max, will return the original curve and
   *  a degenerate curve with the first or last control point
   *
   * \pre Requires \a t in the *interior* of the span of the knots
   */
  void split(T t,
             NURBSCurve<T, NDIMS>& n1,
             NURBSCurve<T, NDIMS>& n2,
             bool normalize = false) const
  {
    SLIC_ASSERT(t > m_knotvec[0] && t < m_knotvec[m_knotvec.getNumKnots() - 1]);

    const bool isCurveRational = this->isRational();
    const int p = getDegree();

    n1 = *this;

    // Will make the multiplicity of the knot equal to p,
    //  even if it is already >= 1
    auto s = n1.insertKnot(t, p);
    auto nkts = n1.getNumKnots();

    // Split the knot vector, add to the returned curves
    KnotVectorType k1, k2;
    n1.getKnots().splitBySpan(s, k1, k2);

    n1.m_knotvec = k1;
    n2.m_knotvec = k2;

    // Copy the control points
    n2.m_controlPoints.resize(nkts - s - 1);
    n2.m_weights.resize(isCurveRational ? nkts - s - 1 : 0);
    for(int i = 0; i < n2.m_controlPoints.size(); ++i)
    {
      n2.m_controlPoints[nkts - s - 2 - i] =
        n1.m_controlPoints[n1.m_controlPoints.size() - 1 - i];

      if(isCurveRational)
      {
        n2.m_weights[nkts - s - 2 - i] =
          n1.m_weights[n1.m_weights.size() - 1 - i];
      }
    }

    n1.m_controlPoints.resize(s - p + 1);
    n1.m_weights.resize(isCurveRational ? s - p + 1 : 0);

    if(normalize)
    {
      n1.normalize();
      n2.normalize();
    }
  }

  /// \brief Normalize the knot vector to the span of [0, 1]
  void normalize() { m_knotvec.normalize(); }

  /*!
   * \brief Rescale the knot vector to the span of [a, b]
   * 
   * \param [in] a The lower bound of the new knot vector
   * \param [in] b The upper bound of the new knot vector
   * 
   * \pre Requires a < b
   */
  void rescale(T a, T b)
  {
    SLIC_ASSERT(a < b);
    m_knotvec.rescale(a, b);
  }

  /*!
   * \brief Splits a NURBS curve (at each internal knot) into several Bezier curves
   *   
   * If the curve is of degree 0, will return a set of disconnected, order 0 Bezier curves
   * 
   * Adapted from Algorithm A5.6 on p. 173 of "The NURBS Book"
   * 
   * \return An array of Bezier curves
   */
  axom::Array<BezierCurve<T, NDIMS>> extractBezier() const
  {
    bool isCurveRational = this->isRational();
    int p = getDegree();
    int n = getNumControlPoints() - 1;
    int ks = m_knotvec.getNumKnotSpans();

    axom::Array<BezierCurve<T, NDIMS>> beziers(ks);
    for(auto& bezier : beziers)
    {
      bezier.setOrder(p);
      if(isCurveRational)
      {
        bezier.makeRational();
      }
    }

    // Handle this special case
    if(p == 0)
    {
      for(int i = 0; i < getNumControlPoints(); ++i)
      {
        beziers[i][0] = m_controlPoints[i];

        if(isCurveRational)
        {
          beziers[i].setWeight(0, m_weights[i]);
        }
      }

      return beziers;
    }

    axom::Array<T> alphas(p - 1);

    int m = n + p + 1;
    int a = p;
    int b = p + 1;

    int nb = 0;

    for(int i = 0; i <= p; ++i)
    {
      beziers[nb][i] = m_controlPoints[i];
      if(isCurveRational)
      {
        beziers[nb].setWeight(i, m_weights[i]);
      }
    }

    while(b < m)
    {
      int i = b;
      while(b < m && m_knotvec[b] == m_knotvec[b + 1])
      {
        b++;
      }

      int mult = b - i + 1;

      if(mult < p)
      {
        T numer = m_knotvec[b] - m_knotvec[a];

        for(int j = p; j > mult; j--)
        {
          alphas[j - mult - 1] = numer / (m_knotvec[a + j] - m_knotvec[a]);
        }

        // Do the knot insertion in-place
        for(int j = 1; j <= p - mult; ++j)
        {
          int save = p - mult - j;
          int s = mult + j;
          for(int k = p; k >= s; k--)
          {
            T weight_k = isCurveRational ? beziers[nb].getWeight(k) : 1.0;
            T weight_km1 = isCurveRational ? beziers[nb].getWeight(k - 1) : 1.0;
            T alpha = alphas[k - s];

            if(isCurveRational)
            {
              beziers[nb].setWeight(
                k,
                alpha * weight_k + (1.0 - alpha) * weight_km1);
            }

            for(int N = 0; N < NDIMS; ++N)
            {
              beziers[nb][k][N] =
                (alpha * beziers[nb][k][N] * weight_k +
                 (1.0 - alpha) * beziers[nb][k - 1][N] * weight_km1) /
                (isCurveRational ? beziers[nb].getWeight(k) : 1.0);
            }
          }

          if(b < m)
          {
            beziers[nb + 1][save] = beziers[nb][p];
            if(isCurveRational)
            {
              beziers[nb + 1].setWeight(save, beziers[nb].getWeight(p));
            }
          }
        }
      }

      ++nb;

      if(b < m)
      {
        for(int i = p - mult; i <= p; ++i)
        {
          beziers[nb][i] = m_controlPoints[b - p + i];
          if(isCurveRational)
          {
            beziers[nb].setWeight(i, m_weights[b - p + i]);
          }
        }
        a = b;
        b++;
      }
    }

    return beziers;
  }

  /*!
   * \brief Reset the degree and resize arrays of points (and weights)
   *
   * \param [in] npts The target number of control points
   * \param [in] degree The target degree
   * 
   * \warning This method will replace existing knot vector with a uniform one.
   */
  void setParameters(int npts, int degree)
  {
    SLIC_ASSERT(npts > degree);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);

    if(isRational())
    {
      m_weights.resize(npts);
    }

    m_knotvec.makeUniform(npts, degree);
  }

  /*!
   * \brief Reset the knot vector
   *
   * \param [in] degree The target degree
   * 
   * \warning This method does NOT change the existing control points, 
   *  i.e. does not perform degree elevation/reduction. 
   *  Will replace existing knot vector with a uniform one.
   *  
   * \pre Requires target degree < npts and degree >= 0
   */
  void setDegree(int degree)
  {
    SLIC_ASSERT(0 <= degree && degree < getNumControlPoints());

    m_knotvec.makeUniform(getNumControlPoints(), degree);
  }

  /// \brief Returns the degree of the NURBS Curve
  int getDegree() const { return static_cast<int>(m_knotvec.getDegree()); }

  /// \brief Returns the order of the NURBS Curve
  int getOrder() const { return static_cast<int>(m_knotvec.getDegree() + 1); }

  /// \brief Returns the number of control poitns in the NURBS Curve
  int getNumControlPoints() const
  {
    return static_cast<int>(m_controlPoints.size());
  }

  /*!
   * \brief Set the number control points
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT change the existing control points, 
   *  i.e. does not perform degree elevation/reduction. 
   *  Will replace existing knot vector with a uniform one.
   */
  void setNumControlPoints(int npts)
  {
    const int deg = m_knotvec.getDegree();
    SLIC_ASSERT(npts > deg);

    m_controlPoints.resize(npts);
    if(isRational())
    {
      m_weights.resize(npts);
    }

    m_knotvec.makeUniform(npts, deg);
  }

  /// \brief Returns the number of knots in the NURBS Curve
  axom::IndexType getNumKnots() const { return m_knotvec.getNumKnots(); }

  /// \brief Make nonrational. If already rational, do nothing
  void makeNonrational() { m_weights.resize(0); }

  /// \brief Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      m_weights.resize(m_controlPoints.size());
      m_weights.fill(1.0);
    }
  }

  /// \brief Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// \brief Clears the list of control points, make nonrational
  void clear()
  {
    m_controlPoints.clear();
    m_knotvec.clear();
    makeNonrational();
  }

  /*!
   * \brief Retrieve the control point at index \a idx
   *
   * \param [in] idx The index of the control point
   * 
   * \return A reference to the control point at index \a idx
   */
  PointType& operator[](int idx) { return m_controlPoints[idx]; }

  /*!
   * \brief Retrieve the control point at index \a idx
   *
   * \param [in] idx The index of the control point
   * 
   * \return A const reference to the control point at index \a idx
   */
  const PointType& operator[](int idx) const { return m_controlPoints[idx]; }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] idx The index of the weight
   * \pre Requires that the curve be rational
   * 
   * \return The weight at index \a idx
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
   * 
   * \pre Requires that the curve be rational, and the weight be rational
   */
  void setWeight(int idx, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights[idx] = weight;
  }

  /*!
   * \brief Get a specific knot value
   *
   * \param [in] idx The index of the knot
   * 
   * \return The knot value at index \a idx
   */
  const T& getKnot(int idx) const { return m_knotvec[idx]; }

  /*!
   * \brief Set the knot value at a specific index
   *
   * \param [in] idx The index of the knot
   * \param [in] knot The updated value of the knot
   */
  void setKnot(int idx, T knot) { m_knotvec[idx] = knot; }

  /*! 
   * \brief Set the knot vector by an axom::Array
   *
   * \param [in] knots The new knot vector
   */
  void setKnots(const axom::Array<T>& knots, int degree)
  {
    m_knotvec = KnotVectorType(knots, degree);
  }

  /*! 
   * \brief Set the knot vector by a KnotVector object
   *
   * \param [in] knotVector The new knot vector
   */
  void setKnots(const KnotVectorType& knotVector) { m_knotvec = knotVector; }

  /*!
   * \brief Equality operator for NURBS Curves
   * 
   * \param [in] lhs The left-hand side NURBS curve
   * \param [in] rhs The right-hand side NURBS curve
   * 
   * \return True if the two curves are equal, false otherwise
   */
  friend inline bool operator==(const NURBSCurve<T, NDIMS>& lhs,
                                const NURBSCurve<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_knotvec == rhs.m_knotvec) && (lhs.m_weights == rhs.m_weights);
  }

  /*!
   * \brief Inequality operator for NURBS Curves
   * 
   * \param [in] lhs The left-hand side NURBS curve
   * \param [in] rhs The right-hand side NURBS curve
   * 
   * \return True if the two curves are not equal, false otherwise
   */
  friend inline bool operator!=(const NURBSCurve<T, NDIMS>& lhs,
                                const NURBSCurve<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /// \brief Returns a copy of the NURBS curve's control points
  CoordsVec getControlPoints() const { return m_controlPoints; }

  /// \brief Returns a copy of the NURBS curve's control points
  WeightsVec getWeights() const { return m_weights; }

  /// \brief Return a copy of the knot vector
  KnotVectorType getKnots() const { return m_knotvec; }

  /// \brief Return a copy of the knot vector as an array
  axom::Array<T> getKnotsArray() const { return m_knotvec.getArray(); }

  /// \brief Reverses the order of the NURBS curve's control points and weights
  void reverseOrientation()
  {
    const int npts = getNumControlPoints();
    const int npts_mid = (npts + 1) / 2;
    for(int i = 0; i < npts_mid; ++i)
    {
      axom::utilities::swap(m_controlPoints[i], m_controlPoints[npts - i - 1]);
    }

    if(isRational())
    {
      for(int i = 0; i < npts_mid; ++i)
      {
        axom::utilities::swap(m_weights[i], m_weights[npts - i - 1]);
      }
    }

    // Reverse the orientation of the knots
    m_knotvec.reverse();
  }

  /// \brief Returns an axis-aligned bounding box containing the NURBS curve
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /// \brief Returns an oriented bounding box containing the NURBS curve
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(),
                                   static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Simple formatted print of a NURBS Curve instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    int npts = getNumControlPoints();
    int nkts = m_knotvec.getNumKnots();

    int deg = m_knotvec.getDegree();

    os << "{ degree " << deg << " NURBS Curve ";
    for(int p = 0; p < npts; ++p)
    {
      os << m_controlPoints[p] << (p < npts - 1 ? "," : "");
    }

    if(isRational())
    {
      os << ", weights [";
      for(int p = 0; p < npts; ++p)
      {
        os << m_weights[p] << (p < npts - 1 ? ", " : "]");
      }
    }

    os << ", knots [";
    for(int i = 0; i < nkts; ++i)
    {
      os << m_knotvec[i] << (i < nkts - 1 ? ", " : "]");
    }
    os << "}";

    return os;
  }

  /// \brief Function to check if the NURBS curve is valid
  bool isValidNURBS() const
  {
    // Check monotonicity, open-ness, continuity
    if(!m_knotvec.isValid())
    {
      return false;
    }

    // Number of knots must match the number of control points
    int p = m_knotvec.getDegree();
    if(m_knotvec.getNumKnots() != m_controlPoints.size() + p + 1)
    {
      return false;
    }

    if(isRational())
    {
      if(m_weights.size() != m_controlPoints.size())
      {
        // Weights must match the number of control points
        return false;
      }

      for(int i = 0; i < m_weights.size(); ++i)
      {
        if(m_weights[i] <= 0)
        {
          // Weights must be positive
          return false;
        }
      }
    }

    return true;
  }

private:
  CoordsVec m_controlPoints;
  WeightsVec m_weights;
  KnotVectorType m_knotvec;
};

//------------------------------------------------------------------------------
/// Free functions related to NURBSCurve
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSCurve<T, NDIMS>& nCurve)
{
  nCurve.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::NURBSCurve using fmt
template <typename T, int NDIMS>
struct axom::fmt::formatter<axom::primal::NURBSCurve<T, NDIMS>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_NURBSCURVE_HPP_