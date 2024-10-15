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

#include "axom/primal/geometry/NumericArray.hpp"
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

/*! \brief Overloaded output operator for Bezier Curves*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSCurve<T, NDIMS>& bCurve);

/*!
 * \class NURBSCurve
 *
 * \brief Represents a NURBS curve defined by an array of control points and knots
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * A NURBS curve has degree `p`, `n+1` control points, optionally `n+1` weights, 
 * and a knot vector of length `k+1`. A valid curve has k+1 = n+p+2
 * and knot values t_{i-1} <= t_i for i=1,...,k.
 * The curve must be c0 continuous and is parametrized from t=0 to t=1.
 * 
 * Nonrational Bezier curves are identified by an empty weights array.
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
  using KnotsVec = axom::Array<T>;
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
   * \brief Constructor for a simple NURBS Curve that reserves space for
   *  the minimum (sensible) number of points for the given degree
   *
   * \param [in] deg the degree of the resulting curve
   * \pre degree is greater than or equal to -1.
   */
  explicit NURBSCurve(int degree = -1)
  {
    SLIC_ASSERT(degree >= -1);
    m_controlPoints.resize(degree + 1);
    m_knots.resize(2 * degree + 2);
    makeNonrational();

    for(int i = 0; i < 2 * (degree + 1); ++i)
    {
      // Using integer division to get a uniform knot vector
      m_knots[i] = static_cast<T>(i / (degree + 1));
    }
  }

  /*!
   * \brief Constructor for a NURBS curve from a Bezier curve
   *
   * \param [in] bezierCurve the Bezier curve to convert to a NURBS curve 
   */
  explicit NURBSCurve(BezierCurve<T, NDIMS> bezierCurve)
  {
    m_controlPoints = bezierCurve.getControlPoints();
    m_weights = bezierCurve.getWeights();

    int degree = bezierCurve.getOrder();
    m_knots.resize(2 * (degree + 1));

    for(int i = 0; i < 2 * (degree + 1); ++i)
    {
      // Using integer division to get a uniform knot vector
      m_knots[i] = static_cast<T>(i / (degree + 1));
    }
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points and degree
   *
   * \param [in] pts the control points of the curve
   * \param [in] npts the number of control points
   * \param [in] degree the degree of the curve
   * 
   * The knot vector is constructed such that the curve is continuous
   * \pre requires npts >= degree+1 and degree >= 0
   */
  NURBSCurve(PointType* pts, int npts, int degree)
  {
    SLIC_ASSERT(npts >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    m_knots.resize(npts + degree + 1);
    makeNonrational();

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
    }

    makeKnotsUniform(degree, npts);
  }

  /*!
   * \brief Constructor for a NURBS Curve with control points, weights, and degree
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] npts the number of control points and weights
   * \param [in] degree the degree of the curve
   * 
   * The knot vector is constructed such that the curve is continuous
   * \pre requires npts >= degree+1 and degree >= 0
   */
  NURBSCurve(PointType* pts, T* weights, int npts, int degree)
  {
    SLIC_ASSERT(npts >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    m_weights.resize(npts);
    m_knots.resize(npts + degree + 1);

    for(int i = 0; i < npts; ++i)
    {
      m_controlPoints[i] = pts[i];
      m_weights[i] = weights[i];
    }

    // Knots for the clamped curve
    for(int i = 0; i < degree + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + degree - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - degree - 1; ++i)
    {
      m_knots[degree + 1 + i] = (i + 1.0) / (npts - degree);
    }
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of data
   *
   * \param [in] controlPoints the control points of the curve
   * \param [in] degree the degree of the curve
   *
   * The knot vector is constructed such that the curve is continuous
   * \pre requires npts >= degree+1 and degree >= 0
   */
  NURBSCurve(const axom::Array<PointType>& pts, int degree)
  {
    SLIC_ASSERT(pts.size() >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints = pts;

    int npts = pts.size();
    m_knots.resize(npts + degree + 1);
    makeNonrational();

    // Knots for the clamped curve
    for(int i = 0; i < degree + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + degree - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - degree - 1; ++i)
    {
      m_knots[degree + 1 + i] = (i + 1.0) / (npts - degree);
    }
  }

  /*!
   * \brief Constructor for a NURBS Curve with axom arrays of data
   *
   * \param [in] pts the control points of the curve
   * \param [in] weights the weights of the control points
   * \param [in] degree the degree of the curve
   *
   * The knot vector is constructed such that the curve is clamped
   * \pre requires npts >= degree+1 and degree >= 0
   */
  NURBSCurve(const axom::Array<PointType>& pts,
             const axom::Array<T>& weights,
             int degree)
  {
    SLIC_ASSERT(pts.size() >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints = pts;
    m_weights = weights;

    int npts = pts.size();
    m_knots.resize(npts + degree + 1);

    // Knots for the clamped curve
    for(int i = 0; i < degree + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + degree - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - degree - 1; ++i)
    {
      m_knots[degree + 1 + i] = (i + 1.0) / (npts - degree);
    }
  }

  /*!
   * \brief Evaluate a NURBS Curve at a particular parameter value \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \return p The point of the NURBS curve at t
   * 
   * Adapted from Algorithm A4.1 on page 124 of "The NURBS Book"
   * 
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  PointType evaluate(T t) const
  {
    const auto span = findSpan(t);
    const auto N_evals = calculateBasisFunctions(span, t);
    const int degree = getDegree();

    PointType P;

    if(isRational())
    {
      Point<T, NDIMS + 1> H;
      for(int j = 0; j <= degree; ++j)
      {
        const int offset = span - degree + j;
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
        const int offset = span - degree + j;
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
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  VectorType dt(T t) const
  {
    PointType eval;
    axom::Array<VectorType> ders;

    evaluate_derivatives(t, 1, eval, ders);
    return ders[0];
  }

  /*!
   * \brief Evaluate the second derivative of NURBS Curve at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \return p The vector of NURBS curve's 2nd derivative at t
   * 
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  VectorType dtdt(T t) const
  {
    PointType eval;
    axom::Array<VectorType> ders;

    evaluate_derivatives(t, 2, eval, ders);
    return ders[1];
  }

  /*!
   * \brief Evaluate the curve and the first derivative at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \param [out] eval The point on the curve at t
   * \param [out] Dt The first derivative of the curve at t
   * 
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  void evaluate_first_derivative(T t, PointType& eval, VectorType& Dt) const
  {
    axom::Array<VectorType> ders;

    evaluate_derivatives(t, 1, eval, ders);
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
   * \note We typically evaluate the curve at \a t between 0 and 1
   */
  void evaluate_second_derivative(T t,
                                  PointType& eval,
                                  VectorType& Dt,
                                  VectorType& DtDt) const
  {
    axom::Array<VectorType> ders;

    evaluate_derivatives(t, 2, eval, ders);
    Dt = ders[0];
    DtDt = ders[1];
  }

  /*!
   * \brief Evaluate the curve and the first \a d derivatives at parameter \a t
   *
   * \param [in] t The parameter value at which to evaluate
   * \param [out] eval The point on the curve at t
   * \param [out] ders An array of the first d derivatives at t
   * 
   * Implementation adapted from Algorithm A3.2 on p. 93 of "The NURBS Book".
   * Rational derivatives from Algorithm A4.2 on p. 127 of "The NURBS Book".
   */
  void evaluate_derivatives(T t,
                            int d,
                            PointType& eval,
                            axom::Array<VectorType>& ders) const
  {
    const int p = getDegree();
    ders.resize(d);

    const bool isRational = this - this->isRational();

    int du = std::min(d, p);
    const auto span = findSpan(t);
    axom::Array<axom::Array<T>> N_evals(d + 1);
    derivativeBasisFunctions(span, t, du, N_evals);

    // Store w(u) in Awders[NDIMS][0], w'(u) in Awders[NDIMS][1], ...
    axom::Array<Point<T, NDIMS + 1>> Awders(d + 1);

    // Compute the homogenous point and its d derivatives
    for(int k = 0; k <= du; k++)
    {
      Point<T, NDIMS + 1> Pw(0.0);
      for(int j = 0; j <= p; j++)
      {
        auto offset = span - p + j;
        const double weight = isRational ? m_weights[offset] : 1.0;

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
   * \param [in] multiplicity The multiplicity of the knot to insert
   * 
   * Algorithm A5.1 on p. 151 of "The NURBS Book"
   * 
   * \note If the knot is already present, it will be inserted
   *  up to the given multiplicity
   */
  void insertKnot(T t, int multiplicity = 1)
  {
    SLIC_ASSERT(t >= 0.0 && t <= 1.0);
    SLIC_ASSERT(multiplicity > 0);

    const bool isRational = this - this->isRational();

    const int n = getNumControlPoints() - 1;
    const int p = getNumKnots() - n - 2;

    // Find the span of the knot
    const auto span = findSpan(t);

    // Find the current multiplicity of the knot
    int s = 0;
    for(auto i = (t == 1.0) ? p + n + 1 : span; i >= 0; --i)
    {
      if(m_knots[i] == t)
      {
        s++;
      }
      else
      {
        break;
      }
    }

    // Fix the maximum multiplicity of the knot
    multiplicity = std::min(multiplicity, p - s);
    if(multiplicity <= 0)
    {
      return; // Early exit if no knots to add
    }

    // Get new vectors of knots and control points
    KnotsVec newKnots(m_knots);
    for(int r = 0; r < multiplicity; ++r)
    {
      newKnots.insert(newKnots.begin() + span + 1, t);
    }

    // Save unaltered control points
    CoordsVec newControlPoints(m_controlPoints.size() + multiplicity);
    WeightsVec newWeights(isRational ? m_weights.size() + multiplicity : 0);
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
      newControlPoints[i + multiplicity] = m_controlPoints[i];
      if(isRational)
      {
        newWeights[i + multiplicity] = m_weights[i];
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

    // Insert the knot multiplicity times
    axom::IndexType L;
    for(int j = 1; j <= multiplicity; ++j)
    {
      L = span - p + j;
      for(int i = 0; i <= p - j - s; ++i)
      {
        T alpha = (t - m_knots[L + i]) / (m_knots[i + span + 1] - m_knots[L + i]);

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
        newControlPoints[span + multiplicity - j - s][N] =
          tempControlPoints[p - j - s][N] /
          (isRational ? tempWeights[p - j - s] : 1.0);
      }

      if(isRational)
      {
        newWeights[L] = tempWeights[0];
        newWeights[span + multiplicity - j - s] = tempWeights[p - j - s];
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
    m_knots = newKnots;
    m_controlPoints = newControlPoints;
    m_weights = newWeights;
  }

  /*!
   * \brief Splits a NURBS curve into two curves at a given parameter value
   *
   * \param [in] t parameter value between 0 and 1 at which to evaluate
   * \param [out] n1 First output NURBS curve
   * \param [out] n2 Second output NURBS curve
   *
   * \pre Parameter \a t must be between 0 and 1
   */
  void split(T t, NURBSCurve<T, NDIMS>& n1, NURBSCurve<T, NDIMS>& n2) const
  {
    SLIC_ASSERT(t >= 0.0 && t <= 1.0);

    const bool isRational = this - this->isRational();
    const int p = getDegree();

    // Handle the special case of splitting at the endpoints
    if(t == 0.0)
    {
      n1 = *this;
      n2.clear();
      n2.setDegree(p);
      for(int i = 0; i <= p; ++i)
      {
        n2[i] = m_controlPoints[0];
        if(isRational)
        {
          n2.makeRational();
          n2.setWeight(i, m_weights[0]);
        }
      }
      return;
    }
    else if(t == 1.0)
    {
      n1.clear();
      n1.setDegree(p);
      for(int i = 0; i <= p; ++i)
      {
        n1[i] = m_controlPoints[getNumControlPoints() - 1];
        if(isRational)
        {
          n1.makeRational();
          n1.setWeight(i, m_weights[getNumControlPoints() - 1]);
        }
      }
      n2 = *this;
      return;
    }


    // Will make the multiplicity of the knot equal to p,
    //  even if it is already >= 1
    n1 = *this;
    n1.insertKnot(t, p);

    int k = n1.getNumKnots();
    auto s = n1.findSpan(t);

    // Copy the (rescaled) knots first
    n2.m_knots.resize(k - s + p);
    n2.m_knots[0] = 0;
    for(int i = 1; i < k - s + p; ++i)
    {
      n2.m_knots[i] = (n1.m_knots[s - p + i] - t) / (1 - t);
    }

    // Copy the control points
    n2.m_controlPoints.resize(k - s - 1);
    n2.m_weights.resize(isRational ? k - s - 1 : 0);
    for(int i = 0; i < n2.m_controlPoints.size(); ++i)
    {
      n2.m_controlPoints[k - s - 2 - i] =
        n1.m_controlPoints[n1.m_controlPoints.size() - 1 - i];

      if(isRational)
      {
        n2.m_weights[k - s - 2 - i] = n1.m_weights[n1.m_weights.size() - 1 - i];
      }
    }

    // Resize the knots and control points of the first curve
    n1.m_knots.resize(s + 2);
    n1.m_knots[s + 1] = 1;
    for(int i = 0; i < s + 1; ++i)
    {
      n1.m_knots[i] = n1.m_knots[i] / t;
    }

    n1.m_controlPoints.resize(s - p + 1);
    n1.m_weights.resize(isRational ? s - p + 1 : 0);
  }

  /*!
   * \brief Splits a NURBS curve (at each internal knot) into several Bezier curves
   *   
   * \return An array of Bezier curves
   */
  axom::Array<BezierCurve<T, NDIMS>> extractBezier() const
  {
    axom::Array<BezierCurve<T, NDIMS>> beziers;

    const bool isRational = this->isRational();

    int p = getDegree();
    if(p == 0)  // Handle this special case
    {
      for(int i = 0; i < getNumControlPoints(); ++i)
      {
        BezierCurve<T, NDIMS> bezier(0);
        bezier[0] = m_controlPoints[i];

        if(isRational)
        {
          bezier.makeRational();
          bezier.setWeight(0, m_weights[i]);
        }

        beziers.push_back(bezier);
      }

      return beziers;
    }

    // Split the curve at each knot value
    int numBeziers = 1;
    NURBSCurve<T, NDIMS> n1(*this);
    for(int i = p + 1; i < getNumKnots() - p - 1; ++i)
    {
      int old_knot_count = n1.getNumKnots();
      n1.insertKnot(m_knots[i], p);
      int new_knot_count = n1.getNumKnots();

      if(new_knot_count != old_knot_count)
      {
        numBeziers++;
      }
    }

    // For each Bezier, copy the control nodes into Bezier curves
    BezierCurve<T, NDIMS> bezier(p);
    if(isRational)
    {
      bezier.makeRational();
    }

    int the_node = 0;
    for(int i = 0; i < n1.getNumControlPoints(); ++i)
    {
      bezier[the_node] = n1[i];
      if(isRational)
      {
        bezier.setWeight(the_node, n1.getWeight(i));
      }

      the_node++;
      if(the_node == (p + 1))
      {
        --i;  // Endpoint nodes are shared between Bezier curves
        the_node = 0;
        beziers.push_back(bezier);
      }
    }

    return beziers;
  }

  /*!
   * \brief Reset the degree and number of points in the curve
   *
   * \param [in] npts The target number of control points
   * \param [in] degree The target degree
   * 
   * \note Will clear any data already in these arrays.
   */
  void setParameters(int npts, int degree)
  {
    SLIC_ASSERT(npts >= degree + 1);
    SLIC_ASSERT(degree >= 0);

    m_controlPoints.resize(npts);
    m_knots.resize(npts + degree + 1);
    makeNonrational();

    makeKnotsUniform(degree, npts);
  }

  /*!
   * \brief Reset the knot vector and increase the number of control points
   *
   * \param [in] degree The target degree
   * 
   * \warning This method does NOT change the existing control points,
   *  i.e. is not performing degree elevation or reduction.
   */
  void setDegree(int degree)
  {
    SLIC_ASSERT(degree >= 0);

    int npts = getNumControlPoints();
    SLIC_ASSERT(degree < npts);

    if(npts < degree + 1)
    {
      npts = degree + 1;
      m_controlPoints.resize(degree + 1);
    }

    makeKnotsUniform(degree, npts);
  }

  /// \brief Returns the degree of the NURBS Curve
  int getDegree() const
  {
    return static_cast<int>(m_knots.size() - m_controlPoints.size() - 1);
  }

  /// \brief Returns the order of the NURBS Curve
  int getOrder() const
  {
    return static_cast<int>(m_knots.size() - m_controlPoints.size());
  }

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
   * \warning This method does NOT maintain the curve shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints(int npts)
  {
    const int old_npts = getNumControlPoints();
    const int old_deg = getDegree();
    SLIC_ASSERT(npts > old_deg);
    m_controlPoints.resize(npts);
    m_knots.resize(npts + old_deg + 1);

    if(npts != old_npts)
    {
      makeKnotsUniform(getDegree(), npts);
    }
  }

  /// \brief Returns the number of knots in the NURBS Curve
  int getNumKnots() const { return m_knots.size(); }

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
    m_knots.clear();
    makeNonrational();
  }

  /// \brief Retrieve the control point at index \a idx
  PointType& operator[](int idx) { return m_controlPoints[idx]; }

  /// \brief Retrieve the control point at index \a idx
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
  }

  /*!
   * \brief Get a specific knot value
   *
   * \param [in] idx The index of the knot
   */
  const T& getKnot(int idx) const { return m_knots[idx]; }

  /*!
   * \brief Set the knot value at a specific index
   *
   * \param [in] idx The index of the knot
   * \param [in] knot The updated value of the knot
   * \pre Requires that the knot value be internal and between its neighbors
   */
  void setKnot(int idx, T knot)
  {
    SLIC_ASSERT(idx > getDegree() && idx < getNumKnots() - getDegree() - 1);
    SLIC_ASSERT(knot >= m_knots[idx - 1] && knot <= m_knots[idx + 1]);

    m_knots[idx] = knot;
  }

  /// \brief Checks equality of two Bezier Curve
  friend inline bool operator==(const NURBSCurve<T, NDIMS>& lhs,
                                const NURBSCurve<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_knots == rhs.m_knots) && (lhs.m_weights == rhs.m_weights);
  }

  /// \brief Checks inequality of two Bezier Curve
  friend inline bool operator!=(const NURBSCurve<T, NDIMS>& lhs,
                                const NURBSCurve<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /// \brief Returns a copy of the Bezier curve's control points
  CoordsVec getControlPoints() const { return m_controlPoints; }

  /// \brief Returns a copy of the Bezier curve's control points
  WeightsVec getWeights() const { return m_weights; }

  KnotsVec getKnots() const { return m_knots; }

  /// \brief Reverses the order of the Bezier curve's control points and weights
  void reverseOrientation()
  {
    const int deg = getDegree();
    const int mid = (deg + 1) / 2;
    for(int i = 0; i < mid; ++i)
    {
      axom::utilities::swap(m_controlPoints[i], m_controlPoints[deg - i]);
    }

    if(isRational())
    {
      for(int i = 0; i < mid; ++i)
      {
        axom::utilities::swap(m_weights[i], m_weights[deg - i]);
      }
    }

    // Reverse the orientation of the knots
    const int num_knots = getNumKnots();
    const int knot_mid = (num_knots + 1) / 2;
    for(int i = 0; i < knot_mid; ++i)
    {
      axom::utilities::swap(m_knots[i], m_knots[num_knots - i]);
    }

    // Replace each knot with 1 - knot_value
    for(int i = 0; i < num_knots; ++i)
    {
      m_knots[i] = 1.0 - m_knots[i];
    }
  }

  /// \brief Returns an axis-aligned bounding box containing the Bezier curve
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /// \brief Returns an oriented bounding box containing the Bezier curve
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
      {
        os << m_weights[p] << (p < ord ? ", " : "]");
      }
    }

    os << ", knots {";
    for(int i = 0; i < m_knots.size(); ++i)
    {
      os << m_knots[i] << (i < m_knots.size() - 1 ? ", " : "");
    }
    os << "}";

    return os;
  }

private:
  /// \brief Private function to make the knots uniform
  void makeKnotsUniform(int p, int npts)
  {
    m_knots.resize(npts + p + 1);

    // Knots for the clamped curve
    for(int i = 0; i < p + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + p - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - p - 1; ++i)
    {
      m_knots[p + 1 + i] = (i + 1.0) / (npts - p);
    }
  }

  /// \brief Private function to check if the NURBS curve is valid
  bool isValidNURBS() const
  {
    int p = getDegree();
    for(int i = 0; i < p + 1; ++i)
    {
      if(m_knots[i] != 0.0)
      {
        // First p+1 knots must be 0
        return false;
      }

      if(m_knots[m_knots.size() - p - 1 + i] != 1.0)
      {
        // Last p+1 knots must be 1
        return false;
      }
    }

    for(int i = 1; i < m_knots.size(); ++i)
    {
      if(m_knots[i] <= m_knots[i - 1])
      {
        // Knots must be non-decreasing
        return false;
      }
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

      return true;
    }

    return true;
  }

  /*!
   * \brief Returns the index of the knot span containing parameter t
   * 
   * Implementation adapted from Algorithm A2.1 on page 68 of "The NURBS Book"
   */
  axom::IndexType findSpan(double t) const
  {
    const auto n = m_controlPoints.size() - 1;

    if(t <= m_knots[0])
    {
      return getDegree();
    }

    if(m_knots[n] <= t)
    {
      return n;
    }

    // perform binary search on the knots,
    auto low = getDegree();
    auto high = n + 1;
    auto mid = (low + high) / 2;
    while(t < m_knots[mid] || t >= m_knots[mid + 1])
    {
      (t < m_knots[mid]) ? high = mid : low = mid;
      mid = (low + high) / 2;
    }

    return mid;
  }

  /*!
   * \brief Evaluates the NURBS basis functions for span at parameter value t
   * 
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   */
  axom::Array<T> calculateBasisFunctions(axom::IndexType span, T u) const
  {
    const int p = getDegree();

    axom::Array<T> N(p + 1);
    axom::Array<T> left(p + 1);
    axom::Array<T> right(p + 1);

    // Note: This implementation avoids division by zero and redundant computation
    // that might arise from a direct implementation of the recurrence relation
    // for basis functions. See "The NURBS Book" for details.
    N[0] = 1.0;
    for(int j = 1; j <= p; ++j)
    {
      left[j] = u - m_knots[span + 1 - j];
      right[j] = m_knots[span + j] - u;
      T saved = 0.0;
      for(int r = 0; r < j; ++r)
      {
        const T tmp = N[r] / (right[r + 1] + left[j - r]);
        N[r] = saved + right[r + 1] * tmp;
        saved = left[j - r] * tmp;
      }
      N[j] = saved;
    }
    return N;
  }
 
   /*!
   * \brief Evaluates the NURBS basis functions and derivatives for span at parameter value t
   * 
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   */
  void derivativeBasisFunctions(axom::IndexType span,
                                T t,
                                int n,
                                axom::Array<axom::Array<T>>& ders) const
  {
    const int p = getDegree();

    axom::Array<axom::Array<T>> ndu(p + 1), a(2);
    axom::Array<T> left(p + 1), right(p + 1);
    for(int j = 0; j <= p; j++)
    {
      ndu[j].resize(p + 1);
    }
    for(int j = 0; j <= n; j++)
    {
      ders[j].resize(p + 1);
    }
    a[0].resize(p + 1);
    a[1].resize(p + 1);

    ndu[0][0] = 1.;
    for(int j = 1; j <= p; j++)
    {
      left[j] = t - m_knots[span + 1 - j];
      right[j] = m_knots[span + j] - t;
      T saved = 0.0;
      for(int r = 0; r < j; r++)
      {
        // lower triangle
        ndu[j][r] = right[r + 1] + left[j - r];
        T temp = ndu[r][j - 1] / ndu[j][r];
        // upper triangle
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    // Load basis functions
    for(int j = 0; j <= p; j++)
    {
      ders[0][j] = ndu[j][p];
    }

    // This section computes the derivatives (Eq. [2.9])

    // Loop over function index.
    for(int r = 0; r <= p; r++)
    {
      int s1 = 0, s2 = 1;  // Alternate rows in array a
      a[0][0] = 1.;
      // Loop to compute kth derivative
      for(int k = 1; k <= n; k++)
      {
        T d = 0.;
        int rk = r - k;
        int pk = p - k;
        if(r >= k)
        {
          a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
          d = a[s2][0] * ndu[rk][pk];
        }
        int j1 = (rk >= -1) ? 1 : -rk;
        int j2 = (r - 1 <= pk) ? (k - 1) : (p - r);
        for(int j = j1; j <= j2; j++)
        {
          a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j];
          d += a[s2][j] * ndu[rk + j][pk];
        }
        if(r <= pk)
        {
          a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r];
          d += a[s2][k] * ndu[r][pk];
        }
        ders[k][r] = d;
        // Switch rows
        std::swap(s1, s2);
      }
    }
    // Multiply through by the correct factors (Eq. [2.9])
    T r = static_cast<T>(p);
    for(int k = 1; k <= n; k++)
    {
      for(int j = 0; j <= p; j++)
      {
        ders[k][j] *= r;
      }
      r *= static_cast<T>(p - k);
    }
  }

  CoordsVec m_controlPoints;
  WeightsVec m_weights;
  KnotsVec m_knots;
};

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_NURBSCURVE_HPP_