// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/quest/readers/C2CReader.hpp"

#ifndef AXOM_USE_C2C
  #error C2CReader should only be included when Axom is configured with C2C
#endif

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/fmt.hpp"

// These macros enable some debugging utilities for linearization.
// #define AXOM_DEBUG_LINEARIZE_VERBOSE
// #define AXOM_DEBUG_SOLVEU_VERBOSE
// #define AXOM_DEBUG_WRITE_ERROR_CURVES
// #define AXOM_DEBUG_WRITE_LINES

/// Overload to format a c2c::Point using fmt
template <>
struct axom::fmt::formatter<c2c::Point> : ostream_formatter
{ };

namespace axom
{
namespace quest
{
//---------------------------------------------------------------------------
// TODO: Use MFEM for quadrature
struct QuadratureElement
{
    double w;
    double x;
};

// Quadrature table from: https://pomax.github.io/bezierinfo/legendre-gauss.html
static const std::vector<QuadratureElement> quad64{
/*1 */{0.0486909570091397,-0.0243502926634244},
/*2 */{0.0486909570091397, 0.0243502926634244},
/*3 */{0.0485754674415034,-0.0729931217877990},
/*4 */{0.0485754674415034, 0.0729931217877990},
/*5 */{0.0483447622348030,-0.1214628192961206},
/*6 */{0.0483447622348030, 0.1214628192961206},
/*7 */{0.0479993885964583,-0.1696444204239928},
/*8 */{0.0479993885964583, 0.1696444204239928},
/*9 */{0.0475401657148303,-0.2174236437400071},
/*10*/{0.0475401657148303, 0.2174236437400071},
/*11*/{0.0469681828162100,-0.2646871622087674},
/*12*/{0.0469681828162100, 0.2646871622087674},
/*13*/{0.0462847965813144,-0.3113228719902110},
/*14*/{0.0462847965813144, 0.3113228719902110},
/*15*/{0.0454916279274181,-0.3572201583376681},
/*16*/{0.0454916279274181, 0.3572201583376681},
/*17*/{0.0445905581637566,-0.4022701579639916},
/*18*/{0.0445905581637566, 0.4022701579639916},
/*19*/{0.0435837245293235,-0.4463660172534641},
/*20*/{0.0435837245293235, 0.4463660172534641},
/*21*/{0.0424735151236536,-0.4894031457070530},
/*22*/{0.0424735151236536, 0.4894031457070530},
/*23*/{0.0412625632426235,-0.5312794640198946},
/*24*/{0.0412625632426235, 0.5312794640198946},
/*25*/{0.0399537411327203,-0.5718956462026340},
/*26*/{0.0399537411327203, 0.5718956462026340},
/*27*/{0.0385501531786156,-0.6111553551723933},
/*28*/{0.0385501531786156, 0.6111553551723933},
/*29*/{0.0370551285402400,-0.6489654712546573},
/*30*/{0.0370551285402400, 0.6489654712546573},
/*31*/{0.0354722132568824,-0.6852363130542333},
/*32*/{0.0354722132568824, 0.6852363130542333},
/*33*/{0.0338051618371416,-0.7198818501716109},
/*34*/{0.0338051618371416, 0.7198818501716109},
/*35*/{0.0320579283548516,-0.7528199072605319},
/*36*/{0.0320579283548516, 0.7528199072605319},
/*37*/{0.0302346570724025,-0.7839723589433414},
/*38*/{0.0302346570724025, 0.7839723589433414},
/*39*/{0.0283396726142595,-0.8132653151227975},
/*40*/{0.0283396726142595, 0.8132653151227975},
/*41*/{0.0263774697150547,-0.8406292962525803},
/*42*/{0.0263774697150547, 0.8406292962525803},
/*43*/{0.0243527025687109,-0.8659993981540928},
/*44*/{0.0243527025687109, 0.8659993981540928},
/*45*/{0.0222701738083833,-0.8893154459951141},
/*46*/{0.0222701738083833, 0.8893154459951141},
/*47*/{0.0201348231535302,-0.9105221370785028},
/*48*/{0.0201348231535302, 0.9105221370785028},
/*49*/{0.0179517157756973,-0.9295691721319396},
/*50*/{0.0179517157756973, 0.9295691721319396},
/*51*/{0.0157260304760247,-0.9464113748584028},
/*52*/{0.0157260304760247, 0.9464113748584028},
/*53*/{0.0134630478967186,-0.9610087996520538},
/*54*/{0.0134630478967186, 0.9610087996520538},
/*55*/{0.0111681394601311,-0.9733268277899110},
/*56*/{0.0111681394601311, 0.9733268277899110},
/*57*/{0.0088467598263639,-0.9833362538846260},
/*58*/{0.0088467598263639, 0.9833362538846260},
/*59*/{0.0065044579689784,-0.9910133714767443},
/*60*/{0.0065044579689784, 0.9910133714767443},
/*61*/{0.0041470332605625,-0.9963401167719553},
/*62*/{0.0041470332605625, 0.9963401167719553},
/*63*/{0.0017832807216964,-0.9993050417357722},
/*64*/{0.0017832807216964, 0.9993050417357722}
};

//---------------------------------------------------------------------------
/*!
 * \brief Transform a 2D point and return a 2D point.
 *
 * \param transform A 4x4 transformation matrix.
 * \param pt The input 2D point.
 *
 * \return The transformed 2D point.
 */
inline primal::Point<double, 2>
transformPoint(const numerics::Matrix<double> &transform,
                 const primal::Point<double, 2> &pt)
{
  // Turn the point2 into a vec4.
  double pt4[] = {0., 0., 0., 1.}, newpt4[] = {0., 0., 0., 1.};
  pt4[0] = pt[0];
  pt4[1] = pt[1];

  // Transform the point.
  axom::numerics::matrix_vector_multiply(transform, pt4, newpt4);

  // Make a new 2D point from the 2D components.
  return primal::Point<double, 2>{newpt4[0], newpt4[1]};
}

//---------------------------------------------------------------------------
/*!
 * \brief Append new segments to the mint mesh, performing any point-welding
 *        that may be needed.
 *
 * \param mesh The input mesh to which we're adding segments/points.
 * \param pts The input set of points.
 * \param EPS_SQ A point-welding tolerance.
 */
/*
static void
appendPoints(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
             std::vector<primal::Point<double, 2>> &pts,
             double EPS_SQ)
{
  // Check for simple vertex welding opportunities at endpoints of newly interpolated points
  {
    int numNodes = mesh->getNumberOfNodes();
    if(numNodes > 0)  // this is not the first Piece
    {
      primal::Point<double, 2> meshPt;
      // Fix start point if necessary; check against most recently added vertex in mesh
      mesh->getNode(numNodes - 1, meshPt.data());
      if(primal::squared_distance(pts[0], meshPt) < EPS_SQ)
      {
        pts[0] = meshPt;
      }

      // Fix end point if necessary; check against 0th vertex in mesh
      const int endIdx = pts.size() - 1;
      mesh->getNode(0, meshPt.data());
      if(primal::squared_distance(pts[endIdx], meshPt) < EPS_SQ)
      {
        pts[endIdx] = meshPt;
      }
    }
    else  // This is the first, and possibly only span, check its endpoint, fix if necessary
    {
      int endIdx = pts.size() - 1;
      if(primal::squared_distance(pts[0], pts[endIdx]) < EPS_SQ)
      {
        pts[endIdx] = pts[0];
      }
    }
  }

  // Add the new points and segments to the mesh, respecting welding checks from previous block
  {
    const int startNode = mesh->getNumberOfNodes();
    const int numNewNodes = pts.size();
    mesh->reserveNodes(startNode + numNewNodes);

    for(int i = 0; i < numNewNodes; ++i)
    {
      mesh->appendNode(pts[i][0], pts[i][1]);
    }

    const int startCell = mesh->getNumberOfCells();
    const int numNewSegments = pts.size() - 1;
    mesh->reserveCells(startCell + numNewSegments);
    for(int i = 0; i < numNewSegments; ++i)
    {
      IndexType seg[2] = {startNode + i, startNode + i + 1};
      mesh->appendCell(seg, mint::SEGMENT);
    }
  }
}
*/

class Segments
{
public:
  using PointType = primal::Point<double, 2>;

  struct segment
  {
    double    u;              // The u parameter for the start of the segment.
    PointType point;          // The point at parameter u.
    double    length;         // Distance from point to next segment's point.

    double    next_u;         // Best place to split this segment.
    double    next_length[2]; // left/right length values at next_u.
  };

  Segments() : segments()
  {
  }

  void reserve(size_t n)
  {
    segments.reserve(n);
  }

  size_t size() const
  {
    return segments.size();
  }

  size_t next_longest() const
  {
    size_t maxIdx = 0;
    double max_length = 0.;
    size_t n = segments.size();
    for(size_t i = 0; i < n; i++)
    {
      double next_length = segments[i].next_length[0] + segments[i].next_length[1];
      if(next_length > max_length)
      {
        max_length = next_length;
        maxIdx = i;
      }
    }
    return maxIdx;
  }

  segment &operator[](size_t idx)
  {
    return segments[idx];
  }

  const segment &operator[](size_t idx) const
  {
    return segments[idx];
  }

  void push_back(const segment &seg)
  {
    segments.push_back(seg);
  }

  void insert(size_t idx, const segment &seg)
  {
    segments.insert(segments.begin() + idx, seg);
  }

private:
  std::vector<segment> segments;
};

static void
appendPoints(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
             Segments &S,
             double EPS_SQ)
{
  // Check for simple vertex welding opportunities at endpoints of newly interpolated points
  {
    int numNodes = mesh->getNumberOfNodes();
    if(numNodes > 0)  // this is not the first Piece
    {
      primal::Point<double, 2> meshPt;
      // Fix start point if necessary; check against most recently added vertex in mesh
      mesh->getNode(numNodes - 1, meshPt.data());
      if(primal::squared_distance(S[0].point, meshPt) < EPS_SQ)
      {
        S[0].point = meshPt;
      }

      // Fix end point if necessary; check against 0th vertex in mesh
      const int endIdx = S.size() - 1;
      mesh->getNode(0, meshPt.data());
      if(primal::squared_distance(S[endIdx].point, meshPt) < EPS_SQ)
      {
        S[endIdx].point = meshPt;
      }
    }
    else  // This is the first, and possibly only span, check its endpoint, fix if necessary
    {
      int endIdx = S.size() - 1;
      if(primal::squared_distance(S[0].point, S[endIdx].point) < EPS_SQ)
      {
        S[endIdx].point = S[0].point;
      }
    }
  }

  // Add the new points and segments to the mesh, respecting welding checks from previous block
  {
    const int startNode = mesh->getNumberOfNodes();
    const int numNewNodes = S.size();
    mesh->reserveNodes(startNode + numNewNodes);

    for(int i = 0; i < numNewNodes; ++i)
    {
      mesh->appendNode(S[i].point[0], S[i].point[1]);
    }

    const int startCell = mesh->getNumberOfCells();
    const int numNewSegments = S.size() - 1;
    mesh->reserveCells(startCell + numNewSegments);
    for(int i = 0; i < numNewSegments; ++i)
    {
      IndexType seg[2] = {startNode + i, startNode + i + 1};
      mesh->appendCell(seg, mint::SEGMENT);
    }
  }
}

//---------------------------------------------------------------------------
/*!
 * \brief Helper class for interpolating points on a NURBS curve
 *
 * This class was adapted from a similar class in the C2C library's testing framework.
 * The algorithms are adapted from: Piegl and Tiller's "The NURBS Book", 2nd Ed,  (Springer, 1997) 
 */
struct NURBSInterpolator
{
  using BasisVector = std::vector<double>;
  using PointType = primal::Point<double, 2>;

  // EPS is used for computing span intervals
  NURBSInterpolator(const c2c::NURBSData& curve, double EPS = 1E-9)
    : m_curve(curve)
  {
    const int p = m_curve.order - 1;
    const int knotSize = m_curve.knots.size();

    AXOM_UNUSED_VAR(p);  // silence warnings in release configs
    AXOM_UNUSED_VAR(knotSize);

    SLIC_ASSERT(p >= 1);
    SLIC_ASSERT(knotSize >= 2 * p);
    computeSpanIntervals(EPS);
  }

  /// Helper function to compute the start and end parametric coordinates of each knot span
  void computeSpanIntervals(double EPS)
  {
    using axom::utilities::isNearlyEqual;
    const auto& U = m_curve.knots;
    const int knotSize = U.size();
    const int p = m_curve.order - 1;
    const int n = knotSize - 1 - p - 1;

    for(int i = p; i < n + 1; ++i)
    {
      const double left = U[i];
      const double right = U[i + 1];
      if(!isNearlyEqual(left, right, EPS))
      {
        m_spanIntervals.push_back(std::make_pair(left, right));
      }
    }
  }

  /*!
   * \brief Checks that the knots of the NURBS curve are closed
   *
   * The knots vector is closed when it begins with \a m_curve.order equal knot values
   * and ends with the same number of equal knot values
   */
  bool areKnotsClosed(double EPS = 1E-9) const
  {
    using axom::utilities::isNearlyEqual;
    const auto& U = m_curve.knots;

    int p = m_curve.order - 1;

    const double startKnot = U[0];
    for(int i = 1; i <= p; ++i)
    {
      if(!isNearlyEqual(startKnot, U[i], EPS))
      {
        return false;
      }
    }

    const int endIndex = U.size() - 1;
    const double endKnot = U[endIndex];
    for(int i = 1; i <= p; ++i)
    {
      const int index = endIndex - p - 1;
      if(!isNearlyEqual(endKnot, U[index], EPS))
      {
        return false;
      }
    }

    return true;
  }

  int numSpans() const { return m_spanIntervals.size(); }

  double startParameter(int span = -1) const
  {
    const bool inRange = span >= 0 && span < numSpans();
    return inRange ? m_spanIntervals[span].first : m_curve.knots[0];
  }

  double endParameter(int span = -1) const
  {
    const bool inRange = span >= 0 && span < numSpans();
    return inRange ? m_spanIntervals[span].second
                   : m_curve.knots[m_curve.knots.size() - 1];
  }

  /*!
   * \brief Finds the index of the knot span containing parameter \a u
   * 
   * Implementation adapted from Algorithm A2.1 on page 68 of "The NURBS Book"
   */
  int findSpan(double u) const
  {
    const auto& U = m_curve.knots;
    const int knotSize = U.size();
    const int p = m_curve.order - 1;
    const int n = knotSize - 1 - p - 1;

    if(U[n] <= u && u <= U[n + 1])
    {
      return n;
    }

    // perform binary search on the knots
    int low = p;
    int high = n + 1;
    int mid = (low + high) / 2;
    while(u < U[mid] || u >= U[mid + 1])
    {
      (u < U[mid]) ? high = mid : low = mid;
      mid = (low + high) / 2;
    }
    return mid;
  }

  /*!
   * \brief Evaluates the B-spline basis functions for span \a span at parameter value \a u 
   * 
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   */
  BasisVector calculateBasisFunctions(int span, double u) const
  {
    const int p = m_curve.order - 1;
    const auto& U = m_curve.knots;

    BasisVector N(p + 1);
    BasisVector left(p + 1);
    BasisVector right(p + 1);

    // Note: This implementation avoids division by zero and redundant computation
    // that might arise from a direct implementation of the recurrence relation
    // for basis functions. See "The NURBS Book" for details.
    N[0] = 1.0;
    for(int j = 1; j <= p; ++j)
    {
      left[j] = u - U[span + 1 - j];
      right[j] = U[span + j] - u;
      double saved = 0.0;
      for(int r = 0; r < j; ++r)
      {
        double temp = N[r] / (right[r + 1] + left[j - r]);
        N[r] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      N[j] = saved;
    }
    return N;
  }

  /*!
   * \brief Finds the point on the curve at parameter \a u
   *
   * Adapted from Algorithm A4.1 on page 124 of "The NURBS Book"
   */
  PointType at(double u) const
  {
    using GrassmanPoint = primal::Point<double, 3>;
    GrassmanPoint cw {0.0};

    const auto span = findSpan(u);
    const auto N = calculateBasisFunctions(span, u);
    const int p = m_curve.order - 1;
    for(int j = 0; j <= p; ++j)
    {
      const int offset = span - p + j;
      const double weight = m_curve.weights[offset];
      const auto& controlPoint = m_curve.controlPoints[offset];

      cw[0] += N[j] * weight * controlPoint.getZ().getValue();
      cw[1] += N[j] * weight * controlPoint.getR().getValue();
      cw[2] += N[j] * weight;
    }

    // Return projected point
    // All units should have been normalized by c2c::toNurbs(piece, units)
    return PointType {cw[0] / cw[2], cw[1] / cw[2]};
  }

  /*!
   * \brief Evaluates the B-spline derivative basis functions for span \a span 
   *        at parameter value \a u
   * 
   * \param span The span of interest.
   * \param u The u value at which to evaluate derivatives.
   * \param n The number of derivatives to compute.
   * \param[out] ders Store the basis functions.
   *
   * \note ders[0] stores the basis functions. ders[1] stores the 1st derivative
   *       basis functions, etc.
   *
   * Implementation adapted from Algorithm A2.3 on pp. 72-73 of "The NURBS Book".
   */
  void derivativeBasisFunctions(int span,
                                double u,
                                int n,
                                std::vector<BasisVector>& ders) const
  {
    const int p = m_curve.order - 1;
    const auto& U = m_curve.knots;

    std::vector<BasisVector> ndu(p + 1), a(2);
    BasisVector left(p + 1), right(p + 1);
    for(int j = 0; j <= p; j++) ndu[j].resize(p + 1);
    for(int j = 0; j <= n; j++) ders[j].resize(p + 1);
    a[0].resize(p + 1);
    a[1].resize(p + 1);

    ndu[0][0] = 1.;
    for(int j = 1; j <= p; j++)
    {
      left[j] = u - U[span + 1 - j];
      right[j] = U[span + j] - u;
      double saved = 0.0;
      for(int r = 0; r < j; r++)
      {
        // lower triangle
        ndu[j][r] = right[r + 1] + left[j - r];
        double temp = ndu[r][j - 1] / ndu[j][r];
        // upper triangle
        ndu[r][j] = saved + right[r + 1] * temp;
        saved = left[j - r] * temp;
      }
      ndu[j][j] = saved;
    }
    // Load basis functions
    for(int j = 0; j <= p; j++) ders[0][j] = ndu[j][p];

    // This section computes the derivatives (Eq. [2.9])

    // Loop over function index.
    for(int r = 0; r <= p; r++)
    {
      int s1 = 0, s2 = 1;  // Alternate rows in array a
      a[0][0] = 1.;
      // Loop to compute kth derivative
      for(int k = 1; k <= n; k++)
      {
        double d = 0.;
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
    double r = static_cast<double>(p);
    for(int k = 1; k <= n; k++)
    {
      for(int j = 0; j <= p; j++)
        ders[k][j] *= r;
      r *= static_cast<double>(p - k);
    }
  }

  /*!
   * \brief Evaluates derivatives at the provided value \a u.
   *
   * \param u The u value at which to evaluate derivatives.
   * \param d The number of derivatives to create (1 for 1st derivative,
   *          2 for both 1st and 2nd derivatives, etc.)
   * \param[out] CK An array of points to contain the output derivatives.
   *
   * Implementation adapted from Algorithm A3.2 on p. 93 of "The NURBS Book".
   */
  void derivativesAt(double u, int d, PointType* CK) const
  {
    const int p = m_curve.order - 1;
    int du = std::min(d, p);

    const auto span = findSpan(u);
    std::vector<BasisVector> N(d + 1);
    derivativeBasisFunctions(span, u, d, N);

    for(int k = 1; k <= du; k++)
    {
      double x = 0., y = 0.;
      for(int j = 0; j <= p; j++)
      {
        int offset = span - p + j;
        // TODO: We may need to include weights.
        x = x + N[k][j] * m_curve.controlPoints[offset].getZ().getValue();
        y = y + N[k][j] * m_curve.controlPoints[offset].getR().getValue();
      }
      CK[k - 1] = PointType {x, y};
    }
  }

  /*!
   * \brief Evaluates the B-spline curvature at parameter value \a u
   * 
   * \param u The parameter value.
   *
   * \return The curvature value at u.
   */
  double curvature(double u) const
  {
    // Evaluate 1st and 2nd derivatives at u.
    PointType derivs[2];
    derivativesAt(u, 2, derivs);
    const PointType& D1 = derivs[0];
    const PointType& D2 = derivs[1];
#if 1
    double xp   = D1.data()[0]; // x'
    double xpp  = D2.data()[0]; // x''

    double yp   = D1.data()[1]; // y'
    double ypp  = D2.data()[1]; // y''

    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - y'x'') / pow(x'x' + y'y', 3./2.)
    double xp2_plus_yp2 = xp * xp + yp * yp;
    return (xp * ypp - yp * xpp) / pow(xp2_plus_yp2, 3. / 2.);
#else
    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - y'x'') / pow(x'x' + y'y', 3./2.)
    double numerator =
      (D1.data()[0] * D2.data()[1]) - (D2.data()[0] * D1.data()[1]);
    double D1mag2 = (D1.data()[0] * D1.data()[0]) + (D1.data()[1] * D1.data()[1]);
    double one_over_denominator = pow(D1mag2, -3. / 2.);
    return numerator * one_over_denominator;
#endif
  }

  /*!
   * \brief Evaluates the B-spline curvature derivatives (up to 2nd) at
   *        parameter value \a u
   * 
   * \param[in] u The parameter value.
   * \param[in] d The number of derivatives to compute (1=1st deriv, 2=1st & 2nd derivs)
   * \param[out] ders An array that will contain the curvature derivatives.
   */
  void curvatureDerivatives(double u, int d, double *ders) const
  {
    // Evaluate 1st, 2nd, 3rd curve derivatives at u.
    PointType derivs[3];
    derivativesAt(u, 3, derivs);
    const PointType& D1 = derivs[0];
    const PointType& D2 = derivs[1];
    const PointType& D3 = derivs[2];

    double xp   = D1.data()[0]; // x'
    double xpp  = D2.data()[0]; // x''
    double xppp = D3.data()[0]; // x'''

    double yp   = D1.data()[1]; // y'
    double ypp  = D2.data()[1]; // y''
    double yppp = D3.data()[1]; // y'''

    // 1st derivative of curvature.
    double xp2_plus_yp2 = xp * xp + yp * yp;
    double A = -3. * (xp * ypp - yp * xpp) * 2 *(xp * xpp + yp * ypp);
    double B = 2. * pow(xp2_plus_yp2, 5. / 2.);
    double C = xp * yppp - yp * xppp;
    double D = pow(xp2_plus_yp2, 3. / 2.);
    ders[0] = A/B + C/D;

    if(d >= 2)
    {
      // 2nd derivative of curvature.
      double E = 15. * (-yp*xpp + xp * ypp) * pow(2. * xp * xpp + 2. * yp * ypp, 2.) / (4. * pow(xp2_plus_yp2, 7. / 2.));
      double F = 3. * (2. * xp * xpp + 2. * yp * ypp) * (-yp * xppp + xp*yppp) / pow(xp2_plus_yp2, 5. / 2.);
      double G = 3. * (-yp*xpp + xp*ypp)* (2.*(xpp*xpp)+2.*(ypp*ypp) + 2.*xp*xppp + 2.*yp*yppp) / (2. * pow(xp2_plus_yp2, 5. / 2.));
      double H = (-ypp*xppp + xpp*yppp) / pow(xp2_plus_yp2, 3. / 2.);

      ders[1] = E - F - G + H;
    }
  }

  /* \brief Looks at the curvature function and returns intervals that
   *        need to be sampled. We'd like to include curvature values
   *        at any of the u values returned here. We try to include
   *        starting endpoints and any u values within that cause
   *        curvature extrema.
   *
   * \brief umin The minimum u value for the starting interval.
   * \brief umax The maximum u value for the starting interval.
   */
  std::vector<double> curvatureIntervals(double umin, double umax) const
  {
    using axom::utilities::lerp;
    constexpr int numSamples = 10;
    constexpr double denom = static_cast<double>(numSamples - 1);
    // Sample curvature 1st derivatives across the interval.
    std::vector<double> uValues(numSamples), curvDeriv(numSamples);
    for(int i = 0; i < numSamples; i++)
    {
       uValues[i] = lerp(umin, umax, i / denom);
       curvatureDerivatives(uValues[i], 1, &curvDeriv[i]);
    }
    // Build intervals.
    std::vector<double> intervals;
    intervals.push_back(umin);
    for(int i = 1; i < numSamples; i++)
    {
      int s0 = (curvDeriv[i - 1] < 0) ? -1 : 1;
      int s1 = (curvDeriv[i] < 0) ? -1 : 1;
      if(s0 != s1)
      {
        // There is a derivative sign change in this interval so there must be
        // a zero crossing, indicating a max or a min in the curvature.
        double newu = 0.;
        if(solveCurvature(uValues[i-1], uValues[i], newu))
          intervals.push_back(newu);
      }
    }
    intervals.push_back(umax);
    return std::move(intervals);
  }

  bool solveCurvature(double umin, double umax, double &u) const
  {
    constexpr double EPS = 1.e-6;
    constexpr int MAX_ITERATIONS = 32;
    bool solved = false;

    // Use the middle of the range as a 1st guess.
    u = (umin + umax) * 0.5;

    // Find a zero-crossing in the curvature derivative using Newton's Method.
    for(int iteration = 0; iteration < MAX_ITERATIONS && !solved; iteration++)
    {
      double ders[2];
      curvatureDerivatives(u, 2, ders);

      if(axom::utilities::isNearlyEqual(ders[0], 0., EPS))
      {
std::cout << iteration << ": Stopping: u=" << u << ", cp=" << ders[0] << ", cpp=" << ders[1] << std::endl;
        solved = true;
        break;
      }
      else
      {
        // Next step.
        u = u - ders[0] / ders[1];
std::cout << iteration << ": next u=" << u << std::endl;

#if 1
        if(u <= umin)
        {
          std::cout << iteration << ": ERROR new u " << u << " is less than umin " << umin << std::endl;
          u = umin;
          break;
        }
        if(u >= umax)
        {
          std::cout << iteration << ": ERROR new u " << u << " is greater than umax " << umax << std::endl;
          u = umax;
          break;
        }
#endif
      }
    }
    return solved;
  }
  

  /* \brief Solve for a zero crossing for the curvature derivative and use the
            determined u value to return the curvature there.
   *
   */
  double curvatureExtreme(double umin, double umax, double &u) const
  {
    constexpr double EPS = 1.e-6;
    double c = 0.;

    // Use the middle of the range as a 1st guess.
    u = umin; //(umin + umax) * 0.5;

    // See whether there is enough curvature to bother with.
    double c0 = curvature(umin);
    double cmid = curvature(u);
    double c1 = curvature(umax);
    double dc0c1 = fabs(c1 - c0);
    double dc0cmid = fabs(cmid - c0);
    double dcmidc1 = fabs(c1 - cmid);
    if(dc0c1 > EPS || dc0cmid > EPS || dcmidc1 > EPS)
    {
std::cout <<"!! Solve for curvature extremes" << std::endl;
std::cout <<"umin=" << umin << std::endl;
std::cout <<"umax=" << umax << std::endl;
std::cout <<"u=" << u << std::endl;

      bool keepGoing = true;
      // Find a zero-crossing in the curvature derivative using Newton's Method.
      int iteration = 0;
      while(keepGoing && iteration < 20)
      {
        double ders[2];
        curvatureDerivatives(u, 2, ders);

std::cout << iteration << ": u=" << u << ", cp=" << ders[0] << ", cpp=" << ders[1] << std::endl;
        if(axom::utilities::isNearlyEqual(ders[0], 0., EPS))
        {
std::cout << iteration << ": Stopping: u=" << u << ", cp=" << ders[0] << ", cpp=" << ders[1] << std::endl;
          keepGoing = false;
        }
        else
        {

          // Next step.
          u = u - ders[0] / ders[1];
std::cout << iteration << ": next u=" << u << std::endl;

#if 1
          if(u <= umin)
          {
            std::cout << iteration << ": ERROR new u " << u << " is less than umin " << umin << std::endl;
            u = umin;
            keepGoing = false;
          }
          if(u >= umax)
          {
            std::cout << iteration << ": ERROR new u " << u << " is greater than umax " << umax << std::endl;
            u = umax;
            keepGoing = false;
          }
#endif
        }
        iteration++;
      }

      // Check to see whether we're close to the endpoints. If so, use the
      // endpoints.
      if(axom::utilities::isNearlyEqual(u, umin, EPS))
      {
        u = umin;
        c = curvature(u);
      }
      else if(axom::utilities::isNearlyEqual(u, umax, EPS))
      {
        u = umax;
        c = curvature(u);
      }
    }
    else
    {
      c = c0;
      u = umin;
    }

    return c;
  }

  /*!
   * \brief Compute the revolved volume of the curve across its entire
   *        parametric interval [0,1] using quadrature.
   *
   * \note We include a transformation of the points so any transforms
   *       can figure into the integration.
   *
   * \return The revolved curve volume.
   */
  double revolvedVolume(const numerics::Matrix<double> &transform) const
  {
    constexpr double a = 0, b = 1;
    constexpr double scale = M_PI * ((b - a) / 2.);

    // Make a transform with no translation. We use this to transform
    // the derivative since we want to permit scaling and rotation but
    // translating it does not make sense.
    numerics::Matrix<double> transform2 = transform;
    transform2(0, 3) = 0.;
    transform2(1, 3) = 0.;
    transform2(2, 3) = 0.;

    // Approximate the integral "Int pi*x'(u)*y(u)^2du" using quadrature.
    double sum = 0.;
    size_t n = quad64.size();
    for(size_t i = 0; i < n; i++)
    {
      // Map quad point x value [-1,1] to [a,b].
      double u = quad64[i].x * ((b - a) / 2.) + ((b + a) / 2.);

      // Compute y(u) to get radius
      PointType p_u = at(u);
      PointType p_uT = transformPoint(transform, p_u);
      double r = p_uT[1];

      // Compute x'(u)
      PointType xprime;
      derivativesAt(u, 1, &xprime);
      PointType xprimeT = transformPoint(transform2, xprime);
      double xp = xprimeT[0];

      // Accumulate weight times dx*r^2.
      sum += quad64[i].w * xp * (r * r);
    }
    double vol = fabs(scale * sum);
    return vol;
  }

private:
  const c2c::NURBSData& m_curve;
  std::vector<std::pair<double, double>> m_spanIntervals;
};

void C2CReader::clear() { m_nurbsData.clear(); }

int C2CReader::read()
{
  SLIC_WARNING_IF(m_fileName.empty(), "Missing a filename in C2CReader::read()");

  using axom::utilities::string::endsWith;

  int ret = 1;

  if(endsWith(m_fileName, ".contour"))
  {
    ret = readContour();
  }
  else if(endsWith(m_fileName, ".assembly"))
  {
    SLIC_WARNING("Input was an assembly! This is not currently supported");
  }
  else
  {
    SLIC_WARNING("Not a valid c2c file");
  }

  return ret;
}

int C2CReader::readContour()
{
  c2c::Contour contour = c2c::parseContour(m_fileName);

  SLIC_INFO(
    fmt::format("Loading contour with {} pieces", contour.getPieces().size()));

  for(auto* piece : contour.getPieces())
  {
    m_nurbsData.emplace_back(c2c::toNurbs(*piece, m_lengthUnit));
  }

  return 0;
}

void C2CReader::log()
{
  std::stringstream sstr;

  sstr << fmt::format("The contour has {} pieces\n", m_nurbsData.size());

  int index = 0;
  for(const auto& nurbs : m_nurbsData)
  {
    sstr << fmt::format("Piece {}\n{{", index);
    sstr << fmt::format("\torder: {}\n", nurbs.order);
    sstr << fmt::format("\tknots: {}\n", fmt::join(nurbs.knots, " "));
    sstr << fmt::format("\tknot spans: {}\n",
                        NURBSInterpolator(nurbs).numSpans());
    sstr << fmt::format("\tweights: {}\n", fmt::join(nurbs.weights, " "));
    sstr << fmt::format("\tcontrol points: {}\n",
                        fmt::join(nurbs.controlPoints, " "));
    sstr << "}\n";
    ++index;
  }

  SLIC_INFO(sstr.str());
}

void C2CReader::getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                              int segmentsPerKnotSpan)
{
  std::vector<double> d1, d2, uvec, curv, sp;
  getLinearMesh(mesh, segmentsPerKnotSpan, d1, d2, uvec, curv, sp);
}

//----------------------------------------------------------------------------
/*static void
write_lines(const std::string &filename, std::vector<primal::Point<double, 2>> &pts)
{
    FILE *f = fopen(filename.c_str(), "wt");
    fprintf(f, "# vtk DataFile Version 4.2\n");
    fprintf(f, "vtk output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET POLYDATA\n");
    fprintf(f, "FIELD FieldData 2\n");
    fprintf(f, "CYCLE 1 1 int\n");
    fprintf(f, "1\n");
    fprintf(f, "TIME 1 1 double\n");
    fprintf(f, "1.0\n");
    // Write points
    int npts = pts.size();
    fprintf(f, "POINTS %d float\n", npts);
    for(int i = 0; i < npts; i += 3)
    {
        fprintf(f, "%1.16lf %1.16lf 0. ", pts[i][0], pts[i][1]);
        if((i+1) < npts)
            fprintf(f, "%1.16lf %1.16lf 0. ", pts[i+1][0], pts[i+1][1]);
        if((i+2) < npts)
            fprintf(f, "%1.16lf %1.16lf 0. ", pts[i+2][0], pts[i+2][1]);
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

    int nspans = npts - 1;

    // Write ncells
    fprintf(f, "LINES %d %d\n", nspans, 3 * nspans);
    for(int ispan = 0; ispan < nspans; ispan++)
    {
        fprintf(f, "2 %d %d\n", ispan, ispan+1);
    }

    fclose(f);
}*/

static void
write_lines(const std::string &filename, Segments &S)
{
    FILE *f = fopen(filename.c_str(), "wt");
    fprintf(f, "# vtk DataFile Version 4.2\n");
    fprintf(f, "vtk output\n");
    fprintf(f, "ASCII\n");
    fprintf(f, "DATASET POLYDATA\n");
    fprintf(f, "FIELD FieldData 2\n");
    fprintf(f, "CYCLE 1 1 int\n");
    fprintf(f, "1\n");
    fprintf(f, "TIME 1 1 double\n");
    fprintf(f, "1.0\n");
    // Write points
    int npts = S.size();
    fprintf(f, "POINTS %d float\n", npts);
    for(int i = 0; i < npts; i += 3)
    {
        fprintf(f, "%1.16lf %1.16lf 0. ", S[i].point[0], S[i].point[1]);
        if((i+1) < npts)
            fprintf(f, "%1.16lf %1.16lf 0. ", S[i+1].point[0], S[i+1].point[1]);
        if((i+2) < npts)
            fprintf(f, "%1.16lf %1.16lf 0. ", S[i+2].point[0], S[i+2].point[1]);
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

    int nspans = npts - 1;

    // Write ncells
    fprintf(f, "LINES %d %d\n", nspans, 3 * nspans);
    for(int ispan = 0; ispan < nspans; ispan++)
    {
        fprintf(f, "2 %d %d\n", ispan, ispan+1);
    }

    fclose(f);
}

//---------------------------------------------------------------------------
template <typename T>
static std::ostream &operator << (std::ostream &os, const std::vector<T> &vec)
{
    os << "{";
    for(const auto &v : vec)
        os << v << ", ";
    os << "}";
    return os;
}

static std::ostream &operator << (std::ostream &os, const primal::Point<double, 2> &pt)
{
    os << "(" << pt[0] << ", " << pt[1] << ")";
    return os;
}

//---------------------------------------------------------------------------
void C2CReader::getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                              const numerics::Matrix<double> &transform,
                              double percentError,
                              double &revolvedVolume)
{
  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT,
                "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(percentError <= 0.,
                "C2C reader: percentError must be greater than zero.");
  SLIC_ERROR_IF(percentError >= 1.,
                "C2C reader: percentError must be less than one.");

  using PointType = primal::Point<double, 2>;
  using PointsArray = std::vector<PointType>;

  /*!
   * \brief Checks curve lengths against tolerances and determines whether more
   *        refinement is needed.
   *
   * \param len The curve length that has been calculated.
   * \param maxlen The upper bound curve length that was calculated with a lot of segments.
   *
   * \note Assume that maxlen is longer than len because hi-res sampling should
   *       result in a longer arc length.
   */
  auto error_percent = [](double len, double maxlen) -> double
  {
      double errPct = 0.;
      if(maxlen > len)
      {
         // pct is in [0,1). Generally, as the curve gets better pct should approach 1.
         double pct = (len / maxlen);
         // If we have a good curve then errPct should be small.
         errPct = 1. - pct;
      }
      return errPct;
    };

  /*!
   * \brief Determines a u value within the interval [u0, u1] that should make
   *        the longest line segments from p0->curve(u)->p1.
   *
   * \param interpolator The NURBSInterpolator we can use to query the curve.
   * \param u0 The parametric value for the start of the interval.
   * \param u1 The parametric value for the end of the interval.
   * \param p0 The point at the start of the interval.
   * \param p1 The point at the end of the interval.
   *
   * \return The u value for the new point.
   */
  auto solveu = [](NURBSInterpolator &interpolator, double u0, double u1,
                   const PointType &p0, const PointType &p1) -> double
  {
    // Set up an initial ui at the midpoint and get the point.
    double ui = (u0 + u1) / 2.;
    double ustart = ui;

    // Interval endpoints.
    double u0x = p0[0];
    double u0y = p0[1];
    double u1x = p1[0];
    double u1y = p1[1];

    constexpr double EPS = 1.e-6;
    constexpr int MAX_ITERATIONS = 32;
    bool solved = false;

    // Solve for derivative == 0.
    for(int iteration = 0; iteration < MAX_ITERATIONS && !solved; iteration++)
    {
      // Get the point at ui
      PointType p = interpolator.at(ui);

      // Deltas from interval endpoints to current point at ui.
      double dx0 = p[0] - u0x;
      double dy0 = p[1] - u0y;
      double dx1 = p[0] - u1x;
      double dy1 = p[1] - u1y;

      // Get the 1st and 2nd derivatives of curve at ui
      PointType d[2];
      interpolator.derivativesAt(ui, 2, d);
      double xp = d[0][0];
      double yp = d[0][1];
      double xpp = d[1][0];
      double ypp = d[1][1];

      // 1st derivative of length L(u), the length of the new line segment pair.
      double D1L = (dx0 * xp + dy0 * yp) / sqrt(dx0 * dx0 + dy0 * dy0) +
                   (dx1 * xp + dy1 * yp) / sqrt(dx1 * dx1 + dy1 * dy1);

      if(axom::utilities::isNearlyEqual(D1L, 0., EPS))
      {
#ifdef AXOM_DEBUG_SOLVEU_VERBOSE
         std::cout << iteration << ": Stopping: Moved u from " << ustart << " to " << ui << std::endl;
#endif
        solved = true;
        break;
      }
      else
      {
        // 2nd derivative of length L(u)
        double D2L = ((-2. * pow(dx0*xp + dy0*yp, 2.) + 2 * (dx0*dx0 + dy0*dy0) * (xp*xp + yp*yp + dx0 * xpp + dy0 * ypp)) /
                      (2 * pow(dx0 * dx0 + dy0 * dy0, 3. / 2.)))
                     +
                     ((-2. * pow(dx1*xp + dy1*yp, 2.) + 2 * (dx1*dx1 + dy1*dy1) * (xp*xp + yp*yp + dx1 * xpp + dy1 * ypp)) /
                      (2 * pow(dx1 * dx1 + dy1 * dy1, 3. / 2.)));

        // Newton step
        ui = ui - D1L / D2L;

        // If we go out of bounds, use the ustart value.
        if(ui <= u0)
        {
#ifdef AXOM_DEBUG_SOLVEU_VERBOSE
          std::cout << iteration << ": ERROR new ui " << ui << " is less than u0 " << u0 << std::endl;
#endif
          ui = ustart;
          break;
        }
        if(ui >= u1)
        {
#ifdef AXOM_DEBUG_SOLVEU_VERBOSE
          std::cout << iteration << ": ERROR new ui " << ui << " is greater than u1 " << u1 << std::endl;
#endif
          ui = ustart;
          break;
        }
      }
    }

    return ui;
  };

  /*!
   * \brief Examines the history values to determine if the deltas between
   *        iterations are sufficiently small that the iterations should
   *        terminate, even though it might not have reached the desired
   *        target error.
   *
   * \param iteration The solve iteration.
   * \param history A circular buffer of the last several history values.
   * \param nhistory The number of history values.
   * \param percentError The percent error to achieve.
   *
   * \note This function assumes that history values increase.
   */
  auto diminishing_returns = [](int iteration,
                                const double *history,
                                int nhistory,
                                double percentError) -> bool
  {
    bool dr = false;
    // We have enough history to decide if there are diminishing returns.
    if(iteration >= nhistory)
    {
      // Compute a series of percents between pairs of history values.
      double sum = 0.;
      for(int i = 0; i < nhistory - 1; i++)
      {
        int cur = (iteration - i) % nhistory;
        int prev = (iteration - i - 1) % nhistory;
        double pct = 1. - history[prev] / history[cur];
        sum += fabs(pct);
      }
      double avg_pct = sum / static_cast<double>(nhistory - 1);
#if 1//def AXOM_DEBUG_LINEARIZE_VERBOSE
      std::cout << "Diminishing returns. iteration = " << iteration
                << ", avg_pct = " << avg_pct
                << ", percentError = " << percentError<< "." << std::endl;
#endif
      dr = avg_pct < percentError;
    }
    return dr;
  };

  constexpr size_t INITIAL_GUESS_NPTS = 100;
  const double EPS_SQ = m_vertexWeldThreshold * m_vertexWeldThreshold;

  // Initialize the revolved volume.
  revolvedVolume = 0.;

  // Clamp the lower error bound so it does not get impractically small.
  percentError = axom::utilities::clampLower(percentError, 1.e-10);

//#define AXOM_DEBUG_WRITE_ERROR_CURVES
#ifdef AXOM_DEBUG_WRITE_ERROR_CURVES
  FILE *ferr = fopen("error.curve", "wt");

  FILE *fhcl = fopen("hicurvelen.curve", "wt");
  fprintf(fhcl, "# hicurvelen\n");

  FILE *fcl = fopen("curvelen.curve", "wt");
  fprintf(fcl, "# curvelen\n");

  FILE *fthresh = fopen("threshold.curve", "wt");
  fprintf(fthresh, "# threshold\n");
#endif
  int contourCount = -1;

  // Iterate over the contours and linearize each of them.
  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs, m_vertexWeldThreshold);
    contourCount++;
#ifdef AXOM_DEBUG_WRITE_ERROR_CURVES
    fprintf(ferr, "# contour%d\n", contourCount);
#endif
    // Add the contour's revolved volume to the total.
    revolvedVolume += interpolator.revolvedVolume(transform);

    // Get the contour start/end parameters.
    const double u0 = interpolator.startParameter(0);
    const double u1 = interpolator.endParameter(interpolator.numSpans() - 1);

    // Compute the start/end points of the whole curve.
    const PointType p0 = interpolator.at(u0);
    const PointType p1 = interpolator.at(u1);

    // This segment represents the whole [0,1] interval.
    Segments::segment first;
    first.u = u0;
    first.point = p0;
    first.length = sqrt(primal::squared_distance(p0, p1));
    first.next_u = solveu(interpolator, u0, u1, p0, p1);
    PointType next_point = interpolator.at(first.next_u);
    first.next_length[0] = sqrt(primal::squared_distance(p0, next_point));
    first.next_length[1] = sqrt(primal::squared_distance(next_point, p1));

    // This segment just contains the end point of the interval.
    Segments::segment last;
    last.u = u1;
    last.point = p1;
    last.length = 0.;
    last.next_u = u1;
    last.next_length[0] = 0.;
    last.next_length[1] = 0.;

    // Store the initial segments.
    Segments S;
    S.reserve(INITIAL_GUESS_NPTS);
    S.push_back(first);
    S.push_back(last);

    // If there is a single span in the contour then we already added its points.
    if(interpolator.numSpans() > 1)
    {
      //---------------------------------------------------------------------
      // Approximate the arc length of the whole curve by using a lot of line segments.
      // Computers are fast, sampling is cheap.
      double hiCurveLen = 0.;
      constexpr int NUMBER_OF_SAMPLES = 2000;
      PointType prev = first.point;
      for(int i = 1; i < NUMBER_OF_SAMPLES; i++)
      {
        double u = i / static_cast<double>(NUMBER_OF_SAMPLES - 1);
        PointType cur = interpolator.at(u);
        hiCurveLen += sqrt(primal::squared_distance(prev, cur));
        prev = cur;
      }
      //---------------------------------------------------------------------

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
      // Print initial iteration.
      SLIC_INFO(fmt::format("hiCurveLen: {}", hiCurveLen);
      SLIC_INFO(fmt::format("curveLength: {}", curveLength);
      SLIC_INFO(fmt::format("percentError: {}", percentError);
#endif

      // The initial curve length.
      double curveLength = first.length;
      constexpr double ADAPTIVE_CUTOFF = 0.;//01;

      // Iterate until the difference between iterations is under the percent error
      // or there are diminishing returns to continuing.
      constexpr int MAX_HISTORY = 10;
      double history[MAX_HISTORY] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
      int iteration = 0;
      while((error_percent(curveLength, hiCurveLen) > percentError)// &&
//            !diminishing_returns(iteration - 1, history, MAX_HISTORY, percentError)
           )
      {
        // Get the index of the segment we'll split.
        size_t splitIndex = S.next_longest();
        size_t nextIndex = splitIndex + 1;

        // Partition segment.
        Segments::segment &left = S[splitIndex];
        Segments::segment right;

        // The old segment length pre-split.
        double oldSegLength = left.length;

        // The left and right segments will be the sum of the parent
        // segment's 2 pieces.
        double newSegLength = left.next_length[0] + left.next_length[1];

        // Split the left segment into left, right. Note that if the interval is
        // small enough, we stop using solveu and just use the midpoint to save time.
        right.u = left.next_u;
        right.point = interpolator.at(right.u);
        right.length = left.next_length[1];

// NOTE: the solveu must be making it worse...

//        if((S[nextIndex].u - right.u) > ADAPTIVE_CUTOFF)
//          right.next_u = solveu(interpolator, right.u, S[nextIndex].u, right.point, S[nextIndex].point);
//        else
          right.next_u = (right.u + S[nextIndex].u) / 2.;
        PointType right_next = interpolator.at(right.next_u);
        right.next_length[0] = sqrt(primal::squared_distance(right.point, right_next));
        right.next_length[1] = sqrt(primal::squared_distance(right_next, S[nextIndex].point));

        left.length = left.next_length[0];
//        if((right.u - left.u) > ADAPTIVE_CUTOFF)
//          left.next_u = solveu(interpolator, left.u, right.u, left.point, right.point);
//        else
          left.next_u = (left.u + right.u) / 2.;
        PointType left_next = interpolator.at(left.next_u);
        left.next_length[0] = sqrt(primal::squared_distance(left.point, left_next));
        left.next_length[1] = sqrt(primal::squared_distance(left_next, right.point));

        // Insert the right segment after the left segment.
        S.insert(nextIndex, right);

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
        // Print initial iteration.
        SLIC_INFO(fmt::format("Iteration: {}", iteration));
        SLIC_INFO(fmt::format("\thiCurveLen: {}", hiCurveLen));
        SLIC_INFO(fmt::format("\tcurveLength: {}", curveLength));
#endif

        // Update the curve length.
        curveLength = curveLength - oldSegLength + newSegLength;

        // Save in history. This is used in diminishing_returns.
        history[iteration % MAX_HISTORY] = curveLength;

#ifdef AXOM_DEBUG_WRITE_LINES
        // Make a filename.
        char filename[512];
        int npts = S.size();
        sprintf(filename, "lines%d_%05d.vtk", contourCount, npts);
        write_lines(filename, S);

        std::cout << "Wrote " << filename << ": hiCurveLen=" << std::setprecision(9) << hiCurveLen
                  << ", curveLength=" << std::setprecision(9) << curveLength
                  << ", dCL=" << (hiCurveLen - curveLength)
                  << std::endl;
#endif

        iteration++;
      } // while(keep_going())

      SLIC_INFO(
        fmt::format("getLinearMesh: "
                "percentError = {}"
                ", hiCurveLength = {}"
                ", curveLength = {}",
                percentError,
                hiCurveLen,
                curveLength));
    } // if(interpolator.numSpans() > 1)

    // Add the points to the mesh.
    appendPoints(mesh, S, EPS_SQ);
  }

#ifdef AXOM_DEBUG_WRITE_ERROR_CURVES
  fclose(ferr);
  fclose(fhcl);
  fclose(fcl);
  fclose(fthresh);
#endif
}

//---------------------------------------------------------------------------
// NOTE: This API change is temporary while I am pulling data out with the curve segments.
void C2CReader::getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                              int segmentsPerKnotSpan,
                              std::vector<double>& d1vec,
                              std::vector<double>& d2vec,
                              std::vector<double>& uvec,
                              std::vector<double>& curvvec,
                              std::vector<double>& sp)
{
  using axom::utilities::lerp;

  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT,
                "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(segmentsPerKnotSpan < 1,
                "C2C reader: Need at least one segment per NURBs span");

  using PointType = primal::Point<double, 2>;
  using PointsArray = std::vector<PointType>;

  const double EPS_SQ = m_vertexWeldThreshold * m_vertexWeldThreshold;

  // Makes a set of uniformly spaced points in the interval.
  auto makeUniformPoints = [](NURBSInterpolator &interpolator,
                              double startParameter,
                              double endParameter,
                              int segmentsPerKnotSpan,
                              std::vector<double> &uValues)//PointsArray &pts)
  {
    double denom = static_cast<double>(segmentsPerKnotSpan);
    for(int i = 0; i <= segmentsPerKnotSpan; ++i)
    {
      double u = lerp(startParameter, endParameter, i / denom);
      //pts.emplace_back(interpolator.at(u));
      uValues.push_back(u);
    }
  };

  // Examine the vector of curvature values for the intervals and see if they
  // seem flat.
  auto looksFlat = [](const std::vector<double> &curvatures) -> bool
  {
    constexpr double EPS = 1.e-4;
    bool flat = true;
    for(size_t j = 0; j < curvatures.size(); j++)
    {
      for(size_t i = 0; i < j; i++)
      {
        flat &= axom::utilities::isNearlyEqual(curvatures[j], curvatures[i], EPS);
      }
    }
    return flat;
  };

  // Determine the u value that produces targetCurv.
  auto solveForU = [](NURBSInterpolator &interpolator,
                      double startParameter,
                      double endParameter,
                      double startCurv,
                      double endCurv,
                      double targetCurv,
                      double curvTolerance) -> double
  {
    double left = startParameter;
    double right = endParameter;
    double u = left;
    while(left <= right)
    {
      double umid = (left + right) / 2;
      double curv_at_umid = interpolator.curvature(umid);
      if(axom::utilities::isNearlyEqual(curv_at_umid, targetCurv, curvTolerance))
      {
        u = umid;
        break;
      }
      else if(startCurv < endCurv)
      {
        if(targetCurv > curv_at_umid)
          left = umid;
        else
          right = umid;
      }
      else
      {
        if(targetCurv > curv_at_umid)
          right = umid;
        else
          left = umid;
      }
    }
    return u;
  };

  // Non-uniform sampling in [0,1].
  auto sample = [](double t) -> double
  {
// sigmoid - highlight the higher curvature areas.
//    double s1 = pow(tanh(4.) * tanh(3.), 0.4);
//    double s = pow(tanh(4 * t * t) * tanh(3 * t * t), 0.4) / s1;

// highlight middle curvature
//    double s = (4. * pow(t - 0.5, 3) + 0.5 + t) / 2;

// highlight lower curvature.
//    double s = pow((t*t + t) / 2., 4./3.);

// highlight lower curvature even more
    double s = (t*t*t + t/50.) / 1.02;
    return s;
  };

  auto area = [](const std::vector<PointType> &pts) -> double
  {
    auto npts = static_cast<int>(pts.size());
    // Compute the areas
    double A = 0.;
    for(int i = 0; i < npts; i++)
    {
        int nexti = (i + 1) % npts;
        A += pts[i][0] * pts[nexti][1] - pts[i][1] * pts[nexti][0];
    }
    return A * 0.5;
  };

  auto filter = [area](NURBSInterpolator &interpolator, std::vector<double> &u)
  {
    // Make the curve points.
    std::vector<PointType> pts(u.size());
    for(size_t i = 0; i < u.size(); i++)
      pts[i] = interpolator.at(u[i]);

    std::vector<PointType> tri(3);

    double m = 0.01;
    for(size_t i = 1; i < u.size() - 1; i++)
    {
      double u0 = u[i-1];
      double u1 = u[i+1];
      double du = (u1 - u0) * m;
      double ustart = u0 + du;
      double uend = u1 - du;

      tri[0] = pts[i+1];
      tri[1] = pts[i-1];
      double a = -1.e6;

      // Sweep the point through the interval and see if we can increase
      // the triangle area. If we can, we move the point.
      while(ustart < uend)
      {
        // Make a new point using ustart
        tri[2] = interpolator.at(ustart);

        // Make the area of this new triangle.
        double newa = fabs(area(tri));

        if(newa > a)
        {
          // The area increased. keep the point.
          a = newa;
          u[i] = ustart;
          pts[i] = tri[2];
        }

        ustart += du;
      }
    }
  };

  auto filter2 = [area](NURBSInterpolator &interpolator, std::vector<double> &u)
  {
    PointType d[2];
    for(size_t i = 1; i < u.size() - 1; i++)
    {
      double u0 = u[i-1];
      double u1 = u[i+1];
      double ui = (u0 + u1) / 2.;

      auto p0 = interpolator.at(u0);
      auto p1 = interpolator.at(u1);

      double u0x = p0[0];
      double u0y = p0[1];
      double u1x = p1[0];
      double u1y = p1[1];

      constexpr double EPS = 1.e-6;
      constexpr int MAX_ITERATIONS = 32;
      bool solved = false;

      // Solve for derivative == 0.
      for(int iteration = 0; iteration < MAX_ITERATIONS && !solved; iteration++)
      {
        // Get the 1st and 2nd derivatives of curve at ui
        interpolator.derivativesAt(ui, 2, d);
        double xp = d[0][0];
        double yp = d[0][1];
        double xpp = d[1][0];
        double ypp = d[1][1];

        // 1st derivative of area A(u)
        double D1A = 0.5 * ((u0x - u1x) * yp + (u1y - u0y) * xp);
        if(axom::utilities::isNearlyEqual(D1A, 0., EPS))
        {
std::cout << iteration << ": Stopping: Moved u["<<i<<"]=" << u[i] << " to " << ui << std::endl;

          u[i] = ui;

          solved = true;
          break;
        }
        else
        {
          // 2nd derivative of area A(u)
          double D2A = 0.5 * ((u0x - u1x) * ypp + (u1y - u0y) * xpp);

          ui = ui - D1A / D2A;

          if(ui <= u0)
          {
            std::cout << iteration << ": ERROR new ui " << ui << " is less than u0 " << u0 << std::endl;
            ui = u0;
            break;
          }
          if(ui >= u1)
          {
            std::cout << iteration << ": ERROR new ui " << ui << " is greater than u1 " << u1 << std::endl;
            ui = u1;
            break;
          }
        }
      }
    }
  };

  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs, m_vertexWeldThreshold);

    // For each knot span
    for(int span = 0; span < interpolator.numSpans(); ++span)
    {
std::cout << "span " << span << " --------------------------------------------------" << std::endl;
      // Generate points on the curve
      PointsArray pts;
      pts.reserve(segmentsPerKnotSpan + 1);

      const double startParameter = interpolator.startParameter(span);
      const double endParameter = interpolator.endParameter(span);

      std::vector<double> uValues;

      std::string method("UNIFORM");
      if(getenv("AXOM_METHOD") != nullptr)
        method = getenv("AXOM_METHOD");
      bool allowFilter = true;
      if(method == "UNIFORM")
      {
        makeUniformPoints(interpolator, startParameter, endParameter,
                          segmentsPerKnotSpan, uValues);
        allowFilter = false;
      }
      else if(method == "BISECT") 
      {
//------------------------------------------------------------------------------

        // Start out with the start/end u value for this span.
        uValues.reserve(segmentsPerKnotSpan + 1);
        uValues.push_back(startParameter);
        uValues.push_back(endParameter);

        int nSegments = 1;
        std::vector<PointType> tri(3);
        while(nSegments < segmentsPerKnotSpan)
        {
          double maxa = -std::numeric_limits<double>::max(), maxu = uValues[0];
          int maxseg = 0;
          // Look at each segment to see which we should bisect. We do this by
          // computing a triangle area and we select the one that makes the
          // largest area.
          for(int seg = 0; seg < nSegments; seg++)
          {
            double umid = (uValues[seg] + uValues[seg + 1]) / 2.;
            tri[0] = interpolator.at(uValues[seg]);
            tri[1] = interpolator.at(umid);
            tri[2] = interpolator.at(uValues[seg + 1]);
            double a = area(tri);

if(a < 0.)
{
//std::cout << "a=" << a << std::endl;
a *= -1.;
}
            if(a > maxa)
            {
              maxu = umid;
              maxa = a;
              maxseg = seg + 1;
            }
          }
          uValues.insert(uValues.begin() + maxseg, maxu);
          nSegments++;

          // For smaller numbers of segments, allow the points to move a bit more
          // so we converge on area slightly faster.
          if(nSegments < 10)
          {
            filter2(interpolator, uValues);
          }
        }
        allowFilter = false;
//------------------------------------------------------------------------------
      }
      else
      {
        // Get the intervals for this span.
        std::vector<double> intervals(interpolator.curvatureIntervals(startParameter, endParameter));
        std::cout << "intervals for (" << startParameter << ", " << endParameter << ") = {";
        for(auto value : intervals)
           std::cout << value << ", ";
        std::cout << "}" << std::endl;
        // Get the curvatures for the interval values.
        std::vector<double> curvatures(intervals.size());
        for(size_t i = 0; i < intervals.size(); i++)
          curvatures[i] = interpolator.curvature(intervals[i]);

        if(segmentsPerKnotSpan == 1 || looksFlat(curvatures))
        {
std::cout << "span " << span << ": looks flat. Make uniform points. " << std::endl;
          makeUniformPoints(interpolator, startParameter, endParameter,
                            segmentsPerKnotSpan, uValues);
        }
        else if(static_cast<int>(intervals.size()-1) == segmentsPerKnotSpan)
        {
std::cout << "span " << span << ": Number of intervals matches segmentsPerKnotSpan." << std::endl;
          // The intervals define the u values for the segment points.
          uValues = intervals;
/*          for(const double u : intervals)
          {
            pts.emplace_back(interpolator.at(u));
std::cout << "span " << span << ": Adding point at u=" << u << std::endl;
          }
*/
        }
        else if(intervals.size() == 2)
        {
std::cout << "span " << span << ": Single interval: segmentsPerKnotSpan=" << segmentsPerKnotSpan << std::endl;
          // We have 2 interval points so the curvature is monotonic.
          double curvStart = curvatures[0];
          double curvEnd = curvatures[1];
          double curvTolerance = fabs(curvEnd - curvStart) / 10000.;
          double denom = static_cast<double>(segmentsPerKnotSpan);
//          pts.emplace_back(interpolator.at(startParameter));
//std::cout << "span " << span << ": Adding point at u=" << startParameter << std::endl;
uValues.push_back(startParameter);
          for(int i = 1; i < segmentsPerKnotSpan; ++i)
          {
            // Divide curvature range uniformly.
            double t = i / denom;
            // Feed the uniform value through some other functions to highlight
            // higher curvature values. Then make the target curvature value
            double s, targetCurv;
            if(curvStart > curvEnd)
            {
              s = sample(1. - t);
              targetCurv = lerp(curvEnd, curvStart, s);
            }
            else
            {
              s = sample(t);
              targetCurv = lerp(curvStart, curvEnd, s);
            }
#if 1
            std::cout << std::setw(2)  << span << ", "
                      << std::setw(12) << intervals[0] << ", "
                      << std::setw(12) << intervals[1] << ", "
                      << std::setw(12) << curvStart << ", "
                      << std::setw(12) << curvEnd << ", "
                      << std::setw(12) << i << ", "
                      << std::setw(12) << t << ", "
                      << std::setw(12) << s << ", "
                      << std::setw(12) << targetCurv
                      << std::endl;
#endif
            // Solve for the u value that has curvature targetCurv.
            double u = solveForU(interpolator,
                                 intervals[0],
                                 intervals[1],
                                 curvStart,
                                 curvEnd,
                                 targetCurv,
                                 curvTolerance);

            // Make the new point at u.
//            pts.emplace_back(interpolator.at(u));
//std::cout << "span " << span << ": Adding point at u=" << u << std::endl;
            uValues.push_back(u);

          }
//          pts.emplace_back(interpolator.at(endParameter));
//std::cout << "span " << span << ": Adding point at u=" << endParameter << std::endl;
          uValues.push_back(endParameter);

        }
        else
        {
std::cout << "span " << span << ": Divide span into intervals." << std::endl;
          int numIntervals = static_cast<int>(intervals.size() - 1);

          // If we get here then these conditions are true:
          // 1. segmentsPerKnotSpan > 1
          // 2. numIntervals > 1
          // 3. segmentsPerKnotSpan != numIntervals

          if(segmentsPerKnotSpan < numIntervals)
          {
std::cout << "span " << span << ": segmentsPerKnotSpan < numIntervals: "
          << segmentsPerKnotSpan << " < " << numIntervals << std::endl;

            // We desire fewer segments than are prescribed for this span.
            // Combine the 2 lowest curvature intervals until we have the
            // desired number of segments.
            while(numIntervals > segmentsPerKnotSpan)
            {
              size_t removeIdx = 0;
              double lowestDC = 0.;
              for(int ii = 1; ii < numIntervals; ii++)
              {
                double dC = fabs(curvatures[ii + 1] - curvatures[ii - 1]);
                if(ii == 1 || dC < lowestDC)
                {
                  lowestDC = dC;
                  removeIdx = ii;
                }
              }

              // We found the index of the interval value we can remove (not
              // endpoints) that we can remove. Doing this combines the intervals.
std::cout << "span " << span << ": Removing index " << removeIdx << std::endl;
              intervals.erase(intervals.begin() + removeIdx);
              curvatures.erase(curvatures.begin() + removeIdx);
              numIntervals = static_cast<int>(intervals.size() - 1);
            }

            // Now we can add all the interval points.
//            for(const double u : intervals)
//            {
//              pts.emplace_back(interpolator.at(u));
//std::cout << "span " << span << ": Adding point at u=" << u << std::endl;
//            }
            uValues = intervals;
          }
          else
          {
std::cout << "span " << span << ": More segments than intervals" << std::endl;

            // numIntervals < segmentsPerKnotSpan. This means we can give
            // each interval 1 segment and then apportion the remaining
            // segments to intervals that have higher curvature differences.
            double totalCurvature = 0.;
            std::vector<double> curvatureRange(numIntervals);
            for(size_t ii = 1; ii < intervals.size(); ii++)
            {
              double dC = fabs(curvatures[ii] - curvatures[ii - 1]);
              curvatureRange[ii - 1] = dC;
              totalCurvature += dC;
            }
std::cout << "span " << span << ": totalCurvature=" << totalCurvature << std::endl;

            // Start each interval out with 1 segment.
            std::vector<int> segments(intervals.size() - 1, 1);
            // Iterate over the ranges and determine the fair share of line segments
            // for each interval based on total curvature difference.
            int availableSegments = segmentsPerKnotSpan - (numIntervals);
            int totalSegments = 0;
            for(size_t ii = 0; ii < curvatureRange.size(); ii++)
            {
              double segFraction = segmentsPerKnotSpan * (curvatureRange[ii] / totalCurvature);
              double dseg = static_cast<double>(segments[ii]);
              if(segFraction > dseg)
              {
                auto additionalSegments = static_cast<int>(std::round(segFraction - dseg));
                int allowableSegments = std::min(additionalSegments, availableSegments);
                segments[ii] += allowableSegments;
                availableSegments -= allowableSegments;
              }
              totalSegments += segments[ii];
            }

            if(totalSegments != segmentsPerKnotSpan)
            {
              std::cout << "The number of segments " << totalSegments
                        << " != " << segmentsPerKnotSpan << std::endl;
            }
            for(const auto &seg : segments)
              std::cout << "\t" << seg << " segments." << std::endl;

            // Iterate over the intervals.
            for(size_t ii = 1; ii < intervals.size(); ii++)
            {
              double curvStart = curvatures[ii - 1];
              double curvEnd = curvatures[ii];
              double curvTolerance = fabs(curvEnd - curvStart) / 10000.;

              int iend = (ii == intervals.size()-1) ? segments[ii-1] : (segments[ii-1] - 1);
std::cout << "span " << span << ": ii=" << ii << ", curvStart=" << curvStart << ", curvEnd=" << curvEnd << ", curvTolerance=" << curvTolerance << ", iend=" << iend << std::endl;

              for(int i = 0; i <= iend; ++i)
              {
                // Divide curvature range uniformly.
                double t = i / static_cast<double>(segments[ii-1]);
std::cout << "\ti=" << i << ", t=" << t << std::endl;

                // Feed the uniform value through some other functions to highlight
                // higher curvature values. Then make the target curvature value
                double s, targetCurv;
                if(curvStart > curvEnd)
                {
                  s = sample(1. - t);
                  targetCurv = lerp(curvEnd, curvStart, s);
                }
                else
                {
                  s = sample(t);
                  targetCurv = lerp(curvStart, curvEnd, s);
                }
std::cout << "\ts=" << s << ", targetCurv=" << targetCurv << std::endl;

                // Solve for the u value that has curvature targetCurv.
                double u = solveForU(interpolator,
                                     intervals[ii - 1],
                                     intervals[ii],
                                     curvStart,
                                     curvEnd,
                                     targetCurv,
                                     curvTolerance);
#if 1
                std::cout << std::setw(2)  << span << ", "
                          << std::setw(12) << intervals[ii - 1] << ", "
                          << std::setw(12) << intervals[ii] << ", "
                          << std::setw(12) << curvStart << ", "
                          << std::setw(12) << curvEnd << ", "
                          << std::setw(12) << u << ", "
                          << std::setw(12) << targetCurv
                          << std::endl;
#endif
                // Make the new point at u.
//                pts.emplace_back(interpolator.at(u));
//std::cout << "span " << span << ": Adding point at u=" << u << std::endl;

                uValues.push_back(u);

              }
            } // for intervals
          }
        }
      }

#if 1
      // Filter the points to see if it improves their distribution.
      if(allowFilter)
        filter(interpolator, uValues);
#endif

      // Make the points.
      std::cout << "span " << span << ": u={";
      for(const double u : uValues)
      {
         pts.emplace_back(interpolator.at(u));
         std::cout << u << ", ";
#if 1
        PointType dpts[2];
        interpolator.derivativesAt(u, 2, dpts);

        d1vec.push_back(dpts[0].data()[0]);
        d1vec.push_back(dpts[0].data()[1]);

        d2vec.push_back(dpts[1].data()[0]);
        d2vec.push_back(dpts[1].data()[1]);

        uvec.push_back(u);

        curvvec.push_back(interpolator.curvature(u));

        sp.push_back(span);
#endif
      }
      std::cout << "}" << std::endl;

      // Check for simple vertex welding opportunities at endpoints of newly interpolated points
      {
        int numNodes = mesh->getNumberOfNodes();
        if(numNodes > 0)  // this is not the first Piece
        {
          PointType meshPt;
          // Fix start point if necessary; check against most recently added vertex in mesh
          mesh->getNode(numNodes - 1, meshPt.data());
          if(primal::squared_distance(pts[0], meshPt) < EPS_SQ)
          {
            pts[0] = meshPt;
          }

          // Fix end point if necessary; check against 0th vertex in mesh
          const int endIdx = pts.size() - 1;
          mesh->getNode(0, meshPt.data());
          if(primal::squared_distance(pts[endIdx], meshPt) < EPS_SQ)
          {
            pts[endIdx] = meshPt;
          }
        }
        else  // This is the first, and possibly only span, check its endpoint, fix if necessary
        {
          int endIdx = pts.size() - 1;
          if(primal::squared_distance(pts[0], pts[endIdx]) < EPS_SQ)
          {
            pts[endIdx] = pts[0];
          }
        }
      }

      // Add the new points and segments to the mesh, respecting welding checks from previous block
      {
        const int startNode = mesh->getNumberOfNodes();
        const int numNewNodes = pts.size();
        mesh->reserveNodes(startNode + numNewNodes);

        for(int i = 0; i < numNewNodes; ++i)
        {
          mesh->appendNode(pts[i][0], pts[i][1]);
        }

        const int startCell = mesh->getNumberOfCells();
        const int numNewSegments = pts.size() - 1;
        mesh->reserveCells(startCell + numNewSegments);
        for(int i = 0; i < numNewSegments; ++i)
        {
          IndexType seg[2] = {startNode + i, startNode + i + 1};
          mesh->appendCell(seg, mint::SEGMENT);
        }
      }

    }  // end for each knot span
  }    // end for each NURBS curve
}

}  // end namespace quest
}  // end namespace axom
