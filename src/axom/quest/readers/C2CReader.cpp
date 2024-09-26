// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
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
/*!
 * \brief Transform a 2D point and return a 2D point.
 *
 * \param transform A 4x4 transformation matrix.
 * \param pt The input 2D point.
 *
 * \return The transformed 2D point.
 */
inline primal::Point<double, 2> transformPoint(
  const numerics::Matrix<double>& transform,
  const primal::Point<double, 2>& pt)
{
  // Turn the point2 into a vec4.
  double pt4[] = {pt[0], pt[1], 0., 1.}, newpt4[] = {0., 0., 0., 1.};

  // Transform the point.
  axom::numerics::matrix_vector_multiply(transform, pt4, newpt4);

  // Make a new 2D point from the 2D components.
  return primal::Point<double, 2> {newpt4[0], newpt4[1]};
}

//---------------------------------------------------------------------------
/*!
 * \brief This data type stores a segment, which is used in refining curves.
 */
struct Segment
{
  using PointType = primal::Point<double, 2>;

  double u;         // The u parameter for the start of the segment.
  PointType point;  // The point at parameter u.
  double length;    // Distance from point to next segment's point.

  double next_u;          // Best place to split this segment.
  double next_length[2];  // left/right length values at next_u.
};

/*!
 * \brief Returns the index of the longest segment in the segments.
 * \param segments A vector of Segments.
 * \return The index of the element with the longest next_length value.
 */
size_t nextLongest(const std::vector<Segment>& segments)
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

//---------------------------------------------------------------------------
/*!
 * \brief Append a set of segments to a mint mesh, filtering out like points.
 * \param mesh The mint mesh to which we're adding points/cells.
 * \param S The segments to write.
 * \param EPS_SQ The squared distance for point matching.
 */
static void appendPoints(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                         std::vector<Segment>& S,
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
#ifdef AXOM_DEBUG_WRITE_LINES
/*!
 * \brief Write line segments to a VTK file for visualization.
 * \param filename The name of the file to write.
 * \param S The segments to write.
 */
static void writeLines(const std::string& filename, const std::vector<Segment>& S)
{
  FILE* f = fopen(filename.c_str(), "wt");
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
    if((i + 1) < npts)
      fprintf(f, "%1.16lf %1.16lf 0. ", S[i + 1].point[0], S[i + 1].point[1]);
    if((i + 2) < npts)
      fprintf(f, "%1.16lf %1.16lf 0. ", S[i + 2].point[0], S[i + 2].point[1]);
    fprintf(f, "\n");
  }
  fprintf(f, "\n");

  int nspans = npts - 1;

  // Write ncells
  fprintf(f, "LINES %d %d\n", nspans, 3 * nspans);
  for(int ispan = 0; ispan < nspans; ispan++)
  {
    fprintf(f, "2 %d %d\n", ispan, ispan + 1);
  }

  fclose(f);
}
#endif

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
      {
        ders[k][j] *= r;
      }
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
   * Rational derivatives from Algorithm A4.2 on p. 127 of "The NURBS Book".
   */
  void derivativesAt(double u, int d, PointType* CK) const
  {
    const int p = m_curve.order - 1;

    // Compute the basis functions for the curve and its d derivatives.
    int du = std::min(d, p);
    const auto span = findSpan(u);
    std::vector<BasisVector> N(d + 1);
    derivativeBasisFunctions(span, u, d, N);

    // Store w(u) in wders[0], w'(u) in wders[1], ...
    std::vector<PointType> Aders(d + 1);
    std::vector<double> wders(d + 1);

    // Compute the point and d derivatives
    for(int k = 0; k <= du; k++)
    {
      double x = 0., y = 0., w = 0.;
      for(int j = 0; j <= p; j++)
      {
        int offset = span - p + j;
        const double weight = m_curve.weights[offset];

        // Compute the weighted point.
        double Pw[3];
        Pw[0] = weight * m_curve.controlPoints[offset].getZ().getValue();
        Pw[1] = weight * m_curve.controlPoints[offset].getR().getValue();
        Pw[2] = weight;

        x = x + N[k][j] * Pw[0];
        y = y + N[k][j] * Pw[1];
        w = w + N[k][j] * Pw[2];
      }
      Aders[k] = PointType {x, y};
      wders[k] = w;
    }

    // Do the rational part from A4.2. Note that Aders[0] is the point A(u)
    // and wders[0] is w(u).
    for(int k = 0; k <= d; k++)
    {
      PointType v = Aders[k];
      for(int i = 1; i <= k; i++)
      {
        auto bin = axom::utilities::binomialCoefficient(k, i);
        v[0] = v[0] - bin * wders[i] * CK[k - 1][0];
        v[1] = v[1] - bin * wders[i] * CK[k - 1][1];
      }
      CK[k][0] = v[0] / wders[0];
      CK[k][1] = v[1] / wders[0];
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
    PointType derivs[3];
    derivativesAt(u, 2, derivs);
    const PointType& D1 = derivs[1];
    const PointType& D2 = derivs[2];

    double xp = D1.data()[0];   // x'
    double xpp = D2.data()[0];  // x''

    double yp = D1.data()[1];   // y'
    double ypp = D2.data()[1];  // y''

    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - y'x'') / pow(x'x' + y'y', 3./2.)
    double xp2_plus_yp2 = xp * xp + yp * yp;
    return (xp * ypp - yp * xpp) / pow(xp2_plus_yp2, 3. / 2.);
  }

  /*!
   * \brief Evaluates the B-spline curvature derivatives (up to 2nd) at
   *        parameter value \a u
   * 
   * \param[in] u The parameter value.
   * \param[in] d The number of derivatives to compute (1=1st deriv, 2=1st & 2nd derivs)
   * \param[out] ders An array that will contain the curvature derivatives.
   */
  void curvatureDerivatives(double u, int d, double* ders) const
  {
    // Evaluate 1st, 2nd, 3rd curve derivatives at u.
    PointType derivs[4];
    derivativesAt(u, 3, derivs);
    const PointType& D1 = derivs[1];
    const PointType& D2 = derivs[2];
    const PointType& D3 = derivs[3];

    double xp = D1.data()[0];    // x'
    double xpp = D2.data()[0];   // x''
    double xppp = D3.data()[0];  // x'''

    double yp = D1.data()[1];    // y'
    double ypp = D2.data()[1];   // y''
    double yppp = D3.data()[1];  // y'''

    // 1st derivative of curvature.
    double xp2_plus_yp2 = xp * xp + yp * yp;
    double A = -3. * (xp * ypp - yp * xpp) * 2 * (xp * xpp + yp * ypp);
    double B = 2. * pow(xp2_plus_yp2, 5. / 2.);
    double C = xp * yppp - yp * xppp;
    double D = pow(xp2_plus_yp2, 3. / 2.);
    ders[0] = A / B + C / D;

    if(d >= 2)
    {
      // 2nd derivative of curvature.
      double E = 15. * (-yp * xpp + xp * ypp) *
        pow(2. * xp * xpp + 2. * yp * ypp, 2.) /
        (4. * pow(xp2_plus_yp2, 7. / 2.));
      double F = 3. * (2. * xp * xpp + 2. * yp * ypp) *
        (-yp * xppp + xp * yppp) / pow(xp2_plus_yp2, 5. / 2.);
      double G = 3. * (-yp * xpp + xp * ypp) *
        (2. * (xpp * xpp) + 2. * (ypp * ypp) + 2. * xp * xppp + 2. * yp * yppp) /
        (2. * pow(xp2_plus_yp2, 5. / 2.));
      double H = (-ypp * xppp + xpp * yppp) / pow(xp2_plus_yp2, 3. / 2.);

      ders[1] = E - F - G + H;
    }
  }

  /*!
   * \brief Compute the revolved volume of the curve across its entire
   *        parametric interval [0,1] using quadrature.
   *
   * \param transform A 4x4 transformation matrix that we use to transform
   *                  the points and the derivative values.
   *
   * \note We include a transformation of the points so any transforms
   *       can figure into the integration.
   *
   * \return The revolved curve volume.
   */
  double revolvedVolume(const numerics::Matrix<double>& transform) const
  {
    // Use 5-point Gauss Quadrature.
    const double X[] = {-0.906179845938664,
                        -0.538469310105683,
                        0,
                        0.538469310105683,
                        0.906179845938664};
    const double W[] = {0.23692688505618908,
                        0.47862867049936647,
                        0.5688888888888889,
                        0.47862867049936647,
                        0.23692688505618908};

    // Make a transform with no translation. We use this to transform
    // the derivative since we want to permit scaling and rotation but
    // translating it does not make sense.
    numerics::Matrix<double> transform2 = transform;
    transform2(0, 3) = 0.;
    transform2(1, 3) = 0.;
    transform2(2, 3) = 0.;

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
    SLIC_INFO(fmt::format("revolvedVolume"));
#endif

    // Break up the [0,1] interval and compute quadrature in the subintervals.
    double vol = 0.;
    for(const auto& interval : m_spanIntervals)
    {
      double ad = interval.first;
      double bd = interval.second;
#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
      SLIC_INFO(fmt::format("interval ({}, {})", ad, bd));
#endif

      // Approximate the integral "Int pi*x'(u)*y(u)^2du" using quadrature.
      double scale = M_PI * ((bd - ad) / 2.);
      double sum = 0.;
      for(size_t i = 0; i < 5; i++)
      {
        // Map quad point x value [-1,1] to [a,b].
        double u = X[i] * ((bd - ad) / 2.) + ((bd + ad) / 2.);

        // Compute y(u) to get radius
        PointType p_u = at(u);
        PointType p_uT = transformPoint(transform, p_u);
        double r = p_uT[1];

        // Compute x'(u)
        PointType xprime[2];
        derivativesAt(u, 1, xprime);
        PointType xprimeT = transformPoint(transform2, xprime[1]);
        double xp = xprimeT[0];

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
        SLIC_INFO(fmt::format("\ti={}, u={}, p(u)=({},{}), xp=({},{})",
                              i,
                              u,
                              p_u[0],
                              p_u[1],
                              xprime[0],
                              xprime[1]));
#endif

        // Accumulate weight times dx*r^2.
        sum += W[i] * xp * (r * r);
      }
      // Guard against volumes being negative (if the curve went the wrong way)
      vol += fabs(scale * sum);
    }
#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
    SLIC_INFO(fmt::format("revolvedVolume={}", vol));
#endif
    return vol;
  }

  /*!
   * \brief Approximate the arc length of the whole curve by using
   *        the segments between a number of points.
   *
   * \param nSamples The number of points to use when making segments.
   *
   * \note Quadrature could probably be used too.
   *
   * \return The arc length of the whole curve when broken into segments.
   */
  double getArcLength(int nSamples) const
  {
    constexpr double u0 = 0.;
    double arcLength = 0.;
    PointType prev = at(u0);
    for(int i = 1; i < nSamples; i++)
    {
      double u = i / static_cast<double>(nSamples - 1);
      PointType cur = at(u);
      arcLength += sqrt(primal::squared_distance(prev, cur));
      prev = cur;
    }
    return arcLength;
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

//---------------------------------------------------------------------------
void C2CReader::getLinearMeshUniform(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                                     int segmentsPerKnotSpan)
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

  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs, m_vertexWeldThreshold);

    // For each knot span
    for(int span = 0; span < interpolator.numSpans(); ++span)
    {
      // Generate points on the curve
      PointsArray pts;
      pts.reserve(segmentsPerKnotSpan + 1);

      const double startParameter = interpolator.startParameter(span);
      const double endParameter = interpolator.endParameter(span);

      double denom = static_cast<double>(segmentsPerKnotSpan);
      for(int i = 0; i <= segmentsPerKnotSpan; ++i)
      {
        double u = lerp(startParameter, endParameter, i / denom);
        pts.emplace_back(interpolator.at(u));
      }

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

//---------------------------------------------------------------------------
void C2CReader::getLinearMeshNonUniform(
  mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
  double percentError)
{
  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT,
                "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(
    percentError <= 0.,
    axom::fmt::format(
      "C2C reader: percentError must be greater than 0. {} supplied.",
      percentError));
  SLIC_ERROR_IF(
    percentError >= 100.,
    axom::fmt::format(
      "C2C reader: percentError must be less than 100. {} supplied.",
      percentError));

  using PointType = primal::Point<double, 2>;

  /*!
   * \brief Checks curve lengths against tolerances and determines whether more
   *        refinement is needed.
   *
   * \param len The curve length that has been calculated.
   * \param maxlen The upper bound curve length that was calculated with a lot
   *               of segments.
   *
   * \note Assume that maxlen is longer than len because hi-res sampling should
   *       result in a longer arc length.
   */
  auto error_percent = [](double len, double maxlen) -> double {
    double errPct = 0.;
    if(maxlen > len)
    {
      errPct = 100. * (1. - (len / maxlen));
    }
    return errPct;
  };

  constexpr size_t INITIAL_GUESS_NPTS = 100;
  const double EPS_SQ = m_vertexWeldThreshold * m_vertexWeldThreshold;

  // Clamp the lower error bound so it does not get impractically small.
  percentError = axom::utilities::clampLower(percentError, 1.e-10);

#ifdef AXOM_DEBUG_WRITE_ERROR_CURVES
  // Make some curves for debugging.
  FILE* ferr = fopen("error.curve", "wt");

  FILE* fhcl = fopen("hicurvelen.curve", "wt");
  fprintf(fhcl, "# hicurvelen\n");

  FILE* fcl = fopen("curvelen.curve", "wt");
  fprintf(fcl, "# curvelen\n");

  FILE* fthresh = fopen("threshold.curve", "wt");
  fprintf(fthresh, "# threshold\n");
#endif

  int contourCount {-1};

  // clang complains about contourCount (-Wunused-but-set-variable); since we want keep it, let's mark it
  AXOM_UNUSED_VAR(contourCount);

  // Iterate over the contours and linearize each of them.
  std::vector<Segment> S;
  S.reserve(INITIAL_GUESS_NPTS);
  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs, m_vertexWeldThreshold);
    ++contourCount;
#ifdef AXOM_DEBUG_WRITE_ERROR_CURVES
    fprintf(ferr, "# contour%d\n", contourCount);
#endif

    // Get the contour start/end parameters.
    const double u0 = interpolator.startParameter(0);
    const double u1 = interpolator.endParameter(interpolator.numSpans() - 1);

    // Compute the start/end points of the whole curve.
    const PointType p0 = interpolator.at(u0);
    const PointType p1 = interpolator.at(u1);

    // This segment represents the whole [0,1] interval.
    Segment first;
    first.u = u0;
    first.point = p0;
    first.length = sqrt(primal::squared_distance(p0, p1));
    first.next_u = (u0 + u1) / 2.;
    PointType next_point = interpolator.at(first.next_u);
    first.next_length[0] = sqrt(primal::squared_distance(p0, next_point));
    first.next_length[1] = sqrt(primal::squared_distance(next_point, p1));

    // This segment just contains the end point of the interval.
    Segment last;
    last.u = u1;
    last.point = p1;
    last.length = 0.;
    last.next_u = u1;
    last.next_length[0] = 0.;
    last.next_length[1] = 0.;

    // Store the initial segments.
    S.clear();
    S.push_back(first);
    S.push_back(last);

    // If there is a single span in the contour then we already added its points.
    if(interpolator.numSpans() > 1)
    {
      // Compute the arc length of the curve.
      constexpr int MAX_NUMBER_OF_SAMPLES = 2000;
      double hiCurveLen = interpolator.getArcLength(MAX_NUMBER_OF_SAMPLES);

      // The initial curve length.
      double curveLength = first.length;

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
      // Print initial iteration.
      SLIC_INFO(fmt::format("hiCurveLen: {}", hiCurveLen));
      SLIC_INFO(fmt::format("curveLength: {}", curveLength));
      SLIC_INFO(fmt::format("percentError: {}", percentError));
#endif

      // Iterate until the difference between iterations is under the percent
      // error or we've reached the max number of samples. Each iteration, we
      // split the longest line segment so the line overall should get longer
      // until it reaches/approaches the upfront arc length. We get to put the
      // new points where they matter most early on so we should ideally get to
      // an acceptable arc length before we reach MAX_NUMBER_OF_SAMPLES.
      int iteration = 0;
      while((error_percent(curveLength, hiCurveLen) > percentError) &&
            (iteration < MAX_NUMBER_OF_SAMPLES))
      {
        // Get the index of the segment we'll split.
        size_t splitIndex = nextLongest(S);
        size_t nextIndex = splitIndex + 1;

        // Partition segment.
        Segment& left = S[splitIndex];
        Segment right;

        // The old segment length pre-split.
        double oldSegLength = left.length;

        // The left and right segments will be the sum of the parent
        // segment's 2 pieces.
        double newSegLength = left.next_length[0] + left.next_length[1];

        // Split the left segment into left, right. We set next_u as the midpoint
        // of the subintervals. It was being set via a Newton solveu() function
        // that yielded the u value with the longest left+right line segments but
        // using the midpoint ended up with better behavior.
        right.u = left.next_u;
        right.point = interpolator.at(right.u);
        right.length = left.next_length[1];

        right.next_u = (right.u + S[nextIndex].u) / 2.;
        PointType right_next = interpolator.at(right.next_u);
        right.next_length[0] =
          sqrt(primal::squared_distance(right.point, right_next));
        right.next_length[1] =
          sqrt(primal::squared_distance(right_next, S[nextIndex].point));

        left.length = left.next_length[0];
        left.next_u = (left.u + right.u) / 2.;
        PointType left_next = interpolator.at(left.next_u);
        left.next_length[0] =
          sqrt(primal::squared_distance(left.point, left_next));
        left.next_length[1] =
          sqrt(primal::squared_distance(left_next, right.point));

        // Insert the right segment after the left segment.
        S.insert(S.begin() + nextIndex, right);

#ifdef AXOM_DEBUG_LINEARIZE_VERBOSE
        // Print initial iteration.
        SLIC_INFO(fmt::format("Iteration: {}", iteration));
        SLIC_INFO(fmt::format("\thiCurveLen: {}", hiCurveLen));
        SLIC_INFO(fmt::format("\tcurveLength: {}", curveLength));
#endif

        // Update the curve length.
        curveLength = curveLength - oldSegLength + newSegLength;

#ifdef AXOM_DEBUG_WRITE_LINES
        // Make a filename.
        char filename[512];
        int npts = S.size();
        sprintf(filename, "lines%d_%05d.vtk", contourCount, npts);
        writeLines(filename, S);
        SLIC_INFO(fmt::format("Wrote {}", filename));
#endif

        iteration++;
      }  // while(error_percent())

      SLIC_INFO(
        fmt::format("getLinearMesh: "
                    "percentError = {}"
                    ", hiCurveLength = {}"
                    ", curveLength = {}",
                    percentError,
                    hiCurveLen,
                    curveLength));
    }  // if(interpolator.numSpans() > 1)

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

double C2CReader::getRevolvedVolume(const numerics::Matrix<double>& transform) const
{
  double revolvedVolume = 0.;
  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs, m_vertexWeldThreshold);

    // Add the contour's revolved volume to the total.
    revolvedVolume += interpolator.revolvedVolume(transform);
  }
  return revolvedVolume;
}

}  // end namespace quest
}  // end namespace axom
