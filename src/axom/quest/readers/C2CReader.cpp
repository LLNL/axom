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

/// Overload to format a c2c::Point using fmt
template <>
struct axom::fmt::formatter<c2c::Point> : ostream_formatter
{ };

namespace axom
{
namespace quest
{
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
  BasisVector calculateBasisFunctions(int span, double u)
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
  PointType at(double u)
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
  void derivativeBasisFunctions(int span, double u, int n,
                                std::vector<BasisVector> &ders) const
  {
    const int p = m_curve.order - 1;
    const auto& U = m_curve.knots;

    std::vector<BasisVector> ndu(p + 1), a(2);
    BasisVector left(p + 1), right(p + 1);
    for(int j = 0; j <= p; j++)
      ndu[j].resize(p + 1);
    for(int j = 0; j <= n; j++)
      ders[j].resize(p + 1);
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
      ders[0][j] = ndu[j][p];

    // This section computes the derivatives (Eq. [2.9])

    // Loop over function index.
    for(int r = 0; r <= p; r++)
    {
      int s1 = 0, s2 = 1; // Alternate rows in array a
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
    int r = p;
    for(int k = 1; k <= n; k++)
    {
      for(int j = 0; j <= p; j++)
        ders[k][j] *= r;
      r *= (p - k);
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
  void derivativesAt(double u, int d, PointType *CK) const
  {
    const int p = m_curve.order - 1;
    int du = std::min(d, p);

    const auto span = findSpan(u);
    std::vector<BasisVector> N(du + 1);
    derivativeBasisFunctions(span, u, d, N);

    for(int k = 1; k <= du; k++)
    {
      double x = 0., y = 0.;
      for(int j = 0; j <= p; j++)
      {
        int offset = span - p + j;
        // TODO: We likely need to include the weights and then compensate.
        x = x + N[k][j] * m_curve.controlPoints[offset].getZ().getValue();
        y = y + N[k][j] * m_curve.controlPoints[offset].getR().getValue();
      }
      CK[k-1] = PointType{x,y};
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
    const PointType &D1 = derivs[0];
    const PointType &D2 = derivs[1];

    // This is signed curvature as formulated at:
    // https://en.wikipedia.org/wiki/Curvature#Curvature_of_a_graph
    // k = (x'y'' - x''y') / pow(x'x' + y'y', 3./2.)
    double numerator = (D1.data()[0] * D2.data()[1]) -
                       (D2.data()[0] * D1.data()[1]);
    double D1mag2 = (D1.data()[0] * D1.data()[0]) +
                    (D1.data()[1] * D1.data()[1]);
    double one_over_denominator = pow(D1mag2, -3. / 2.);
    return numerator * one_over_denominator;
  }

private:
  const c2c::NURBSData& m_curve;
  std::vector<std::pair<double, double>> m_spanIntervals;
};

// NOTE: We would eventually like to be able to pass an error term to the
//       c2c reader that lets it figure out how many segements it needs to make
//       to get a linearized curve that is precise enough (when integrated) that
//       analytic_solution - this_linearization < error_tolerance.

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

// NOTE: This API change is temporary while I am pulling data out with the curve segments.
void C2CReader::getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                              int segmentsPerKnotSpan,
                              std::vector<double> &d1vec,
                              std::vector<double> &d2vec,
                              std::vector<double> &uvec,
                              std::vector<double> &curvvec,
                              std::vector<double> &sp)
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
#if 1
  FILE *d1, *d2, *cf;
  d1 = fopen("d1.curve", "wt");
  fprintf(d1, "# d1\n");
  d2 = fopen("d2.curve", "wt");
  fprintf(d2, "# d2\n");
  cf = fopen("cf.curve", "wt");
  fprintf(cf, "# curvature\n");
#endif
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

#if 0
      double denom = static_cast<double>(segmentsPerKnotSpan);
      for(int i = 0; i <= segmentsPerKnotSpan; ++i)
      {
        double u = lerp(startParameter, endParameter, i / denom);
        pts.emplace_back(interpolator.at(u));
      }
#else
      std::vector<double> uvalues;
// NOTE: This is temporary. It's nice to be able to revert to the old algorithm
//       when making data for comparisons.
if(getenv("AXOM_UNIFORM") != nullptr)
{
      double denom = static_cast<double>(segmentsPerKnotSpan);
      for(int i = 0; i <= segmentsPerKnotSpan; ++i)
      {
        double u = lerp(startParameter, endParameter, i / denom);
        pts.emplace_back(interpolator.at(u));

        uvalues.push_back(u);
      }
}
else
{
      double denom = static_cast<double>(segmentsPerKnotSpan);
      double curvStart = interpolator.curvature(startParameter);
      double curvEnd = interpolator.curvature(endParameter);
      constexpr double EPS = 1.e-6;
      if(axom::utilities::isNearlyEqual(curvStart, curvEnd, EPS))
      {
        // There is essentially no curvature in the span. Break up the 
        // span uniformly.
        for(int i = 0; i <= segmentsPerKnotSpan; ++i)
        {
          double u = lerp(startParameter, endParameter, i / denom);
          pts.emplace_back(interpolator.at(u));
          uvalues.push_back(u);
        }
      }
      else
      {
        // There is a range of curvature. We break up that range of curvature
        // in the span uniformly and then figure out which u values produce
        // those curvature values. We assume for the time being that within
        // a span that there are no curvature values outside the range determined
        // by the span endpoints.
        pts.emplace_back(interpolator.at(startParameter));
        uvalues.push_back(startParameter);
std::cout << "span,curveStart,curveEnd,targetCurv,i,t,s,u" << std::endl;

        double CURV_EPS = fabs(curvEnd - curvStart) / 1000.;
        for(int i = 1; i < segmentsPerKnotSpan; ++i)
        {
          // Divide curvature range uniformly.
          double t = i / denom;
          // Feed the uniform value through some other functions to highlight
          // higher curvature values. Then make the target curvature value
          double s, targetCurv;
          if(curvStart > curvEnd)
          {
              t = 1. - t;
              s = pow(tanh(4 * t * t) * tanh(3 * t * t), 0.4);
              targetCurv = lerp(curvEnd, curvStart, s);
          }
          else
          {
              s = pow(tanh(4 * t * t) * tanh(3 * t * t), 0.4);
              targetCurv = lerp(curvStart, curvEnd, s);
          }

          // Determine the u value that produces targetCurv.
          double left = startParameter;
          double right = endParameter;
          double u = 0.;
          while(left <= right)
          {
            double umid = (left + right) / 2;
            double curv_at_umid = interpolator.curvature(umid);
            if(axom::utilities::isNearlyEqual(curv_at_umid, targetCurv, CURV_EPS))
            {
              u = umid;
              break;
            }
            else if(curvStart < curvEnd)
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
std::cout << span << ", "
          << curvStart << ", "
          << curvEnd << ", "
          << targetCurv << ", "
          << i << ", "
          << t << ", "
          << s << ", "
          << u << std::endl;

          pts.emplace_back(interpolator.at(u));
          uvalues.push_back(u);
        }
        pts.emplace_back(interpolator.at(endParameter));
        uvalues.push_back(endParameter);
      }
}
      // Now that we know all the u values for this span, compute some quantities of interest.
      for(auto u : uvalues)
      {
        PointType dpts[2];
        interpolator.derivativesAt(u, 2, dpts);
        fprintf(d1, "%lg %lg\n", dpts[0].data()[0], dpts[0].data()[1]);
        fprintf(d2, "%lg %lg\n", dpts[1].data()[0], dpts[1].data()[1]);
        fprintf(cf, "%lg %lg\n", u, interpolator.curvature(u));

        d1vec.push_back(dpts[0].data()[0]);
        d1vec.push_back(dpts[0].data()[1]);

        d2vec.push_back(dpts[1].data()[0]);
        d2vec.push_back(dpts[1].data()[1]);

        uvec.push_back(u);

        curvvec.push_back(interpolator.curvature(u));

        sp.push_back(span);
      }
#endif

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
#if 1
  fclose(d1);
  fclose(d2);
  fclose(cf);
#endif
}

}  // end namespace quest
}  // end namespace axom
