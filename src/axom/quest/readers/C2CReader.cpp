// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
#include "fmt/fmt.hpp"

namespace axom
{
namespace quest
{
/// returns the linear interpolation of \a A and \a B at \a t. i.e. (1-t)A+tB
inline double lerp(double A, double B, double t) { return (1 - t) * A + t * B; }

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

  NURBSInterpolator(const c2c::NURBSData& curve) : m_curve(curve)
  {
    int p = m_curve.order - 1;
    int knotSize = m_curve.knots.size();
    SLIC_ASSERT(p >= 1);
    SLIC_ASSERT(knotSize >= 2 * p);
  }

  /*!
   * \brief Checks that the knots of the NURBS curve are closed
   *
   * The knots vector is closed when it begins with \a m_curve.order equal knot values
   * and ends with the same number of equal knot values
   */
  bool isClosed(double EPS = 1E-8) const
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

  double startParameter() const { return m_curve.knots[0]; }
  double endParameter() const
  {
    return m_curve.knots[m_curve.knots.size() - 1];
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

      cw[0] += N[j] * weight * controlPoint.getR().getValue();
      cw[1] += N[j] * weight * controlPoint.getZ().getValue();
      cw[2] += N[j] * weight;
    }

    // Return projected point
    // All units should have been normalized by c2c::toNurbs(piece, units)
    return PointType {cw[0] / cw[2], cw[1] / cw[2]};
  }

private:
  const c2c::NURBSData& m_curve;
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
    fmt::format("Loading contour with {} pieces", contour.getNumEntries()));

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
    sstr << fmt::format("\tweights: {}\n", fmt::join(nurbs.weights, " "));
    sstr << fmt::format("\tcontrol points: {}\n",
                        fmt::join(nurbs.controlPoints, " "));
    sstr << "}\n";
    ++index;
  }

  SLIC_INFO(sstr.str());
}

void C2CReader::getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                              int segmentsPerPiece)
{
  // Sanity checks
  SLIC_ERROR_IF(mesh == nullptr, "supplied mesh is null!");
  SLIC_ERROR_IF(mesh->getDimension() != 2, "C2C reader expects a 2D mesh!");
  SLIC_ERROR_IF(mesh->getCellType() != mint::SEGMENT,
                "C2C reader expects a segment mesh!");
  SLIC_ERROR_IF(segmentsPerPiece < 1,
                "C2C reader: Need at least one sample per NURBs span");

  using PointType = primal::Point<double, 2>;
  using PointsArray = std::vector<PointType>;

  for(const auto& nurbs : m_nurbsData)
  {
    NURBSInterpolator interpolator(nurbs);

    PointsArray pts;

    // Generate points on the curve
    {
      pts.reserve(segmentsPerPiece + 1);

      const double startParameter = interpolator.startParameter();
      const double endParameter = interpolator.endParameter();

      double denom = static_cast<double>(segmentsPerPiece);
      for(int i = 0; i <= segmentsPerPiece; ++i)
      {
        double u = lerp(startParameter, endParameter, i / denom);
        pts.emplace_back(interpolator.at(u));
      }
    }

    // Add the new points and segments to the mesh
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
}

}  // end namespace quest
}  // end namespace axom
