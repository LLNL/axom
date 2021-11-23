// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_
#define AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_

/*!
 * \file compute_moments.hpp
 *
 * \brief Consists of a set of methods to compute areas/volumes and centroids 
 * for Polygons and CurvedPolygons composed of BezierCurves
 */

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/geometry/CurvedPolygon.hpp"

#include <vector>
#include <map>

namespace axom
{
namespace primal
{
namespace internal
{
template <typename T>
class MemoizedSectorAreaWeights;
template <typename T>
class MemoizedSectorCentroidWeights;
}  // namespace internal

/*!
   * \brief Calculates the sector area of a planar Bezier Curve
   *
   * The sector area is the area between the curve and the origin.
   * The equation and derivation is described in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
template <typename T>
T sector_area(const primal::BezierCurve<T, 2>& curve)
{
  // Weights for each polynomial order are precomputed and memoized
  static internal::MemoizedSectorAreaWeights<T> s_weights;

  T A = 0;
  const int ord = curve.getOrder();
  const auto& weights = s_weights.getWeights(ord);

  for(int p = 0; p <= ord; ++p)
  {
    for(int q = 0; q <= ord; ++q)
    {
      A += weights(p, q) * curve[p][1] * curve[q][0];
    }
  }
  return A;
}

/*!
   * \brief Calculates the sector centroid of a planar Bezier Curve
   *
   * This is the centroid of the region between the curve and the origin.
   * The equation and derivation are generalizations of:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
template <typename T>
primal::Point<T, 2> sector_centroid(const primal::BezierCurve<T, 2>& curve)
{
  // Weights for each polynomial order's centroid are precomputed and memoized
  static internal::MemoizedSectorCentroidWeights<T> s_weights;

  T Mx = 0;
  T My = 0;
  const int ord = curve.getOrder();
  for(int r = 0; r <= ord; ++r)
  {
    const auto& weights_r = s_weights.getWeights(ord, r);
    for(int p = 0; p <= ord; ++p)
    {
      for(int q = 0; q <= ord; ++q)
      {
        Mx += weights_r(p, q) * curve[p][1] * curve[q][0] * curve[r][0];
        My += weights_r(p, q) * curve[p][1] * curve[q][0] * curve[r][1];
      }
    }
  }
  return primal::Point<T, 2> {Mx, My};
}

/// \brief Returns the area enclosed by the CurvedPolygon
template <typename T>
T area(const primal::CurvedPolygon<T, 2>& poly, double tol = 1e-8)
{
  const int ngon = poly.numEdges();
  T A = 0.0;
  if(!poly.isClosed(1e3 * tol))
  {
    SLIC_DEBUG(
      "Warning! The area is 0 because the curved polygon is not closed.");
    return A;
  }
  else
  {
    for(int ed = 0; ed < ngon; ++ed)
    {
      A += primal::sector_area(poly[ed]);
    }
    return A;
  }
}

/// \brief Returns the centroid of the CurvedPolygon
template <typename T>
primal::Point<T, 2> centroid(const primal::CurvedPolygon<T, 2>& poly,
                             double tol = 1e-8)
{
  using PointType = primal::Point<T, 2>;

  const int ngon = poly.numEdges();
  PointType M;

  if(!poly.isClosed(1e3 * tol))
  {
    SLIC_DEBUG(
      "Warning! The moments are 0 because the curved polygon is not closed.");
    return M;
  }
  else
  {
    const T A = area(poly, tol);
    if(A != 0.)
    {
      for(int ed = 0; ed < ngon; ++ed)
      {
        M.array() += primal::sector_centroid(poly[ed]).array();
      }
      M.array() /= A;
    }
    return M;
  }
}

namespace internal
{
/// Utility class that caches precomputed coefficient matrices for sectorArea computation
template <typename T>
class MemoizedSectorAreaWeights
{
public:
  using SectorWeights = numerics::Matrix<T>;

  MemoizedSectorAreaWeights() = default;

  ~MemoizedSectorAreaWeights()
  {
    for(auto& p : m_sectorWeightsMap)
    {
      delete[] p.second->data();  // delete the matrix's data
      delete p.second;            // delete the matrix
      p.second = nullptr;
    }
    m_sectorWeightsMap.clear();
  }

  /// Returns a memoized matrix of coeficients for sector area computation
  const SectorWeights& getWeights(int order) const
  {
    // Compute and cache the weights if they are not already available
    if(m_sectorWeightsMap.find(order) == m_sectorWeightsMap.end())
    {
      SectorWeights* weights = generateBezierCurveSectorWeights(order);
      m_sectorWeightsMap[order] = weights;
    }

    return *(m_sectorWeightsMap[order]);
  }

private:
  /*!
   * \brief Computes the weights for BezierCurve's sectorArea() function
   *
   * \param order The polynomial order of the curve
   * \return An anti-symmetric matrix with (order+1)*{order+1) entries
   * containing the integration weights for entry (i,j)
   *
   * The derivation is provided in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
  SectorWeights* generateBezierCurveSectorWeights(int ord) const
  {
    const bool memoryIsExternal = true;
    const int SZ = ord + 1;
    SectorWeights* weights =
      new SectorWeights(SZ, SZ, new T[SZ * SZ], memoryIsExternal);

    T binom_2n_n = static_cast<T>(utilities::binomialCoefficient(2 * ord, ord));
    for(int i = 0; i <= ord; ++i)
    {
      (*weights)(i, i) = 0.;  // zero on the diagonal
      for(int j = i + 1; j <= ord; ++j)
      {
        double val = 0.;
        if(i != j)
        {
          T binom_ij_i = static_cast<T>(utilities::binomialCoefficient(i + j, i));
          T binom_2nij_nj = static_cast<T>(
            utilities::binomialCoefficient(2 * ord - i - j, ord - j));

          val = ((j - i) * ord) / binom_2n_n *
            (binom_ij_i / static_cast<T>(i + j)) *
            (binom_2nij_nj / (2. * ord - j - i));
        }
        (*weights)(i, j) = val;  // antisymmetric
        (*weights)(j, i) = -val;
      }
    }
    return weights;
  }

private:
  mutable std::map<int, SectorWeights*> m_sectorWeightsMap;
};

/// Utility class that caches precomputed coefficient matrices for sector_centroid() computation
template <typename T>
class MemoizedSectorCentroidWeights
{
public:
  using SectorWeights = numerics::Matrix<T>;

  MemoizedSectorCentroidWeights() = default;

  ~MemoizedSectorCentroidWeights()
  {
    for(auto& p : m_sectorWeightsMap)  // for each matrix of weights
    {
      delete[] p.second->data();  // delete the matrix's data
      delete p.second;            // delete the matrix
      p.second = nullptr;
    }
    m_sectorWeightsMap.clear();
  }

  /// Returns a memoized matrix of sector moment coeficients for component \a dim or order \a order
  const SectorWeights& getWeights(int order, int dim) const
  {
    // Compute and cache the weights if they are not already available
    if(m_sectorWeightsMap.find(std::make_pair(order, dim)) ==
       m_sectorWeightsMap.end())
    {
      auto vec = generateBezierCurveSectorCentroidWeights(order);
      for(int d = 0; d <= order; ++d)
      {
        m_sectorWeightsMap[std::make_pair(order, d)] = vec[d];
      }
    }

    return *(m_sectorWeightsMap[std::make_pair(order, dim)]);
  }

  /*!
   * \brief Computes the weights for BezierCurve's sectorCentroid() function
   *
   * \param order The polynomial order of the curve
   * \return An anti-symmetric matrix with (order+1)*{order+1) entries
   * containing the integration weights for entry (i,j)
   *
   * The derivation is provided in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
  std::vector<SectorWeights*> generateBezierCurveSectorCentroidWeights(int ord) const
  {
    const bool memoryIsExternal = true;
    const int SZ = ord + 1;

    std::vector<SectorWeights*> weights;
    weights.resize(SZ);
    for(int k = 0; k <= ord; ++k)
    {
      SectorWeights* weights_k =
        new SectorWeights(SZ, SZ, new T[SZ * SZ], memoryIsExternal);
      for(int i = 0; i <= ord; ++i)
      {
        (*weights_k)(i, i) = 0.;  // zero on the diagonal
        for(int j = i + 1; j <= ord; ++j)
        {
          double val = 0.;
          if(i != j)
          {
            T binom_n_i = static_cast<T>(utilities::binomialCoefficient(ord, i));
            T binom_n_j = static_cast<T>(utilities::binomialCoefficient(ord, j));
            T binom_n_k = static_cast<T>(utilities::binomialCoefficient(ord, k));
            T binom_3n2_ijk1 = static_cast<T>(
              utilities::binomialCoefficient(3 * ord - 2, i + j + k - 1));

            val = (1. * (j - i)) / (3. * (3 * ord - 1)) *
              (1. * binom_n_i * binom_n_j * binom_n_k / (1. * binom_3n2_ijk1));
          }
          (*weights_k)(i, j) = val;  // antisymmetric
          (*weights_k)(j, i) = -val;
        }
      }
      weights[k] = weights_k;
    }
    return weights;
  }

private:
  mutable std::map<std::pair<int, int>, SectorWeights*> m_sectorWeightsMap;
};

}  // namespace internal

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_COMPUTE_MOMENTS_HPP_
