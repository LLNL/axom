// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
#include <map>
#include <ostream>

namespace axom
{
namespace primal
{
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

/// Utility class that caches precomputed coefficient matrices for sectorMoment computation
template <typename T>
class MemoizedSectorMomentWeights
{
public:
  using SectorWeights = numerics::Matrix<T>;

  MemoizedSectorMomentWeights() = default;

  ~MemoizedSectorMomentWeights()
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
      auto vec = generateBezierCurveSectorMomentsWeights(order);
      for(int d = 0; d <= order; ++d)
      {
        m_sectorWeightsMap[std::make_pair(order, d)] = vec[d];
      }
    }

    return *(m_sectorWeightsMap[std::make_pair(order, dim)]);
  }

  /*!
   * \brief Computes the weights for BezierCurve's sectorMoment() function
   *
   * \param order The polynomial order of the curve
   * \return An anti-symmetric matrix with (order+1)*{order+1) entries
   * containing the integration weights for entry (i,j)
   *
   * The derivation is provided in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
  std::vector<SectorWeights*> generateBezierCurveSectorMomentsWeights(int ord) const
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
 */

template <typename T, int NDIMS>
class BezierCurve
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using NumArrayType = NumericArray<T, NDIMS>;
  using SegmentType = Segment<T, NDIMS>;
  using CoordsVec = std::vector<PointType>;
  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;

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
    AXOM_STATIC_ASSERT_MSG(
      (NDIMS == 2) || (NDIMS == 3),
      "A Bezier Curve object may be defined in 2-D or 3-D");
    AXOM_STATIC_ASSERT_MSG(
      std::is_arithmetic<T>::value,
      "A Bezier Curve must be defined using an arithmetic type");
    SLIC_ASSERT(ord >= -1);
    const int sz = utilities::max(-1, ord + 1);
    m_controlPoints.resize(sz);
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
    AXOM_STATIC_ASSERT_MSG(
      (NDIMS == 2) || (NDIMS == 3),
      "A Bezier Curve object may be defined in 2-D or 3-D");
    AXOM_STATIC_ASSERT_MSG(
      std::is_arithmetic<T>::value,
      "A Bezier Curve must be defined using an arithmetic type");
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
    AXOM_STATIC_ASSERT_MSG(
      (NDIMS == 2) || (NDIMS == 3),
      "A Bezier Curve object may be defined in 2-D or 3-D");
    AXOM_STATIC_ASSERT_MSG(
      std::is_arithmetic<T>::value,
      "A Bezier Curve must be defined using an arithmetic type");
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(ord >= 0);

    const int sz = utilities::max(0, ord + 1);
    m_controlPoints.resize(sz);

    for(int p = 0; p <= ord; ++p)
    {
      m_controlPoints[p] = pts[p];
    }
  }

  /*! Sets the order of the Bezier Curve*/
  void setOrder(int ord) { m_controlPoints.resize(ord + 1); }

  /*! Returns the order of the Bezier Curve*/
  int getOrder() const { return static_cast<int>(m_controlPoints.size()) - 1; }

  /*! Clears the list of control points*/
  void clear()
  {
    const int ord = getOrder();
    for(int p = 0; p <= ord; ++p)
    {
      m_controlPoints[p] = PointType();
    }
  }

  /*! Retrieves the control point at index \a idx */
  PointType& operator[](int idx) { return m_controlPoints[idx]; }

  /*! Retrieves the control point at index \a idx */
  const PointType& operator[](int idx) const { return m_controlPoints[idx]; }

  /* Checks equality of two Bezier Curve */
  friend inline bool operator==(const BezierCurve<T, NDIMS>& lhs,
                                const BezierCurve<T, NDIMS>& rhs)
  {
    return lhs.m_controlPoints == rhs.m_controlPoints;
  }

  friend inline bool operator!=(const BezierCurve<T, NDIMS>& lhs,
                                const BezierCurve<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /*! Returns a copy of the Bezier curve's control points */
  CoordsVec getControlPoints() const { return m_controlPoints; }

  /*! Reverses the order of the Bezier curve's control points */
  void reverseOrientation()
  {
    const int ord = getOrder();
    CoordsVec old_controlPoints = m_controlPoints;
    for(int i = 0; i <= ord; ++i)
    {
      m_controlPoints[i] = old_controlPoints[ord - i];
    }
  }

  /*! Returns an axis-aligned bounding box containing the Bezier curve */
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /*! Returns an oriented bounding box containing the Bezier curve */
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
    PointType ptval;

    const int ord = getOrder();
    std::vector<T> dCarray(ord + 1);

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
          dCarray[k] = (1 - t) * dCarray[k] + t * dCarray[k + 1];
        }
      }
      ptval[i] = dCarray[0];
    }

    return ptval;
  }

  /*!
   * \brief Computes the tangent of  a Bezier curve at a particular parameter value \a t
   *
   * \param [in] t parameter value at which to compute tangent 
   * \return p the tangent vector of the Bezier curve at t
   *
   * \note We typically find the tangent of the curve at \a t between 0 and 1
   */

  PointType dt(T t) const
  {
    PointType ptval;

    const int ord = getOrder();
    std::vector<T> dCarray(ord + 1);

    // Run de Casteljau algorithm on each dimension
    for(int i = 0; i < NDIMS; ++i)
    {
      for(int p = 0; p <= ord; ++p)
      {
        dCarray[p] = m_controlPoints[p][i];
      }

      for(int p = 1; p <= ord - 1; ++p)
      {
        const int end = ord - p;
        for(int k = 0; k <= end; ++k)
        {
          dCarray[k] = (1 - t) * dCarray[k] + t * dCarray[k + 1];
        }
      }
      ptval[i] = ord * (dCarray[1] - dCarray[0]);
    }

    return ptval;
  }

  /*!
   * \brief Splits a Bezier curve into two Bezier curves at particular parameter
   * value between 0 and 1
   *
   * \param [in] t parameter value between 0 and 1 at which to evaluate
   * \param [out] c1, c2 Bezier curves that split the original
   */
  void split(T t, BezierCurve& c1, BezierCurve& c2) const
  {
    int ord = getOrder();
    SLIC_ASSERT(ord >= 0);

    // Note: the second curve's control points are computed inline
    //       as we find the first curve's control points
    c2 = *this;

    c1.setOrder(ord);
    c1[0] = c2[0];

    // Run de Casteljau algorithm
    // After each iteration, save the first control point into c1
    for(int p = 1; p <= ord; ++p)
    {
      const int end = ord - p;
      for(int k = 0; k <= end; ++k)
      {
        PointType& pt1 = c2[k];
        const PointType& pt2 = c2[k + 1];
        for(int i = 0; i < NDIMS; ++i)
        {
          pt1[i] = (1 - t) * pt1[i] + t * pt2[i];
        }
      }
      c1[p] = c2[0];
    }

    return;
  }

  /*!
   * \brief Calculates the sector moment of a Bezier Curve
   *
   * The sector moment is the moment between the curve and the origin.
   * The equation and derivation are generalizations of:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
  PointType sectorMoment() const
  {
    // Weights for each polynomial order's moments are precomputed and memoized
    static internal::MemoizedSectorMomentWeights<T> s_weights;

    T Mx = 0;
    T My = 0;
    const int ord = getOrder();
    for(int r = 0; r <= ord; ++r)
    {
      const auto& weights_r = s_weights.getWeights(ord, r);
      for(int p = 0; p <= ord; ++p)
      {
        for(int q = 0; q <= ord; ++q)
        {
          Mx += weights_r(p, q) * m_controlPoints[p][1] *
            m_controlPoints[q][0] * m_controlPoints[r][0];
          My += weights_r(p, q) * m_controlPoints[p][1] *
            m_controlPoints[q][0] * m_controlPoints[r][1];
        }
      }
    }
    PointType M = PointType::make_point(Mx, My);
    return M;
  }

  /*!
   * \brief Calculates the sector area of a Bezier Curve
   *
   * The sector area is the area between the curve and the origin.
   * The equation and derivation is described in:
   *  Ueda, K. "Signed area of sectors between spline curves and the origin"
   *  IEEE International Conference on Information Visualization, 1999.
   */
  T sectorArea() const
  {
    // Weights for each polynomial order are precomputed and memoized
    static internal::MemoizedSectorAreaWeights<T> s_weights;

    T A = 0;
    const int ord = getOrder();
    const auto& weights = s_weights.getWeights(ord);

    for(int p = 0; p <= ord; ++p)
    {
      for(int q = 0; q <= ord; ++q)
      {
        A += weights(p, q) * m_controlPoints[p][1] * m_controlPoints[q][0];
      }
    }
    return A;
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
    for(int p = 1; p < ord && sqDist < tol; ++p)  // check interior control points
    {
      sqDist += squared_distance(m_controlPoints[p], seg);
    }
    return (sqDist < tol);
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
    os << "}";

    return os;
  }

private:
  CoordsVec m_controlPoints;
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
