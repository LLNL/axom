// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSPatch.hpp
 *
 * \brief A (trimmed) NURBSPatch primitive
 */

#ifndef AXOM_PRIMAL_NURBSPATCH_HPP_
#define AXOM_PRIMAL_NURBSPATCH_HPP_

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include "axom/primal/geometry/NumericArray.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/NURBSCurve.hpp"
#include "axom/primal/geometry/BezierPatch.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"

#include "axom/primal/operators/squared_distance.hpp"

#include <ostream>

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T, int NDIMS>
class NURBSPatch;

/*! \brief Overloaded output operator for NURBS Patches*/
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSPatch<T, NDIMS>& nPatch);

/*!
 * \class NURBSPatch
 *
 * \brief Represents a 3D NURBS patch defined by a 2D array of control points
 * \tparam T the coordinate type, e.g., double, float, etc.
 *
 * A NURBS patch has degrees `p` and `q` with knot vectors of length 
 * `r+1` and `s+1` respectively. There is a control net of (n + 1) * (m + 1) points
 * with r+1 = n+p+2 and s+1 = m+q+2. Optionally has weights for rational patches.
 * The curve must be c0 continuous and is parametrized from u=0 to u=1 and v=0 to v=1.
 * 
 * Nonrational NURBS patches are identified by an empty weights array.
 */
template <typename T, int NDIMS>
class NURBSPatch
{
public:
  using PointType = Point<T, NDIMS>;
  using VectorType = Vector<T, NDIMS>;
  using PlaneType = Plane<T, NDIMS>;

  using CoordsVec = axom::Array<PointType, 1>;
  using CoordsMat = axom::Array<PointType, 2>;
  using WeightsVec = axom::Array<T, 1>;
  using WeightsMat = axom::Array<T, 2>;
  using KnotVectorType = KnotVector<T>;

  using BoundingBoxType = BoundingBox<T, NDIMS>;
  using OrientedBoundingBoxType = OrientedBoundingBox<T, NDIMS>;
  using NURBSCurveType = primal::NURBSCurve<T, NDIMS>;

  AXOM_STATIC_ASSERT_MSG(
    (NDIMS == 1) || (NDIMS == 2) || (NDIMS == 3),
    "A NURBS Patch object may be defined in 1-, 2-, or 3-D");

  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A NURBS Patch must be defined using an arithmetic type");

public:
  /*!
   * \brief Constructor for a simple NURBS surface that reserves space for
   *  the minimum (sensible) number of points for the given degrees
   *
   * Constructs an empty patch by default (no nodes/weights on either axis)
   * 
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre deg_u, deg_v greater than or equal to -1.
   */
  NURBSPatch(int deg_u = -1, int deg_v = -1)
  {
    SLIC_ASSERT(deg_u >= -1 && deg_v >= -1);

    m_controlPoints.resize(deg_u + 1, deg_v + 1);
    m_knotvec_u = KnotVectorType(deg_u + 1, deg_u);
    m_knotvec_v = KnotVectorType(deg_v + 1, deg_v);

    makeNonrational();
  }

  /*!
   * \brief Constructor for a NURBS surface from a Bezier surface
   *
   * \param [in] bezierPatch the Bezier patch to convert to a NURBS patch 
   */
  explicit NURBSPatch(BezierPatch<T, NDIMS> bezierPatch)
  {
    m_controlPoints = bezierPatch.getControlPoints();
    m_weights = bezierPatch.getWeights();

    int deg_u = bezierPatch.getOrder_u();
    int deg_v = bezierPatch.getOrder_v();

    m_knotvec_u = KnotVectorType(deg_u + 1, deg_u);
    m_knotvec_v = KnotVectorType(deg_v + 1, deg_v);
  }

  /*!
   * \brief Constructor for a NURBS Patch from an array of coordinates and degrees
   *
   * \param [in] pts A 1D C-style array of npts_u*npts_v control points
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(PointType* pts, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(npts_u >= deg_u + 1 && npts_v >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    makeNonrational();

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);
  }

  /*!
   * \brief Constructor for a NURBS Patch from arrays of coordinates and weights
   *
   * \param [in] pts A 1D C-style array of (ord_u+1)*(ord_v+1) control points
   * \param [in] weights A 1D C-style array of (ord_u+1)*(ord_v+1) positive weights
   * \param [in] npts_u The number of control points on the first axis
   * \param [in] npts_v The number of control points on the second axis
   * \param [in] deg_u The patch's degree on the first axis
   * \param [in] deg_v The patch's degree on the second axis
   * \pre Requires that npts_d >= deg_d + 1 and deg_d >= 0 for d = u, v
   *
   * The knot vectors are constructed such that the patch is uniform
   * 
   * Elements of pts[k] are mapped to control nodes (p, q) lexicographically, i.e.
   * pts[k] = nodes[ k // (npts_u + 1), k % npts_v ]
   */
  NURBSPatch(PointType* pts, T* weights, int npts_u, int npts_v, int deg_u, int deg_v)
  {
    SLIC_ASSERT(pts != nullptr);
    SLIC_ASSERT(weights != nullptr);
    SLIC_ASSERT(npts_u >= deg_u + 1 && npts_v >= deg_v + 1);
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    m_controlPoints.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_controlPoints.flatIndex(t) = pts[t];
    }

    m_weights.resize(npts_u, npts_v);
    for(int t = 0; t < npts_u * npts_v; ++t)
    {
      m_weights.flatIndex(t) = weights[t];
    }

    SLIC_ASSERT(isValidRational());

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);
  }

  /*!
   * \brief Reset the knot vector in u
   *
   * \param [in] deg The target degree
   *
   * \warning This method does NOT change the existing control points,
   *  i.e. is not performing degree elevation or reduction.
   * \pre Requires deg_u < npts_u
   */
  void setDegree_u(int deg)
  {
    SLIC_ASSERT(deg >= 0);

    int npts_u = m_controlPoints.shape()[0];
    SLIC_ASSERT(deg < npts_u);

    if(npts_u < deg + 1)
    {
      npts_u = deg + 1;
    }

    m_controlPoints.resize(npts_u, getNumControlPoints_v());
    if(isRational())
    {
      m_weights.resize(npts_u, getNumControlPoints_v());
    }

    m_knotvec_u = KnotVectorType(npts_u, deg);
  }

  /*!
   * \brief Reset the knot vector in v
   *
   * \param [in] deg The target degree
   *
   * \warning This method does NOT change the existing control points,
   *  i.e. is not performing degree elevation or reduction.
   * \pre Requires deg_v < npts_v
   */
  void setDegree_v(int deg)
  {
    SLIC_ASSERT(deg >= 0);

    int npts_v = m_controlPoints.shape()[1];
    SLIC_ASSERT(deg < npts_v);

    if(npts_v < deg + 1)
    {
      npts_v = deg + 1;
    }

    m_controlPoints.resize(getNumControlPoints_u(), npts_v);
    if(isRational())
    {
      m_weights.resize(getNumControlPoints_u(), npts_v);
    }

    m_knotvec_v = KnotVectorType(npts_v, deg);
  }

  /*!
   * \brief Reset the knot vector and increase the number of control points
   *
   * \param [in] deg_u The target degree in u
   * \param [in] deg_v The target degree in v
   *
   * \warning This method does NOT change the existing control points,
   *  i.e. is not performing degree elevation or reduction.
   * \pre Requires deg_u < npts_u and deg_v < npts_v
   */
  void setDegree(int deg_u, int deg_v)
  {
    SLIC_ASSERT(deg_u >= 0 && deg_v >= 0);

    int npts_u = m_controlPoints.shape()[0];
    SLIC_ASSERT(deg_u < npts_u);

    int npts_v = m_controlPoints.shape()[1];
    SLIC_ASSERT(deg_v < npts_v);

    if(npts_u < deg_u + 1)
    {
      npts_u = deg_u + 1;
    }

    if(npts_v < deg_v + 1)
    {
      npts_v = deg_v + 1;
    }

    m_controlPoints.resize(npts_u, npts_v);
    if(isRational())
    {
      m_weights.resize(npts_u, npts_v);
    }

    m_knotvec_u = KnotVectorType(npts_u, deg_u);
    m_knotvec_v = KnotVectorType(npts_v, deg_v);
  }

  /*!
   * \brief Set the number control points in u
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints(int npts_u, int npts_v)
  {
    SLIC_ASSERT(npts_u > getDegree_u());
    SLIC_ASSERT(npts_v > getDegree_v());

    m_controlPoints.resize(npts_u, npts_v);
    m_knotvec_u.makeUniform(npts_u, getDegree_u());
    m_knotvec_v.makeUniform(npts_v, getDegree_v());
  }

  /*!
   * \brief Set the number control points in u
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints_u(int npts)
  {
    SLIC_ASSERT(npts > getDegree_u());

    m_controlPoints.resize(npts, getNumControlPoints_v());
    m_knotvec_u.makeUniform(npts, getDegree_u());
  }

  /*!
   * \brief Set the number control points in v
   *
   * \param [in] npts The target number of control points
   * 
   * \warning This method does NOT maintain the patch shape,
   *  i.e. is not performing knot insertion/removal.
   */
  void setNumControlPoints_v(int npts)
  {
    SLIC_ASSERT(npts > getDegree_v());

    m_controlPoints.resize(getNumControlPoints_u(), npts);
    m_knotvec_v.makeUniform(npts, getDegree_v());
  }

  /// \brief Returns the degree of the NURBS Patch on the first axis
  int getDegree_u() const { return m_knotvec_u.getDegree(); }

  /// \brief Returns the degree of the NURBS Patch on the second axis
  int getDegree_v() const { return m_knotvec_v.getDegree(); }

  /// \brief Returns the number of control points in the NURBS Patch on the first axis
  int getNumControlPoints_u() const
  {
    return static_cast<int>(m_controlPoints.shape()[0]);
  }

  /// \brief Returns the number of control points in the NURBS Patch on the second axis
  int getNumControlPoints_v() const
  {
    return static_cast<int>(m_controlPoints.shape()[1]);
  }

  /// \brief Return the length of the knot vector on the first axis
  int getNumKnots_u() const { return m_knotvec_u.getNumKnots(); }

  /// \brief Return the length of the knot vector on the second axis
  int getNumKnots_v() const { return m_knotvec_v.getNumKnots(); }

  /// Make trivially rational. If already rational, do nothing
  void makeRational()
  {
    if(!isRational())
    {
      auto patch_shape = m_controlPoints.shape();
      m_weights.resize(patch_shape[0], patch_shape[1]);
      m_weights.fill(1.0);
    }
  }

  /// Make nonrational by shrinking array of weights
  void makeNonrational() { m_weights.clear(); }

  /// Use array size as flag for rationality
  bool isRational() const { return !m_weights.empty(); }

  /// Clears the list of control points, make nonrational
  void clear()
  {
    m_controlPoints.clear();
    m_knotvec_u.clear();
    m_knotvec_v.clear();
    makeNonrational();
  }

  /// Retrieves the control point at index \a (idx_p, idx_q)
  PointType& operator()(int ui, int vi) { return m_controlPoints(ui, vi); }

  /// Retrieves the vector of control points at index \a idx
  const PointType& operator()(int ui, int vi) const
  {
    return m_controlPoints(ui, vi);
  }

  /*!
   * \brief Get a specific weight
   *
   * \param [in] ui The index of the weight on the first axis
   * \param [in] vi The index of the weight on the second axis
   * \pre Requires that the surface be rational
   */
  const T& getWeight(int ui, int vi) const
  {
    SLIC_ASSERT(isRational());
    return m_weights(ui, vi);
  }

  /*!
   * \brief Set the weight at a specific index
   *
   * \param [in] ui The index of the weight in on the first axis
   * \param [in] vi The index of the weight in on the second axis
   * \param [in] weight The updated value of the weight
   * \pre Requires that the surface be rational
   * \pre Requires that the weight be positive
   */
  void setWeight(int ui, int vi, T weight)
  {
    SLIC_ASSERT(isRational());
    SLIC_ASSERT(weight > 0);

    m_weights(ui, vi) = weight;
  }

  /// Checks equality of two Bezier Patches
  friend inline bool operator==(const NURBSPatch<T, NDIMS>& lhs,
                                const NURBSPatch<T, NDIMS>& rhs)
  {
    return (lhs.m_controlPoints == rhs.m_controlPoints) &&
      (lhs.m_weights == rhs.m_weights) && (lhs.m_knotvec_u == rhs.m_knotvec_u) &&
      (lhs.m_knotvec_v == rhs.m_knotvec_v);
  }

  friend inline bool operator!=(const NURBSPatch<T, NDIMS>& lhs,
                                const NURBSPatch<T, NDIMS>& rhs)
  {
    return !(lhs == rhs);
  }

  /// Returns a copy of the Bezier patch's control points
  CoordsMat getControlPoints() const { return m_controlPoints; }

  /// Returns a copy of the Bezier patch's weights
  WeightsMat getWeights() const { return m_weights; }

  /*!
   * \brief Reverses the order of one direction of the Bezier patch's control points and weights
   *
   * \param [in] axis orientation of patch. 0 to reverse in u, 1 for reverse in v
   */
  void reverseOrientation(int axis)
  {
    if(axis == 0)
    {
      reverseOrientation_u();
    }
    else
    {
      reverseOrientation_v();
    }
  }

  /// \brief Reverses the order of the control points, weights, and knots on the first axis
  void reverseOrientation_u()
  {
    auto patch_shape = m_controlPoints.shape();
    const int npts_u_mid = (patch_shape[0] + 1) / 2;

    for(int q = 0; q < patch_shape[1]; ++q)
    {
      for(int i = 0; i < npts_u_mid; ++i)
      {
        axom::utilities::swap(m_controlPoints(i, q),
                              m_controlPoints(patch_shape[0] - i - 1, q));
      }

      if(isRational())
      {
        for(int i = 0; i < npts_u_mid; ++i)
        {
          axom::utilities::swap(m_weights(i, q),
                                m_weights(patch_shape[0] - i, q));
        }
      }
    }

    m_knotvec_u.reverse();
  }

  /// \brief Reverses the order of the control points, weights, and knots on the second axis
  void reverseOrientation_v()
  {
    auto patch_shape = m_controlPoints.shape();
    const int npts_v_mid = (patch_shape[1] + 1) / 2;

    for(int p = 0; p < patch_shape[0]; ++p)
    {
      for(int i = 0; i < npts_v_mid; ++i)
      {
        axom::utilities::swap(m_controlPoints(p, i),
                              m_controlPoints(p, patch_shape[1] - i - 1));
      }

      if(isRational())
      {
        for(int i = 0; i < npts_v_mid; ++i)
        {
          axom::utilities::swap(m_weights(p, i),
                                m_weights(p, patch_shape[1] - i - 1));
        }
      }
    }

    m_knotvec_v.reverse();
  }

  /// Swap the axes such that s(u, v) becomes s(v, u)
  void swapAxes()
  {
    auto patch_shape = m_controlPoints.shape();

    CoordsMat new_controlPoints(patch_shape[1], patch_shape[0]);

    for(int p = 0; p < patch_shape[0]; ++p)
    {
      for(int q = 0; q < patch_shape[1]; ++q)
      {
        new_controlPoints(q, p) = m_controlPoints(p, q);
      }
    }

    m_controlPoints = new_controlPoints;

    if(isRational())
    {
      WeightsMat new_weights(patch_shape[1], patch_shape[0]);

      for(int p = 0; p < patch_shape[0]; ++p)
      {
        for(int q = 0; q < patch_shape[1]; ++q)
        {
          new_weights(q, p) = m_weights(p, q);
        }
      }

      m_weights = new_weights;
    }

    std::swap(m_knotvec_u, m_knotvec_v);
  }

  /// Returns an axis-aligned bounding box containing the patch
  BoundingBoxType boundingBox() const
  {
    return BoundingBoxType(m_controlPoints.data(),
                           static_cast<int>(m_controlPoints.size()));
  }

  /// Returns an oriented bounding box containing the patch
  OrientedBoundingBoxType orientedBoundingBox() const
  {
    return OrientedBoundingBoxType(m_controlPoints.data(),
                                   static_cast<int>(m_controlPoints.size()));
  }

  /*!
   * \brief Evaluates a slice Bezier patch for a fixed parameter value of \a u or \a v
   *
   * \param [in] u parameter value at which to evaluate the first axis
   * \param [in] v parameter value at which to evaluate the second axis
   * \param [in] axis orientation of curve. 0 for fixed u, 1 for fixed v
   * \return p the value of the Bezier patch at (u, v)
   *
   * \note We typically evaluate the patch at \a u or \a v between 0 and 1
   */
  NURBSCurveType isocurve(T uv, int axis) const
  {
    SLIC_ASSERT((axis == 0) || (axis == 1));

    if(axis == 0)
    {
      return isocurve_u(uv);
    }
    // else
    // {
    //   return isocurve_v(uv);
    // }
  }

  /*!
   * \brief Returns an isocurve with a fixed value of u
   *
   * \param [in] u Parameter value fixed in the isocurve
   * \return c The isocurve C(v) = S(u, v) for fixed u
   */
  NURBSCurveType isocurve_u(T u) const
  {
    using axom::utilities::lerp;

    bool isRationalPatch = isRational();

    auto patch_shape = m_controlPoints.shape();
    const int deg_u = m_knotvec_u.getDegree();
    const int deg_v = m_knotvec_v.getDegree();

    NURBSCurveType c(patch_shape[1], deg_v);

    // Find the control points by evaluating each column of the patch
    const auto span_u = m_knotvec_u.findSpan(u);
    const auto N_evals_u = m_knotvec_u.evalNonZero(span_u);
    for(int q = 0; q < patch_shape[1]; ++q)
    {
      Point<T, NDIMS + 1> H;
      for(int j = 0; j <= deg_u; ++j)
      {
        const auto offset = span_u - deg_u + j;
        const T weight = isRationalPatch ? m_weights(offset, q) : 1.0;
        const auto& controlPoint = m_controlPoints(offset, q);

        for(int i = 0; i < NDIMS; ++i)
        {
          H[i] += N_evals_u[j] * weight * controlPoint[i];
        }
        H[NDIMS] += N_evals_u[j] * weight;
      }

      for(int i = 0; i < NDIMS; ++i)
      {
        c(q)[i] = H[i] / H[NDIMS];
      }
    }

    c.m_knotvec = m_knotvec_v;

    return c;
  }

  //   /*!
  //    * \brief Returns an isocurve with a fixed value of v
  //    *
  //    * \param [in] v Parameter value fixed in the isocurve
  //    * \return c The isocurve C(u) = S(u, v) for fixed v
  //    */
  //   BezierCurveType isocurve_v(T v) const
  //   {
  //     using axom::utilities::lerp;

  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     BezierCurveType c(ord_u);
  //     axom::Array<T> dCarray(ord_v + 1);

  //     if(!isRational())
  //     {
  //       // If looking for an edge, don't need to use De Casteljau
  //       if(v == 0.0)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           c[p] = m_controlPoints(p, 0);
  //         }
  //         return c;
  //       }
  //       if(v == 1.0)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           c[p] = m_controlPoints(p, ord_v);
  //         }
  //         return c;
  //       }

  //       // Otherwise, do De Casteljau along the u-axis
  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int i = 0; i < NDIMS; ++i)
  //         {
  //           for(int q = 0; q <= ord_v; ++q)
  //           {
  //             dCarray[q] = m_controlPoints(p, q)[i];
  //           }

  //           for(int q = 1; q <= ord_v; ++q)
  //           {
  //             const int end = ord_v - q;
  //             for(int k = 0; k <= end; ++k)
  //             {
  //               dCarray[k] = lerp(dCarray[k], dCarray[k + 1], v);
  //             }
  //           }

  //           c[p][i] = dCarray[0];
  //         }
  //       }

  //       return c;
  //     }
  //     else
  //     {
  //       c.makeRational();

  //       // If looking for an edge, don't need to use De Casteljau
  //       if(v == 0.0)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           c[p] = m_controlPoints(p, 0);
  //           c.setWeight(p, m_weights(p, 0));
  //         }
  //         return c;
  //       }
  //       if(v == 1.0)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           c[p] = m_controlPoints(p, ord_v);
  //           c.setWeight(p, m_weights(p, ord_v));
  //         }
  //         return c;
  //       }

  //       // Otherwise, do de Casteljau on the homogeneous control points

  //       // Store BezierPatch of projective weights, (wx, wy, wz)
  //       //  and BezierPatch of weights (w)
  //       BezierPatch<T, NDIMS> projective(ord_u, ord_v);
  //       BezierPatch<T, 1> weights(ord_u, ord_v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           weights(p, q)[0] = m_weights(p, q);

  //           for(int i = 0; i < NDIMS; ++i)
  //           {
  //             projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
  //           }
  //         }
  //       }

  //       BezierCurve<T, NDIMS> P = projective.isocurve_v(v);
  //       BezierCurve<T, 1> W = weights.isocurve_v(v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int i = 0; i < NDIMS; ++i)
  //         {
  //           c[p][i] = P[p][i] / W[p][0];
  //           c.setWeight(p, W[p][0]);
  //         }
  //       }

  //       return c;
  //     }
  //   }

  //   /*!
  //    * \brief Evaluates a Bezier patch at a particular parameter value (\a u, \a v)
  //    *
  //    * \param [in] u parameter value at which to evaluate on the first axis
  //    * \param [in] v parameter value at which to evaluate on the second axis
  //    * \return p the value of the Bezier patch at (u, v)
  //    *
  //    * \note We typically evaluate the patch at \a u and \a v between 0 and 1
  //    */
  //   PointType evaluate(T u, T v) const
  //   {
  //     if(getOrder_u() >= getOrder_v())
  //     {
  //       return isocurve_u(u).evaluate(v);
  //     }
  //     else
  //     {
  //       return isocurve_v(v).evaluate(u);
  //     }
  //   }

  //   /*!
  //    * \brief Evaluates all first derivatives Bezier patch at (\a u, \a v)
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \param [out] eval The point value of the Bezier patch at (u, v)
  //    * \param [out] Du The vector value of S_u(u, v)
  //    * \param [out] Dv The vector value of S_v(u, v)
  //    *
  //    * \note We typically evaluate the patch at \a u and \a v between 0 and 1
  //    */
  //   void evaluate_first_derivatives(T u,
  //                                   T v,
  //                                   Point<T, NDIMS>& eval,
  //                                   Vector<T, NDIMS>& Du,
  //                                   Vector<T, NDIMS>& Dv) const
  //   {
  //     using axom::utilities::lerp;
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

  //     if(!isRational())
  //     {
  //       // Do de Casteljau until we get a 2x2
  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           for(int q = 0; q <= ord_v; ++q)
  //           {
  //             dCmat(p, q) = m_controlPoints(p, q)[i];
  //           }
  //         }

  //         // Store the size after each de Casteljau reduction is made
  //         int end_u = ord_u;
  //         int end_v = ord_v;

  //         // Do de Casteljau over the longer direction first
  //         if(ord_u >= ord_v)
  //         {
  //           // Stop 1 steps early in u
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }

  //           // Stop 1 steps early in v
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }
  //         }
  //         else
  //         {
  //           // Stop 1 steps early in v
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }

  //           // Stop 2 steps early in v over the three columns
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }
  //         }

  //         // By now, the first 2x2 submatrix of dCmat should have
  //         //  all the intermediate values needed

  //         // Compute first order derivatives
  //         // clang-format off
  //         if( end_u == 1 && end_v == 1 )
  //         {
  //           Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) +
  //                                   v * (dCmat(1, 1) - dCmat(0, 1)) );
  //           Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) +
  //                                   u * (dCmat(1, 1) - dCmat(1, 0)) );
  //           eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
  //         }
  //         else if (end_u == 1 && end_v == 0)
  //         {
  //           Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
  //           Dv[i] = 0.0;
  //           eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
  //         }
  //         else if (end_u == 0 && end_v == 1)
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
  //           eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
  //         }
  //         else
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = 0.0;
  //           eval[i] = dCmat(0, 0);
  //         }
  //         // clang-format on
  //       }
  //     }
  //     else
  //     {
  //       // Store BezierPatch of projective weights, (wx, wy, wz)
  //       //  and BezierPatch of weights (w)
  //       BezierPatch<T, NDIMS> projective(ord_u, ord_v);
  //       BezierPatch<T, 1> weights(ord_u, ord_v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           weights(p, q)[0] = m_weights(p, q);

  //           for(int i = 0; i < NDIMS; ++i)
  //           {
  //             projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
  //           }
  //         }
  //       }

  //       Point<T, NDIMS> P;
  //       Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

  //       Point<T, 1> W;
  //       Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

  //       projective.evaluate_second_derivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
  //       weights.evaluate_second_derivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         eval[i] = P[i] / W[0];
  //         Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
  //         Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
  //       }
  //     }
  //   }

  //   /*!
  //    * \brief Evaluates all linear derivatives Bezier patch at (\a u, \a v)
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \param [out] eval The point value of the Bezier patch at (u, v)
  //    * \param [out] Du The vector value of S_u(u, v)
  //    * \param [out] Dv The vector value of S_v(u, v)
  //    * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
  //    *
  //    * \note We typically evaluate the patch at \a u and \a v between 0 and 1
  //    */
  //   void evaluate_linear_derivatives(T u,
  //                                    T v,
  //                                    Point<T, NDIMS>& eval,
  //                                    Vector<T, NDIMS>& Du,
  //                                    Vector<T, NDIMS>& Dv,
  //                                    Vector<T, NDIMS>& DuDv) const
  //   {
  //     using axom::utilities::lerp;
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

  //     if(!isRational())
  //     {
  //       // Do de Casteljau until we get a 2x2
  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           for(int q = 0; q <= ord_v; ++q)
  //           {
  //             dCmat(p, q) = m_controlPoints(p, q)[i];
  //           }
  //         }

  //         // Store the size after each de Casteljau reduction is made
  //         int end_u = ord_u;
  //         int end_v = ord_v;

  //         // Do de Casteljau over the longer direction first
  //         if(ord_u >= ord_v)
  //         {
  //           // Stop 1 steps early in u
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }

  //           // Stop 1 steps early in v
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }
  //         }
  //         else
  //         {
  //           // Stop 1 steps early in v
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }

  //           // Stop 2 steps early in v over the three columns
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }
  //         }

  //         // By now, the first 2x2 submatrix of dCmat should have
  //         //  all the intermediate values needed

  //         // Compute first order derivatives
  //         // clang-format off
  //         if( end_u == 1 && end_v == 1 )
  //         {
  //           Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) +
  //                                   v * (dCmat(1, 1) - dCmat(0, 1)) );
  //           Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) +
  //                                   u * (dCmat(1, 1) - dCmat(1, 0)) );
  //           DuDv[i] = ord_u * ord_v * (dCmat(1, 1) - dCmat(1, 0) - dCmat(0, 1) + dCmat(0, 0));
  //           eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
  //         }
  //         else if (end_u == 1 && end_v == 0)
  //         {
  //           Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
  //           Dv[i] = 0.0;
  //           DuDv[i] = 0.0;
  //           eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
  //         }
  //         else if (end_u == 0 && end_v == 1)
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
  //           DuDv[i] = 0.0;
  //           eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
  //         }
  //         else
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = 0.0;
  //           DuDv[i] = 0.0;
  //           eval[i] = dCmat(0, 0);
  //         }
  //         // clang-format on
  //       }
  //     }
  //     else
  //     {
  //       // Store BezierPatch of projective weights, (wx, wy, wz)
  //       //  and BezierPatch of weights (w)
  //       BezierPatch<T, NDIMS> projective(ord_u, ord_v);
  //       BezierPatch<T, 1> weights(ord_u, ord_v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           weights(p, q)[0] = m_weights(p, q);

  //           for(int i = 0; i < NDIMS; ++i)
  //           {
  //             projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
  //           }
  //         }
  //       }

  //       Point<T, NDIMS> P;
  //       Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

  //       Point<T, 1> W;
  //       Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

  //       projective.evaluate_second_derivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
  //       weights.evaluate_second_derivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         eval[i] = P[i] / W[0];
  //         Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
  //         Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
  //         DuDv[i] =
  //           (P_uv[i] - Du[i] * W_v[0] - Dv[i] * W_u[0] - eval[i] * W_uv[0]) / W[0];
  //       }
  //     }
  //   }

  //   /*!
  //    * \brief Evaluates all second derivatives Bezier patch at (\a u, \a v)
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \param [out] eval The point value of the Bezier patch at (u, v)
  //    * \param [out] Du The vector value of S_u(u, v)
  //    * \param [out] Dv The vector value of S_v(u, v)
  //    * \param [out] DuDu The vector value of S_uu(u, v)
  //    * \param [out] DvDv The vector value of S_vv(u, v)
  //    * \param [out] DuDv The vector value of S_uv(u, v) == S_vu(u, v)
  //    *
  //    * \note We typically evaluate the patch at \a u and \a v between 0 and 1
  //    */
  //   void evaluate_second_derivatives(T u,
  //                                    T v,
  //                                    Point<T, NDIMS>& eval,
  //                                    Vector<T, NDIMS>& Du,
  //                                    Vector<T, NDIMS>& Dv,
  //                                    Vector<T, NDIMS>& DuDu,
  //                                    Vector<T, NDIMS>& DvDv,
  //                                    Vector<T, NDIMS>& DuDv) const
  //   {
  //     using axom::utilities::lerp;
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);

  //     if(!isRational())
  //     {
  //       // Do de Casteljau until we get a 3x3
  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           for(int q = 0; q <= ord_v; ++q)
  //           {
  //             dCmat(p, q) = m_controlPoints(p, q)[i];
  //           }
  //         }

  //         // Store the size after each de Casteljau reduction is made
  //         int end_u = ord_u;
  //         int end_v = ord_v;

  //         // Do de Casteljau over the longer direction first
  //         if(ord_u >= ord_v)
  //         {
  //           // Stop 2 steps early in u
  //           while(end_u > 2)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }

  //           // Stop 2 steps early in v
  //           while(end_v > 2)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }
  //         }
  //         else
  //         {
  //           // Stop 2 steps early in v
  //           while(end_v > 2)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }

  //           // Stop 2 steps early in v over the three columns
  //           while(end_u > 2)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= end_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }
  //         }

  //         // By now, the first 3x3 submatrix of dCmat should have
  //         //  all the intermediate values needed

  //         // clang-format off
  //         // Get second order derivatives
  //         if( ord_u >= 2 )
  //         {
  //           if (ord_v == 0)
  //           {
  //             DuDu[i] = (ord_u - 1) * ord_u *
  //                                     (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0));
  //           }
  //           else if (ord_v == 1)
  //           {
  //             DuDu[i] = (ord_u - 1) * ord_u *
  //                         ( (1 - v) * (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0)) +
  //                                 v * (dCmat(2, 1) - 2 * dCmat(1, 1) + dCmat(0, 1)) );
  //           }
  //           else
  //           {
  //             DuDu[i] = (ord_u - 1) * ord_u *
  //               ( (1 - v) * (1 - v) * (dCmat(2, 0) - 2 * dCmat(1, 0) + dCmat(0, 0)) +
  //                   2 * v * (1 - v) * (dCmat(2, 1) - 2 * dCmat(1, 1) + dCmat(0, 1)) +
  //                             v * v * (dCmat(2, 2) - 2 * dCmat(1, 2) + dCmat(0, 2)) );
  //           }
  //         }
  //         else
  //         {
  //           DuDu[i] = 0.0;
  //         }

  //         if( ord_v >= 2 )
  //         {
  //           if( ord_u == 0 )
  //           {
  //             DvDv[i] = (ord_v - 1) * ord_v *
  //                                     (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0));
  //           }
  //           else if (ord_u == 1)
  //           {
  //             DvDv[i] = (ord_v - 1) * ord_v *
  //                         ( (1 - u) * (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0)) +
  //                                 u * (dCmat(1, 2) - 2 * dCmat(1, 1) + dCmat(1, 0)) );
  //           }
  //           else
  //           {
  //             DvDv[i] = (ord_v - 1) * ord_v *
  //               ( (1 - u) * (1 - u) * (dCmat(0, 2) - 2 * dCmat(0, 1) + dCmat(0, 0)) +
  //                   2 * u * (1 - u) * (dCmat(1, 2) - 2 * dCmat(1, 1) + dCmat(1, 0)) +
  //                             u * u * (dCmat(2, 2) - 2 * dCmat(2, 1) + dCmat(2, 0)) );
  //           }
  //         }
  //         else
  //         {
  //           DvDv[i] = 0.0;
  //         }

  //         // Do one de Casteljau iteration in u
  //         while( end_u > 1 )
  //         {
  //           end_u -= 1;
  //           for(int k = 0; k <= end_u; ++k)
  //           {
  //             for(int q = 0; q <= end_v; ++q)
  //             {
  //               dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //             }
  //           }
  //         }

  //         // Do (at most) one de Casteljau iteration in v
  //         while(end_v > 1)
  //         {
  //           end_v -= 1;
  //           for(int k = 0; k <= end_v; ++k)
  //           {
  //             for(int p = 0; p <= end_u; ++p)
  //             {
  //               dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //             }
  //           }
  //         }

  //         // Compute first order derivatives
  //         // clang-format off
  //         if( end_u == 1 && end_v == 1 )
  //         {
  //           Du[i] = ord_u * ( (1 - v) * (dCmat(1, 0) - dCmat(0, 0)) +
  //                                   v * (dCmat(1, 1) - dCmat(0, 1)) );
  //           Dv[i] = ord_v * ( (1 - u) * (dCmat(0, 1) - dCmat(0, 0)) +
  //                                   u * (dCmat(1, 1) - dCmat(1, 0)) );
  //           DuDv[i] = ord_u * ord_v * (dCmat(1, 1) - dCmat(1, 0) - dCmat(0, 1) + dCmat(0, 0));
  //           eval[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) + (1 - u) * v * dCmat(0, 1) + u * v * dCmat(1, 1);
  //         }
  //         else if (end_u == 1 && end_v == 0)
  //         {
  //           Du[i] = ord_u * (dCmat(1, 0) - dCmat(0, 0));
  //           Dv[i] = 0.0;
  //           DuDv[i] = 0.0;
  //           eval[i] = (1 - u) * dCmat(0, 0) + u * dCmat(1, 0);
  //         }
  //         else if (end_u == 0 && end_v == 1)
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = ord_v * (dCmat(0, 1) - dCmat(0, 0));
  //           DuDv[i] = 0.0;
  //           eval[i] = (1 - v) * dCmat(0, 0) + v * dCmat(0, 1);
  //         }
  //         else
  //         {
  //           Du[i] = 0.0;
  //           Dv[i] = 0.0;
  //           DuDv[i] = 0.0;
  //           eval[i] = dCmat(0, 0);
  //         }
  //         // clang-format on
  //       }
  //     }
  //     else
  //     {
  //       // Store BezierPatch of projective weights, (wx, wy, wz)
  //       //  and BezierPatch of weights (w)
  //       BezierPatch<T, NDIMS> projective(ord_u, ord_v);
  //       BezierPatch<T, 1> weights(ord_u, ord_v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           weights(p, q)[0] = m_weights(p, q);

  //           for(int i = 0; i < NDIMS; ++i)
  //           {
  //             projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
  //           }
  //         }
  //       }

  //       Point<T, NDIMS> P;
  //       Vector<T, NDIMS> P_u, P_v, P_uu, P_vv, P_uv;

  //       Point<T, 1> W;
  //       Vector<T, 1> W_u, W_v, W_uu, W_vv, W_uv;

  //       projective.evaluate_second_derivatives(u, v, P, P_u, P_v, P_uu, P_vv, P_uv);
  //       weights.evaluate_second_derivatives(u, v, W, W_u, W_v, W_uu, W_vv, W_uv);

  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         eval[i] = P[i] / W[0];
  //         Du[i] = (P_u[i] - eval[i] * W_u[0]) / W[0];
  //         Dv[i] = (P_v[i] - eval[i] * W_v[0]) / W[0];
  //         DuDu[i] = (P_uu[i] - 2 * W_u[0] * Du[i] - eval[i] * W_uu[0]) / W[0];
  //         DvDv[i] = (P_vv[i] - 2 * W_v[0] * Dv[i] - eval[i] * W_vv[0]) / W[0];
  //         DuDv[i] =
  //           (P_uv[i] - Du[i] * W_v[0] - Dv[i] * W_u[0] - eval[i] * W_uv[0]) / W[0];
  //       }
  //     }
  //   }

  //   /*!
  //    * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the u axis
  //    *
  //    * \param [in] u parameter value at which to evaluate on the first axis
  //    * \param [in] v parameter value at which to evaluate on the second axis
  //    * \return vec a tangent vector in u of the Bezier patch at (u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType du(T u, T v) const { return isocurve_v(v).dt(u); }

  //   /*!
  //    * \brief Computes the second derivative of a Bezier patch at (\a u, \a v) along the u axis
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \return vec The vector value of S_uu(u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType dudu(T u, T v) const { return isocurve_v(v).dtdt(u); }

  //   /*!
  //    * \brief Computes a tangent of a Bezier patch at a particular parameter value (\a u, \a v) along the v axis
  //    *
  //    * \param [in] u parameter value at which to evaluate on the first axis
  //    * \param [in] v parameter value at which to evaluate on the second axis
  //    * \return vec a tangent vector in v of the Bezier patch at (u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType dv(T u, T v) const { return isocurve_u(u).dt(v); }

  //   /*!
  //    * \brief Computes the second derivative of a Bezier patch at (\a u, \a v) along the v axis
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \return vec The vector value of S_vv(u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType dvdv(T u, T v) const { return isocurve_u(u).dtdt(v); }

  //   /*!
  //    * \brief Computes the mixed second derivative of a Bezier patch at (\a u, \a v)
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \return vec The vector value of S_uv(u, v) == S_vu(u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType dudv(T u, T v) const
  //   {
  //     using axom::utilities::lerp;

  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     // If the curve is nonrational, we can use standard de Casteljau
  //     if(!isRational())
  //     {
  //       Vector<T, NDIMS> val;

  //       axom::Array<T, 2> dCmat(ord_u + 1, ord_v + 1);
  //       axom::Array<T, 2> dWmat(ord_u + 1, ord_v + 1);

  //       // Do de Casteljau until we get a 2x2
  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         for(int p = 0; p <= ord_u; ++p)
  //         {
  //           for(int q = 0; q <= ord_v; ++q)
  //           {
  //             dCmat(p, q) = m_controlPoints(p, q)[i];
  //           }
  //         }

  //         // Store the size after each de Casteljau reduction is made
  //         int end_u = ord_u;
  //         int end_v = ord_v;

  //         // Do de Casteljau over the longer direction first
  //         if(ord_u >= ord_v)
  //         {
  //           // Stop 1 steps early in u
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= ord_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }

  //           // Stop 1 steps early in u
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }
  //         }
  //         else
  //         {
  //           // Stop 1 step early in v
  //           while(end_v > 1)
  //           {
  //             end_v -= 1;
  //             for(int k = 0; k <= end_v; ++k)
  //             {
  //               for(int p = 0; p <= end_u; ++p)
  //               {
  //                 dCmat(p, k) = lerp(dCmat(p, k), dCmat(p, k + 1), v);
  //               }
  //             }
  //           }

  //           // Stop 1 step early in v over the two columns
  //           while(end_u > 1)
  //           {
  //             end_u -= 1;
  //             for(int k = 0; k <= end_u; ++k)
  //             {
  //               for(int q = 0; q <= ord_v; ++q)
  //               {
  //                 dCmat(k, q) = lerp(dCmat(k, q), dCmat(k + 1, q), u);
  //               }
  //             }
  //           }
  //         }

  //         // Do the evaluation, accounting for insufficient control points
  //         if(end_u == 0 || end_v == 0)
  //         {
  //           val[i] = 0.0;
  //         }
  //         else
  //         {
  //           // clang-format off
  //         val[i] = (1 - u) * (1 - v) * dCmat(0, 0) + u * (1 - v) * dCmat(1, 0) +
  //                        (1 - u) * v * dCmat(0, 1) +       u * v * dCmat(1, 1);
  //           // clang-format on
  //         }
  //       }

  //       return val;
  //     }
  //     // If rational, construct the 4D homogeneous surface,
  //     //  which requires all first derivatives
  //     else
  //     {
  //       Vector<T, NDIMS> val;

  //       // Store BezierPatch of projective weights, (wx, wy, wz)
  //       //  and BezierPatch of weights (w)
  //       BezierPatch<T, NDIMS> projective(ord_u, ord_v);
  //       BezierPatch<T, 1> weights(ord_u, ord_v);

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           weights(p, q)[0] = m_weights(p, q);

  //           for(int i = 0; i < NDIMS; ++i)
  //           {
  //             projective(p, q)[i] = m_controlPoints(p, q)[i] * m_weights(p, q);
  //           }
  //         }
  //       }

  //       Point<T, NDIMS> P;
  //       Vector<T, NDIMS> P_u, P_v, P_uv;

  //       Point<T, 1> W;
  //       Vector<T, 1> W_u, W_v, W_uv;

  //       projective.evaluate_linear_derivatives(u, v, P, P_u, P_v, P_uv);
  //       weights.evaluate_linear_derivatives(u, v, W, W_u, W_v, W_uv);

  //       // Store values used in each coordinate computation
  //       double weight_prod = 2 * W_u[0] * W_v[0] - W[0] * W_uv[0];
  //       double weight_cubed = W[0] * W[0] * W[0];
  //       for(int i = 0; i < NDIMS; ++i)
  //       {
  //         val[i] = W[0] * W[0] * P_uv[i] -
  //           W[0] * (P_u[i] * W_v[0] + P_v[i] * W_u[0]) + P[i] * weight_prod;
  //         val[i] /= weight_cubed;
  //       }

  //       return val;
  //     }
  //   }

  //   /*!
  //    * \brief Computes the mixed second derivative of a Bezier patch at (\a u, \a v)
  //    *
  //    * \param [in] u Parameter value at which to evaluate on the first axis
  //    * \param [in] v Parameter value at which to evaluate on the second axis
  //    * \return vec The vector value of S_uv(u, v) == S_vu(u, v)
  //    *
  //    * \note We typically find the tangent of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType dvdu(T u, T v) const { return dudv(u, v); }

  //   /*!
  //    * \brief Computes the normal vector of a Bezier patch at a particular parameter value (\a u, \a v)
  //    *
  //    * \param [in] u parameter value at which to evaluate on the first axis
  //    * \param [in] v parameter value at which to evaluate on the second axis
  //    * \return vec the normal vector of the Bezier patch at (u, v)
  //    *
  //    * \note We typically find the normal of the patch at \a u and \a v between 0 and 1
  //    */
  //   VectorType normal(T u, T v) const
  //   {
  //     Point<T, NDIMS> eval;
  //     Vector<T, NDIMS> Du, Dv;
  //     evaluate_first_derivatives(u, v, eval, Du, Dv);
  //     return VectorType::cross_product(Du, Dv);
  //   }

  //   /*!
  //    * \brief Splits a Bezier patch into two Bezier patches
  //    *
  //    * \param [in] uv parameter value between 0 and 1 at which to bisect the patch
  //    * \param [in] axis orientation of split. 0 for fixed u, 1 for fixed v
  //    * \param [out] p1 First output Bezier patch
  //    * \param [out] p2 Second output Bezier patch
  //    *
  //    * \pre Parameter \a uv must be between 0 and 1
  //    */
  //   void split(T uv, int axis, BezierPatch& p1, BezierPatch& p2) const
  //   {
  //     SLIC_ASSERT((axis == 0) || (axis == 1));

  //     if(axis == 0)
  //     {
  //       split_u(uv, p1, p2);
  //     }
  //     else
  //     {
  //       split_v(uv, p1, p2);
  //     }
  //   }

  //   /// Split the patch along a fixed value of u
  //   void split_u(T u, BezierPatch& p1, BezierPatch& p2) const
  //   {
  //     using axom::utilities::lerp;

  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     p1.setOrder(ord_u, ord_v);

  //     // Note: The second patch's control points are computed inline
  //     //       as we find the first patch's control points
  //     p2 = *this;

  //     if(isRational())
  //     {
  //       p1.makeRational();  // p2 already rational

  //       for(int q = 0; q <= ord_v; ++q)
  //       {
  //         // Do the rational de Casteljau algorithm
  //         p1(0, q) = p2(0, q);
  //         p1.setWeight(0, q, p2.getWeight(0, q));

  //         for(int p = 1; p <= ord_u; ++p)
  //         {
  //           const int end = ord_u - p;

  //           for(int k = 0; k <= end; ++k)
  //           {
  //             double temp_weight =
  //               lerp(p2.getWeight(k, q), p2.getWeight(k + 1, q), u);

  //             for(int i = 0; i < NDIMS; ++i)
  //             {
  //               p2(k, q)[i] = lerp(p2.getWeight(k, q) * p2(k, q)[i],
  //                                  p2.getWeight(k + 1, q) * p2(k + 1, q)[i],
  //                                  u) /
  //                 temp_weight;
  //             }

  //             p2.setWeight(k, q, temp_weight);
  //           }

  //           p1(p, q) = p2(0, q);
  //           p1.setWeight(p, q, p2.getWeight(0, q));
  //         }
  //       }
  //     }
  //     else
  //     {
  //       for(int q = 0; q <= ord_v; ++q)
  //       {
  //         p1(0, q) = m_controlPoints(0, q);

  //         for(int i = 0; i < NDIMS; ++i)
  //         {
  //           for(int p = 1; p <= ord_u; ++p)
  //           {
  //             const int end = ord_u - p;
  //             for(int k = 0; k <= end; ++k)
  //             {
  //               p2(k, q)[i] = lerp(p2(k, q)[i], p2(k + 1, q)[i], u);
  //             }

  //             p1(p, q)[i] = p2(0, q)[i];
  //           }
  //         }
  //       }
  //     }
  //   }

  //   void split_v(T v, BezierPatch& p1, BezierPatch& p2) const
  //   {
  //     using axom::utilities::lerp;

  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     p1.setOrder(ord_u, ord_v);

  //     // Note: The second patch's control points are computed inline
  //     //       as we find the first patch's control points
  //     p2 = *this;

  //     if(isRational())
  //     {
  //       p1.makeRational();  // p2 already rational

  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         // Do the rational de Casteljau algorithm
  //         p1(p, 0) = p2(p, 0);
  //         p1.setWeight(p, 0, p2.getWeight(p, 0));

  //         for(int q = 1; q <= ord_v; ++q)
  //         {
  //           const int end = ord_v - q;

  //           for(int k = 0; k <= end; ++k)
  //           {
  //             double temp_weight =
  //               lerp(p2.getWeight(p, k), p2.getWeight(p, k + 1), v);

  //             for(int i = 0; i < NDIMS; ++i)
  //             {
  //               p2(p, k)[i] = lerp(p2.getWeight(p, k) * p2(p, k)[i],
  //                                  p2.getWeight(p, k + 1) * p2(p, k + 1)[i],
  //                                  v) /
  //                 temp_weight;
  //             }

  //             p2.setWeight(p, k, temp_weight);
  //           }

  //           p1(p, q) = p2(p, 0);
  //           p1.setWeight(p, q, p2.getWeight(p, 0));
  //         }
  //       }
  //     }
  //     else
  //     {
  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         p1(p, 0) = m_controlPoints(p, 0);

  //         for(int i = 0; i < NDIMS; ++i)
  //         {
  //           for(int q = 1; q <= ord_v; ++q)
  //           {
  //             const int end = ord_v - q;
  //             for(int k = 0; k <= end; ++k)
  //             {
  //               p2(p, k)[i] = lerp(p2(p, k)[i], p2(p, k + 1)[i], v);
  //             }

  //             p1(p, q)[i] = p2(p, 0)[i];
  //           }
  //         }
  //       }
  //     }
  //   }

  //   /*!
  //    * \brief Splits a Bezier patch into four Bezier patches
  //    *
  //    * \param [in] u parameter value between 0 and 1 at which to bisect on the first axis
  //    * \param [in] v parameter value between 0 and 1 at which to bisect on the second axis
  //    * \param [out] p1 First output Bezier patch
  //    * \param [out] p2 Second output Bezier patch
  //    * \param [out] p3 Third output Bezier patch
  //    * \param [out] p4 Fourth output Bezier patch
  //    *
  //    *   v = 1
  //    *   ----------------------
  //    *   |         |          |
  //    *   |   p3    |    p4    |
  //    *   |         |          |
  //    *   --------(u,v)---------
  //    *   |         |          |
  //    *   |   p1    |    p2    |
  //    *   |         |          |
  //    *   ---------------------- u = 1
  //    *
  //    * \pre Parameter \a u and \a v must be between 0 and 1
  //    */
  //   void split(T u,
  //              T v,
  //              BezierPatch& p1,
  //              BezierPatch& p2,
  //              BezierPatch& p3,
  //              BezierPatch& p4) const
  //   {
  //     // Bisect the patch along the u direction
  //     split_u(u, p1, p2);

  //     // Temporarily store the result in each half and split again
  //     BezierPatch p0(p1);
  //     p0.split_v(v, p1, p3);

  //     p0 = p2;
  //     p0.split_v(v, p2, p4);
  //   }

  //   /*!
  //    * \brief Predicate to check if the Bezier patch is approximately planar
  //    *
  //    * This function checks if all control points of the BezierPatch
  //    * are approximately on the plane defined by its four corners
  //    *
  //    * \param [in] tol Threshold for sum of squared distances
  //    * \return True if c1 is near-planar
  //    */
  //   bool isPlanar(double tol = 1E-8) const
  //   {
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     if(ord_u <= 0 && ord_v <= 0)
  //     {
  //       return true;
  //     }
  //     if(ord_u == 1 && ord_v == 0)
  //     {
  //       return true;
  //     }
  //     if(ord_u == 0 && ord_v == 1)
  //     {
  //       return true;
  //     }

  //     // Check that the four corners aren't coplanar
  //     VectorType v1(m_controlPoints(0, 0), m_controlPoints(0, ord_v));
  //     VectorType v2(m_controlPoints(0, 0), m_controlPoints(ord_u, 0));
  //     VectorType v3(m_controlPoints(0, 0), m_controlPoints(ord_u, ord_v));
  //     if(!axom::utilities::isNearlyEqual(
  //          VectorType::scalar_triple_product(v1, v2, v3),
  //          0.0,
  //          tol))
  //     {
  //       return false;
  //     }

  //     // Find three points that produce a nonzero normal
  //     Vector3D plane_normal = VectorType::cross_product(v1, v2);
  //     if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, tol))
  //     {
  //       plane_normal = VectorType::cross_product(v1, v3);
  //     }
  //     if(axom::utilities::isNearlyEqual(plane_normal.norm(), 0.0, tol))
  //     {
  //       plane_normal = VectorType::cross_product(v2, v3);
  //     }
  //     plane_normal = plane_normal.unitVector();

  //     double sqDist = 0.0;

  //     // Check all control points for simplicity
  //     for(int p = 0; p <= ord_u && sqDist <= tol; ++p)
  //     {
  //       for(int q = ((p == 0) ? 1 : 0); q <= ord_v && sqDist <= tol; ++q)
  //       {
  //         const double signedDist =
  //           plane_normal.dot(m_controlPoints(p, q) - m_controlPoints(0, 0));
  //         sqDist += signedDist * signedDist;
  //       }
  //     }

  //     return (sqDist <= tol);
  //   }

  //   /*!
  //    * \brief Predicate to check if the patch can be approximated by a polygon
  //    *
  //    * This function checks if a BezierPatch lies in a plane
  //    *  and that the edged are linear up to tolerance `tol`
  //    *
  //    * \param [in] tol Threshold for sum of squared distances
  //    * \return True if c1 is near-planar-polygonal
  //    */
  //   bool isPolygonal(double tol = 1E-8) const
  //   {
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     if(ord_u <= 0 && ord_v <= 0)
  //     {
  //       return true;
  //     }
  //     if(ord_u == 1 && ord_v == 0)
  //     {
  //       return true;
  //     }
  //     if(ord_u == 0 && ord_v == 1)
  //     {
  //       return true;
  //     }

  //     // Check if the patch is planar
  //     if(!isPlanar(tol))
  //     {
  //       return false;
  //     }

  //     // Check if each bounding curve is linear
  //     if(!isocurve_u(0).isLinear(tol))
  //     {
  //       return false;
  //     }
  //     if(!isocurve_v(0).isLinear(tol))
  //     {
  //       return false;
  //     }
  //     if(!isocurve_u(1).isLinear(tol))
  //     {
  //       return false;
  //     }
  //     if(!isocurve_v(1).isLinear(tol))
  //     {
  //       return false;
  //     }

  //     return true;
  //   }

  //   /*!
  //    * \brief Simple formatted print of a Bezier Patch instance
  //    *
  //    * \param os The output stream to write to
  //    * \return A reference to the modified ostream
  //    */
  //   std::ostream& print(std::ostream& os) const
  //   {
  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     os << "{ order (" << ord_u << ',' << ord_v << ") Bezier Patch ";

  //     for(int p = 0; p <= ord_u; ++p)
  //     {
  //       for(int q = 0; q <= ord_v; ++q)
  //       {
  //         os << m_controlPoints(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
  //       }
  //     }

  //     if(isRational())
  //     {
  //       os << ", weights [";
  //       for(int p = 0; p <= ord_u; ++p)
  //       {
  //         for(int q = 0; q <= ord_v; ++q)
  //         {
  //           os << m_weights(p, q) << ((p < ord_u || q < ord_v) ? "," : "");
  //         }
  //       }
  //     }
  //     os << "}";

  //     return os;
  //   }

  // private:
  //   /// \brief Private function to make the given knot vector uniform
  //   void makeKnotsUniform(KnotsVec& knots, int npts, int deg)
  //   {
  //     knots.resize(npts + deg + 1);

  //     // Knots for the clamped curve
  //     for(int i = 0; i < deg + 1; ++i)
  //     {
  //       knots[i] = 0.0;
  //       knots[npts + deg - i] = 1.0;
  //     }

  //     // Interior knots (if any)
  //     for(int i = 0; i < npts - deg - 1; ++i)
  //     {
  //       knots[deg + 1 + i] = (i + 1.0) / (npts - deg);
  //     }
  //   }

  //   /// Check that the weights used are positive, and
  //   ///  that there is one for each control node
  //   bool isValidRational() const
  //   {
  //     if(!isRational())
  //     {
  //       return true;
  //     }

  //     const int ord_u = getOrder_u();
  //     const int ord_v = getOrder_v();

  //     if(m_weights.shape()[0] != (ord_u + 1) || m_weights.shape()[1] != (ord_v + 1))
  //     {
  //       return false;
  //     }

  //     for(int p = 0; p <= ord_u; ++p)
  //     {
  //       for(int q = 0; q <= ord_v; ++q)
  //       {
  //         if(m_weights(p, q) <= 0)
  //         {
  //           return false;
  //         }
  //       }
  //     }

  //     return true;
  //   }

  CoordsMat m_controlPoints;
  WeightsMat m_weights;
  KnotVectorType m_knotvec_u, m_knotvec_v;
};

//------------------------------------------------------------------------------
/// Free functions related to BezierPatch
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const NURBSPatch<T, NDIMS>& nPatch)
{
  nPatch.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_NURBSPATCH_HPP_
