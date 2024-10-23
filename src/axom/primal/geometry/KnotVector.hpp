// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file NURBSCurve.hpp
 *
 * \brief A NURBS curve primitive
 */

#ifndef AXOM_PRIMAL_KNOTVECTOR_HPP
#define AXOM_PRIMAL_KNOTVECTOR_HPP

#include "axom/core.hpp"
#include "axom/slic.hpp"

#include <ostream>
#include "axom/fmt.hpp"

namespace axom
{
namespace primal
{
// Forward declare the templated classes and operator functions
template <typename T>
class KnotVector;

/*! \brief Overloaded output operator for Bezier Curves*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const KnotVector<T>& kVec);

/*!
 * \class KnotVector
 *
 * \brief Represents the knot vector for a B-Spline/NURBS Curve/Surface
 * \tparam T the coordinate type, e.g., double, float, etc.
 * 
 * Contains methods for knot vector manipulation shared between NURBS curves and surfaces.
 * This class maintains the degree and requirements on a knot vector, namely
 *  > Monotonicity (u_i <= u_i+1)
 *  > Open/Clamped-ness (degree + 1 knots at beginning and end)
 *  > Continuity (maximum internal knot multiplicity = degree)
 */
template <typename T>
class KnotVector
{
public:
  AXOM_STATIC_ASSERT_MSG(
    std::is_arithmetic<T>::value,
    "A knot vector must be defined using an arithmetic type");

public:
  /// \brief Default constructor
  KnotVector() : m_deg(-1) { m_knots.resize(0); }

  /*!
   * \brief Constructor for a normalized, uniform knot vector
   *
   * \param [in] npts the number of control points
   * \param [in] deg the degree
   * \pre degree is greater than or equal to -1, and deg < npts.
   */
  KnotVector(axom::IndexType npts, int deg) { makeUniform(npts, deg); }

  /// \brief Constructor from a user-supplied knot vector (C-style)
  KnotVector(const T* knots, axom::IndexType nkts, int degree)
  {
    m_knots.resize(nkts);
    for(int i = 0; i < nkts; ++i)
    {
      m_knots[i] = knots[i];
    }

    m_deg = degree;
  }

  /// \brief Constructor from a user-supplied knot vector (axom::Array)
  KnotVector(const axom::Array<T>& knots, int degree)
    : m_knots(knots)
    , m_deg(degree)
  { }

  /*!
   * \brief Make the knot vector uniform
   */
  void makeUniform(axom::IndexType npts, int deg)
  {
    // Set the degree
    m_deg = deg;

    // Initialize the knot vector
    m_knots.resize(npts + deg + 1);

    // Knots for clamped-ness
    for(int i = 0; i < deg + 1; ++i)
    {
      m_knots[i] = 0.0;
      m_knots[npts + deg - i] = 1.0;
    }

    // Interior knots (if any)
    for(int i = 0; i < npts - deg - 1; ++i)
    {
      m_knots[deg + 1 + i] = (i + 1.0) / static_cast<T>(npts - deg);
    }
  }

  /// \brief Accessor for the knot vector
  const axom::Array<T>& getArray() const { return m_knots; }

  /// \brief Getter for the knot vector
  const T& operator[](int i) const { return m_knots[i]; }

  /// \brief Setter for the knot vector
  T& operator[](int i) { return m_knots[i]; }

  /// \brief Normalize the knot vector to the span of [0, 1]
  void normalize()
  {
    T min_knot = m_knots[0];
    T max_knot = m_knots[m_knots.size() - 1];
    T span = max_knot - min_knot;

    for(int i = 0; i < m_knots.size(); ++i)
    {
      m_knots[i] = (m_knots[i] - min_knot) / span;
    }
  }

  /// \brief Insert a knot into the knot vector r times
  void insertKnot(axom::IndexType span, T t, int r = 1)
  {
    for(int i = 0; i < r; ++i)
    {
      m_knots.insert(m_knots.begin() + span + 1, t);
    }
  }

  /*!
   * \brief Split a knot vector at the value at index `span`
   * 
   * \pre Assumes that the multiplicity of the specific knot is equal
   *    to the degree of the knot vector
   */
  void split(axom::IndexType span, KnotVector& k1, KnotVector& k2) const
  {
    const auto nkts = getNumKnots();

    // Copy the degree
    k1.m_deg = m_deg;
    k2.m_deg = m_deg;

    // Copy the knots
    k1.m_knots.resize(span + 2);
    k2.m_knots.resize(nkts - span + m_deg);

    k1.m_knots[span + 1] = m_knots[span];
    for(int i = 0; i < span + 1; ++i)
    {
      k1.m_knots[i] = m_knots[i];
    }

    k2.m_knots[0] = m_knots[span];
    for(int i = 0; i < nkts - span + m_deg - 1; ++i)
    {
      k2.m_knots[i + 1] = m_knots[span + i - m_deg + 1];
    }
  }

  /*!
   * \brief Evaluates the NURBS basis functions for span at parameter value t
   * 
   * Implementation adapted from Algorithm A2.2 on page 70 of "The NURBS Book".
   */
  axom::Array<T> calculateBasisFunctions(axom::IndexType span, T t) const
  {
    axom::Array<T> N(m_deg + 1);
    axom::Array<T> left(m_deg + 1);
    axom::Array<T> right(m_deg + 1);

    // Note: This implementation avoids division by zero and redundant computation
    // that might arise from a direct implementation of the recurrence relation
    // for basis functions. See "The NURBS Book" for details.
    N[0] = 1.0;
    for(int j = 1; j <= m_deg; ++j)
    {
      left[j] = t - m_knots[span + 1 - j];
      right[j] = m_knots[span + j] - t;
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
    const int m_deg = getDegree();

    axom::Array<axom::Array<T>> ndu(m_deg + 1), a(2);
    axom::Array<T> left(m_deg + 1), right(m_deg + 1);
    for(int j = 0; j <= m_deg; j++)
    {
      ndu[j].resize(m_deg + 1);
    }
    for(int j = 0; j <= n; j++)
    {
      ders[j].resize(m_deg + 1);
    }
    a[0].resize(m_deg + 1);
    a[1].resize(m_deg + 1);

    ndu[0][0] = 1.;
    for(int j = 1; j <= m_deg; j++)
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
    for(int j = 0; j <= m_deg; j++)
    {
      ders[0][j] = ndu[j][m_deg];
    }

    // This section computes the derivatives (Eq. [2.9])

    // Loop over function index.
    for(int r = 0; r <= m_deg; r++)
    {
      int s1 = 0, s2 = 1;  // Alternate rows in array a
      a[0][0] = 1.;
      // Loop to compute kth derivative
      for(int k = 1; k <= n; k++)
      {
        T d = 0.;
        int rk = r - k;
        int pk = m_deg - k;
        if(r >= k)
        {
          a[s2][0] = a[s1][0] / ndu[pk + 1][rk];
          d = a[s2][0] * ndu[rk][pk];
        }
        int j1 = (rk >= -1) ? 1 : -rk;
        int j2 = (r - 1 <= pk) ? (k - 1) : (m_deg - r);
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
    T r = static_cast<T>(m_deg);
    for(int k = 1; k <= n; k++)
    {
      for(int j = 0; j <= m_deg; j++)
      {
        ders[k][j] *= r;
      }
      r *= static_cast<T>(m_deg - k);
    }
  }

  /// \brief Return the degree of the knot vector
  int getDegree() const { return m_deg; }

  /// \brief Return the number of knots in the knot vector
  axom::IndexType getNumKnots() const { return m_knots.size(); }

  /// \brief Return the number of control points implied by the knot vector
  axom::IndexType getNumControlPoints() const
  {
    return m_knots.size() - m_deg - 1;
  }

  /// \brief Clear the list of knots
  void clear()
  {
    m_deg = -1;
    m_knots.clear();
  }

  /*!
   * \brief Returns the index of the knot span containing parameter t
   * 
   * Implementation adapted from Algorithm A2.1 on page 68 of "The NURBS Book"
   */
  axom::IndexType findSpan(double t) const
  {
    const axom::IndexType nkts = m_knots.size();

    if(t <= m_knots[0])
    {
      return m_deg;
    }

    if(t >= m_knots[nkts - 1])
    {
      return nkts - m_deg - 2;
    }

    // perform binary search on the knots,
    axom::IndexType low = m_deg;
    axom::IndexType high = nkts - m_deg - 1;
    axom::IndexType mid = (low + high) / 2;
    while(t < m_knots[mid] || t >= m_knots[mid + 1])
    {
      (t < m_knots[mid]) ? high = mid : low = mid;
      mid = (low + high) / 2;
    }

    return mid;
  }

  /// \brief Returns the index of the knot span and the multiplicity
  axom::IndexType findSpan(double t, int& multiplicity)
  {
    const auto nkts = m_knots.size();
    const auto span = findSpan(t);

    multiplicity = 0;
    for(auto i = (t == m_knots[nkts - 1]) ? nkts - 1 : span; i >= 0; --i)
    {
      if(m_knots[i] == t)
      {
        multiplicity++;
      }
      else
      {
        break;
      }
    }

    return span;
  }

  /// \brief Reverse the knot vector
  void reverse()
  {
    const axom::IndexType nkts = m_knots.size();
    const axom::IndexType knot_mid = (nkts + 1) / 2;

    // Flip the vector
    for(int i = 0; i < knot_mid; ++i)
    {
      axom::utilities::swap(m_knots[i], m_knots[nkts - 1 - i]);
    }

    // Replace each knot with sum - knot_value
    const T the_sum = m_knots[0] + m_knots[nkts - 1];
    for(int i = 0; i < nkts; ++i)
    {
      m_knots[i] = the_sum - m_knots[i];
    }
  }

  /// \brief Return if the knot vector is valid
  bool isValid() const
  {
    // Check degree
    if(m_deg < 0)
    {
      return false;
    }

    // Check for monotonicity
    for(int i = 0; i < m_knots.size() - 1; ++i)
    {
      if(m_knots[i] > m_knots[i + 1])
      {
        return false;
      }
    }

    // Check for clamped-ness
    auto nkts = m_knots.size();
    for(int i = 0; i < m_deg + 1; ++i)
    {
      if(m_knots[i] != m_knots[0])
      {
        return false;
      }
      if(m_knots[nkts - 1 - i] != m_knots[nkts - 1])
      {
        return false;
      }
    }

    // Check for continuity
    T this_knot = m_knots[m_deg];
    int this_multiplicity = 1;
    for(int i = m_deg + 1; i < nkts - m_deg - 1; ++i)
    {
      if(m_knots[i] != this_knot)
      {
        this_knot = m_knots[i];
        this_multiplicity = 1;
      }
      else
      {
        this_multiplicity++;
        if(this_multiplicity > m_deg)
        {
          return false;
        }
      }
    }

    return true;
  }

  /// \brief Check equality of two knot vectors
  friend inline bool operator==(const KnotVector<T>& lhs, const KnotVector<T>& rhs)
  {
    return lhs.m_deg == rhs.m_deg && lhs.m_knots == rhs.m_knots;
  }

  /// \brief Check inequality of two knot vectors
  friend inline bool operator!=(const KnotVector<T>& lhs, const KnotVector<T>& rhs)
  {
    return !(lhs == rhs);
  }

  /*!
   * \brief Simple formatted print of a knot vector instance
   *
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const
  {
    int nkts = m_knots.size();

    os << "{ degree " << m_deg << " knot vector ";
    for(int i = 0; i < nkts; ++i)
    {
      os << m_knots[i] << (i < nkts - 1 ? ", " : "");
    }
    os << "}";

    return os;
  }

private:
  int m_deg;
  axom::Array<T> m_knots;
};

//------------------------------------------------------------------------------
/// Free functions related to KnotVector
//------------------------------------------------------------------------------
template <typename T>
std::ostream& operator<<(std::ostream& os, const KnotVector<T>& kvector)
{
  kvector.print(os);
  return os;
}

}  // namespace primal
}  // namespace axom

/// Overload to format a primal::NURBSCurve using fmt
template <typename T>
struct axom::fmt::formatter<axom::primal::KnotVector<T>> : ostream_formatter
{ };

#endif  // AXOM_PRIMAL_KNOTVECTOR_HPP