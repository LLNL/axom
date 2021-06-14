// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_SPHERE_HPP_
#define AXOM_PRIMAL_SPHERE_HPP_

#include "axom/core/Macros.hpp"  // for AXOM macros

#include "axom/core/utilities/Utilities.hpp"  // utilities::isNearlyEqual()
#include "axom/core/numerics/matvecops.hpp"   // for dot_product()

#include "axom/primal/geometry/OrientationResult.hpp"  // OrientationResult

#include "axom/slic/interface/slic.hpp"  // for SLIC macros

namespace axom
{
namespace primal
{
/// \name Forward Declared Overloaded Operators
/// @{

template <typename T, int NDIMS>
class Sphere;

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Sphere<T, NDIMS>& s);

/// @}

/*!
 * \class Sphere
 *
 * \brief Defines an oriented Sphere in 2-D (i.e., a circle) or 3-D given by
 *  its center, \f$ \mathcal{X} \f$ and radius \f$ \mathcal{R} \f$. The Sphere
 *  object provides associated operations on a sphere, such as, signed distance
 *  and orientation.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Sphere
{
public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Constructs a Sphere centered at origin with the given radius.
   * \param [in] radius the radius of the Sphere (optional).
   * \note If a radius is not supplied, the default radius is 1.0.
   */
  explicit Sphere(T radius = 1.0);

  /*!
   * \brief Constructs a Sphere with the given center and radius.
   *
   * \param [in] center user-supplied center.
   * \param [in] radius the radius of the Sphere (optional).
   *
   * \note If a radius is not supplied, the default radius is 1.0.
   *
   * \pre center != nullptr
   */
  explicit Sphere(const T* center, T radius = 1.0);

  /// @}

  /*!
   * \brief Destructor.
   */
  ~Sphere();

  /*!
   * \brief Returns the radius of the Sphere.
   * \return r the radius of the Sphere.
   */
  inline T getRadius() const { return m_radius; };

  /*!
   * \brief Returns the center of the Sphere.
   *
   * \return c pointer to array that holds the center of the Sphere.
   * \note c points to an array that is NDIMS long.
   * \post c != nullptr
   */
  inline const T* getCenter() const { return m_center; };

  /*!
   * \brief Computes the signed distance of a point to the Sphere's boundary.
   *
   * \param [in] q pointer to buffer consisting of query point coordinates.
   * \return d the computed signed distance of the point, q, to the sphere.
   *
   * \note The signed distance of a point, q, is:
   *  <ul>
   *   <li> negative inside the sphere </li>
   *   <li> positive outside the sphere </li>
   *   <li> zero on the boundary </li>
   *  </ul>
   *
   * \pre q != nullptr
   * \pre q must be a pointer to an array that is at least NDIMS long
   */
  inline T computeSignedDistance(const T* q) const;

  /*!
   * \brief Computes the orientation of a point with respect to the Sphere.
   *
   * \param [in] q pointer to user-supplied point q.
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orient the orientation of q with respect to the sphere.
   *
   *  \note This method returns one of the following values:
   *   <ul>
   *    <li> <b>ON_BOUNDARY</b>      : if `q` is on the sphere's boundary </li>
   *    <li> <b>ON_POSITIVE_SIDE</b> : if `q` is outside the sphere </li>
   *    <li> <b>ON_NEGATIVE_SIDE</b> : if `q` is inside the sphere </li>
   *  </ul>
   *
   * \see OrientedSide for the list of possible return values.
   *
   * \pre q != nullptr
   * \pre q must be a pointer to an array that is at least NDIMS long
   *
   */
  inline int getOrientation(const T* q, double TOL = 1.e-9) const;

  /*!
   * \brief Tests if this sphere instance intersects with another sphere.
   *
   * \param [in] sphere the sphere object to check for intersection
   * \param [in] TOL tolerance for intersection test. Optional. If not specified
   *  the default tolerance is set to 1.e-9.
   *
   * \return status true if the sphere intersects, false otherwise.
   */
  inline bool intersectsWith(const Sphere<T, NDIMS>& sphere,
                             double TOL = 1.e-9) const;

  /*!
   * \brief Prints the Sphere information in the given output stream.
   * \param [in,out] os the output stream to write to.
   * \note This method is primarily used for debugging.
   * \return s the modified output stream object.
   */
  std::ostream& print(std::ostream& os) const;

private:
  T m_center[NDIMS]; /*!< sphere center */
  T m_radius;        /*!< sphere radius */
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
// Sphere Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Sphere<T, NDIMS>::Sphere(T radius) : m_radius(radius)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A Sphere object may be defined in 2-D or 3-D");

  for(int i = 0; i < NDIMS; ++i)
  {
    m_center[i] = static_cast<T>(0.0);
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Sphere<T, NDIMS>::Sphere(const T* center, T radius) : m_radius(radius)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A Sphere object may be defined in 2-D or 3-D");

  SLIC_ASSERT(center != nullptr);
  for(int i = 0; i < NDIMS; ++i)
  {
    m_center[i] = center[i];
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Sphere<T, NDIMS>::~Sphere()
{ }

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline T Sphere<T, NDIMS>::computeSignedDistance(const T* q) const
{
  SLIC_ASSERT(q != nullptr);

  T d = 0.0;
  for(int i = 0; i < NDIMS; ++i)
  {
    const T dx = q[i] - m_center[i];
    d += (dx * dx);
  }

  return (std::sqrt(d) - m_radius);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline int Sphere<T, NDIMS>::getOrientation(const T* q, double TOL) const
{
  SLIC_ASSERT(q != nullptr);

  T signed_distance = this->computeSignedDistance(q);

  int orient = -1;

  if(axom::utilities::isNearlyEqual(signed_distance, 0.0, TOL))
  {
    orient = ON_BOUNDARY;
  }
  else if(signed_distance < 0.0f)
  {
    // inside
    orient = ON_NEGATIVE_SIDE;
  }
  else
  {
    // outside
    orient = ON_POSITIVE_SIDE;
  }

  return orient;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Sphere<T, NDIMS>::print(std::ostream& os) const
{
  os << " center: [ ";
  for(int i = 0; i < NDIMS; ++i)
  {
    os << m_center[i] << " ";
  }
  os << "] radius: " << m_radius;

  return os;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline bool Sphere<T, NDIMS>::intersectsWith(const Sphere<T, NDIMS>& sphere,
                                             double TOL) const
{
  double r[NDIMS];
  for(int i = 0; i < NDIMS; ++i)
  {
    r[i] = m_center[i] - sphere.getCenter()[i];
  }

  const T distance_squared = numerics::dot_product(r, r, NDIMS);
  const T sum_of_radii = m_radius + sphere.getRadius();
  const T sum_of_radii_2 = sum_of_radii * sum_of_radii;

  return (distance_squared < sum_of_radii_2 ||
          utilities::isNearlyEqual(distance_squared, sum_of_radii_2, TOL));
}

//------------------------------------------------------------------------------
//  implementation of free functions
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Sphere<T, NDIMS>& s)
{
  return (s.print(os));
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_SPHERE_HPP_
