// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_PRIMAL_PLANE_HPP_
#define AXOM_PRIMAL_PLANE_HPP_

#include "axom/core/Macros.hpp"                        // for Axom macros
#include "axom/core/numerics/matvecops.hpp"            // for vector operators
#include "axom/primal/geometry/OrientationResult.hpp"  // for OrientedSide
#include "axom/primal/geometry/Vector.hpp"             // for primal::Vector
#include "axom/primal/geometry/Point.hpp"              // for primal::Point

#include "axom/slic/interface/slic.hpp"  // for SLIC macros

namespace axom
{
namespace primal
{
/// \name Forward Declared Overloaded Operators
/// @{

template <typename T, int NDIMS>
class Plane;

template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Plane<T, NDIMS>& p);

/// @}

/*!
 * \class Plane
 *
 * \brief Defines an oriented plane within a 2-D or 3-D Euclidean space
 *  and provides associated operators, such as, projection, signed distance,
 *  orientation, etc.
 *
 *  The Plane object defines an oriented plane in Hessian Normal form:
 *
 *           \f$ \mathcal{N} \cdot x - d = 0 \f$
 *
 *  where, \f$ \mathcal{N} \f$ is a unit normal and \f$ d \f$ is the offset of
 *  the plane (in the direction of the specified normal) to the origin.
 *
 *  The Plane class defines a co-dimension-one plane that splits the ambient
 *  space into two halfspaces, such that, the signed distance for any point
 *  is: (a) negative below the plance, (b) positive above the plane and (c)
 *  zero on the plane.
 *
 *  A Plane object may be constructed in three ways:
 *  1. Specifying a normal and a point that the plane passes through
 *  2. Specifying a normal and an offset from the origin, or
 *  3. Specifying the points that the plane passes through. A minimum of two
 *     points are required to define a plane in 2D, or, three points for 3D.
 *
 * \tparam T the underlying data type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 *
 * \pre NDIMS==2 || NDIMS==3
 */
template <typename T, int NDIMS>
class Plane
{
public:
  using VectorType = primal::Vector<T, NDIMS>;
  using PointType = primal::Point<T, NDIMS>;

public:
  /// \name Constructors
  /// @{

  /*!
   * \brief Constructs a Plane with a given normal, \f$ \mathcal{N} \f$, that
   *  passes through the specified point, \f$  x \f$
   *
   * \param [in] normal the supplied plane normal, \f$ \mathcal{N} \f$
   * \param [in] x a specified point that the plane passes through, \f$ x \f$
   *
   * \note The supplied normal will be normalized, such that, the Plane's normal
   *  will always be a unit normal.
   */
  Plane(const VectorType& normal, const PointType& x);

  /*!
   * \brief Constructs a plane with a specified normal, \f$ \mathcal{N} \f$,
   *  located at the given offset from origin.
   *
   * \param [in] normal the supplied plane normal, \f$ \mathcal{N} \f$
   * \param [in] offset the plane offset from origin
   *
   * \note The supplied normal will be normalized, such that, the Plane's normal
   *  will always be a unit normal.
   */
  Plane(const VectorType& normal, T offset);

  /*!
   * \brief Constructs a Plane that goes through the specified points.
   *
   * \param [in] x1 coordinates of the first point.
   * \param [in] x2 coordinates of the second point.
   * \param [in] x3 coordinates of the third point, or nullptr for 2D.
   *
   * \note Two points are required to define a plane in 2D, in which case, the
   *  caller must pass nullptr as the third argument to this constructor.
   *
   * \note The specified points, x1, x2, x3, must point to buffers that are
   *  at least NDIMS long
   *
   * \pre x1 != nullptr
   * \pre x2 != nullptr
   * \pre x3 != nullptr \f$ \iff \f$ NDIMS==3
   * \pre In 2D, x1 should not be equal to x2
   * \pre In 3D, the user-supplied points, x1, x2, x3 should not be collinear.
   */
  AXOM_HOST_DEVICE Plane(const T* x1, const T* x2, const T* x3);

  /// @}

  /*!
   * \brief Returns the dimension of the Plane
   * \return dim the dimension of the Plane.
   * \post (dim==2) || (dim==3)
   */
  inline int getDimension() const { return NDIMS; };

  /*!
   * \brief Returns a const pointer to the plane's unit normal.
   * \return N the unit normal.
   * \post N != nullptr
   */
  AXOM_HOST_DEVICE const VectorType& getNormal() const { return m_normal; }

  /*!
   * \brief Returns the offset of the plane from origin.
   * \return offSet the offset of the plane from origin.
   */
  AXOM_HOST_DEVICE T getOffset() const { return m_offset; }

  /*!
   * \brief Computes the signed distance of the specified point to this Plane.
   *
   * \param [in] x buffer consisting of the query point coordinates.
   * \return d the signed distance of the point x to the plane.
   */
  AXOM_HOST_DEVICE T signedDistance(const PointType& x) const
  {
    return m_normal.dot(VectorType(x.array())) - m_offset;
  }

  /*!
   * \brief Computes the projection of a given point, x, onto this Plane.
   *
   * \param [in] x buffer consisting of the coordinates of the point to project.
   * \return projx the coordinates of the projected point.
   *
   * \post this->getOrientedSide( projx ) == ON_BOUNDARY
   */
  PointType projectPoint(const PointType& x) const;

  /*!
   * \brief Flips the orientation of the plane.
   */
  inline void flip();

  /*!
   * \brief Computes the orientation of the point with respect to this Plane.
   *
   * \param [in] x buffer consisting of the query point coordinates.
   * \param [in] TOL user-supplied tolerance. Optional. Default is 1.e-9.
   * \return orientation the orientation of the point
   *
   * \note This method returns one of the following values:
   *  <ul>
   *   <li> <b>ON_BOUNDARY</b>      : if `x` lies on the plane  </li>
   *   <li> <b>ON_POSITIVE_SIDE</b> : if `x` is above the plane </li>
   *   <li> <b>ON_NEGATIVE_SIDE</b> : if `x` is below the plane </li>
   *  </ul>
   *
   * \see OrientationResult
   */
  AXOM_HOST_DEVICE inline int getOrientation(const PointType& x,
                                             double TOL = 1.e-9) const;

  /*!
   * \brief Prints the Plane information in the given output stream.
   * \param [in] os The output stream to write to.
   * \note This method is primarily used for debugging.
   * \return s the modified output stream object.
   */
  std::ostream& print(std::ostream& os) const;

private:
  /*!
   * \brief Sets the normal of this Plane instance.
   *
   * \param [in] normal pointer to a buffer consisting of the normal
   * \note The supplied buffer must be at least NDIMS long.
   * \pre normal != nullptr
   */
  AXOM_HOST_DEVICE inline void setNormal(const VectorType& normal);

  VectorType m_normal; /*!< plane unit-normal  */
  T m_offset;          /*!< offset from origin */
};

} /* namespace primal */
} /* namespace axom */

//------------------------------------------------------------------------------
//  IMPLEMENTATION
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{
template <typename T, int NDIMS>
Plane<T, NDIMS>::Plane(const VectorType& normal, const PointType& x)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A plane object may be defined in 2-D or 3-D");
  SLIC_ASSERT_MSG(!normal.is_zero(),
                  "Normal vector of a plane should be non-zero");

  this->setNormal(normal);

  m_offset = m_normal.dot(VectorType(x.array()));
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Plane<T, NDIMS>::Plane(const VectorType& normal, T offset) : m_offset(offset)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A plane object may be defined in 2-D or 3-D");
  SLIC_ASSERT_MSG(!normal.is_zero(),
                  "Normal vector of a plane should be non-zero");

  this->setNormal(normal);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Plane<T, NDIMS>::Plane(const T* x1, const T* x2, const T* x3)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A plane object may be defined in 2-D or 3-D");

  SLIC_ASSERT(x1 != nullptr);
  SLIC_ASSERT(x2 != nullptr);

  VectorType normal;

  if(x3 == nullptr)
  {
    SLIC_ASSERT(NDIMS == 2);
    normal[0] = x1[1] - x2[1];  // -dy
    normal[1] = x2[0] - x1[0];  //  dx

  }  // END if
  else
  {
    SLIC_ASSERT(NDIMS == 3);

    VectorType r1(x1, x2);
    VectorType r2(x1, x3);

    normal = VectorType(VectorType::cross_product(r1, r2).data(), NDIMS);

  }  // END else

  // check for degenerate line or triangle
  bool degenerate = normal.is_zero();

  SLIC_CHECK_MSG(!degenerate,
                 "Supplied points form a degenerate "
                   << ((NDIMS == 2) ? "line" : "triangle"));

  this->setNormal(normal);
  m_offset = m_normal.dot(x1);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline typename Plane<T, NDIMS>::PointType Plane<T, NDIMS>::projectPoint(
  const PointType& x) const
{
  const T signed_distance = this->signedDistance(x);
  PointType projx =
    PointType((VectorType(x) - (m_normal * signed_distance)).array());

  return projx;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline void Plane<T, NDIMS>::flip()
{
  m_normal *= static_cast<T>(-1.0);
  m_offset *= static_cast<T>(-1.0);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline int Plane<T, NDIMS>::getOrientation(const PointType& x, double TOL) const
{
  int oriented_side = -1;

  const T signed_distance = this->signedDistance(x);

  if(utilities::isNearlyEqual(signed_distance, 0.0, TOL))
  {
    // on the plane
    oriented_side = ON_BOUNDARY;
  }
  else if(signed_distance < 0.0)
  {
    // below the plane
    oriented_side = ON_NEGATIVE_SIDE;
  }
  else
  {
    // above  the plane
    oriented_side = ON_POSITIVE_SIDE;
  }

  return oriented_side;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& Plane<T, NDIMS>::print(std::ostream& os) const
{
  os << "unit normal: [ ";
  for(int i = 0; i < NDIMS; ++i)
  {
    os << m_normal[i] << " ";
  }
  os << "] offset: " << m_offset;
  return (os);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline void Plane<T, NDIMS>::setNormal(const VectorType& normal)
{
  m_normal = normal.unitVector();
}

//------------------------------------------------------------------------------
//  implementation of free functions
//------------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Plane<T, NDIMS>& p)
{
  return (p.print(os));
}

}  // namespace primal
}  // namespace axom

#endif  // AXOM_PRIMAL_PLANE_HPP_
