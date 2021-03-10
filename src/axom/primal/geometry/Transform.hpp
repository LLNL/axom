// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef PRIMAL_TRANSFORM_HPP_
#define PRIMAL_TRANSFORM_HPP_

#include "NumericArray.hpp"
#include "Octahedron.hpp"
#include "Point.hpp"
#include "Vector.hpp"

#include "axom/core/Macros.hpp"                        // for Axom macros
#include "axom/core/numerics/Matrix.hpp"               // for Matrix class
#include "axom/core/numerics/matvecops.hpp"            // for vector operators

namespace axom
{
namespace primal
{

enum RotationAxis
{
   X = 0,
   Y,
   Z
};

/*!
 * \class Transform
 *
 * \brief The Transform class represents the operation of transforming
 * another geometric primitive.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template <typename T, int NDIMS>
class Transform
{
public:
/// \name Create transforms
///@{
   /// Default constructor: identity transform
   Transform();
   /// Create a transform from a matrix
   Transform(const numerics::Matrix<T> & m);
   /// Create a translation from a vector
   Transform(const Vector<T, NDIMS> & v);
   /// Create a rotation from the X-, Y-, or Z-axis and an angle
   Transform(int axis, double theta);
   /// Map an axis to the canonical X axis
   static Transform<T, NDIMS> RotateToXAxisFrom(const Vector<T, NDIMS> & u);
///@}


/// \name Compose transforms
///@{
   Transform<T, NDIMS> operator*(const Transform<T, NDIMS>& rh) const;
///@}


/// \name Apply a transform to a geometric primitive
///@{
   NumericArray<T, NDIMS> operator*(const NumericArray<T, NDIMS>& rh) const;
   Point<T, NDIMS> operator*(const Point<T, NDIMS>& rh) const;
   Octahedron<T, NDIMS> operator*(const Octahedron<T, NDIMS>& rh) const;
///@}

   /// Retrieve this Transform's matrix
   const numerics::Matrix<T> & getTransform() const { return m_transform; }

private:
   /// The current transform.
   numerics::Matrix<T> m_transform;
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
Transform<T, NDIMS>::Transform() : m_transform(NDIMS+1, NDIMS+1)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A transform object may be defined in 2-D or 3-D");

  m_transform = numerics::Matrix<T>::identity(NDIMS+1);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Transform<T, NDIMS>::Transform(const numerics::Matrix<T> & m) : m_transform(m)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A transform object may be defined in 2-D or 3-D");
  SLIC_ASSERT_MSG((m.isSquare() && m.getNumColumns() == NDIMS),
                  "A homogeneous transform object must be defined by a square "
                  "3x3 or 4x4 matrix");
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Transform<T, NDIMS>::Transform(const Vector<T, NDIMS> & v) :
   m_transform(NDIMS+1, NDIMS+1)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A transform object may be defined in 2-D or 3-D");

  m_transform = numerics::Matrix<T>::identity(NDIMS+1);

  // Get a reference to the last column
  T* lastcol = m_transform.getColumn(NDIMS-1);
  // Set v into the last column
  for (int row = 0; row < NDIMS; ++row)
  {
     lastcol[row] = v[row];
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Transform<T, NDIMS>::Transform(int axis, double theta) :
   m_transform(NDIMS+1, NDIMS+1)
{
  AXOM_STATIC_ASSERT_MSG((NDIMS == 2) || (NDIMS == 3),
                         "A transform object may be defined in 2-D or 3-D");

  m_transform = numerics::Matrix<T>::identity(NDIMS+1);

  double s = sin(theta);
  double c = cos(theta);
  
  if (NDIMS == 2)
  {
     m_transform(0, 0) = c;
     m_transform(0, 1) = s;
     m_transform(1, 0) = -s;
     m_transform(1, 1) = c;
  }
  else
  {
     int c1 = (axis + 1) % 3;
     int c2 = (axis + 2) % 3;
     m_transform(c1, c1) = c;
     m_transform(c1, c2) = s;
     m_transform(c2, c1) = -s;
     m_transform(c2, c2) = c;
  }
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Transform<T, NDIMS>
Transform<T, NDIMS>::RotateToXAxisFrom(const Vector<T, NDIMS> & u)
{
   using VectorType = Vector<T, NDIMS>;
   // The constraints for this transform are that it must map v to the X-axis,
   // and it must preserve length.
   numerics::Matrix<T> t = numerics::Matrix<T>::identity(NDIMS+1);

   VectorType uhat = u.unitVector();

   VectorType yhat{0, 1, 0};
   double udoty = uhat.dot(yhat);
   VectorType zhat{0, 0, 1};
   double udotz = uhat.dot(zhat);

   VectorType m;
   double udotm;
   if (udoty < udotz)
   {
      m = yhat;
      udotm = udoty;
   }
   else
   {
      m = zhat;
      udotm = udotz;
   }
   VectorType vhat = (m - udotm*uhat).unitVector();
   VectorType what = VectorType::cross_product(uhat, vhat);

   for (int d = 0; d < 3; ++d)
   {
      t(d, 0) = uhat[d];
      t(d, 1) = vhat[d];
      t(d, 2) = what[d];
   }

   numerics::Matrix<T> tt(NDIMS+1, NDIMS+1);
   matrix_transpose(t, tt);
   
   return Transform<T, NDIMS>(tt);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Transform<T, NDIMS>
Transform<T, NDIMS>::operator*(const Transform<T, NDIMS>& rh) const
{
   numerics::Matrix<T> result(NDIMS+1, NDIMS+1);
   numerics::matrix_multiply(getTransform(), rh.getTransform(), result);
   return Transform<T, NDIMS>(result);
}


//------------------------------------------------------------------------------
template <typename T, int NDIMS>
NumericArray<T, NDIMS>
Transform<T, NDIMS>::operator*(const NumericArray<T, NDIMS>& rh) const
{
   NumericArray<T, NDIMS+1> tx, tb;
   for (int d = 0; d < NDIMS; ++d) { tx[d] = rh[d]; }

   numerics::matrix_vector_multiply(getTransform(), tx.data(), tb.data());

   NumericArray<T, NDIMS> result;
   for (int d = 0; d < NDIMS; ++d) { result[d] = tb[d]; }
   return result;
}


//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Point<T, NDIMS>
Transform<T, NDIMS>::operator*(const Point<T, NDIMS>& rh) const
{
   return Point<T, NDIMS>((*this)*(rh.array()));
}


//------------------------------------------------------------------------------
template <typename T, int NDIMS>
Octahedron<T, NDIMS>
Transform<T, NDIMS>::operator*(const Octahedron<T, NDIMS>& rh) const
{
   Octahedron<T, NDIMS> out;
   for (int p = 0; p < Octahedron<T, NDIMS>::NUM_OCT_VERTS; ++p)
   {
      out[p] = (*this)*rh[p];
   }
   return out;
}

//-----------------------------------------------------------------------------
template <typename T, int NDIMS>
std::ostream& operator<<(std::ostream& os, const Transform<T, NDIMS>& t)
{
  const numerics::Matrix<T>& tmat = t.getTransform();
  const int nrows = tmat.getNumRows();
  const int ncols = tmat.getNumColumns();

  os << "Homogeneous " << nrows << "D transform:\n";
  
  for(typename numerics::Matrix<T>::IndexType i = 0; i < nrows; ++i)
  {
    os << "[ ";
    for(int j = 0; j < ncols; ++j)
    {
      os << tmat(i, j) << " ";
    }  // END for all j

    os << "]\n";

  }  // END for all i

  return (os);
}

} /* namespace primal */
} /* namespace axom */

#endif /* PRIMAL_TRANSFORM_HPP_ */
