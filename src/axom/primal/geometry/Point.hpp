// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef POINT_HXX_
#define POINT_HXX_

#include "axom/core/Macros.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/primal/geometry/NumericArray.hpp"

// C/C++ includes
#include <cstring> // For memcpy()
#include <ostream> // For print() and operator <<

namespace axom
{
namespace primal
{

// Forward declare the templated classes and operator functions
template < typename T,int NDIMS >
class Point;

/// \name Forward Declared Overloaded Operators
///@{

/*!
 * \brief Equality comparison operator for points
 */
template < typename T,int NDIMS >
bool operator==(const Point< T, NDIMS > & lhs, const Point< T, NDIMS >& rhs);

/*!
 * \brief Inequality comparison operator for points
 */
template < typename T,int NDIMS >
bool operator!=(const Point< T, NDIMS > & lhs, const Point< T, NDIMS >& rhs);

/*!
 * \brief Overloaded output operator for points
 */
template < typename T,int NDIMS >
std::ostream& operator<<(std::ostream & os, const Point< T,NDIMS > & pt);

///@}

/*!
 * \class Point
 *
 * \brief The point class represents a point, \f$ p \in \mathcal{R}^d \f$ . It
 *  provides access methods to set and query the point coordinates.
 *
 * \tparam T the coordinate type, e.g., double, float, etc.
 * \tparam NDIMS the number of dimensions
 */
template < typename T, int NDIMS >
class Point
{
public:
  enum
  {
    DIMENSION = NDIMS,
    NBYTES    = NDIMS * sizeof(T)
  };

  typedef Point< T,NDIMS > PointType;
  typedef T CoordType;

public:

  /*!
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to.  Defaults to zero
   * \param [in] sz The number of coordinates to set to val.
   * The rest will be set to zero.  Defaults is NDIMS.
   * If sz is greater than NDIMS, we set all coordinates to val
   */
  AXOM_HOST_DEVICE
  explicit Point(T val = T(), int sz = NDIMS) : m_components(val,sz) { }

  /*!
   * \brief Constructor from a numeric array
   * \param [in] arr The numeric array to copy from
   */
  AXOM_HOST_DEVICE
  Point(const NumericArray< T,NDIMS >& arr) : m_components(arr) { }

  /*!
   * \brief Creates a point from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz num values to copy from the vals array. Defaults to NDIMS.
   * \note If sz is greater than NDIMS, we only take the first NDIMS values.
   */
  AXOM_HOST_DEVICE
  Point(const T* vals, int sz = NDIMS) : m_components(vals,sz) { }

  /*!
   * \brief Copy constructor.
   * \param [in] other The point to copy
   */
  AXOM_HOST_DEVICE
  Point( const Point& other) : m_components( other.m_components) { }

  /*!
   * \brief Destructor.
   */
  AXOM_HOST_DEVICE
  ~Point() { }

  /*!
   * \brief Returns the dimension of this point instance.
   * \return d the dimension of the point.
   * \post d >= 1.
   */
  static int dimension() { return NDIMS; };

  /*!
   * \brief Assignment operator.
   * \param [in] rhs a point instance on the right hand side.
   */
  AXOM_HOST_DEVICE
  Point& operator=(const Point& rhs)
  { m_components=rhs.m_components; return *this; }

  /// \name Overloaded [] operator methods
  ///@{

  /*!
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   */
  AXOM_HOST_DEVICE
  const T& operator[](int i) const { return m_components[i]; }

  AXOM_HOST_DEVICE
  T& operator[](int i)             { return m_components[i]; }

  ///@}

  /// \name Raw data member access methods
  ///@{

  /*!
   * \brief Returns a pointer to the underlying data.
   */
  const T* data() const { return m_components.data(); }
  T* data()             { return m_components.data(); }

  ///@}

  /*!
   * \brief Returns a reference to the underlying NumericArray.
   */
  AXOM_HOST_DEVICE
  const NumericArray< T,NDIMS >& array() const { return m_components; }

  AXOM_HOST_DEVICE
  NumericArray< T,NDIMS >& array()             { return m_components; }

  /*!
   * \brief Output the point's coordinates to the array
   * \param arr The array that we are outputting to.
   * \pre The user needs to make sure that the array has been allocated
   * and has sufficient space for NDIMS coordinates.
   */
  void to_array(T* arr) const { m_components.to_array(arr); }

  /*!
   * \brief Equality comparison operator for points
   */
  friend inline bool operator==(const Point& lhs, const Point& rhs)
  { return lhs.m_components == rhs.m_components; }

  /*!
   * \brief Inequality operator for points
   */
  friend inline bool operator!=(const Point& lhs, const Point& rhs)
  { return !(lhs == rhs); }

  /*!
   * \brief Simple formatted print of a point instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   */
  std::ostream& print(std::ostream& os) const;

  /*!
   * \brief Utility function to constructs a Point with the given coordinates.
   * \param [in] x the x--coordinate of the point.
   * \param [in] y the y--coordinate of the point.
   * \param [in] z the z--coordinate of the point. Default is 0.0.
   * \return p a Point instance with the given coordinates.
   */
  AXOM_HOST_DEVICE
  static Point make_point( const T& x, const T& y, const T& z=0.0 );

  /*!
   * \brief Returns the midpoint between point A and point B.
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \return p point at the midpoint A and B.
   */
  static Point midpoint( const Point& A, const Point& B );

  /*!
   * \brief Performs linear interpolation/extrapolation between two points.
   *
   * Given points A, B and a weight \f$ \alpha \f$, this method computes
   * point \f$ P\f$ defined by \f$ P = (1-\alpha) \cdot A + \alpha \cdot B\f$
   *
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \param [in] alpha user-supplied scalar weight \f$ \alpha\f$
   *
   * \return P the computed point.
   *
   * \note \f$ \forall\alpha \in [0,1] \f$ this method linearly interpolates
   *  between the two points A, B.
   * \note \f$ \forall\alpha \not\in [0,1] \f$ the method extrapolates.
   *
   * \post \f$ P==A\f$ when \f$ \alpha=0.0\f$
   * \post \f$ P==B\f$ when \f$ \alpha=1.0\f$
   * \post The return point, P, and the user-supplied points A, B are collinear.
   */
  static Point lerp( const Point& A, const Point& B, T alpha);

  /*!
   * \brief Helper function to return a point whose coordinates are all 0
   */
  static Point zero() { return Point(); }

  /*!
   * \brief Helper function to return a point whose coordinates are all 1
   * This is equivalent to using the single value constructor: Point(1)
   * (with the appropriate casting) and is only valid for Points with
   * a numerical type (i.e. where static_cast<T>(1) is valid.
   */
  static Point ones() { return Point(static_cast< T >(1)); }

private:
  NumericArray< T,NDIMS > m_components;

};

/// \name Pre-defined point types
/// @{

typedef Point< double,2 > Point2D;
typedef Point< double,3 > Point3D;

/// @}

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  Point implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace primal
{

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Point< T, NDIMS > Point< T,NDIMS >::make_point( const T& x,
                                                       const T& y,
                                                       const T& z )
{
  T tmp_array[3] = { x, y, z};
  return Point(tmp_array, NDIMS);
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Point< T,NDIMS > Point< T,NDIMS >::midpoint(
  const Point< T,NDIMS >& A,
  const Point< T,NDIMS >& B )
{
  Point< T,NDIMS > mid_point;

  for ( int i=0 ; i < NDIMS ; ++i )
  {

    mid_point[ i ] = static_cast<T>(0.5*( A[i]+B[i] ) );
  }

  return mid_point;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Point< T,NDIMS > Point< T,NDIMS >::lerp(
  const Point< T,NDIMS >& A,
  const Point< T,NDIMS >& B,
  T alpha)
{
  PointType res;
  const T beta = 1.-alpha;
  for ( int i=0 ; i < NDIMS ; ++i )
  {
    res[ i ] = beta*A[ i ] + alpha*B[ i ];
  }
  return res;
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
std::ostream& Point< T, NDIMS >::print(std::ostream& os) const
{
  os <<"(";
  for (int dim=0 ; dim < NDIMS -1 ; ++dim)
  {
    os << static_cast< typename NonChar< T >::type >(m_components[dim]) << ",";
  }
  os << static_cast< typename NonChar< T >::type >(m_components[NDIMS-1])
     << ")";

  return os;
}

//------------------------------------------------------------------------------
/// Free functions implementing Point's operators
//------------------------------------------------------------------------------

template < typename T, int NDIMS >
std::ostream& operator<<(std::ostream & os, const Point< T,NDIMS > & pt)
{
  pt.print(os);
  return os;
}

} /* namespace primal*/
} /* namespace axom */

#endif /* POINT_HXX_ */
