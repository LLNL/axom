/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Aug 28, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef POINT_HXX_
#define POINT_HXX_

#include "slic/slic.hpp"


namespace {
    /*!
     *****************************************************************************
     * \brief Utility function that clamps an input val to a given range.
     * \param [in] val  The value to clamp
     * \param [in] lower The lower range
     * \param [in] upper The upper range
     * \return The clamped value.
     * \post lower <= returned value <= upper.
     *****************************************************************************
     */
    template<typename T>
    T clamp(T val, T lower, T upper)
    {
        SLIC_ASSERT( lower <= upper);
        return std::min( std::max( val, lower), upper);
    }
}


namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Point;

/*!
 * \brief Equality comparison operator for points
 */
template<typename T, int DIM> bool operator==(const Point<T, DIM> & lhs, const Point<T, DIM>& rhs);

/*!
 * \brief Inequality comparison operator for points
 */
template<typename T, int DIM> bool operator!=(const Point<T, DIM> & lhs, const Point<T, DIM>& rhs);


/*!
 *******************************************************************************
 * \class Point
 *
 * \brief The point class represents a point, \f$ p \in \mathcal{R}^d \f$ . It
 *  provides access methods to set and query the point coordinates.
 *******************************************************************************
 */
template < typename T, int DIM >
class Point
{
public:
    enum { NDIMS = DIM
         , NBYTES = DIM * sizeof(T)
         };

public:

  /*!
   *****************************************************************************
   * \brief Fill the first sz coordinates with val and zeros the rest
   * \param [in] val The value to set the coordinates to.  Defaults to zero
   * \param [in] sz The number of coordinates to set to val.
   * The rest will be set to zero.  Defaults is DIM.
   * If sz is greater than DIM, we set all coordinates to val
   *****************************************************************************
   */
  explicit Point(T val = T(), int sz = DIM);

  /*!
   *****************************************************************************
   * \brief Creates a point from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz The number of coordinates to take from the array.  Defaults
   * to DIM.
   * It sz is greater than DIM, we only take the first DIM values.
   *****************************************************************************
   */
  Point(T* vals, int sz = DIM);

  /*!
   *****************************************************************************
   * \brief Copy constructor.
   * \param [in] rhs
   *****************************************************************************
   */
  Point( const Point& rhs ) { *this = rhs; };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
   ~Point() {}

  /*!
   *****************************************************************************
   * \brief Returns the dimension of this point instance.
   * \return d the dimension of the point.
   * \post d >= 1.
   *****************************************************************************
   */
  int dimension() const { return DIM; };

  /*!
   *****************************************************************************
   * \brief Assignment operator.
   * \param [in] rhs a point instance on the right hand side.
   *****************************************************************************
   */
  Point& operator=(const Point& rhs);

  /*!
   *****************************************************************************
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   *****************************************************************************
   */
  const T& operator[](int i) const;
  T& operator[](int i);


  /*!
   *****************************************************************************
   * \brief Returns a pointer to the underlying data.
   *****************************************************************************
   */
  const T* data() const;
  T* data();

  /*!
   *****************************************************************************
   * \brief Output the point's coordinates to the array
   * \param arr The array that we are outputting to.
   * \pre The user needs to make sure that the array has been allocated
   * and has sufficient space for DIM coordinates.
   *****************************************************************************
   */
  void to_array(T* arr) const;

  /*!
   *****************************************************************************
   * \brief Utility function to constructs a Point with the given coordinates.
   * \param [in] x the x--coordinate of the point.
   * \param [in] y the y--coordinate of the point.
   * \param [in] z the z--coordinate of the point. Default is 0.0.
   * \return p a Point instance with the given coordinates.
   *****************************************************************************
   */
  static Point make_point( const T& x, const T& y, const T& z=0.0 );

  /*!
   *****************************************************************************
   * \brief Returns the midpoint between point A and point B.
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \return p point at the midpoint A and B.
   *****************************************************************************
   */
  static Point midpoint( const Point& A, const Point& B );


  /*!
   *****************************************************************************
   * \brief Helper function to return a point whose coordinates are all 0
   *****************************************************************************
   */
  static Point zero() { return Point(); }

  /*!
   *****************************************************************************
   * \brief Helper function to return a point whose coordinates are all 1
   * This is equivalent to using the single value constructor: Point(1)
   * (with the appropriate casting) and is only valid for Points with
   * a numerical type (i.e. where static_cast<T>(1) is valid.
   *****************************************************************************
   */
  static Point ones() { return Point(static_cast<T>(1)); }

private:
  void verifyIndex(int idx) const { SLIC_ASSERT(idx >= 0 && idx < DIM); }

protected:
  T m_components[ DIM ];
};

/// \name Pre-defined point types
/// @{

typedef Point<double,2> Point2D;
typedef Point<double,3> Point3D;

/// @}

} /* namespace quest */

//------------------------------------------------------------------------------
//  Point implementation
//------------------------------------------------------------------------------
namespace quest {

//------------------------------------------------------------------------------
template < typename T, int DIM >
Point< T, DIM >::Point(T val, int sz)
{
  // NOTE (KW): This should be a static assert in the class
  SLIC_ASSERT( DIM >= 1 );

  const int nvals = clamp(sz, 0, DIM);

  std::fill( m_components, m_components+nvals, val );

  // Fill any remaining coordinates with zero
  if(nvals < DIM)
  {
      std::fill( m_components+nvals, m_components+DIM, T() );
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
Point< T, DIM >::Point(T* vals, int sz)
{
  SLIC_ASSERT( DIM >= 1 );

  const int nvals = clamp(sz, 0, DIM);

  std::copy( vals, vals+nvals, m_components);

  // Fill any remaining coordinates with zero
  if(nvals < DIM)
  {
      std::fill( m_components+nvals, m_components+DIM, T());
  }
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Point<T,DIM>& Point< T,DIM >::operator=(const Point<T,DIM>& rhs )
{

  if( this == &rhs ) {
    return *this;
  }

  // copy all the data
  memcpy( m_components, rhs.m_components, NBYTES);
  return *this;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Point< T, DIM > Point< T,DIM >::make_point( const T& x,
                                                   const T& y,
                                                   const T& z )
{
  T tmp_array[3];
  tmp_array[0] = x;
  tmp_array[1] = y;
  tmp_array[2] = z;

  return Point(tmp_array, DIM);
}

//------------------------------------------------------------------------------
template< typename T, int DIM >
inline Point< T,DIM > Point< T,DIM >::midpoint(
        const Point<T,DIM>& A,
        const Point<T,DIM>& B )
{
  Point< T,DIM > mid_point;

  for ( int i=0; i < DIM; ++i ) {

     mid_point[ i ] = 0.5*( A[i]+B[i] );
  }

  return mid_point;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T& Point< T, DIM >::operator[](int i)
{
    verifyIndex(i);
    return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T& Point< T, DIM >::operator[](int i) const
{
  verifyIndex(i);
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T* Point< T, DIM >::data() const
{
  return m_components;
}
//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T* Point< T, DIM >::data()
{
    return m_components;
}

template < typename T, int DIM >
void Point< T, DIM >::to_array(T* arr) const
{
    SLIC_ASSERT( arr != ATK_NULLPTR);
    memcpy( arr, m_components, NBYTES );
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
bool operator==(const Point<T, DIM>& lhs, const Point<T, DIM>& rhs)
{
    for(int dim=0;dim<DIM;++dim)
    {
        if( lhs[dim] != rhs[dim])
            return false;
    }
    return true;
}

//------------------------------------------------------------------------------

template<typename T, int DIM>
bool operator!=(const Point<T, DIM>& lhs, const Point<T, DIM>& rhs)
{
    return !(lhs == rhs);
}


} /* namespace quest*/

#endif /* POINT_HXX_ */
