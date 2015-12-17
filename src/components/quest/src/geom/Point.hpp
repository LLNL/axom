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

namespace quest
{

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
    enum {
        DIMENSION = DIM
        , BYTES = DIM * sizeof(T)
    };

public:

  /*!
   *****************************************************************************
   * \brief Default Constructor. Creates a point at the origin.
   *****************************************************************************
   */
  Point();

  /*!
   *****************************************************************************
   * \brief Copy constructor.
   * \param [in] rhs
   *****************************************************************************
   */
  Point( const Point<T,DIM>& rhs ) { *this = rhs; };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
   ~Point();

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
  Point<T,DIM>& operator=(const Point<T,DIM>& rhs);

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
  const T& data() const;
  T& data();

  /*!
   *****************************************************************************
   * \brief Constructs a Point instance with the given coordinates.
   * \param [in] x the x--coordinate of the point.
   * \param [in] y the y--coordinate of the point.
   * \param [in] z the z--coordinate of the point. Default is 0.0.
   * \return p a Point instance with the given coordinates.
   *****************************************************************************
   */
  static Point<T,DIM> make_point( const T& x, const T& y, const T& z=0.0 );

  /*!
   *****************************************************************************
   * \brief Returns the midpoint between point A and point B.
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \return p point at the midpoint A and B.
   *****************************************************************************
   */
  static Point<T,DIM> midpoint( const Point<T,DIM>& A,
                                  const Point<T,DIM>& B );


  static Point zero() { return Point(); }
  static Point ones() { return Point::fromValue( static_cast<T>(1)); }
  static Point fromValue(T val);

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
Point< T, DIM >::Point()
{
  SLIC_ASSERT( DIM >= 1 );
  std::fill( m_components, m_components+DIM, T() );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
Point< T, DIM >::~Point()
{

}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Point< T, DIM > Point< T, DIM >::fromValue(T val)
{
  Point pt;
  std::fill( pt.data(), pt.data()+DIM, val );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline Point<T,DIM>& Point< T,DIM >::operator=(const Point<T,DIM>& rhs )
{

  if( this == &rhs ) {
    return *this;
  }

  // copy all the data
  memcpy( m_components, rhs.m_components, DIM*sizeof( T ) );
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

  Point< T, DIM > p;

  for ( int i=0; i < DIM; ++i ) {
     p[ i ] = tmp_array[ i ];
  }

  return( p );
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
  SLIC_ASSERT( (i >= 0) && (i < DIM) );
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T& Point< T, DIM >::operator[](int i) const
{
  SLIC_ASSERT( (i >= 0) && (i < DIM) );
  return m_components[ i ];
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
inline const T& Point< T, DIM >::data() const
{
  return &m_components[ 0 ];
}
//------------------------------------------------------------------------------
template < typename T, int DIM >
inline T& Point< T, DIM >::data()
{
    return &m_components[ 0 ];
}


} /* namespace quest*/

#endif /* POINT_HXX_ */
