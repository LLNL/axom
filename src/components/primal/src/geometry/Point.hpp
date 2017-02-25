/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef POINT_HXX_
#define POINT_HXX_

#include "slic/slic.hpp"

#include "primal/NumericArray.hpp"

// C/C++ includes
#include <cstring> // For memcpy()
#include <ostream> // For print() and operator <<


namespace quest
{

// Forward declare the templated classes and operator functions
template<typename T, int DIM> class Point;

/*!
 * \brief Equality comparison operator for points
 */
template<typename T, int DIM>
bool operator==(const Point<T, DIM> & lhs, const Point<T, DIM>& rhs);

/*!
 * \brief Inequality comparison operator for points
 */
template<typename T, int DIM>
bool operator!=(const Point<T, DIM> & lhs, const Point<T, DIM>& rhs);

/*!
 * \brief Overloaded output operator for points
 */
template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Point<T,DIM> & pt);

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

    typedef Point<T,DIM> PointType;
    typedef T CoordType;

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
  explicit Point(T val = T(), int sz = DIM) : m_components(val,sz) {}


  /*!
   *****************************************************************************
   * \brief Constructor from a numeric array
   * \param [in] arr The numeric array to copy from
   *****************************************************************************
   */
  Point(const NumericArray<T,DIM>& arr) : m_components(arr) {}

  /*!
   *****************************************************************************
   * \brief Creates a point from the first sz values of the input array.
   * \param [in] vals An array containing at least sz values
   * \param [in] sz The number of coordinates to take from the array.  Defaults
   * to DIM.
   * It sz is greater than DIM, we only take the first DIM values.
   *****************************************************************************
   */
  Point(T* vals, int sz = DIM) : m_components(vals,sz) {}

  /*!
   *****************************************************************************
   * \brief Copy constructor.
   * \param [in] other The point to copy
   *****************************************************************************
   */
  Point( const Point& other) : m_components( other.m_components) {}

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
  static int dimension() { return DIM; };

  /*!
   *****************************************************************************
   * \brief Assignment operator.
   * \param [in] rhs a point instance on the right hand side.
   *****************************************************************************
   */
  Point& operator=(const Point& rhs) { m_components = rhs.m_components; return *this;}

  /*!
   *****************************************************************************
   * \brief Access operator for individual components.
   * \param [in] i the component index to access
   * \return p[i] the value at the given component index.
   * \pre (i >= 0) && (i < ndims)
   *****************************************************************************
   */
  const T& operator[](int i) const { return m_components[i]; }
  T& operator[](int i)             { return m_components[i]; }


  /*!
   *****************************************************************************
   * \brief Returns a pointer to the underlying data.
   *****************************************************************************
   */
  const T* data() const             { return m_components.data(); }
  T* data()                         { return m_components.data(); }

  /*!
   *****************************************************************************
   * \brief Returns a reference to the underlying NumericArray.
   *****************************************************************************
   */
  const NumericArray<T,DIM>& array() const  { return m_components; }
  NumericArray<T,DIM>& array()              { return m_components; }


  /*!
   *****************************************************************************
   * \brief Output the point's coordinates to the array
   * \param arr The array that we are outputting to.
   * \pre The user needs to make sure that the array has been allocated
   * and has sufficient space for DIM coordinates.
   *****************************************************************************
   */
  void to_array(T* arr) const       { m_components.to_array(arr); }



  /*!
   * \brief Equality comparison operator for points
   */
  friend bool operator==(const Point& lhs, const Point& rhs)
  { return lhs.m_components == rhs.m_components; }

  /*!
   * \brief Inequality operator for points
   */
  friend bool operator!=(const Point& lhs, const Point& rhs)
  { return !(lhs == rhs); }

  /*!
   *****************************************************************************
   * \brief Simple formatted print of a point instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const;

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
   * \brief Linearly interpolates two points.
   * \param [in] A user-supplied point
   * \param [in] B user-supplied point
   * \param [in] alpha weight with which to interpolate
   * \return p Linearly interpolated point: p = (1-alpha)A + alpha*B.
   *****************************************************************************
   */
  static Point lerp( const Point& A, const Point& B, double alpha);
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
  NumericArray<T,DIM> m_components;

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
inline Point< T, DIM > Point< T,DIM >::make_point( const T& x,
                                                   const T& y,
                                                   const T& z )
{
  T tmp_array[3] = { x, y, z};
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
inline Point< T,DIM > Point< T,DIM >::lerp(
        const Point<T,DIM>& A,
        const Point<T,DIM>& B,
        double alpha)
{
  return PointType((1.-alpha)*A.array() + alpha*B.array());
}


//------------------------------------------------------------------------------
template < typename T, int DIM >
std::ostream& Point< T, DIM >::print(std::ostream& os) const
{
    os <<"(";
    for(int dim=0; dim < DIM -1; ++ dim)
        os << static_cast<typename NonChar<T>::type>(m_components[dim]) << ",";
    os << static_cast<typename NonChar<T>::type>(m_components[DIM-1]) << ")";

    return os;
}

//------------------------------------------------------------------------------
/// Free functions implementing Point's operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const Point<T,DIM> & pt)
{
    pt.print(os);
    return os;
}



} /* namespace quest*/

#endif /* POINT_HXX_ */
