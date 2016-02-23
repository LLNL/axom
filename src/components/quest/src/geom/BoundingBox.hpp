/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef BOUNDINGBOX_HPP_
#define BOUNDINGBOX_HPP_

#include <limits>

#include "quest/Point.hpp"
#include "quest/Vector.hpp"


namespace quest
{


// Forward declare the templated classes and operator functions
template<typename T, int DIM> class BoundingBox;

/*!
 * \brief Equality comparison operator for bounding boxes.
 * Two bounding boxes are equal when they have the same bounds
 */
template<typename T, int DIM>
bool operator==(const BoundingBox<T, DIM> & lhs,const BoundingBox<T, DIM>& rhs);

/*!
 * \brief Inequality comparison operator for bounding boxes.
 * Two bounding boxes are unequal when they have different bounds
 */
template<typename T, int DIM>
bool operator!=(const BoundingBox<T, DIM> & lhs,const BoundingBox<T, DIM>& rhs);

/*!
 * \brief Overloaded output operator for bounding boxes
 */
template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const BoundingBox<T,DIM> & pt);


/**
 * \brief Type trait to find the highest and lowest values for a given numeric type
 *
 * \note numeric_limits::max() always provides the highest possible value for all numeric type.
 * \note For integral values, numeric_limits ::min() provides the lowest value,
 *       but for float and double, it provides the smallest positive number.
 *       This was fixed in cxx11 with the function numeric_limits::lowest()
 */
template<typename T>
struct ValueRange
{
    /** \brief Returns the highest representable value of type T */
    static const T highest() { return std::numeric_limits<T>::max(); }

    /** \brief Returns the lowest representable value of type T */
    static const T lowest() {
        #ifdef USE_CXX11
            return std::numeric_limits<T>::lowest();
        #else
            return std::numeric_limits<T>::min();
        #endif
    }
};

#ifndef USE_CXX11

/**
 * \brief Template specialization of ValueRange for float types
 * \note Only necessary for pre-CXX11
 */
template<> struct ValueRange<float>
{
    typedef float T;
    static const T highest() { return std::numeric_limits<T>::max(); }
    static const T lowest()  { return -std::numeric_limits<T>::max(); }
};

/**
 * \brief Template specialization of ValueRange for double types
 * \note Only necessary for pre-CXX11
 */
template<> struct ValueRange<double>
{
    typedef double T;
    static const T highest() { return std::numeric_limits<T>::max(); }
    static const T lowest()  { return -std::numeric_limits<T>::max(); }
};
#endif


/*!
 *******************************************************************************
 * \class
 *
 * \brief The bounding box represents and axis-aligned bounding box.
 * It is defined by two points: its lowermost point and its uppermost point.
 * A bounding box is considered to be closed on all sides
 * -- it contains its min and max points.
 * \note Do we need to consider half-open bounding boxes
 * -- where the upper boundaries are not contained?
 * \note Does user code have to set a tolerance, or should the bounding box support this?
 * This would make it easier for user code to set this at the beginning and not have to worry
 * about maintaining a fuzzy tolerance while adding the points.
 *******************************************************************************
 */
template<typename CoordType, int DIM>
class BoundingBox
{
public:
    /*! \brief The underlying Point type of the bounding box */
    typedef Point<CoordType, DIM> PointType;

    /*! \brief The underlying Vector type of the bounding box */
    typedef Vector<CoordType, DIM> VectorType;

    /*! \brief The corresponding BoxType */
    typedef BoundingBox<CoordType,DIM> BoxType;
public:

  /*!
   *****************************************************************************
   * \brief Constructor. Creates a bounding box with an invalid bound
   * The lower bound is set to the greatest possible point and the upper bound
   * is set to the smallest possible point.  This way adding any point resets
   * the bounds to a valid range.
   *****************************************************************************
   */
  BoundingBox()
    : m_min( PointType( ValueRange< CoordType>::highest() ) )
    , m_max( PointType( ValueRange< CoordType>::lowest() ) ) {}


  /*!
   *****************************************************************************
   * \brief Constructor. Creates a bounding box containing a single point
   *****************************************************************************
   */
  BoundingBox(const PointType& pt)
      : m_min( pt), m_max( pt)
  {
  }

  /*!
   *****************************************************************************
   * \brief Constructor. Creates a bounding box with a given min and max point
   * The code ensures that the bounds are valid.
   *****************************************************************************
   */
  BoundingBox(const PointType& lowerPt, const PointType& upperPt)
      : m_min( lowerPt), m_max( upperPt)
  {
      checkAndFixBounds();
  }

  /*!
   *****************************************************************************
   * \brief Copy Constructor.
   * \param [in] rhs
   *****************************************************************************
   */
  BoundingBox( const BoundingBox& rhs ) { *this = rhs; };

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~BoundingBox() {}


  /*!
   *****************************************************************************
   * \brief Resets the bounds to those of the default constructor
   * \note This invalidates the bounding box (i.e. isValid() will be false)
   *****************************************************************************
   */
  void clear();

  /*!
   *****************************************************************************
   * \brief Returns const reference to the min corner of the bounding box.
   * \return const reference to the min corner of the bounding box.
   *****************************************************************************
   */
  const PointType& getMin() const { return m_min; };

  /*!
   *****************************************************************************
   * \brief Returns const reference to the max corner of the bounding box.
   * \return const reference to the max corner of the bounding box.
   *****************************************************************************
   */
  const PointType& getMax() const { return m_max; };

  /*!
   *****************************************************************************
   * \brief Returns a vector from the min to the max points of the bounding box
   * \return Vector from min point to max point of bounding box.
   *****************************************************************************
   */
  VectorType range() const { return VectorType(getMin(), getMax()); };

  /*!
   *****************************************************************************
   * \brief Updates bounds to include the provided point.
   * \param [in] pt to include.
   *****************************************************************************
   */
  template<typename OtherType>
  void addPoint(const Point<OtherType,DIM>& pt);

  /*!
   *****************************************************************************
   * \brief Updates bounds to include the provided bounding box.
   * Convenience function -- equivalent to adding the min and max point of bbox
   * \param [in] bbox to include.
   *****************************************************************************
   */
  template<typename OtherType>
  void addBox(const BoundingBox<OtherType,DIM>& bbox);

  /*!
   *****************************************************************************
   * \brief Computes the centroid of this bounding box instance.
   * \return pt point at the bounding box centroid.
   *****************************************************************************
   */
  PointType centroid() const;

  /*!
   *****************************************************************************
   * \brief Returns the dimension of the ambient space for this bounding box.
   * \return d the dimension of this bounding box instance.
   * \post d >= 1.
   *****************************************************************************
   */
  int dimension() const { return DIM; };

  /*!
   *****************************************************************************
   * \brief Finds the longest dimension of the bounding box
   * \return idx the index of the longest dimension.
   * \post idx >= 0 < DIM
   * \note In the case of ties, where the bounding box has more than one side
   *  with the same length, the code picks the first dimension as the longest
   *  dimension.
   *****************************************************************************
   */
  int getLongestDimension() const;

  /*!
   *****************************************************************************
   * \brief Expands the lower and upper bounds by the given amount.
   * \param [in] expansionAmount an absolute amount to expand
   * Moves min point expansionAmount away from center and
   * max point expansionAmount away from center (component-wise).
   * This function checks to ensure that the bounding box is valid after expansion.
   * \note If expansionAmount is negative, the bounding box will contract
   *****************************************************************************
   */
  void expand(CoordType expansionAmount);

  /*!
   *****************************************************************************
   * \brief Scales the bounding box towards its center by a given amount.
   * \param [in] scaleFactor the multiplicative factor by which to scale
   * \note Checks to ensure that the bounding box is valid after inflation.
   * \note If scaleFactor is less than 1, the bounding box will shrink.
   * \note If scaleFactor is 0, the bounding box will shrink to its midpoint
   * \note The sign of the shrinkFactor has no effect since we are shrinking
   *  towards the center, and we fix the bounds after shrinking
   *****************************************************************************
   */
  void scale(double scaleFactor);

  /*!
   *****************************************************************************
   * \brief Shifts the bounding box by a fixed displacement.
   * \param [in] displacement the amount with which to move the bounding box
   *****************************************************************************
   */
  void shift(const VectorType& displacement);

  /*!
   *****************************************************************************
   * \brief Overloaded assignment operator.
   * \param [in] rhs bounding box instance on the right-hand side
   * \return
   *****************************************************************************
   */
  BoundingBox& operator=(const BoundingBox& rhs );

  /*!
   *****************************************************************************
   * \brief Checks whether the box contains the point
   * \param [in] otherPt the point that we are checking
   * \return status true if point inside the box, else false.
   *****************************************************************************
   */
  template<typename OtherType>
  bool contains( const Point<OtherType, DIM>& otherPt) const;

  /*!
   *****************************************************************************
   * \brief Checks whether the box contains another bounding box
   * \param [in] otherBB the bounding box that we are checking
   * \return status true if bb is inside the box, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two CoordTypes are comparable with
   *  operator<().
   *****************************************************************************
   */
  template<typename OtherType>
  bool contains( const BoundingBox<OtherType, DIM>& otherBB) const;

  /*!
   *****************************************************************************
   * \param [in] otherBB the bounding box that we are checking.
   * \return status true if bb intersects, else false.
   * \note We are allowing the other bounding box to have a different coordinate
   *  type. This should work as long as the two CoordTypes are comparable with
   *  operator<().
   *****************************************************************************
   */
  template < typename OtherType >
  bool intersects( const BoundingBox< OtherType, DIM >& otherBB ) const;

  /*!
   *****************************************************************************
   * \brief Checks that we have a valid bounding box.
   * \note A bounding box is valid when the length of each dimension is greater
   *  than or equal to zero.
   * \return status true if point inside the box, else false.
   *****************************************************************************
   */
  bool isValid() const;

  /*!
   *****************************************************************************
   * \brief Subdivides this bounding box instance into two sub-boxes by
   *  splitting along the given dimension. If a dimension is not provided, this
   *  method will split the bounding box along the longest dimension.
   * \param [in/out] right the right sub-box.
   * \param [in/out] left  the left sub-box.
   * \param [in] dimension the dimension to split along (optional)
   * \pre dimension >= -1 && dimension < DIM
   * \note if dimension==-1, the bounding box is split along its longest edge.
   *****************************************************************************
   */
  void bisect( BoxType& right, BoxType& left, int dimension=-1);

  /*!
   *****************************************************************************
   * \brief Simple formatted print of a bounding box instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const;

  /// \name Static methods
  /// @{

  /*!
   *****************************************************************************
   * \brief Returns the list of points of a 2-D BoundingBox instance.
   * \param [in]  bb user-supplied instance of the bounding box.
   * \param [out] pnts the list of points
   * \post pnts.size() == 4
   * \note The ordering of the points has as follows:
   * \verbatim
   *
   *     3      2
   *     +------+
   *     |      |
   *     |      |
   *     +------+
   *     0      1
   *
   * \endverbatim
   *****************************************************************************
   */
  static void getPoints( const BoundingBox< CoordType,2 >& bb,
                         std::vector< Point< CoordType,2 > >& pnts );

  /*!
   *****************************************************************************
   * \brief Returns the list of points of a 3-D BoundingBox instance.
   * \param [in]  bb user-supplied instance of the bounding box.
   * \param [out] pnts the list of points
   * \post pnts.size() == 8
   * \note The ordering of the points has as follows:
   * \verbatim
   *
   *         4          7
   *         +----------+
   *        /|         /|
   *     5 / |      6 / |
   *      +--|-------+  |
   *      |  +-------|--+
   *      | / 0      | / 3
   *      |/         |/
   *      +----------+
   *      1          2
   *
   * \endverbatim
   *****************************************************************************
   */
  static void getPoints( const BoundingBox< CoordType,3 >& bb,
                         std::vector< Point< CoordType,3 > >& pnts );

  /// @}

private:

  /*!
   *****************************************************************************
   * \brief Sets the min point for this bounding box instance.
   * \param [in] newMin the new min point.
   *****************************************************************************
   */
  inline void setMin( const PointType& newMin ) { m_min = newMin; };

  /*!
   *****************************************************************************
   * \brief Sets the max point for this bounding box instance.
   * \param [in] newMax the new max point.
   *****************************************************************************
   */
  inline void setMax( const PointType& newMax ) { m_max = newMax; };

  /*!
   *****************************************************************************
   * \brief Ensures that the bounds are valid.
   * A bounding box is valid when its extent in each dimension
   * (max coordinate minus min coordinate) is greater than or equal to zero
   *****************************************************************************
   */
  void checkAndFixBounds();

private:
  PointType m_min;
  PointType m_max;
};

} /* namespace quest */



//------------------------------------------------------------------------------
//  BoundingBox implementation
//------------------------------------------------------------------------------
namespace quest {


//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
BoundingBox<CoordType, DIM>&
BoundingBox<CoordType, DIM>::operator=(const BoundingBox& rhs )
{

  if ( this != &rhs ) {
    m_min = rhs.m_min;
    m_max = rhs.m_max;
  }

  return *this;
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
template<typename OtherCoordType>
bool BoundingBox<CoordType, DIM>::contains(
        const Point<OtherCoordType,DIM>& otherPt) const
{
    for(int dim = 0; dim < DIM; ++dim)
    {
        if( otherPt[dim] < m_min[dim] || otherPt[dim] >  m_max[dim])
            return false;
    }
    return true;
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
template<typename OtherCoordType>
bool BoundingBox<CoordType, DIM>::contains(
        const BoundingBox<OtherCoordType,DIM>& otherBB) const
{
    return this->contains(otherBB.getMin()) && this->contains(otherBB.getMax());
}

//------------------------------------------------------------------------------
template < typename CoordType, int DIM >
template < typename OtherType  >
bool BoundingBox< CoordType,DIM >::intersects(
        const BoundingBox< OtherType, DIM >& otherBB ) const
{
  for ( int i=0; i < DIM; ++i ) {

     if ( (m_max[ i ] < otherBB.m_min[ i ]) ||
          (m_min[ i ] > otherBB.m_max[ i ]) ) {

         return false;

     } // END if

  } // END for all dimensions

  return true;
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
bool BoundingBox<CoordType, DIM>::isValid() const
{
  for(int dim = 0; dim < DIM; ++dim)
  {
      if( m_min[dim] >  m_max[dim])
          return false;
  }
  return true;

}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
template<typename OtherCoordType>
void BoundingBox<CoordType, DIM>::addPoint(
        const Point<OtherCoordType,DIM>& pt)
{
    for (int dim=0; dim < DIM; ++dim ) {

        CoordType coord = static_cast<CoordType>(pt[dim]);

        if ( coord < m_min[dim] ) {
            m_min[dim] = coord;
        }
        if ( coord > m_max[dim] ) {
            m_max[dim] = coord;
        }
    }
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
template<typename OtherCoordType>
void BoundingBox<CoordType, DIM>::addBox(
        const BoundingBox<OtherCoordType,DIM>& bbox)
{
    addPoint(bbox.getMin());
    addPoint(bbox.getMax());
}

//------------------------------------------------------------------------------
template< typename CoordType, int DIM >
quest::Point< CoordType, DIM > BoundingBox< CoordType,DIM >::centroid() const
{
  SLIC_ASSERT( this->isValid() );

  PointType pt;
  for ( int i=0; i < DIM; ++i ) {

      pt[ i ] = 0.5 * ( m_min[ i ] + m_max[ i ] );

  }

  return( pt );
}

//------------------------------------------------------------------------------
template< typename CoordType, int DIM >
int BoundingBox<CoordType,DIM>::getLongestDimension() const
{
  SLIC_ASSERT( this->isValid() );

  int maxDim = 0;
  CoordType max = std::numeric_limits< CoordType >::min();
  for ( int i=0; i < DIM; ++i ) {

     CoordType dx = m_max[ i ] - m_min[ i ];
     if ( dx > max ) {
        max    = dx;
        maxDim = i;
     }

  } // END for all dimensions

  return maxDim;
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
void BoundingBox<CoordType, DIM>::expand(CoordType expansionAmount)
{
    for (int dim=0; dim < DIM; ++dim ) {
        m_min[dim] -= expansionAmount;
        m_max[dim] += expansionAmount;
    }

    checkAndFixBounds();
}


//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
void BoundingBox<CoordType, DIM>::scale(double scaleFactor)
{
    const PointType midpoint = PointType::midpoint(getMin(), getMax());
    const VectorType r = scaleFactor * 0.5 * range();

    m_min = PointType( midpoint.array() - r.array());
    m_max = PointType( midpoint.array() + r.array());

    checkAndFixBounds();
}


//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
void BoundingBox<CoordType, DIM>::shift(const VectorType& displacement)
{
    m_min.array() += displacement.array();
    m_max.array() += displacement.array();
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
void BoundingBox<CoordType, DIM>::checkAndFixBounds ()
{
    for (int dim=0; dim < DIM; ++dim ) {
        if ( m_min[dim] > m_max[dim] ) {
            std::swap( m_min[dim], m_max[dim] );
        }
    }
}

//------------------------------------------------------------------------------
template<typename CoordType, int DIM>
void BoundingBox<CoordType, DIM>::clear()
{
    m_min = PointType( ValueRange< CoordType>::highest() );
    m_max = PointType( ValueRange< CoordType>::lowest() );
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
std::ostream& BoundingBox< T, DIM >::print(std::ostream& os) const
{
    os <<"{ min:"<<m_min <<"; max:"<< m_max <<"; range:"<< range() << " }";
    return os;
}

//------------------------------------------------------------------------------
template < typename T, int DIM >
void BoundingBox< T,DIM >::bisect( BoxType& right, BoxType& left, int dim )
{
   SLIC_ASSERT( this->isValid() );

   if ( dim < 0 ) {

      dim = this->getLongestDimension();

   }
   SLIC_ASSERT( dim >=0 && dim < DIM );

   // calculate mid along the given dimension
   T mid = 0.5*( m_max[dim] + m_min[dim] );

   // update right
   right.setMin( this->getMin() );
   PointType new_right_max = this->getMax();
   new_right_max[ dim ]    = mid;
   right.setMax( new_right_max );
   SLIC_ASSERT( right.isValid() );

   // update left
   left.setMax( this->getMax() );
   PointType new_left_min = this->getMin();
   new_left_min[ dim ] = mid;
   left.setMin( new_left_min );
   SLIC_ASSERT( left.isValid() );
}

//------------------------------------------------------------------------------
//    Implementation of static methods
//------------------------------------------------------------------------------
template < typename CoordType, int DIM >
inline void BoundingBox< CoordType,DIM >::getPoints(
        const BoundingBox< CoordType,2 >& bb,
        std::vector< Point< CoordType,2 > >& pnts )
{
  pnts.resize( 4 );
  const Point< CoordType,2 >& min = bb.getMin();
  const Point< CoordType,2 >& max = bb.getMax();

  pnts[ 0 ] = Point< CoordType,2 >::make_point( min[0], min[1] );
  pnts[ 1 ] = Point< CoordType,2 >::make_point( max[0], min[1] );
  pnts[ 2 ] = Point< CoordType,2 >::make_point( max[0], max[1] );
  pnts[ 3 ] = Point< CoordType,2 >::make_point( min[0], max[1] );
}

//------------------------------------------------------------------------------
template < typename CoordType, int DIM >
inline void BoundingBox< CoordType,DIM >::getPoints(
        const BoundingBox< CoordType,3 >& bb,
        std::vector< Point< CoordType,3 > >& pnts )
{
  pnts.resize( 8 );
  const Point< CoordType,3 >& min = bb.getMin();
  const Point< CoordType,3 >& max = bb.getMax();

  pnts[ 0 ] = Point< CoordType,3 >::make_point( min[0], min[1], min[2] );
  pnts[ 1 ] = Point< CoordType,3 >::make_point( max[0], min[1], min[2] );
  pnts[ 2 ] = Point< CoordType,3 >::make_point( max[0], max[1], min[2] );
  pnts[ 3 ] = Point< CoordType,3 >::make_point( min[0], max[1], min[2] );

  pnts[ 4 ] = Point< CoordType,3 >::make_point( min[0], min[1], max[2] );
  pnts[ 5 ] = Point< CoordType,3 >::make_point( max[0], min[1], max[2] );
  pnts[ 6 ] = Point< CoordType,3 >::make_point( max[0], max[1], max[2] );
  pnts[ 7 ] = Point< CoordType,3 >::make_point( min[0], max[1], max[2] );
}

//------------------------------------------------------------------------------
/// Free functions implementing comparison and arithmetic operators
//------------------------------------------------------------------------------

template<typename T, int DIM>
bool operator==(const BoundingBox<T, DIM>& lhs, const BoundingBox<T, DIM>& rhs)
{
    return lhs.getMin() == rhs.getMin() && lhs.getMax() == rhs.getMax();
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
bool operator!=(const BoundingBox<T, DIM>& lhs, const BoundingBox<T, DIM>& rhs)
{
    return !(lhs == rhs);
}

//------------------------------------------------------------------------------
template<typename T, int DIM>
std::ostream& operator<<(std::ostream & os, const BoundingBox<T,DIM> & bb)
{
    bb.print(os);
    return os;
}


} // end namespace quest


#endif /* BOUNDINGBOX_HPP_ */
