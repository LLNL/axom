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


namespace quest
{

template<typename CoordType, int DIM>
class BoundingBox
{
public:
    typedef Point<CoordType, DIM> PointType;
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
    : m_min( PointType( std::numeric_limits< CoordType>::max() ) )
    , m_max( PointType( std::numeric_limits< CoordType>::min() ) ) {}


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
   * \brief Sets the min corner of this bounding box instance.
   * \param [in] min user-supplied min coordinates.
   *****************************************************************************
   */
  void setMin(const PointType& pt) { m_min = pt; } ;

  /*!
   *****************************************************************************
   * \brief Sets the max corner of this bounding box instance.
   * \param [in] max user-supplied max coordinates.
   *****************************************************************************
   */
  void setMax(const PointType& pt) { m_max = pt; } ;

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
   * \return sonct reference to the max corner of the bounding box.
   *****************************************************************************
   */
  const PointType& getMax() const { return m_max; };


  /*!
   *****************************************************************************
   * \brief Updates bounds to include the provided point.
   * \param [in] point to include.
   *****************************************************************************
   */
  void addPoint(const PointType& pt);

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
   * \param [in] x the x--coordinate
   * \param [in] y the y--coordinate
   * \param [in] z the z--coordinate
   * \return status true if point inside the box, else false.
   *****************************************************************************
   */
  template<typename OtherType>
  bool contains( const Point<OtherType, DIM>& otherPt) const;


  /*!
   *****************************************************************************
   * \brief Checks that we have a valid bounding box.
   * A bounding box is valid when its extent in each dimension
   * (max coordinate minus min coordinate) is greater than or equal to zero
   * \return status true if point inside the box, else false.
   *****************************************************************************
   */
  bool isValid() const;

static BoundingBox box_union( const BoundingBox& b1,
                                const BoundingBox& b2 );
private:
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
namespace quest{


    //------------------------------------------------------------------------------
    template<typename CoordType, int DIM>
    BoundingBox<CoordType, DIM>& BoundingBox<CoordType, DIM>::operator=(const BoundingBox& rhs )
    {

      if ( this != &rhs ) {
        m_min = rhs.m_min;
        m_max = rhs.m_max;
      }

      return *this;
    }

    //------------------------------------------------------------------------------
    template<typename CoordType, int DIM>
    template<typename OtherType>
    bool BoundingBox<CoordType, DIM>::contains(const Point<OtherType,DIM>& otherPt) const
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
    void BoundingBox<CoordType, DIM>::addPoint (const PointType& pt)
    {
        for (int dim=0; dim < DIM; ++dim ) {
            if ( pt[dim] < m_min[dim] ) {
                m_min[dim] = pt[dim];
            }
            if ( pt[dim] > m_max[dim] ) {
                m_max[dim] = pt[dim];
            }
        }
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

template<int DIM>
inline BoundingBox<DIM> BoundingBox<DIM>::box_union( const BoundingBox& b1,
                                           const BoundingBox& b2 )
{
   BoundingBox bb;

   double min[3];
   double max[3];
   for (int i=0; i < 3; ++i ) {
      min[ i ] = std::min( b1.getMin()[i], b2.getMin()[i] );
      max[ i ] = std::max( b1.getMax()[i], b2.getMax()[i] );
   }

   bb.setMin( min );
   bb.setMax( max );
   return( bb );
}

} // end namespace quest


#endif /* BOUNDINGBOX_HPP_ */
