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

#include "quest/Point.hpp"


namespace quest
{

template<int DIM>
class BoundingBox
{
public:
    typedef double CoordType;
    typedef Point<CoordType, DIM> PointType;
public:

  /*!
   *****************************************************************************
   * \brief Constructor. Creates a bounding box instance within [0,1]
   *****************************************************************************
   */
  BoundingBox() : m_min( PointType::zero()), m_max( PointType::ones() ) {};

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
   ~BoundingBox();

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
   * \brief Returns const pointer to the min corner of the bounding box.
   * \return ptr pointer to the min corner of the bounding box.
   * \post ptr != ATK_NULLPTR
   *****************************************************************************
   */
  const PointType& getMin() const { return m_min; };

  /*!
   *****************************************************************************
   * \brief Returns const pointer to the max corner of the bounding box.
   * \return ptr pointer to the max corner of the bounding box.
   * \post ptr != ATK_NULLPTR
   *****************************************************************************
   */
  const PointType& getMax() const { return m_max; };

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
   * @param [in] x the x--coordinate
   * @param [in] y the y--coordinate
   * @param [in] z the z--coordinate
   * @return status true if point inside the box, else false.
   *****************************************************************************
   */
  template<typename T>
  bool hasPoint( const Point<T, DIM>& otherPt) const;


  /*!
   *****************************************************************************
   * \brief Checks that we have a valid bounding box.
   * A bounding box is valid its extent in each dimension
   * (max point minus min point) is greater than or equal to zero
   * @return status true if point inside the box, else false.
   *****************************************************************************
   */
  bool isValid() const;

static BoundingBox box_union( const BoundingBox& b1,
                                const BoundingBox& b2 );
private:
  PointType m_min;
  PointType m_max;
};

} /* namespace quest */


//------------------------------------------------------------------------------
//              BoundingBox Implementation
//------------------------------------------------------------------------------
namespace quest {

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

} /* namespace quest */


#endif /* BOUNDINGBOX_HPP_ */
