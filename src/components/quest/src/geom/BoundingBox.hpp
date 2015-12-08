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

#include <cstring> // for memcpy

namespace quest
{

class BoundingBox
{
public:

  /*!
   *****************************************************************************
   * \brief Constructor. Creates a bounding box instance within [0,1]
   *****************************************************************************
   */
  BoundingBox();

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
  virtual ~BoundingBox();

  /*!
   *****************************************************************************
   * \brief Sets the min corner of this bounding box instance.
   * \param [in] min user-supplied min coordinates.
   *****************************************************************************
   */
  void setMin( double min[3] ) { memcpy( m_min, min, 3*sizeof(double) ); } ;

  /*!
   *****************************************************************************
   * \brief Sets the max corner of this bounding box instance.
   * \param [in] max user-supplied max coordinates.
   *****************************************************************************
   */
  void setMax( double max[3] ) { memcpy( m_max, max, 3*sizeof(double) ); } ;

  /*!
   *****************************************************************************
   * \brief Returns const pointer to the min corner of the bounding box.
   * \return ptr pointer to the min corner of the bounding box.
   * \post ptr != ATK_NULLPTR
   *****************************************************************************
   */
  const double* getMin() const { return &m_min[0]; };

  /*!
   *****************************************************************************
   * \brief Returns const pointer to the max corner of the bounding box.
   * \return ptr pointer to the max corner of the bounding box.
   * \post ptr != ATK_NULLPTR
   *****************************************************************************
   */
  const double* getMax() const { return &m_max[0]; };

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
  bool hasPoint( double x, double y, double z=0.0 ) const;

  /*!
   *****************************************************************************
   * \brief Computes the union of two bounding boxes, b1, b2.
   * \param [in] b1 user-supplied bounding box.
   * \param [in] b2 user-supplied bounding box.
   * \return BB the bounding box that contains both b1 and b2.
   *****************************************************************************
   */
  static BoundingBox box_union( const BoundingBox& b1,
                                const BoundingBox& b2 );

private:
  double m_min[3];
  double m_max[3];
};

} /* namespace quest */


//------------------------------------------------------------------------------
//              BoundingBox Implementation
//------------------------------------------------------------------------------
namespace quest {

inline BoundingBox BoundingBox::box_union( const BoundingBox& b1,
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
