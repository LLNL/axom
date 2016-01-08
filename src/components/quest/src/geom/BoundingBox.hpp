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
template<typename T, int DIM> bool operator==(const BoundingBox<T, DIM> & lhs, const BoundingBox<T, DIM>& rhs);

/*!
 * \brief Inequality comparison operator for bounding boxes.
 * Two bounding boxes are unequal when they have different bounds
 */
template<typename T, int DIM> bool operator!=(const BoundingBox<T, DIM> & lhs, const BoundingBox<T, DIM>& rhs);

/*!
 * \brief Overloaded output operator for bounding boxes
 */
template<typename T, int DIM> std::ostream& operator<<(std::ostream & os, const BoundingBox<T,DIM> & pt);



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
    typedef Point<CoordType, DIM> PointType;
    typedef Vector<CoordType, DIM> VectorType;
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
   * \param [in] point to include.
   *****************************************************************************
   */
  void addPoint(const PointType& pt);

  /*!
   *****************************************************************************
   * \brief Updates bounds to include the provided bounding box.
   * Convenience function -- equivalent to adding the min and max point of bbox
   * \param [in] bbox to include.
   *****************************************************************************
   */
  void addBox(const BoundingBox& bbox);



  /*!
   *****************************************************************************
   * \brief Expands the lower and upper bounds by the given amount.
   * \param [in] expansionAmount an absolute amount to expand
   * Moves min point expansionAmount away from center and
   * max point expansionAmount away from center (component-wise).
   * If expansionAmount is negative, the bounding box will contract
   * This function checks to ensure that the bounding box is valid after expansion.
   *****************************************************************************
   */
  void expand(CoordType expansionAmount);

  /*!
   *****************************************************************************
   * \brief Scales the bounding box by a fixed amount.
   * \param [in] scaleFactor the multiplicative factor in which to scale
   * If scaleFactor is less than 1, the bounding box will shrink
   * \pre scaleFactor must be greater than or equal to 0
   * This function checks to ensure that the bounding box is valid after inflation.
   *****************************************************************************
   */
  void scale(double scaleFactor);

  /*!
   *****************************************************************************
   * \brief Shifts the bounding box by a fixed amount.
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


  /*!
   *****************************************************************************
   * \brief Simple formatted print of a bounding box instance
   * \param os The output stream to write to
   * \return A reference to the modified ostream
   *****************************************************************************
   */
  std::ostream& print(std::ostream& os) const;

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
    void BoundingBox<CoordType, DIM>::addBox (const BoundingBox& bbox)
    {
        addPoint(bbox.getMin());
        addPoint(bbox.getMax());
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
        SLIC_ASSERT(scaleFactor >= 0.);

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
        m_min = PointType( std::numeric_limits< CoordType>::max() );
        m_max = PointType( std::numeric_limits< CoordType>::min() );
    }


    //------------------------------------------------------------------------------
    template < typename T, int DIM >
    std::ostream& BoundingBox< T, DIM >::print(std::ostream& os) const
    {
        os <<"{ min:"<<m_min <<"; max:"<< m_max <<"; range:"<< range() << " }";
        return os;
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
