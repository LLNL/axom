/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \brief Consists of a set of templated (overloaded) routines used to calculate
 *  the orientation of a given point to another geometric entity.
 *
 * \date Dec 17, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef ORIENTATION_HPP_
#define ORIENTATION_HPP_

#include "quest/Determinants.hpp"
#include "quest/Point.hpp"
#include "quest/Triangle.hpp"
#include "quest/fuzzy_compare.hpp"
#include "slic/slic.hpp"

namespace quest
{


/*!
 *******************************************************************************
 * \brief Enumerates possible return values for orientation methods
 *******************************************************************************
 */
enum OrientedSide {
    ON_BOUNDARY,       /*!< point is on boundary of the given primitive      */
    ON_POSITIVE_SIDE,  /*!< point is on positive side of the given primitive */
    ON_NEGATIVE_SIDE   /*!< point is on negative side of the given primitive */
};

/*!
 *******************************************************************************
 * \brief Computes the orientation of the given point, p, with respect to a
 *  supplied oriented triangle.
 * \param [in] p the query point.
 * \param [in] tri oriented triangle.
 * \return the
 * \note The return value can be one of the following:
 * <ul>
 *  <li> ON_BOUNDARY      </li>
 *  <li> ON_POSITIVE_SIDE </li>
 *  <li> ON_NEGATIVE_SIDE </li>
 * </ul>
 *******************************************************************************
 */
template < typename T, int NDIMS >
inline
int orientation( const Point< T,NDIMS >& p, const Triangle< T,NDIMS >& tri )
{

   double det = math::determinant( tri.A()[0], tri.A()[1], tri.A()[2], 1.0,
                                   tri.B()[0], tri.B()[1], tri.B()[2], 1.0,
                                   tri.C()[0], tri.C()[1], tri.C()[2], 1.0,
                                         p[0],       p[1],       p[2], 1.0  );

   int orient = -1;

   if ( math::fuzzy_compare( det, 0.0, 1.0e-9 ) ) {

       orient = ON_BOUNDARY;

   } else if ( det < 0.0f ) {

       // outside
       orient = ON_POSITIVE_SIDE;

   } else {

       // inside
       orient = ON_NEGATIVE_SIDE;

   }

   return orient;
}


}
#endif /* ORIENTATION_HPP_ */
