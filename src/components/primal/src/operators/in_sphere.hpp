/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/*!
 * \file
 *
 * \brief Consists of methods that tests whether a point is inside a sphere
 * defined by points.
 *
 * This is a well known computational geometry technique, and can be found in
 * Section 3.1.6.4 in "Real-time collision detection" by C. Ericson
 *
 */

#ifndef IN_SPHERE_H_
#define IN_SPHERE_H_

#include "primal/Point.hpp"
#include "axom_utils/Determinants.hpp"

namespace axom
{
namespace primal
{

/*!
 * \brief Computes whether a query point is inside a circle or not.
 *
 * \param [in] q the query point
 * \param [in] p0 first point that defines a circle
 * \param [in] p1 second point that defines a circle
 * \param [in] p2 third point that defines a circle
 * \return true if the point is inside the circle, false otherwise
 */
template < typename T >
inline bool in_sphere( const Point< T, 2 >& q,
                       const Point< T, 2 >& p0,
                       const Point< T, 2 >& p1,
                       const Point< T, 2 >& p2)
{
  double det = axom::numerics::determinant(
    1.0, p0[0], p0[1], p0[0]*p0[0] + p0[1]*p0[1],
    1.0, p1[0], p1[1], p1[0]*p1[0] + p1[1]*p1[1],
    1.0, p2[0], p2[1], p2[0]*p2[0] + p2[1]*p2[1],
    1.0,  q[0],  q[1],  q[0]* q[0] +  q[1]* q[1]);

  return det < 0;
}


/*!
 * \brief Computes whether a query point is inside a sphere or not.
 *
 * \param [in] q the query point
 * \param [in] p0 first point that defines a sphere
 * \param [in] p1 second point that defines a sphere
 * \param [in] p2 third point that defines a sphere
 * \param [in] p3 fourth point that defines a sphere
 * \return true if the point is inside the sphere, false otherwise
 */
template < typename T >
inline bool in_sphere( const Point< T, 3 >& q,
                       const Point< T, 3 >& p0,
                       const Point< T, 3 >& p1,
                       const Point< T, 3 >& p2,
                       const Point< T, 3 >& p3)
{
  double mat_val[] = {
    1.0, p0[0], p0[1], p0[2], p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2],
    1.0, p1[0], p1[1], p1[2], p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2],
    1.0, p2[0], p2[1], p2[2], p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2],
    1.0, p3[0], p3[1], p3[2], p3[0]*p3[0] + p3[1]*p3[1] + p3[2]*p3[2],
    1.0,  q[0],  q[1],  q[2],  q[0]* q[0] +  q[1]* q[1] +  q[2]* q[2]
  };

  axom::numerics::Matrix<double> mat( 5, 5, mat_val, true);

  double det = axom::numerics::determinant(mat);
  return det < 0;
}

} /* namespace primal */
} /* namespace axom */

#endif /* IN_SPHERE_H_ */
