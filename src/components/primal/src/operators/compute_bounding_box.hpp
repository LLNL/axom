/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file compute_bounding_box.hpp
 *
 * \brief Consists of functions to create bounding boxes.
 *******************************************************************************
 */

#ifndef COMPUTE_BOUNDING_BOX_HPP_
#define COMPUTE_BOUNDING_BOX_HPP_

#include "primal/NumericArray.hpp" // for numeric arrays
#include "axom_utils/Utilities.hpp" // for nearly equal
#include "axom_utils/Matrix.hpp" // for Matrix
#include "axom_utils/eigen_solve.hpp" // for eigen_solve
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/OrientedBoundingBox.hpp"

namespace axom {
namespace primal {

/*!
 *****************************************************************************
 * \brief Constructor. Creates an oriented bounding box which contains
 * the collection of passed in points.
 * \param [in] pts array of points
 * \param [in] n number of points
 * \note if n <= 0, invokes default constructor
 *****************************************************************************
 */
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS > OBB_from_points(const Point< T, NDIMS > *pts,
  int n)
{
  if (n <= 0) {
    OrientedBoundingBox< T, NDIMS > res;
    return res;
  }

  const int DEPTH = 100;  // controls depth of power method (i.e. accuracy)

  numerics::Matrix< T > covar = numerics::Matrix< T >(NDIMS, NDIMS);
  NumericArray< T, NDIMS > curr;
  NumericArray< T, NDIMS > c;

  for (int i = 0; i < n; i++) {
    curr = pts[i].array();
    for (int j = 0; j < NDIMS; j++) {
      c[j] += curr[j];
    }
  }

  // average
  for (int j = 0; j < NDIMS; j++) {
    c[j] /= static_cast< T >(n);
  }

  for (int i = 0; i < n; i++) {
    curr = pts[i].array();
    for (int j = 0; j < NDIMS; j++) {
      for (int k = 0; k < NDIMS; k++) {
        covar(j,k) += (curr[j] - c[j])*(curr[k] - c[k]);
      }
    }
  }

  // average the covariance matrix
  for (int j = 0; j < NDIMS; j++) {
    for (int k = 0; k < NDIMS; k++) {
      covar(j,k) /= static_cast< T >(n);
    }
  }

  T u[NDIMS*NDIMS];
  T lambdas[NDIMS];

  SLIC_ASSERT(numerics::eigen_solve< T >(covar, NDIMS, DEPTH, u, lambdas) == 0);

  T max;
  T dot;

  Vector< T, NDIMS > w[NDIMS];
  Vector< T, NDIMS > e;

  for (int i = 0; i < NDIMS; i++) {
    w[i] = Vector< T, NDIMS >(u + NDIMS*i);

    // compute extent in this direction
    max = T();
    for (int j = 0; j < n; j++) {
      dot = T();
      curr = pts[j].array();
      for (int k = 0; k < NDIMS; k++) {
        dot += u[NDIMS*i + k]*curr[k];
      }
      if (dot < T()) dot = -dot;

      if (max < dot) max = dot;
    }
    e[i] = max;
  }

  OrientedBoundingBox< T, NDIMS > res(Point< T, NDIMS >(c), w, e);
  return res;
}

}  /* namespace axom */
}  /* namespace primal */


#endif  /* COMPUTE_BOUNDING_BOX_HPP_ */
