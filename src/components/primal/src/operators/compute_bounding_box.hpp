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
#include "axom_utils/Matrix.hpp" // for Matrix
#include "axom_utils/eigen_solve.hpp" // for eigen_solve
#include "primal/Point.hpp"
#include "primal/Vector.hpp"
#include "primal/OrientedBoundingBox.hpp"

namespace axom {
namespace primal {

template < typename T, int NDIMS >
class OrientedBoundingBox;

/*!
 *****************************************************************************
 * \brief Creates a bounding box which contains the collection of passed in
 * points.
 *
 * \param [in] pts array of points
 * \param [in] n number of points
 * \note if n <= 0, invokes default constructor
 *****************************************************************************
 */
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS > compute_oriented_bounding_box(
  const Point< T, NDIMS > *pts, int n)
{
  if (n <= 0) {
    OrientedBoundingBox< T, NDIMS > res;
    return res;
  }

  numerics::Matrix< T > covar = numerics::Matrix< T >(NDIMS, NDIMS);
  NumericArray< T, NDIMS > curr;
  NumericArray< T, NDIMS > c;  // centroid

  for (int i = 0; i < n; i++) {
    c +=  pts[i].array();
  }
  c /= static_cast< T >(n);

  // save space for pts minus the centroid
  NumericArray< T, NDIMS > *ptsMinusCentroid = new NumericArray< T, NDIMS >[n];

  for (int i = 0; i < n; i++) {
    ptsMinusCentroid[i] = pts[i].array() - c;
    for (int j = 0; j < NDIMS; j++) {
      for(int k = 0; k < NDIMS; k++) {
        covar(j, k) += ptsMinusCentroid[i][j]*ptsMinusCentroid[i][k];
      }
    }
  }

  // average the covariance matrix
  numerics::scalar_multiply(covar, static_cast< T >(1./n));

  // make room for the eigenvectors and eigenvalues
  T u[NDIMS*NDIMS];
  T lambdas[NDIMS];

  int eigen_res = numerics::eigen_solve< T >(covar, NDIMS, u, lambdas);
  SLIC_ASSERT(eigen_res);

  T maxima;
  T dot;

  Vector< T, NDIMS > w[NDIMS];
  Vector< T, NDIMS > e;

  for (int i = 0; i < NDIMS; i++) {
    w[i] = Vector< T, NDIMS >(u + NDIMS*i);

    // compute extent in this direction
    maxima = T();
    for (int j = 0; j < n; j++) {
      dot = T();
      for (int k = 0; k < NDIMS; k++) {
        dot += u[NDIMS*i + k]*((ptsMinusCentroid[j])[k]);
      }

      if (dot < T()) dot = -dot;

      if (maxima < dot) maxima = dot;
    }
    e[i] = maxima;
  }

  OrientedBoundingBox< T, NDIMS > res(Point< T, NDIMS >(c), w, e);

  // free up allocated memory
  delete [] ptsMinusCentroid;

  return res;
}

/*!
 *****************************************************************************
 * \brief Creates an oriented bounding box which contains the passed in OBBs.
 *
 * \param [in] l left obb
 * \param [in] r right obb
 *****************************************************************************
 */
template < typename T, int NDIMS >
OrientedBoundingBox< T, NDIMS > merge_boxes( const OrientedBoundingBox< T, NDIMS >
  &l, const OrientedBoundingBox< T, NDIMS > &r)
{
  if (l.contains(r)) {
    OrientedBoundingBox< T, NDIMS > res = l;
    return res;
  }
  if (r.contains(l)) {
    OrientedBoundingBox< T, NDIMS > res = r;
    return res;
  }


  std::vector< Point< T, NDIMS > > lv = l.vertices();
  std::vector< Point< T, NDIMS > > rv = r.vertices();

  int size = lv.size();

  Point< T, NDIMS > pts[2*size];

  for (int i = 0; i < size; i++) {
    pts[i] = lv[i];
    pts[i + size] = rv[i];
  }

  return compute_oriented_bounding_box< T, NDIMS >(pts, 2*size);
}

/*!
 *****************************************************************************
 * \brief Constructor. Creates a bounding box which contains the passed in
 * bounding boxes.
 *
 * \param [in] l left bb
 * \param [in] r right bb
 *****************************************************************************
 */
template < typename T, int NDIMS >
BoundingBox< T, NDIMS > merge_boxes( const BoundingBox< T, NDIMS >
  &l, const BoundingBox< T, NDIMS > &r)
{
  BoundingBox< T, NDIMS > res(l);
  res.addBox(r);
  return res;
}

}  /* namespace axom */
}  /* namespace primal */


#endif  /* COMPUTE_BOUNDING_BOX_HPP_ */
