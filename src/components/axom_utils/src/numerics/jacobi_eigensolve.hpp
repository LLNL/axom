/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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

#ifndef AXOM_JACOBI_EIGENSOLVE_HPP_
#define AXOM_JACOBI_EIGENSOLVE_HPP_

#include "axom/Types.hpp"              // for AXOM_NULLPTR
#include "axom/Macros.hpp"             // for AXOM_STATIC_ASSERT

#include "axom_utils/Matrix.hpp"       // for numerics::Matrix
#include "axom_utils/eigen_sort.hpp"   // for numerics::eigen_sort()
#include "axom_utils/Utilities.hpp"    // for abs(), isNearlyEqual()

namespace axom
{
namespace numerics
{

constexpr double JACOBI_DEFAULT_TOLERANCE   = 1.e-18;
constexpr int JACOBI_DEFAULT_MAX_ITERATIONS = 20;
constexpr int JACOBI_EIGENSOLVE_SUCCESS     = 0;
constexpr int JACOBI_EIGENSOLVE_FAILURE     = -1;


/*!
 * \brief Computes the eigenvalues and eigenvectors of a real symmetric matrix
 *  using the Jacobi iteration method.
 *
 * \param [in]  A the input matrix whose eigenpairs will be computed
 * \param [out] V matrix to store the eigenvectors of the input matrix
 * \param [out] lambdas buffer where the computed eigenvalues will be stored.
 * \param [in] maxIterations the maximum number of iterations (optional)
 * \param [in] numIterations the number of actual iterations (optional)
 * \param [in] TOL convergence tolerance. Default is set to 1.e-18 (optional)
 *
 * \return status returns JACOBI_EIGENSOLVE_SUCCESS on success, otherwise,
 *  JACOBI_EIGENSOLVE_FAILURE is returned.
 *
 * \note The supplied matrix is expected to be a real symmetric matrix.
 *
 * \note The Jacobi method is absolutely foolproof for all real symmetric
 *  matrices. It is efficient for moderate size matrices and generally can
 *  evaluate the smaller eigenvalues with better relative accuracy than other
 *  methods that convert the matrix to tridiagonal form, e.g., QR. However,
 *  for matrices of order greater than, e.g., 20, the Jacobi algorithm is
 *  generally slower by a significant constant factor.
 *
 * \note The Jacobi iteration will generally converge within 6-15 iterations.
 *  The default max number of iterations is set to 20, but, the caller may
 *  specify a different value if necessary.
 *
 * \pre A.isSquare() == true
 * \pre V.isSquare() == true
 * \pre lambdas != AXOM_NULLPTR
 */
template < typename T >
int jacobi_eigensolve( Matrix < T > A,
                       Matrix< T >& V,
                       T* lambdas,
                       int maxIterations=JACOBI_DEFAULT_MAX_ITERATIONS,
                       int* numIterations=AXOM_NULLPTR,
                       T TOL=JACOBI_DEFAULT_TOLERANCE );

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------
template < typename T >
int jacobi_eigensolve( Matrix < T > A,
                       Matrix< T >& V,
                       T* lambdas,
                       int maxIterations,
                       int* numIterations,
                       T TOL )
{
  int n = A.getNumRows();

  AXOM_STATIC_ASSERT_MSG(std::is_floating_point< T >::value,
                         "pre: T is a floating point type");
  assert("pre: input matrix must be square" && A.isSquare());
  assert("pre: can't have more eigenvectors than rows" &&
         (n <= A.getNumRows()));
  assert("pre: lambdas vector is null" && (lambdas != AXOM_NULLPTR));


  if ( !A.isSquare() || !V.isSquare() )
  {
    return JACOBI_EIGENSOLVE_FAILURE;
  }


  T* bw = new T[n]; // temp array
  T* zw = new T[n]; // temp array

  // initialize
  for (int i = 0 ; i < n ; ++i)
  {
    lambdas[ i ] = A(i,i);
    bw[ i ]      = lambdas[i];
    zw[ i ]      = 0.0;

    for (int j = 0 ; j < n ; ++j)
    {
      V( i,j ) = ( i==j ) ? 1.0 : 0.0;
    }

  }

  // Jacobi solve
  for (int iter = 0 ; iter < maxIterations ; ++iter)
  {
    // compute the sum of all elements in the upper triangular portion of A.
    // The convergence criterion (thresh) is based on the absolute value of
    // this sum.
    T sum = 0.0;

    for (int i = 0 ; i < n ; ++i)
    {
      for (int j = i+1 ; j < n ; ++j)
      {
        sum += A(i,j) * A(i,j);
      }
    }

    const T thresh = sqrt(sum) / (4.0*static_cast< T >( n ) );

    if (utilities::isNearlyEqual( thresh, 0.0, TOL ))
    {
      break;
    }

    // if not converged, then perform next iteration.
    for (int p = 0 ; p < n ; ++p)
    {
      for (int q = p+1 ; q < n ; ++q)
      {

        T gapq  = 10.0 * utilities::abs(A(p,q));
        T termp = gapq + utilities::abs(lambdas[p]);
        T termq = gapq + utilities::abs(lambdas[q]);

        // the Jacobi iteration ignores off diagonal elements close to zero
        if ( 4 < iter &&
             termp == utilities::abs(lambdas[p]) &&
             termq == utilities::abs(lambdas[q]) )
        {
          A(p,q) = 0.0;
        }
        else if ( thresh <= utilities::abs(A(p,q)) )
        {

          T h = lambdas[q] - lambdas[p];
          T term = utilities::abs(h) + gapq;

          T t;

          if (term == utilities::abs(h))
          {
            t = A(p,q)/h;
          }
          else
          {
            T theta = 0.5 * h / A(p,q);
            t = 1.0 / ( utilities::abs(theta) + sqrt(1.0 + theta*theta) );

            if (theta < 0.0)
            {
              t = -t;
            }

          }

          // compute Givens rotation terms: c = cos(theta), s = sin(theta)
          T c = 1.0 / sqrt(1.0 + t*t);
          T s = t * c;
          T tau = s / (1.0 + c);
          h = t * A(p,q);

          // accumulate corrections to diagonals
          zw[p] += -h;
          zw[q] += h;
          lambdas[p] += -h;
          lambdas[q] += h;

          A(p,q) = 0.0;

          // perform the rotation using information from upper triangle of A
          for (int j = 0 ; j < p ; ++j)
          {
            T g1   = A(j,p);
            T g2   = A(j,q);
            A(j,p) = g1 - s * (g2 + g1 * tau);
            A(j,q) = g2 + s * (g1 - g2 * tau);
          }

          for (int j = p+1 ; j < q ; ++j)
          {
            T g1   = A(p,j);
            T g2   = A(j,q);
            A(p,j) = g1 - s * (g2 + g1 * tau);
            A(j,q) = g2 + s * (g1 - g2 * tau);
          }

          for (int j = q+1 ; j < n ; ++j)
          {
            T g1   = A(p,j);
            T g2   = A(q,j);
            A(p,j) = g1 - s * (g2 + g1 * tau);
            A(q,j) = g2 + s * (g1 - g2 * tau);
          }

          // accumulate results into eigenvector matrix
          for (int j = 0 ; j < n ; ++j)
          {
            T g1   = V(j,p);
            T g2   = V(j,q);
            V(j,p) = g1 - s * (g2 + g1 * tau);
            V(j,q) = g2 + s * (g1 - g2 * tau);
          }

        } // END else if
      } // END for all q
    } // END for all p

    for (int i = 0 ; i < n ; ++i)
    {
      bw[i] += zw[i];
      lambdas[i] = bw[i];
      zw[i] = 0.0;
    }

    // update number of actual iterations (if specified)
    if ( numIterations != AXOM_NULLPTR )
    {
      (*numIterations)++;
    }

  } // END for all iterations

  // sort eigenvalues in ascending order
  eigen_sort( lambdas, V );

  delete[] bw;
  delete[] zw;

  return JACOBI_EIGENSOLVE_SUCCESS;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif /* AXOM_JACOBI_EIGENSOLVE_HPP_ */
