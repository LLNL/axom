/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef AXOM_NUMERICS_EIGEN_SOLVE_HPP_
#define AXOM_NUMERICS_EIGEN_SOLVE_HPP_

#include "axom/Types.hpp" // for AXOM_NULLPTR

#include "axom_utils/Determinants.hpp" // for Determinants
#include "axom_utils/LU.hpp"           // for lu_decompose()/lu_solve()
#include "axom_utils/Matrix.hpp"       // for Matrix

// C/C++ includes
#include <cassert> // for assert()
#include <cmath> // for sqrt()
#include <random> // for std::srand
#include <ctime> // for time(0)

namespace axom {
namespace numerics {

/*!
 *******************************************************************************
 * \brief Approximates first k eigenvectors and eigenvalues of the passed in 
 * square PSD matrix.
 *
 * \param [in] A a square input matrix (must be PSD to be valid)
 * \param [in] k number of eigenvectors to find
 * \param [in] depth number of iterations for the power method
 * \param [out] u pointer to k eigenvectors in order by eigenvalue size
 * \param [out] lambdas pointer to k eigenvales in order by size
 * \return rc return value, 0 if the solve is successful.
 *
 * \pre A.isSquare() == true && A is PSD
 * \pre u != AXOM_NULLPTR
 * \pre lambdas != AXOM_NULLPTR
 *******************************************************************************
 */
template < typename T >
int eigen_solve( Matrix< T >& A, int k, int depth, T* u, T* lambdas );

} /* end namespace numerics */
} /* end namespace axom */

//------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------
namespace axom {
namespace numerics {

template < typename T >
int eigen_solve( Matrix< T >& A, int k, int depth, T* u, T* lambdas )
{
  // TODO: assert for PSD-ness?
  assert( "pre: input matrix must be square" && A.isSquare() );
  assert( "pre: eigenvectors pointer is null" && (u != AXOM_NULLPTR) );
  assert( "pre: lambdas vector is null" && (lambdas != AXOM_NULLPTR) );

  if (!A.isSquare()) {
    return LU_NONSQUARE_MATRIX;
  }

  if (k <= 0) return 0;

  static const double EPS = 1E-8;

  int N = A.getNumColumns();

  std::default_random_engine gen(time(0));
  std::normal_distribution< T > dist(static_cast< T >(0.0),
    static_cast< T >(1.0));


  for (int i = 0; i < k; i++) {
    // COMPUTING EIGENVECTORS AND VALUES VIA POWER METHOD:
    //   1: generate random vec
    //   2: make ortho to previous eigenvecs
    //   3: run depth iterations of power on it
    //   4: store

    // 1: generate random vec
    for (int j = 0; j < N; j++) {
      u[i*N + j] = dist(gen);
    }

    // 2: make ortho to previous eigenvecs
    for (int j = 0; j < i; j++) {
      // first compute projection onto this other eigenvec
      T dot = T();
      for (int l = 0; l < N; l++) dot += u[i*N + l]*u[j*N + l];

      // then subtract off the projection
      for (int l = 0; l < N; l++) u[i*N + l] -= dot*u[j*N + l];
    }

    // 3: run depth iterations of power method
    T temp[N];
    for (int j = 0; j < depth; j++) {
      // multiply
      vector_multiply(A, u + i*N, temp);

      // normalize
      double norm = 0.0;
      for (int l = 0; l < N; l++) norm += temp[l]*temp[l];
      T tnorm = static_cast< T >(std::sqrt(norm));

      if (norm < EPS) {
        for (int l = 0; l < N; l++) {
          u[i*N + l] = T();
        }
        break;
      }

      for (int l = 0; l < N; l++) u[i*N + l]= temp[l]/tnorm;
    }

    // 4: store the eigenval (already stored the eigenvec)
    vector_multiply(A, u + i*N, temp);

    lambdas[i] = T();
    for (int l = 0; l < N; l++) lambdas[i] += u[i*N + l]*temp[l];
  }


  return 0;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif
