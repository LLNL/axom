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

#include "axom_utils/vector_utilities.hpp" // for Determinants
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
 * \brief Approximates k eigenvectors and eigenvalues of the passed in square
 * matrix using the power method.
 * 
 * We compute eigenvectors and eigenvalues from the power method as follows:
 *    1. Generate a random vector
 *    2. Make ortho to previous eigenvectors if any
 *    3. Run depth iterations of power method
 *    4. Store the result
 *
 * \param [in] A a square input matrix
 * \param [in] k number of eigenvalue-eigenvectors to find
 * \param [in] depth number of iterations for the power method
 * \param [out] u pointer to k eigenvectors in order by magnitude of eigenvalue
 * \param [out] lambdas pointer to k eigenvales in order by size
 * \return rc return value, nonzero if the solve is successful.
 *
 * \pre A.isSquare() == true
 * \pre u != AXOM_NULLPTR
 * \pre lambdas != AXOM_NULLPTR
 *******************************************************************************
 */
template < typename T >
int eigen_solve(Matrix< T >& A, int k, T* u, T* lambdas, int depth=20);


} /* end namespace numerics */
} /* end namespace axom */


//------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------
namespace axom {
namespace numerics {

namespace { /* anonymous namespace */

  // TODO: make this work with non-CXX11 compilers
  template < typename T >
  T getRandom() {
    static std::default_random_engine gen(time(0));
    static std::normal_distribution< T > dist(static_cast< T >(0.0),
      static_cast< T >(1.0));
    return dist(gen);
  }

} /* end namespace anonymous */

template < typename T >
int eigen_solve(Matrix< T >& A, int k, T* u, T* lambdas, int depth)
{
  assert("pre: input matrix must be square" && A.isSquare());
  assert("pre: eigenvectors pointer is null" && (u != AXOM_NULLPTR));
  assert("pre: lambdas vector is null" && (lambdas != AXOM_NULLPTR));

  if (!A.isSquare()) {
    return 0;
  }

  if (k <= 0) return 1;

  int N = A.getNumColumns();

  // allocate memory for a temp var
  T *temp = new T[N];

  for (int i = 0; i < k; i++) {

    T *vec = &u[i*N];

    // 1: generate random vec
    for (int j = 0; j < N; j++) vec[j] = getRandom< T >();

    // 2: make ortho to previous eigenvecs then normalize
    for (int j = 0; j < i; j++) make_orthogonal< T >(vec, u + j*N, N);

    int res = normalize< T >(vec, N);

    if (!res) {  // something went wrong
      return 0;
    }

    // 3: run depth iterations of power method; note that a loop invariant
    // of this is that vec is normalized
    for (int j = 0; j < depth; j++) {
      // multiply
      vector_multiply(A, vec, temp);

      // make ortho to previous (for stability)
      for (int k = 0; k < i; k++) make_orthogonal< T >(temp, u + k*N, N);

      res = normalize< T >(temp, N);

      if (!res) {  // must be 0 eigenvalue; done in that case, since vec
        // is guaranteed to be orthogonal to previous eigenvecs and normal
        break;
      } else {  // else copy it over
        for (int l = 0; l < N; l++) vec[l] = temp[l];
      }
    }

    // 4: store the eigenval (already stored the eigenvec)
    vector_multiply(A, vec, temp);

    lambdas[i] = dot_product(vec, temp, N);
  }

  // free up allocated memory
  delete [] temp;

  return 1;
}

} /* end namespace numerics */
} /* end namespace axom */

#endif
