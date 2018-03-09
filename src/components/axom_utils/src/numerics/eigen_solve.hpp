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

#ifndef AXOM_NUMERICS_EIGEN_SOLVE_HPP_
#define AXOM_NUMERICS_EIGEN_SOLVE_HPP_

#include "axom/Types.hpp" // for AXOM_NULLPTR

#include "axom_utils/vector_utilities.hpp" // for Determinants
#include "axom_utils/Determinants.hpp" // for Determinants
#include "axom_utils/LU.hpp"           // for lu_decompose()/lu_solve()
#include "axom_utils/Matrix.hpp"       // for Matrix

// C/C++ includes
#include <cassert> // for assert()
#include <cmath> // for sqrt() and rand()
#include <ctime> // for time(0)

namespace axom
{
namespace numerics
{

/*!
 * \brief Approximates k eigenvectors and eigenvalues of the passed in square
 * matrix using the power method.
 *
 * We compute eigenvectors and eigenvalues from the power method as follows:
 *    1. Generate a random vector
 *    2. Make orthonormal to previous eigenvectors if any
 *    3. Run depth iterations of power method
 *    4. Store the result
 *
 * Because step 1 and step 3 require random numbers, this function requires
 * the caller to properly seed the randomizer (for example, by calling
 * srand()).
 *
 * \param [in] A a square input matrix
 * \param [in] k number of eigenvalue-eigenvectors to find
 * \param [out] u pointer to k eigenvectors in order by magnitude of eigenvalue
 * \param [out] lambdas pointer to k eigenvales in order by size
 * \param [in] numIterations optional number of iterations for the power method
 * \note if k <= 0, the solve is declared successful
 * \return rc return value, nonzero if the solve is successful.
 *
 * \pre A.isSquare() == true
 * \pre u != AXOM_NULLPTR
 * \pre lambdas != AXOM_NULLPTR
 * \pre k <= A.getNumRows()
 * \pre T is a floating point type
 * \pre System randomizer has been initialized (for example, with srand())
 */
template < typename T >
int eigen_solve(Matrix< T >& A, int k, T* u, T* lambdas,
                int numIterations=125);

} /* end namespace numerics */
} /* end namespace axom */


//------------------------------------------------------------------------------
// Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace numerics
{

namespace
{

/*!
 * \brief Returns a uniformly distributed random between 0 and 1.
 *
 * Caller must seed the random number generator (for example, with srand()).
 */
template < typename T >
T getRandom()
{
  return static_cast< T >(((double) rand())/RAND_MAX);
}
} /* end anonymous namespace */

template < typename T >
int eigen_solve(Matrix< T >& A, int k, T* u, T* lambdas, int numIterations)
{
  AXOM_STATIC_ASSERT_MSG(std::is_floating_point< T >::value,
                         "pre: T is a floating point type");
  assert("pre: input matrix must be square" && A.isSquare());
  assert("pre: can't have more eigenvectors than rows" &&
         (k <= A.getNumRows()));
  assert("pre: eigenvectors pointer is null" && (u != AXOM_NULLPTR));
  assert("pre: lambdas vector is null" && (lambdas != AXOM_NULLPTR));

  if (!A.isSquare())
  {
    return 0;
  }

  if (k <= 0)
  {
    return 1;
  }

  int N = A.getNumColumns();

  // allocate memory for a temp var
  T* temp = new T[N];

  for (int i = 0 ; i < k ; i++)
  {

    T* vec = &u[i*N];
    // 1: generate random vec
    for (int j = 0 ; j < N ; j++)
    {
      vec[j] = getRandom< T >();
    }

    // 2: make ortho to previous eigenvecs then normalize
    for (int j = 0 ; j < i ; j++)
    {
      make_orthogonal< T >(vec, u + j*N, N);
    }

    bool res = normalize< T >(vec, N);

    if (!res)  // something went wrong
    {
      return 0;
    }

    // 3: run depth iterations of power method; note that a loop invariant
    // of this is that vec is normalized
    for (int j = 0 ; j < numIterations ; j++)
    {
      // multiply
      vector_multiply(A, vec, temp);

      // make ortho to previous (for stability)
      for (int k = 0 ; k < i ; k++)
      {
        make_orthogonal< T >(temp, u + k*N, N);
      }

      res = normalize< T >(temp, N);

      if (!res)  // must be 0 eigenvalue; done in that case, since vec
      {// is guaranteed to be orthogonal to previous eigenvecs and normal
        break;
      }
      else    // else copy it over
      {
        for (int l = 0 ; l < N ; l++)
        {
          vec[l] = temp[l];
        }
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
