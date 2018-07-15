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

#ifndef AXOM_EIGEN_SORT_HPP_
#define AXOM_EIGEN_SORT_HPP_

#include "axom/Types.hpp"  // for AXOM_NULLPTR

namespace axom
{
namespace numerics
{

/*! Sorts eigenvalues in ascending order.
 */
template < typename T >
void eigen_sort( T* lamdas, Matrix< T >& eigen_vectors );

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------
template < typename T >
void eigen_sort( T* lambdas, Matrix< T >& eigen_vectors )
{
  assert( "pre: lambdas vector is null" && (lambdas != AXOM_NULLPTR) );

  int n = eigen_vectors.getNumRows();

  for (int i = 0; i < n-1; ++i)
  {
    int m = i;

    for (int j = i+1; j < n; ++j)
    {
      if (lambdas[j] < lambdas[m])
      {
        m = j;
      }
    }

    if (m != i)
    {
      double t   = lambdas[m];
      lambdas[m] = lambdas[i];
      lambdas[i] = t;

      eigen_vectors.swapColumns(m, i);
    }

  }

}

} /* end namespace numerics */
} /* end namespace axom */

#endif /* AXOM_EIGEN_SORT_HPP_ */
