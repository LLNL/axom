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

#include "gtest/gtest.h"

#include "axom_utils/vector_utilities.hpp"
#include "axom_utils/Utilities.hpp"

namespace numerics  = axom::numerics;
namespace utilities = axom::utilities;

namespace
{

/*!
 * \brief Checks if the given two vectors are equal.
 *
 * \param [in] u rhs vector to compare
 * \param [in] v lhs vector to compare
 * \param [in] N the size of the vector
 */
void expect_vector_eq( const double* u, const double* v, int N)
{
  for ( int i=0 ; i < N ; ++i )
  {
    EXPECT_DOUBLE_EQ( u[i], v[i] );
  }
}

} /* end anonymous namespace */

//------------------------------------------------------------------------------
// UNIT TESTS
//------------------------------------------------------------------------------

TEST( numerics_vector_utilities, linspace_test )
{
  double vc[ 6 ]; // storage for the computed vector

  const double v1[ 6 ] = { -2, -1, 0, 1, 2, 3   };
  const double v2[ 6 ] = {  3,  2, 1, 0, -1, -2 };
  const double v3[ 3 ] = { 0.0, 0.5, 1.0 };

  const double x0 = -2.0;
  const double x1 =  3.0;

  EXPECT_FALSE( numerics::linspace( x0, x1, vc, 0) );
  EXPECT_FALSE( numerics::linspace( x0, x1, vc, 1) );

  EXPECT_TRUE( numerics::linspace( x0, x1, vc, 6 ) );
  expect_vector_eq( vc, v1, 6 );

  EXPECT_TRUE( numerics::linspace( x1, x0, vc, 6) );
  expect_vector_eq( vc, v2, 6 );

  EXPECT_TRUE( numerics::linspace( 0.0, 1.0, vc, 3 ) );
  expect_vector_eq( vc, v3, 3 );
}

//------------------------------------------------------------------------------
TEST( numerics_vector_utilities, cross_product_test )
{
  const int NDIMS = 3;

  const double e1[ 3  ] = { 1.0, 0.0, 0.0  };
  const double e2[ 3  ] = { 0.0, 1.0, 0.0  };
  const double e3[ 3  ] = { 0.0, 0.0, 1.0  };
  const double me3[ 3 ] = { 0.0, 0.0, -1.0 };

  const double u[ 3 ]     = { 2.0,  1.0, -1.0  };
  const double v[ 3 ]     = {-3.0,  4.0,  1.0  };
  const double u_x_v[ 3 ] = { 5.0,  1.0,  11.0 };

  double w[3];

  numerics::cross_product( e1, e2, w);
  expect_vector_eq( w, e3, NDIMS );

  numerics::cross_product( e2, e1, w );
  expect_vector_eq( w, me3, NDIMS );

  numerics::cross_product( e3, e1, w );
  expect_vector_eq( w, e2, NDIMS );

  numerics::cross_product( e2, e3, w );
  expect_vector_eq( w, e1, NDIMS );

  numerics::cross_product( u, v, w );
  expect_vector_eq( w, u_x_v, NDIMS );
}

//------------------------------------------------------------------------------
TEST( numerics_vector_utilities, dot_product_test )
{
  const int dim = 10;
  const double EPS = 1E-4;

  double* u = new double[dim]();   // 0 initialization
  double* v = new double[dim]();   // 0 initialization

  u[0] = 1.;
  v[1] = 1.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), 0., EPS);

  u[1] = 1.;
  v[0] = 1.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), 2., EPS);

  u[5] = 10.;
  v[5] = -10.;

  EXPECT_NEAR(numerics::dot_product< double >(u, v, dim), -98., EPS);

  // free up allocated memory
  delete [] u;
  delete [] v;
}

//------------------------------------------------------------------------------
TEST( numerics_vector_utilities, make_orthogonal_test )
{
  const int dim = 3;
  const double ONE_THIRD = 0.3333333333;
  const double EPS = 1E-4;

  double* u = new double[dim]();   // 0 initialization
  double* v = new double[dim]();   // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = -1.;

  v[0] = 1.;
  v[1] = 1.;
  v[2] = 1.;

  numerics::make_orthogonal< double >(u, v, dim);

  EXPECT_NEAR(u[0], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[1], 1. - ONE_THIRD, EPS);
  EXPECT_NEAR(u[2], -1. - ONE_THIRD, EPS);

  numerics::make_orthogonal(u, u, dim);

  EXPECT_NEAR(u[0], 0., EPS);
  EXPECT_NEAR(u[1], 0., EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  // free up allocated memory
  delete [] u;
  delete [] v;
}

//------------------------------------------------------------------------------
TEST( numerics_vector_utilities, orthonormalize_test )
{
  const int dim = 3;
  const int size = 2;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_TWO = 0.70710678;

  double* basis = new double[dim*dim]();   // 0 initialization

  double* u = basis;
  double* v = basis + dim;

  u[0] = 1.;
  u[1] = 1.;
  v[0] = 2.;

  EXPECT_TRUE(numerics::orthonormalize< double >(basis, size, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(u[2], 0., EPS);

  EXPECT_NEAR(v[0], ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[1], -ONE_OVER_SQRT_TWO, EPS);
  EXPECT_NEAR(v[2], 0., EPS);

  // doesn't span
  EXPECT_FALSE(numerics::orthonormalize< double >(basis, dim, dim));

  // free up allocated memory
  delete [] basis;
}

//------------------------------------------------------------------------------
TEST( numerics_vector_utilities, normalize_test )
{
  const int dim = 3;
  const double EPS = 1E-4;
  const double ONE_OVER_SQRT_THREE = 0.57735026;
  const double ONE_OVER_SQRT_TEN = 0.316227766;

  double* u = new double[dim]();   // 0 initialization
  double* v = new double[dim]();   // 0 initialization

  u[0] = 1.;
  u[1] = 1.;
  u[2] = 1.;

  v[0] = 0.;
  v[1] = 3.;
  v[2] = 1.;

  EXPECT_TRUE(numerics::normalize< double >(u, dim));

  EXPECT_NEAR(u[0], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[1], ONE_OVER_SQRT_THREE, EPS);
  EXPECT_NEAR(u[2], ONE_OVER_SQRT_THREE, EPS);

  EXPECT_TRUE(numerics::normalize< double >(v, dim));

  EXPECT_NEAR(v[0], 0., EPS);
  EXPECT_NEAR(v[1], 3.*ONE_OVER_SQRT_TEN, EPS);
  EXPECT_NEAR(v[2], ONE_OVER_SQRT_TEN, EPS);

  // now make u = (0, 0, 0)
  numerics::make_orthogonal(u, u, dim);

  // normalize should fail
  EXPECT_FALSE(numerics::normalize< double >(u, dim));
}
