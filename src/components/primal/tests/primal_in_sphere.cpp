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

#include "axom/config.hpp"

#include "slic/slic.hpp"

#include "primal/Point.hpp"
#include "primal/in_sphere.hpp"

using namespace axom;

TEST(primal_in_sphere, test_in_sphere_2d)
{
  const int DIM = 2;
  typedef primal::Point< double,DIM >     PointType;
  PointType p0 = PointType::make_point(0,0);
  PointType p1 = PointType::make_point(1,0);
  PointType p2 = PointType::make_point(0,1);

  PointType q_inside =   PointType::make_point(0.1,0.1);
  PointType q_on_circle =  PointType::make_point(1,1);
  PointType q_outside =   PointType::make_point(1.1, 0);

  EXPECT_TRUE( in_sphere( q_inside, p0, p1, p2) );
  EXPECT_FALSE( in_sphere( q_on_circle, p0, p1, p2) );
  EXPECT_FALSE( in_sphere( q_outside, p0, p1, p2) );
}


TEST(primal_in_sphere, test_in_sphere_3d)
{
  const int DIM = 3;
  typedef primal::Point< double,DIM >     PointType;
  PointType p0 = PointType::make_point(-1,-1,1);
  PointType p1 = PointType::make_point(1,-1,-1);
  PointType p2 = PointType::make_point(-1,1,-1);
  PointType p3 = PointType::make_point(1,1,-1);

  PointType q_inside =   PointType::make_point(0.5,0.5,0.5);
  PointType q_origin =   PointType::make_point(0.,0.,0.);
  PointType q_on_circle =  PointType::make_point(1,1,1);
  PointType q_outside =   PointType::make_point(1.1,1,1);
  PointType q_outside2 =   PointType::make_point(-1.1,1,1);

  EXPECT_TRUE( in_sphere( q_inside, p0, p1, p2, p3) );
  EXPECT_TRUE( in_sphere( q_origin, p0, p1, p2, p3) );
  EXPECT_FALSE( in_sphere( q_on_circle, p0, p1, p2, p3) );
  EXPECT_FALSE( in_sphere( q_outside, p0, p1, p2, p3) );
  EXPECT_FALSE( in_sphere( q_outside2, p0, p1, p2, p3) );
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  axom::slic::setLoggingMsgLevel( axom::slic::message::Warning);

  int result = RUN_ALL_TESTS();
  return result;
}
