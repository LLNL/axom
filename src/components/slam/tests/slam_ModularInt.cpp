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


/*
 * \file
 *
 * Unit tests for the modular arithmetic class ModularInt
 */

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "slam/SizePolicies.hpp"
#include "slam/ModularInt.hpp"

TEST(slam_modInt,runtime_modular_int_unitialized_and_full)
{
  typedef axom::slam::ModularInt<axom::slam::policies::RuntimeSize<int> >
    ModularIntType;

#ifdef AXOM_DEBUG
  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  SLIC_INFO("Checking that modular int with modulus zero fails.");
  SLIC_INFO("Note: Expecting a SLIC Failure: ");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntType(0,0),"")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntType(1,0),"")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntType(),   "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntType(1),  "")
    << " SIZE of Modular int not allowed to be zero";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("Checking modular int with value set to "
            <<" modulus equals (i.e. is equivalent to) 0 (runtime)");
  volatile int sz = 5;
  for(int i = 1 ; i< sz ; ++i)
  {
    ModularIntType modIntFull(i,i);
    EXPECT_EQ( modIntFull, 0);
  }
}

TEST(slam_modInt,compile_modular_int_unitialized_and_full)
{
  using namespace axom::slam;

#ifdef AXOM_DEBUG
  typedef ModularInt<policies::CompileTimeSize<int, 0> > ModularIntZero;

  // NOTE: AXOM_DEBUG is disabled in release mode, so this test will only fail
  // in debug mode
  SLIC_INFO("Checking that modular int with modulus "
            <<" zero fails.\nNote: Expecting a SLIC Failure: ");

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntZero(0,0),"")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntZero(1,0),"")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntZero(),   "")
    << " SIZE of Modular int not allowed to be zero";
  EXPECT_DEATH_IF_SUPPORTED(  ModularIntZero(1),  "")
    << " SIZE of Modular int not allowed to be zero";
#else
  SLIC_INFO("Skipped assertion failure check in release mode.");
#endif

  SLIC_INFO("Checking modular int with value set to modulus"
            <<" equals (i.e. is equivalent to) 0 for (compile time)");
  ModularInt<policies::CompileTimeSize<int, 1> > m1(1);
  EXPECT_EQ(  m1, 0);

  ModularInt<policies::CompileTimeSize<int, 2> > m2(2);
  EXPECT_EQ(  m2, 0);

  ModularInt<policies::CompileTimeSize<int, 3> > m3(3);
  EXPECT_EQ(  m3, 0);

  ModularInt<policies::CompileTimeSize<int, 4> > m4(4);
  EXPECT_EQ(  m4, 0);

}


TEST(slam_modInt,runtime_modular_int)
{
  SLIC_INFO("Checking modular int addition and subtraction"
            <<" when supplying the max value at runtime");

  typedef axom::slam::ModularInt<axom::slam::policies::RuntimeSize<int> >
    ModularIntType;

  volatile int sz = 937;

  for(int i = 0 ; i< sz ; ++i)
  {
    ModularIntType modInt(i,sz);
    EXPECT_EQ(  modInt, i);
    EXPECT_EQ(  modInt, modInt + sz);
    EXPECT_EQ(  modInt, modInt + 2 * sz);
    EXPECT_EQ(  modInt, modInt - sz);
  }

  ModularIntType modIntUp(0,sz);
  ModularIntType modIntDn(0,sz);
  const int loopEnd = 3 * modIntUp.modulus();
  for(int i = 0 ; i< loopEnd ; ++i)
  {
    EXPECT_EQ(  modIntUp, modIntUp + sz);
    EXPECT_EQ(  modIntUp, modIntUp + 2 * sz);
    EXPECT_EQ(  modIntUp, modIntUp - sz);

    EXPECT_EQ(  modIntDn, modIntDn + sz);
    EXPECT_EQ(  modIntDn, modIntDn + 2 * sz);
    EXPECT_EQ(  modIntDn, modIntDn - sz);

    ++modIntUp;
    --modIntDn;
  }
}

TEST(slam_modInt,runtime_modular_int_mult)
{
  typedef axom::slam::ModularInt<axom::slam::policies::RuntimeSize<int> >
    ModularIntType;

  volatile int sz = 10;

  SLIC_INFO("Checking modular int multiplication ");


  ModularIntType modInt5(5,sz);
  EXPECT_EQ(  modInt5 * 2,  0);

  ModularIntType modInt2(2,sz);
  EXPECT_EQ(  modInt2 * 5,  0);
  EXPECT_EQ(  modInt2 * 4,  8);
  EXPECT_EQ(  modInt2 * 6,  2);

  ModularIntType modInt3(3,sz);
  EXPECT_EQ(  modInt3 * 0,  0);
  EXPECT_EQ(  modInt3 * 1,  3);
  EXPECT_EQ(  modInt3 * 2,  6);
  EXPECT_EQ(  modInt3 * 3,  9);
  EXPECT_EQ(  modInt3 * 4,  2);

  ModularIntType modInt13(13,sz);
  EXPECT_EQ(  modInt13,     3);
  EXPECT_EQ(  modInt13 * 2, 6);

}

TEST(slam_modInt,compiletime_modular_int)
{
  SLIC_INFO("Checking modular int addition and subtraction"
            << " when supplying the max value at compile time");

  const int SZ = 937;

  typedef
    axom::slam::ModularInt<axom::slam::policies::CompileTimeSize<int,SZ> >
    ModularIntType;

  int sz = SZ;

  ModularIntType modIntZero(sz,sz);
  EXPECT_EQ( modIntZero, 0);


  for(int i = 0 ; i< sz ; ++i)
  {
    ModularIntType modInt(i,sz);
    EXPECT_EQ(  modInt, i);
    EXPECT_EQ(  modInt, modInt + sz);
    EXPECT_EQ(  modInt, modInt + 2 * sz);
    EXPECT_EQ(  modInt, modInt - sz);
  }

  ModularIntType modIntUp(0,sz);
  ModularIntType modIntDn(0,sz);
  const int loopEnd = 3 * modIntUp.modulus();
  for(int i = 0 ; i< loopEnd ; ++i)
  {
    EXPECT_EQ(  modIntUp, modIntUp + sz);
    EXPECT_EQ(  modIntUp, modIntUp + 2 * sz);
    EXPECT_EQ(  modIntUp, modIntUp - sz);

    EXPECT_EQ(  modIntDn, modIntDn + sz);
    EXPECT_EQ(  modIntDn, modIntDn + 2 * sz);
    EXPECT_EQ(  modIntDn, modIntDn - sz);

    ++modIntUp;
    --modIntDn;
  }
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
