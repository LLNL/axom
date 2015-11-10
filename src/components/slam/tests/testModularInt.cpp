/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * \file
 *
 * Unit tests for the modular arithmetic class ModularInt
 */

#include "gtest/gtest.h"

#include "meshapi/SizePolicies.hpp"
#include "meshapi/ModularInt.hpp"

TEST(gtest_meshapi_modInt,runtime_modular_int_unitialized_and_full)
{
    typedef asctoolkit::meshapi::ModularInt<asctoolkit::meshapi::policies::RuntimeSizeHolder<int> > ModularIntType;

#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode
  std::cout << "\n -- Checking that modular int with modulus zero fails.\nNote: Expecting a SLIC Failure: " << std::endl;

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH( ModularIntType(0,0),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntType(1,0),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntType(),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntType(1),"") << " SIZE of Modular int not allowed to be zero";
#else
  std::cout << "Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
#endif

  std::cout << "\n -- Checking modular int with value set to modulus equals (i.e. is equivalent to) 0 (runtime)" << std::endl;
  volatile int sz = 5;
  for(int i= 1; i< sz; ++i)
  {
      ModularIntType modIntFull(i,i);
      EXPECT_EQ( modIntFull, 0);
  }
}

TEST(gtest_meshapi_modInt,compile_modular_int_unitialized_and_full)
{
    using namespace asctoolkit::meshapi;
    typedef ModularInt<policies::CompileTimeSizeHolder<int, 0> > ModularIntZero;

#ifdef ATK_DEBUG
  // NOTE: ATK_ASSSERT is disabled in release mode, so this test will only fail in debug mode
  std::cout << "\n -- Checking that modular int with modulus zero fails.\nNote: Expecting a SLIC Failure: " << std::endl;

  // add this line to avoid a warning in the output about thread safety
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH( ModularIntZero(0,0),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntZero(1,0),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntZero(),"") << " SIZE of Modular int not allowed to be zero";
  ASSERT_DEATH( ModularIntZero(1),"") << " SIZE of Modular int not allowed to be zero";
#else
  std::cout << "Did not check for assertion failure since assertions are compiled out in release mode." << std::endl;
#endif

  std::cout << "\n -- Checking modular int with value set to modulus equals (i.e. is equivalent to) 0 for (compile time)" << std::endl;
  ModularInt<policies::CompileTimeSizeHolder<int, 1> > m1(1);
  EXPECT_EQ( m1, 0);

  ModularInt<policies::CompileTimeSizeHolder<int, 2> > m2(2);
  EXPECT_EQ( m2, 0);

  ModularInt<policies::CompileTimeSizeHolder<int, 3> > m3(3);
  EXPECT_EQ( m3, 0);

  ModularInt<policies::CompileTimeSizeHolder<int, 4> > m4(4);
  EXPECT_EQ( m4, 0);

}


TEST(gtest_meshapi_modInt,runtime_modular_int)
{
  std::cout << "\n -- Checking modular int addition and subtraction when supplying the max value at runtime" << std::endl;

  typedef asctoolkit::meshapi::ModularInt<asctoolkit::meshapi::policies::RuntimeSizeHolder<int> > ModularIntType;

  volatile int sz = 937;

  for(int i= 0; i< sz; ++i)
  {
      ModularIntType modInt(i,sz);
      EXPECT_EQ( modInt, i);
      EXPECT_EQ( modInt, modInt + sz);
      EXPECT_EQ( modInt, modInt + 2*sz);
      EXPECT_EQ( modInt, modInt - sz);
  }

  ModularIntType modIntUp(0,sz);
  ModularIntType modIntDn(0,sz);
  const int loopEnd = 3 * modIntUp.modulus();
  for(int i= 0; i< loopEnd; ++i)
  {
      EXPECT_EQ( modIntUp, modIntUp + sz);
      EXPECT_EQ( modIntUp, modIntUp + 2*sz);
      EXPECT_EQ( modIntUp, modIntUp - sz);

      EXPECT_EQ( modIntDn, modIntDn + sz);
      EXPECT_EQ( modIntDn, modIntDn + 2*sz);
      EXPECT_EQ( modIntDn, modIntDn - sz);

      ++modIntUp;
      --modIntDn;
  }
}

TEST(gtest_meshapi_modInt,runtime_modular_int_mult)
{
  typedef asctoolkit::meshapi::ModularInt<asctoolkit::meshapi::policies::RuntimeSizeHolder<int> > ModularIntType;

  volatile int sz = 10;

  std::cout << "\n -- Checking modular int multiplication " << std::endl;


  ModularIntType modInt5(5,sz);
  EXPECT_EQ( modInt5 * 2, 0);

  ModularIntType modInt2(2,sz);
  EXPECT_EQ( modInt2 * 5, 0);
  EXPECT_EQ( modInt2 * 4, 8);
  EXPECT_EQ( modInt2 * 6, 2);

  ModularIntType modInt3(3,sz);
  EXPECT_EQ( modInt3 * 0, 0);
  EXPECT_EQ( modInt3 * 1, 3);
  EXPECT_EQ( modInt3 * 2, 6);
  EXPECT_EQ( modInt3 * 3, 9);
  EXPECT_EQ( modInt3 * 4, 2);

  ModularIntType modInt13(13,sz);
  EXPECT_EQ( modInt13, 3);
  EXPECT_EQ( modInt13*2, 6);

}

TEST(gtest_meshapi_modInt,compiletime_modular_int)
{
    std::cout << "\n -- Checking modular int addition and subtraction when supplying the max value at compile time" << std::endl;

    const int SZ = 937;
    typedef asctoolkit::meshapi::ModularInt<asctoolkit::meshapi::policies::CompileTimeSizeHolder<int, SZ> > ModularIntType;

    int sz = SZ;

    ModularIntType modIntZero(sz,sz);
    EXPECT_EQ( modIntZero, 0);


    for(int i= 0; i< sz; ++i)
    {
        ModularIntType modInt(i,sz);
        EXPECT_EQ( modInt, i);
        EXPECT_EQ( modInt, modInt + sz);
        EXPECT_EQ( modInt, modInt + 2*sz);
        EXPECT_EQ( modInt, modInt - sz);
    }

    ModularIntType modIntUp(0,sz);
    ModularIntType modIntDn(0,sz);
    const int loopEnd = 3 * modIntUp.modulus();
    for(int i= 0; i< loopEnd; ++i)
    {
        EXPECT_EQ( modIntUp, modIntUp + sz);
        EXPECT_EQ( modIntUp, modIntUp + 2*sz);
        EXPECT_EQ( modIntUp, modIntUp - sz);

        EXPECT_EQ( modIntDn, modIntDn + sz);
        EXPECT_EQ( modIntDn, modIntDn + 2*sz);
        EXPECT_EQ( modIntDn, modIntDn - sz);

        ++modIntUp;
        --modIntDn;
    }
}



//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
 int result = 0;

 ::testing::InitGoogleTest(&argc, argv);

 UnitTestLogger logger;   // create & initialize test logger,
 // finalized when exiting main scope

 result = RUN_ALL_TESTS();

 return result;
}
