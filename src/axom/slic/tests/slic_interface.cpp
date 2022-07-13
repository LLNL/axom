// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

TEST(slic_interface, initialize_finalize)
{
  EXPECT_FALSE(axom::slic::isInitialized());

  // Initialize and finalize slic
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());

  // Check that we can Initialize and finalize slic a second time
  axom::slic::initialize();
  EXPECT_TRUE(axom::slic::isInitialized());

  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

TEST(slic_interface, logging_level)
{
  axom::slic::initialize();

  // test get and set for each level
  for(int i = 0; i < axom::slic::message::Num_Levels; ++i)
  {
    axom::slic::message::Level lev = static_cast<axom::slic::message::Level>(i);

    axom::slic::setLoggingMsgLevel(lev);
    EXPECT_EQ(lev, axom::slic::getLoggingMsgLevel());
  }

  axom::slic::finalize();
}

TEST(slic_interface, aborting)
{
  axom::slic::initialize();
  axom::slic::setLoggingMsgLevel(axom::slic::message::Debug);
  axom::slic::addStreamToAllMsgLevels(
    new axom::slic::GenericOutputStream(&std::cout));

  // test abort on info/debug (no-ops)
  axom::slic::setAbortFlag(true, axom::slic::message::Info);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Info));

  axom::slic::setAbortFlag(true, axom::slic::message::Debug);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Debug));

  // test abort on errors
  EXPECT_TRUE(axom::slic::isAbortOnErrorsEnabled());

  axom::slic::disableAbortOnError();
  EXPECT_FALSE(axom::slic::isAbortOnErrorsEnabled());

  axom::slic::setAbortFlag(true, axom::slic::message::Error);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Error));

  axom::slic::enableAbortOnError();
  EXPECT_TRUE(axom::slic::isAbortOnErrorsEnabled());

  axom::slic::setAbortFlag(true, axom::slic::message::Error);
  EXPECT_TRUE(axom::slic::getAbortFlag(axom::slic::message::Error));

  axom::slic::setAbortFlag(false, axom::slic::message::Error);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Error));

  axom::slic::enableAbortOnError();
  axom::slic::setAbortFlag(true, axom::slic::message::Error);
  axom::slic::disableAbortOnError();
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Error));

  axom::slic::enableAbortOnError();
  axom::slic::setAbortFlag(true, axom::slic::message::Error);
  axom::slic::setAbortOnError(false);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Error));

  // test abort on warnings
  EXPECT_FALSE(axom::slic::isAbortOnWarningsEnabled());

  axom::slic::disableAbortOnWarning();
  axom::slic::setAbortFlag(true, axom::slic::message::Warning);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Warning));

  axom::slic::enableAbortOnWarning();
  EXPECT_TRUE(axom::slic::isAbortOnWarningsEnabled());

  axom::slic::setAbortFlag(true, axom::slic::message::Warning);
  EXPECT_TRUE(axom::slic::getAbortFlag(axom::slic::message::Warning));

  axom::slic::setAbortFlag(false, axom::slic::message::Warning);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Warning));

  axom::slic::enableAbortOnWarning();
  axom::slic::setAbortFlag(true, axom::slic::message::Warning);
  axom::slic::disableAbortOnWarning();
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Warning));

  axom::slic::enableAbortOnWarning();
  axom::slic::setAbortFlag(true, axom::slic::message::Warning);
  axom::slic::setAbortOnWarning(false);
  EXPECT_FALSE(axom::slic::getAbortFlag(axom::slic::message::Warning));

  axom::slic::finalize();
}