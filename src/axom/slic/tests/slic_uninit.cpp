// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "axom/slic/interface/slic.hpp"
#include "axom/slic/interface/slic_macros.hpp"
#include "axom/slic/streams/GenericOutputStream.hpp"

void testInit(const std::string& label, std::function<void(void)> f)
{
  SCOPED_TRACE(label);
  EXPECT_FALSE(axom::slic::isInitialized());
  f();
  EXPECT_TRUE(axom::slic::isInitialized());
  axom::slic::finalize();
  EXPECT_FALSE(axom::slic::isInitialized());
}

TEST(slic_uninit, log_macro)
{
#ifdef AXOM_DEBUG
  testInit("debug message", []() { SLIC_DEBUG("test debug message"); });
  testInit("debug_if",
           []() { SLIC_DEBUG_IF(true, "debug message should print"); });
#endif
  testInit("info message", []() { SLIC_INFO("test info message"); });
  testInit("info_if", []() { SLIC_INFO_IF(true, "info message should print"); });
  testInit("warning message", []() { SLIC_WARNING("test warning message"); });
  testInit("warning_if", []() { SLIC_WARNING_IF(true, "test warning message"); });
  testInit("error message", []() { SLIC_ERROR("test error message"); });
  testInit("error_if", []() { SLIC_ERROR_IF(true, "test error message"); });
}

TEST(slic_uninit, tests_and_asserts)
{
#ifdef AXOM_DEBUG
  testInit("check", []() { SLIC_CHECK(false); });
  testInit("check_msg",
           []() { SLIC_CHECK_MSG(false, "checked false; this should print"); });
  testInit("assert", []() { SLIC_ASSERT(false); });
  testInit("assert_msg", []() {
    SLIC_ASSERT_MSG(false, "asserted false; this should print");
  });
#else
  EXPECT_TRUE(true);
#endif
}

TEST(slic_uninit, manage_loggers)
{
  testInit("createLogger", []() { axom::slic::createLogger("foo"); });
  testInit("activateLogger", []() { (void)axom::slic::activateLogger("foo"); });
  testInit("getActiveLoggerName",
           []() { (void)axom::slic::getActiveLoggerName(); });
  testInit("setLoggingMsgLevel",
           []() { axom::slic::setLoggingMsgLevel(axom::slic::message::Info); });
  testInit("getLoggingMsgLevel", []() { axom::slic::getLoggingMsgLevel(); });

  testInit("setAbortOnError", []() { axom::slic::setAbortOnError(true); });
  testInit("enableAbortOnError", []() { axom::slic::enableAbortOnError(); });
  testInit("disableAbortOnError", []() { axom::slic::disableAbortOnError(); });
  testInit("isAbortOnErrorsEnabled",
           []() { (void)axom::slic::isAbortOnErrorsEnabled(); });

  testInit("setAbortOnWarning", []() { axom::slic::setAbortOnWarning(true); });
  testInit("enableAbortOnWarning", []() { axom::slic::enableAbortOnWarning(); });
  testInit("disableAbortOnWarning",
           []() { axom::slic::disableAbortOnWarning(); });
  testInit("isAbortOnWarningsEnabled",
           []() { (void)axom::slic::isAbortOnWarningsEnabled(); });

  testInit("addStreamToMsgLevel", []() {
    std::string format("[<LEVEL>]: <MESSAGE> \n");
    axom::slic::addStreamToMsgLevel(
      new axom::slic::GenericOutputStream(&std::cout, format),
      axom::slic::message::Info);
  });
  testInit("addStreamToAllMsgLevels", []() {
    std::string format("[<LEVEL>]: <MESSAGE> \n");
    axom::slic::addStreamToAllMsgLevels(
      new axom::slic::GenericOutputStream(&std::cout, format));
  });
}

TEST(slic_uninit, log_methods)
{
  const axom::slic::message::Level lvl {axom::slic::message::Info};
  const std::string fname("no_file");
  const int line = 42;

  testInit("logMessage_A",
           [lvl]() { axom::slic::logMessage(lvl, "an info message"); });
  testInit("logMessage_B", [lvl]() {
    axom::slic::logMessage(lvl, "another info message", "tagB");
  });
  testInit("logMessage_C", [lvl, fname, line]() {
    axom::slic::logMessage(lvl, "info msg", fname, line);
  });
  testInit("logMessage_D", [lvl, fname, line]() {
    axom::slic::logMessage(lvl, "info msg", "tagD", fname, line);
  });

  testInit("logErrorMessage", [fname, line]() {
    axom::slic::logErrorMessage("an error message", fname, line);
  });
  testInit("logWarningMessage", [fname, line]() {
    axom::slic::logWarningMessage("a warning message", fname, line);
  });
}

TEST(slic_uninit, streams)
{
  testInit("flushStreams", []() { axom::slic::flushStreams(); });
  testInit("pushStreams", []() { axom::slic::pushStreams(); });
}